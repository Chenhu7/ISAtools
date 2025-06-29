import numpy as np
import pandas as pd
from numba import jit
from joblib import Parallel, delayed
import multiprocessing as mp
from typing import List, Dict, Any

class ExonChainCollapse:
    """
    Efficient exon chain collapse with parallel processing.
    Collapses similar transcripts within a specified base pair difference.
    """
    
    def __init__(self, collapse_max_diff: int = 10, num_processes: int = 10):
        """
        Args:
            collapse_max_diff: Maximum allowed base pair difference for collapsing (default: 10)
            num_processes: Number of parallel processes to use (default: 10)
        """
        self.collapse_max_diff = collapse_max_diff
        self.num_processes = num_processes
    
    @staticmethod
    @jit(nopython=True)
    def _calc_diff_matrix(coords_list: np.ndarray, max_diff: int) -> np.ndarray:
        n = len(coords_list)
        diff_matrix = np.zeros((n, n), dtype=np.bool_)
        for i in range(n):
            for j in range(i+1, n):
                if len(coords_list[i]) == len(coords_list[j]):
                    total_diff = np.sum(np.abs(coords_list[i] - coords_list[j]))
                    diff_matrix[i, j] = total_diff <= max_diff  # 关键修复：用 <=
                    diff_matrix[j, i] = diff_matrix[i, j]
        return diff_matrix
    
    def _cluster_group(self, coords_array: np.ndarray) -> List[List[int]]:
        if len(coords_array) == 0:
            return []

        merge_mask = self._calc_diff_matrix(coords_array, self.collapse_max_diff)
        n = len(coords_array)
        visited = [False] * n
        clusters = []

        for i in range(n):
            if not visited[i]:
                # 使用BFS找到所有连通节点
                queue = [i]
                cluster = []
                while queue:
                    node = queue.pop(0)
                    if not visited[node]:
                        visited[node] = True
                        cluster.append(node)
                        # 添加所有未访问且可合并的邻居
                        neighbors = np.where(merge_mask[node])[0].tolist()
                        queue.extend([nbr for nbr in neighbors if not visited[nbr]])
                clusters.append(cluster)

        return clusters

    def _cluster_transcripts(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Cluster similar transcripts across the entire DataFrame using TOTAL difference.
        
        Args:
            df: Input DataFrame with exon chain data
            
        Returns:
            DataFrame with added 'cluster' column
        """
        df['sites_num'] = df['exonChain'].str.count('-') + 1
        
        grouped = df.groupby(['Chr', 'Strand', 'Group', 'sites_num'])
        
        def process_group(name: Any, group: pd.DataFrame) -> pd.DataFrame:
            coords_list = [np.array(list(map(int, s.split('-')))) for s in group['exonChain']]
            clusters = self._cluster_group(np.array(coords_list))
            
            cluster_ids = np.zeros(len(group), dtype=int)
            for clus_id, cluster in enumerate(clusters, 1):
                cluster_ids[cluster] = clus_id
                
            return group.assign(cluster=cluster_ids)
        
        results = Parallel(n_jobs=self.num_processes)(
            delayed(process_group)(name, group) for name, group in grouped
        )
        
        return pd.concat(results)
    
    @staticmethod
    def _contains_non_gt_ag(junction_str: str) -> bool:
        """
        Check if junction contains non GT-AG splice sites.
        
        Args:
            junction_str: String of junction pairs
            
        Returns:
            True if any non GT-AG pair is found
        """
        if pd.isna(junction_str):
            return False
            
        return any(pair.upper() != "GT-AG" for pair in junction_str.split(','))
    
    def _collapse_for_chr(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Collapse transcripts for a single chromosome.
        
        Args:
            df: DataFrame for one chromosome
            
        Returns:
            Collapsed DataFrame
        """
        index_collapsed = []
        exonChain_collapsed = []
        junction_collapsed = []
        
        # Group by cluster and process
        for _, df_group in df.groupby(['Chr', 'Strand', 'Group', 'sites_num', 'cluster']):
            if len(df_group) > 1:
                ref_row = df_group.iloc[0]
                for index, row in df_group.iloc[1:].iterrows():
                    if row['frequency'] <= 1 or self._contains_non_gt_ag(row['junction']):
                        index_collapsed.append(index)
                        exonChain_collapsed.append(ref_row['exonChain'])
                        junction_collapsed.append(ref_row['junction'])
        
        # Apply collapses
        df.loc[index_collapsed, 'exonChain'] = exonChain_collapsed
        df.loc[index_collapsed, 'junction'] = junction_collapsed
        
        # Explode and regroup
        df = df.explode(['TrStart', 'TrEnd'], ignore_index=True)
        df = df.groupby(['Chr', 'Strand', 'Group', 'exonChain', 'junction'], as_index=False).agg({
            'TrStart': list,
            'TrEnd': list
        })
        df['frequency'] = df['TrStart'].str.len()
        
        return df
    
    
    def collapse(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Main collapse method with parallel processing by chromosome.
        
        Args:
            df: Input DataFrame with transcript data
            
        Returns:
            Collapsed DataFrame with similar transcripts merged
        """
        # Cluster similar transcripts
        df = self._cluster_transcripts(df)
        
        # Sort by key features (most frequent first within clusters)
        df = df.sort_values(
            ['Chr', 'Strand', 'Group', 'sites_num', 'cluster', 'frequency'],
            ascending=[True, True, True, True, True, False]
        )
        
        # Split by chromosome for parallel processing
        grouped = [group for _, group in df.groupby('Chr')]
        
        with mp.Pool(self.num_processes) as pool:
            results = pool.map(self._collapse_for_chr, grouped)
            
        return pd.concat(results, ignore_index=True)



# class ExonChaincollapse:
#     """
#     ExonChain collapse
#     """
#     def __init__(self, num_processes=10):
#         self.num_processes = num_processes

#     @staticmethod
#     def contains_non_gt_ag(junction_str):
#         """检查 junction 是否包含非 GT-AG 的碱基对（不区分大小写）"""
#         if pd.isna(junction_str):  # 处理可能的 NaN 值
#             return False
#         pairs = junction_str.split(',')
#         for pair in pairs:
#             if pair.upper() != "GT-AG":
#                 return True
#         return False
    
    
#     def collapse_for_chr(self,df):
#         """
#         处理单条 Chr：对正负链分别collapse
#         """
#         # 预处理
#         df['exonChain2'] = df['exonChain'].apply(lambda x : list(map(int,x.split('-'))))
#         df['exonChain_length'] = df['exonChain2'].apply(lambda x: sum([x[i+1] - x[i] for i in range(len(x)-1)]))

#         index_collapsed = []
#         exonChain_collapsed = []
#         junction_collapsed = []
#         df = df.sort_values(['Chr', 'Strand','Group','exonChain_length','frequency'], ascending=[True,True,True,False,False])
#         for _, df_group in df.groupby(['Chr', 'Strand', 'Group']):
#             for _,df_exonChainGroup in df_group.groupby('exonChain_length'):
#                 if len(df_exonChainGroup) > 1:
#                     ref_exonChain = df_exonChainGroup.iloc[0]['exonChain']
#                     ref_junction = df_exonChainGroup.iloc[0]['junction']
#                     for index,row in df_exonChainGroup[1:].iterrows():
#                         if row['frequency'] <= 2 or self.contains_non_gt_ag(row['junction']):
#                             index_collapsed.append(index)
#                             exonChain_collapsed.append(ref_exonChain)
#                             junction_collapsed.append(ref_junction)
                            
#         df.loc[index_collapsed, 'exonChain'] = exonChain_collapsed
#         df.loc[index_collapsed, 'junction'] = junction_collapsed
        
#         df = df.explode(['TrStart', 'TrEnd'], ignore_index=True)
#         df = df.groupby(['Chr', 'Strand', 'Group', 'exonChain','junction'],as_index=False).agg({'TrStart': list,'TrEnd': list}).reset_index(drop=True)
#         df['frequency'] = df['TrStart'].apply(len)

#         return df

    
#     def collapse(self, df):
#         """
#         按 Chr 多进程collapse
#         :param df: 输入的 DataFrame
#         """
#         grouped = [group for _, group in df.groupby('Chr')]

#         with mp.Pool(self.num_processes) as pool:
#             results = pool.map(self.collapse_for_chr, grouped)

#         collapsed_df = pd.concat(results, ignore_index=True)
#         return collapsed_df.reset_index(drop=True)

