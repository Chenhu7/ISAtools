import pandas as pd
import multiprocessing as mp
from functools import partial


class MutuallyExclusiveExonFilter:
    def __init__(self, mutually_exclusive_ref_ratio=0.5, num_processes=10):
        """
        初始化类
        :param mutually_exclusive_ref_ratio: 用于筛选参考的阈值
        :param num_processes: 并行处理的进程数
        """
        self.mutually_exclusive_ref_ratio = mutually_exclusive_ref_ratio
        self.num_processes = num_processes

    def filter_mutually_exclusive_exons_forChr(self, df):
        """
        针对单个Chr分组的数据，筛选互斥exon/junction
        """
        df = df.copy()
        
        # 计算每组的总频率和频率比例
        df['group_freq'] = df.groupby(['Chr', 'Strand', 'Group'])['frequency'].transform('sum')
        df['freq_ratio'] = df['frequency'] / df['group_freq']
        
        # 构建 junction_tuple
        df['junction_tuple'] = df['exonChain2'].apply(lambda x: [(x[i], x[i+1]) for i in range(0, len(x)-1, 2)])
        
        # 分为参考集和潜在互斥集
        df_ref = df[df['freq_ratio'] > self.mutually_exclusive_ref_ratio].copy()
        df_potential = df[df['freq_ratio'] <= self.mutually_exclusive_ref_ratio].copy()
        
        # 如果潜在互斥集为空，直接返回原始数据
        if df_potential.empty:
            return df

        # 构建参考集中的junction组合
        junction_compose = {}
        for _, df_group in df_ref.groupby(['Chr', 'Strand', 'Group']):
            for jt_list in df_group['junction_tuple']:
                for jt in jt_list:
                    other_jt_list = jt_list.copy()
                    other_jt_list.remove(jt)
                    
                    if jt not in junction_compose:
                        junction_compose[jt] = other_jt_list
                    else:
                        junction_compose[jt].extend(other_jt_list)

        # 去重每个 junction 的组合
        junction_compose = {key: list(set(value)) for key, value in junction_compose.items()}

        # 筛选潜在互斥的行
        def filter_row(row):
            """
            判断当前行是否包含互斥的junction
            """
            row_jt_list = row['junction_tuple']
            for row_jt in row_jt_list:
                ref_jt_list = junction_compose.get(row_jt, [])
                row_other_jt_list = [jt for jt in row_jt_list if jt != row_jt]
                if len(set(row_other_jt_list) - set(ref_jt_list)) > 0:
                    return True
            return False

        # 标记需要移除的行
        df_potential['to_remove'] = df_potential.apply(filter_row, axis=1)

        # 删除互斥的行
        df_potential = df_potential[~df_potential['to_remove']].drop(columns=['to_remove'])
        
        # 合并参考集和过滤后的潜在互斥集
        df_result = pd.concat([df_ref, df_potential], ignore_index=True)

        return df_result

    def filter_mxes(self, df):
        """
        主过滤函数，按Chr并行处理
        """
        # 按Chr分组
        dfChr_list = [dfChr for _, dfChr in df.groupby('Chr')]

        # 并行处理
        with mp.Pool(self.num_processes) as pool:
            partial_func = partial(self.filter_mutually_exclusive_exons_forChr)
            results = pool.map(partial_func, dfChr_list)

        # 合并所有分组结果
        df_result = pd.concat(results, ignore_index=True).reset_index(drop=True)

        return df_result