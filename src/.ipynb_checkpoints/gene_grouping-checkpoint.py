import pandas as pd
import multiprocessing as mp


class GeneClustering:
    """
    Gene聚类
    """
    def __init__(self, num_processes=10):
        self.num_processes = num_processes

    @staticmethod
    def _single_strand_interval_clustering(df):
        """
        单条 Chr 单条 Strand 单链区间聚类
        """
        df = df.copy()
        df['start'] = df['exonChain'].str.split('-').str[0].astype(int)
        df['end'] = df['exonChain'].str.split('-').str[-1].astype(int)
        df = df.sort_values(by='start').reset_index(drop=True)

        group = 0
        group_list = []
        end = df.loc[0, 'end']

        for _, row in df.iterrows():
            if row['start'] <= end:
                group_list.append(group)
                end = max(end, row['end'])
            else:
                group += 1
                end = row['end']
                group_list.append(group)

        df['Group'] = group_list
        return df.drop(columns=['start', 'end'])

    @staticmethod
    def _cluster_for_chr(df):
        """
        处理单条 Chr：对正负链分别聚类
        """
        results = []
        for strand, strand_group in df.groupby('Strand'):
            clustered = GeneClustering._single_strand_interval_clustering(strand_group)
            results.append(clustered)

        return pd.concat(results, ignore_index=True)

    def cluster(self, df):
        """
        按 Chr 多进程聚类
        :param df: 输入的 DataFrame
        """
        grouped = [group for _, group in df.groupby('Chr')]

        with mp.Pool(self.num_processes) as pool:
            results = pool.map(self._cluster_for_chr, grouped)

        clustered_df = pd.concat(results, ignore_index=True)
        return clustered_df.reset_index(drop=True)
