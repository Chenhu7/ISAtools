import pandas as pd
import numpy as np
import multiprocessing as mp
import logging

logger = logging.getLogger(__name__)

class GeneClustering:
    """
    gene-level grouping
    """

    def __init__(self, num_processes=None):
        self.num_processes = num_processes

    @staticmethod
    def _single_strand_interval_clustering(df):
        df = df.copy()
        df['start'] = df['SSC'].str.split('-').str[0].astype(np.int32)
        df['end'] = df['SSC'].str.split('-').str[-1].astype(np.int32)

        # 区间聚类（向量化）
        df = df.sort_values(by='start')
        ends = df['end'].cummax().shift(1, fill_value=df['end'].iloc[0])
        group_ids = (df['start'] > ends).cumsum()
        df['Group'] = group_ids

        return df.drop(columns=['start', 'end'])

    @staticmethod
    def _cluster_for_chr(df):
        results = []
        for strand, strand_group in df.groupby('Strand', observed=True):
            clustered = GeneClustering._single_strand_interval_clustering(strand_group)
            results.append(clustered)

        return pd.concat(results, ignore_index=True)

    def cluster(self, df):
        if df.empty:
            return df.assign(Group=pd.Series(dtype=np.int32))

        df = df.astype({
            'Chr': 'category',
            'Strand': 'category',
            'SSC': str
        })

        grouped = [group for _, group in df.groupby('Chr', observed=True)]
        num_processes = min(self.num_processes, len(grouped))

        ctx = mp.get_context("spawn")
        with ctx.Pool(num_processes) as pool:
            results = pool.map(self._cluster_for_chr, grouped)

        clustered_df = pd.concat(results, ignore_index=True)
        
        return clustered_df.reset_index(drop=True)
