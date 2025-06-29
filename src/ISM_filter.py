import pandas as pd
import numpy as np
from multiprocessing import Pool
from functools import partial

class TruncationProcessor:
    def __init__(self, threshold_logFC_truncation_source_freq=0.5, threshold_logFC_truncation_group_freq=0.5,num_processes = 10):
        self.threshold_logFC_truncation_source_freq = threshold_logFC_truncation_source_freq
        self.threshold_logFC_truncation_group_freq = threshold_logFC_truncation_group_freq
        self.num_processes = num_processes

    def _get_truncation_for_Chr(self, df_clustered_Chr):
        df = df_clustered_Chr.copy()
        df['sourceIso_num'] = 0
        df['trun_source_freq'] = 0

        grouped_dict = {}
        for _, row in df.iterrows():
            key = (row['Chr'], row['Strand'], row['Group'])
            grouped_dict.setdefault(key, []).append(row)

        for index, row in df.iterrows():
            key = (row['Chr'], row['Strand'], row['Group'])
            sourceIso_num = 0
            trun_source_freq = 0
            truncation_source = []
            for other_row in grouped_dict[key]:
                if row['SSC'] != other_row['SSC'] and row['SSC'] in other_row['SSC']:
                    sourceIso_num += 1
                    trun_source_freq += other_row['frequency']
                    truncation_source.append(other_row['SSC'])

            df.at[index, 'sourceIso_num'] = sourceIso_num
            df.at[index, 'trun_source_freq'] = trun_source_freq
            df.at[index, 'truncation_source'] = ','.join(truncation_source) if truncation_source else 'n'
        return df

    def get_truncation(self, df_clustered):
        Chr_groups = df_clustered.groupby(['Chr','Strand'],observed=True)
        Chr_list = [group for _, group in Chr_groups]

        with Pool(self.num_processes) as pool:
            results = pool.map(self._get_truncation_for_Chr, Chr_list)

        return pd.concat(results)

    def filter_truncation_logFC(self, df1):
        df = df1.copy()
        df['logFC_source_freq'] = np.where(
            df['truncation_source'] != 'n',
            np.log10(df['frequency']) - np.log10(df['trun_source_freq'] + 1e-10),
            np.inf
        )
        df['group_freq'] = df.groupby(['Chr', 'Strand', 'Group'],observed=True)['frequency'].transform('sum')
        df['logFC_group_freq'] = np.where(
            df['truncation_source'] != 'n',
            np.log10(df['frequency']) - np.log10(df['group_freq'] + 1e-10),
            np.inf
        )

        def truncation_classify(row):
            if row['logFC_source_freq'] >= self.threshold_logFC_truncation_source_freq and \
               row['logFC_group_freq'] >= self.threshold_logFC_truncation_group_freq:
                return 'no'
            else:
                return 'yes'

        df['truncation'] = df.apply(truncation_classify, axis=1)
        df = df[df['truncation'] == 'no']
        return df.drop(columns=['group_freq', 'sourceIso_num', 'trun_source_freq', 'logFC_source_freq','logFC_group_freq','truncation_source','truncation']).reset_index(drop=True)