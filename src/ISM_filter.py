import pandas as pd
import numpy as np
from multiprocessing import Pool
from functools import partial

class TruncationProcessor:
    def __init__(self, delta_Notrun_bp=100, threshold_logFC_truncation_source_freq=0.5, threshold_logFC_truncation_group_freq=0.5,num_processes = 10):
        """
        初始化 TruncationProcessor
        :param delta_Notrun_bp: 用于过滤的bp阈值
        :param threshold_logFC_truncation_source_freq: logFC 截断/来源freq的阈值
        :param threshold_logFC_truncation_group_freq: logFC 截断/group平均freq的阈值
        """
        self.delta_Notrun_bp = delta_Notrun_bp
        self.threshold_logFC_truncation_source_freq = threshold_logFC_truncation_source_freq
        self.threshold_logFC_truncation_group_freq = threshold_logFC_truncation_group_freq
        self.num_processes = num_processes

    def _get_truncation_for_Chr(self, df_clustered_Chr):
        """
        处理单个Chr的数据并添加截断信息
        """
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
                if row['exonChain'] != other_row['exonChain'] and row['exonChain'] in other_row['exonChain']:
                    sourceIso_num += 1
                    trun_source_freq += other_row['frequency']
                    truncation_source.append(other_row['exonChain'])

            df.at[index, 'sourceIso_num'] = sourceIso_num
            df.at[index, 'trun_source_freq'] = trun_source_freq
            df.at[index, 'truncation_source'] = ','.join(truncation_source) if truncation_source else 'n'
        return df

    def get_truncation(self, df_clustered):
        """
        添加截断信息，按Chr多进程运行
        """
        Chr_groups = df_clustered.groupby('Chr')
        Chr_list = [group for _, group in Chr_groups]

        with Pool(self.num_processes) as pool:
            results = pool.map(self._get_truncation_for_Chr, Chr_list)

        return pd.concat(results)

    def filter_truncation_logFC(self, df1):
        """
        根据 logFC 过滤截断信息
        """
        df = df1.copy()
        df['logFC_source_freq'] = np.where(
            df['truncation_source'] != 'n',
            np.log10(df['frequency']) - np.log10(df['trun_source_freq'] + 1e-10),
            np.inf
        )
        df['group_freq'] = df.groupby(['Chr', 'Strand', 'Group'])['frequency'].transform('sum')
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
        return df.drop(columns=['group_freq', 'sourceIso_num', 'trun_source_freq']).reset_index(drop=True)

    def _get_truncation_type(self, row):
        """
        获取截断类型
        """
        exonChain = row['exonChain'].split('-')
        source = row['truncation_source'].split(',')[0].split('-')
        if exonChain[-1] == source[-1]:
            return 'y5'
        elif exonChain[0] == source[0]:
            return 'y3'
        else:
            return 'ym'

    def _split_truncation(self, row):
        """
        拆分截断的行数据
        """
        truncation_type = row['truncation_type']
        sites, ends = np.array(row['TrStart']), np.array(row['TrEnd'])
        unique_sorted, unique_indices = np.unique(np.sort(sites if truncation_type != 'y3' else ends), return_index=True)
        diff = np.diff(unique_sorted)
        split_indices = np.where(diff > self.delta_Notrun_bp)[0] + 1
        indices = np.insert(split_indices, 0, 0)
        indices = np.append(indices, len(unique_sorted))
        groups = [unique_sorted[indices[i]:indices[i + 1]] for i in range(len(indices) - 1)]
        original_groups = [[] for _ in groups]
        original_ends = [[] for _ in groups]
        value_to_group = {}
        for idx, group in enumerate(groups):
            for value in group:
                value_to_group[value] = idx
        for idx, value in enumerate(sites if truncation_type != 'y3' else ends):
            group_idx = value_to_group[value]
            original_groups[group_idx].append(value)
            original_ends[group_idx].append(ends[idx] if truncation_type != 'y3' else sites[idx])
        return original_groups, original_ends

    def process_group(self, group):
        if len(group[group['truncation_source'] != 'n']) == 0:
            group['truncation_type'] = np.nan
            return group
        df1 = group[group['truncation_source'] == 'n'].copy()
        trun = group[group['truncation_source'] != 'n'].copy()
        trun['truncation_type'] = trun.apply(self._get_truncation_type, axis=1)
        expanded_rows = []
        for _, row in trun.iterrows():
            original_groups, original_ends = self._split_truncation(row)
            for group, end_group in zip(original_groups, original_ends):
                new_row = row.copy()
                new_row['TrStart'] = group
                new_row['TrEnd'] = end_group
                expanded_rows.append(new_row)
        expanded_df = pd.DataFrame(expanded_rows)
        expanded_df['frequency'] = expanded_df['TrStart'].apply(len)
        expanded_df = expanded_df.reset_index(drop=True)
        return pd.concat([df1, expanded_df], ignore_index=True).reset_index(drop=True)
    
    def filter_truncation_delta(self, df):
        Chr_groups = df.groupby('Chr')
        df_chrs = [group for _, group in Chr_groups if not group.empty]

        with Pool(10) as pool:
            results = pool.map(self.process_group, df_chrs)

        # Ensure all items in results are DataFrames
        return pd.concat(results, ignore_index=True)