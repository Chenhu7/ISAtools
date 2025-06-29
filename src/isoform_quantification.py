from collections import defaultdict
import pandas as pd
import numpy as np
import os
from src.SSC_assign import SSCAssign
from src.gene_grouping import GeneClustering
import multiprocessing as mp

class IsoformQuantifier:
    def __init__(self, assign_delta=0, assign_SSC_weight=0, truncation_weight=1, num_processes=10):
        self.assign_delta = assign_delta
        self.assign_SSC_weight = assign_SSC_weight
        self.truncation_weight = truncation_weight
        self.num_processes = num_processes

    def assign_quantification(self, result, SSC_count_path):
        SSC_count = pd.read_parquet(SSC_count_path, columns=['Chr', 'Strand', 'SSC', 'frequency'])

        if len(SSC_count) > len(result) * 5:
            return self._assign_method(result, SSC_count)
        else:
            return self._Trun_method(result, SSC_count)

    def _assign_method(self, result, SSC_count):
        result = result.copy()

        ssc_dict = {
            (row.Chr, row.Strand, row.SSC): row.frequency
            for row in SSC_count.itertuples(index=False)
        }
        result['frequency_orig'] = result['frequency']

        result['SSC_frequency'] = result.apply(
            lambda x: ssc_dict.get((x.Chr, x.Strand, x.SSC), np.nan),
            axis=1
        )
        result['frequency_ratio'] = result['frequency'] / result.groupby(['Chr', 'Strand', 'SSC'])['frequency'].transform('sum')
        result['frequency'] = np.where(
            result['SSC_frequency'].isna(),
            result['frequency_orig'],
            result['frequency_ratio'] * result['SSC_frequency']
        )

        result = result.drop(columns=['frequency_ratio', 'SSC_frequency', 'frequency_orig'])

        groups = list(result.groupby(['Chr', 'Strand'], observed=True))
        task_args = [
            (group_df.copy(), SSC_count, self.assign_delta, self.assign_SSC_weight, self.truncation_weight)
            for (_, _), group_df in groups
        ]
        with mp.Pool(processes=self.num_processes) as pool:
            results = pool.map(self.assign_wrapper, task_args)

        return pd.concat(results, ignore_index=True)

    def _Trun_method(self, df, SSC_count):
        df = df.copy()
        SSC_count = SSC_count[~SSC_count['SSC'].isin(df['SSC'])].copy()

        ssc_dict = {
            (row.Chr, row.Strand, row.SSC): row.frequency
            for row in SSC_count.itertuples(index=False)
        }
        df['frequency_orig'] = df['frequency']
        df['SSC_frequency'] = df.apply(
            lambda x: ssc_dict.get((x.Chr, x.Strand, x.SSC), np.nan),
            axis=1
            )
        df['frequency_ratio'] = df['frequency'] / df.groupby(['Chr', 'Strand', 'SSC'])['frequency'].transform('sum')
        df['frequency'] = np.where(
            df['SSC_frequency'].isna(),
            df['frequency_orig'],
            df['frequency_ratio'] * df['SSC_frequency']
            )

        df = df.drop(columns=['frequency_ratio', 'SSC_frequency', 'frequency_orig'])

        df['type'] = 'F'
        SSC_count['type'] = 'D'
        df['type'] = df['type'].astype('category')
        SSC_count['type'] = SSC_count['type'].astype('category')

        df_combined = pd.concat([df, SSC_count])
        gene_clustering = GeneClustering(num_processes=self.num_processes)
        df_combined = gene_clustering.cluster(df_combined)

        grouped = df_combined.groupby(['Chr', 'Strand'], observed=True)
        task_args = [(chr_, strand, group.copy(), df.copy()) for (chr_, strand), group in grouped]

        with mp.Pool(processes=self.num_processes) as pool:
            results = pool.map(self._assign_wrapper, task_args)

        source_SSC_ISMdict = defaultdict(float)
        for result in results:
            for k, v in result.items():
                source_SSC_ISMdict[k] += v

        df['key'] = df['Chr'].astype(str) + df['Strand'].astype(str) + df['SSC'] + df['TrStart'].astype(str) + df['TrEnd'].astype(str)
        df['quantification'] = df.apply(lambda row: row['frequency'] + source_SSC_ISMdict.get(row['key'], 0), axis=1)
        df['quantification'] = df['quantification'].round().astype(int)
        return df.drop(columns=['frequency', 'type', 'key'])

    @staticmethod
    def assign_wrapper(args):
        group_df, SSC_count, assign_delta, assign_SSC_weight, truncation_weight = args

        assigner = SSCAssign(delta=assign_delta)
        
        SSC_count_sub = SSC_count[
            (~SSC_count['SSC'].isin(group_df['SSC'])) &
            (SSC_count['Chr'] == group_df.iloc[0]['Chr']) &
            (SSC_count['Strand'] == group_df.iloc[0]['Strand'])
        ]
        SSC_count_dict = dict(zip(SSC_count_sub['SSC'], SSC_count_sub['frequency']))

        if assign_SSC_weight > 0:
            assigner.build_index(SSC_count_sub['SSC'])
            group_df['assign'] = group_df['SSC'].apply(lambda x: assigner.assign_SSC(x))
            group_df['quantification_math'] = group_df['assign'].apply(lambda x: sum([SSC_count_dict.get(i, 0) for i in x]))
        else:
            group_df['assign'] = 'none'
            group_df['quantification_math'] = 0

        assigner.build_index_fortrun(SSC_count_sub['SSC'])
        group_df['assign_trun'] = group_df['SSC'].apply(lambda x: assigner.assign_truncation(x))

        dict_source_freq = defaultdict(dict)
        dict_truncation_freq = {}
        for _, row in group_df.iterrows():
            for trun in row['assign_trun']:
                dict_source_freq[trun][row['SSC']] = row['frequency']
                if trun not in dict_truncation_freq and trun in SSC_count_dict:
                    dict_truncation_freq[trun] = SSC_count_dict[trun]

        dict_source_ratio = {
            trun: {src: freq / sum(sources.values()) for src, freq in sources.items()}
            for trun, sources in dict_source_freq.items()
        }

        group_df['frequency_trun'] = group_df.apply(
            lambda row: [round(dict_source_ratio[trun].get(row['SSC'], 0) * dict_truncation_freq.get(trun, 0), 2)
                        for trun in row['assign_trun']], axis=1)
        group_df['quantification_trun'] = group_df['frequency_trun'].apply(sum)

        group_df['quantification'] = (
            group_df['frequency'] +
            group_df['quantification_math'] * assign_SSC_weight +
            group_df['quantification_trun'] * truncation_weight
        ).round(0).astype('int32')

        return group_df.drop(columns=[
            'assign', 'quantification_math',
            'assign_trun', 'frequency_trun', 'quantification_trun', 'frequency'
        ])

    @staticmethod
    def _assign_wrapper(args):
        chr_, strand, df_combined_chunk, df_filtered = args
        f_df = df_combined_chunk[df_combined_chunk['type'] == 'F']
        d_df = df_combined_chunk[df_combined_chunk['type'] == 'D']
        f_grouped = {}

        for _, row in f_df.iterrows():
            key = row['Group']
            f_grouped.setdefault(key, []).append(row)

        result_dict = {}
        for _, row in d_df.iterrows():
            group = row['Group']
            candidates = f_grouped.get(group, [])
            ssc = row['SSC']
            sources = [r['SSC'] for r in candidates if ssc != r['SSC'] and ssc in r['SSC']]

            if not sources:
                continue

            df_source = df_filtered[
                (df_filtered['Chr'] == chr_) &
                (df_filtered['Strand'] == strand) &
                (df_filtered['SSC'].isin(sources))
            ]

            if len(sources) == 1 and len(df_source) == 1:
                key = chr_ + strand + sources[0] + str(df_source.iloc[0]['TrStart']) + str(df_source.iloc[0]['TrEnd'])
                result_dict[key] = result_dict.get(key, 0) + row['frequency']
        return result_dict
