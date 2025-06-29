from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
import pandas as pd
from src.exonChain_assign import ExonChainAssign

class Quantification:
    def __init__(self, assign_delta=5, assign_exonChain_weight=1, truncation_weight=0.7, num_processes = 10):
        self.assign_delta = assign_delta
        self.assign_exonChain_weight = assign_exonChain_weight
        self.truncation_weight = truncation_weight
        self.num_processes = num_processes

    def assign_quantify_exonChain_single(self, group_df, exonChain_count):
        # 全匹配过滤
        exonChain_count_sub = exonChain_count[
            (~exonChain_count['exonChain'].isin(group_df['exonChain'])) &
            (exonChain_count['Chr'] == group_df.iloc[0]['Chr']) &
            (exonChain_count['Strand'] == group_df.iloc[0]['Strand'])
        ].copy()
        exonChain_count_dict = dict(zip(exonChain_count_sub['exonChain'], exonChain_count_sub['frequency']))

        assigner = ExonChainAssign(delta=self.assign_delta)
        
        group_df = group_df.copy()
        
        if self.assign_exonChain_weight > 0:
            assigner.build_index(exonChain_count_sub['exonChain'])
            group_df['assign'] = group_df['exonChain'].apply(lambda x: assigner.assign_exonChain(x))
            group_df['quantification_math'] = group_df['assign'].apply(lambda x: sum([exonChain_count_dict.get(i, 0) for i in x]))
        
        else:
            group_df['assign'] = 'none'
            group_df['quantification_math'] = 0
            pass

        # 截断匹配
        assigner.build_index_fortrun(exonChain_count_sub['exonChain'])
        group_df['assign_trun'] = group_df['exonChain'].apply(lambda x: assigner.assign_truncation(x))

        # 构建 dict_source_freq 和 dict_truncation_freq
        dict_source_freq = defaultdict(dict)
        dict_truncation_freq = {}

        for _, row in group_df.iterrows():
            for trun in row['assign_trun']:
                dict_source_freq[trun][row['exonChain']] = row['frequency']
                if trun not in dict_truncation_freq and trun in exonChain_count_dict:
                    dict_truncation_freq[trun] = exonChain_count_dict[trun]

        # 截断占比
        dict_source_ratio = {
            trun: {src: freq / sum(sources.values()) for src, freq in sources.items()}
            for trun, sources in dict_source_freq.items()
        }

        # 计算 frequency_trun
        group_df['frequency_trun'] = group_df.apply(
            lambda row: [round(dict_source_ratio[trun].get(row['exonChain'], 0) * dict_truncation_freq.get(trun, 0), 2)
                         for trun in row['assign_trun']], axis=1)
        group_df['quantification_trun'] = group_df['frequency_trun'].apply(sum)

        group_df['quantification'] = (
            group_df['frequency'] +
            group_df['quantification_math'] * self.assign_exonChain_weight +
            group_df['quantification_trun'] * self.truncation_weight
        )

        return group_df.drop(columns=['assign', 'quantification_math', 'assign_trun', 'frequency_trun', 'quantification_trun','frequency'])

    def assign_quantify_exonChain(self, result, exonChain_count):
        groups = list(result.groupby(['Chr', 'Strand']))
        results = []

        with ThreadPoolExecutor(max_workers=self.num_processes) as executor:
            futures = {
                executor.submit(self.assign_quantify_exonChain_single, group_df, exonChain_count): (chr_, strand_)
                for (chr_, strand_), group_df in groups
            }

            for future in as_completed(futures):
                group_result = future.result()
                results.append(group_result)

        return pd.concat(results, ignore_index=True)
