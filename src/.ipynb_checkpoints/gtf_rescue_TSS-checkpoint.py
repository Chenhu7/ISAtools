import pandas as pd
import numpy as np
from src.isoform_classify import IsoformClassifier
from src.gene_grouping import GeneClustering
from src.ISM_filter import TruncationProcessor
from src.get_terminal_sites import TerminalSitesProcessor

class GTFRescueTSS:
    def __init__(self, 
                 nnc_fake_exon_max_diff = 100, nnc_mismatch_max_diff = 10,
                 nic_little_exon_diff = 10,
                 two_exon_FreqRatio=1,
                 
                 ism_freqRatio_notrun=4,
                 nic_freqratio_mean=0.1, nic_freqratio_group=0.25,
                 nnc_freqratio_mean=0.1, nnc_freqratio_group=0.25,
                 
                 cluster_group_size=1500, eps=10, min_samples=20,
                 num_processes = 10):
        """
        初始化参数
        """
        self.nnc_fake_exon_max_diff = nnc_fake_exon_max_diff
        self.nnc_mismatch_max_diff = nnc_mismatch_max_diff
        self.nic_little_exon_diff = nic_little_exon_diff
        self.two_exon_FreqRatio = two_exon_FreqRatio
        
        self.ism_freqRatio_notrun = ism_freqRatio_notrun
        self.nic_freqratio_mean = nic_freqratio_mean
        self.nic_freqratio_group = nic_freqratio_group
        self.nnc_freqratio_mean = nnc_freqratio_mean
        self.nnc_freqratio_group = nnc_freqratio_group
        
        self.cluster_group_size = cluster_group_size
        self.eps = eps
        self.min_samples = min_samples
        self.num_processes = num_processes
    
    
    def rescue_fsm(self, df_raw, df, ref_anno):
        """
        FSM rescue
        """
        ref_anno["key"] = ref_anno["Chr"] + ref_anno["Strand"] + ref_anno["exonChain"]
        df_raw["key"] = df_raw["Chr"] + df_raw["Strand"] + df_raw["exonChain"]
        df["key"] = df["Chr"] + df["Strand"] + df["exonChain"]

        df_rescue = df_raw[
            df_raw["key"].isin(ref_anno["key"]) & ~df_raw["key"].isin(df["key"])
        ].copy()
        df_rescue["category"] = "FSM"

        df = pd.concat([df[['Chr','Strand','exonChain','category','TrStart','TrEnd','key']],
                    df_rescue[['Chr','Strand','exonChain','category','TrStart','TrEnd','key']]],ignore_index = True).reset_index(drop = True)  
        return df
    
    
    def nnc_correct(self,df,ref_anno):
        """
        NNC correct
        """
        # 初始化并预处理数据
        df['frequency'] = df['TrStart'].apply(len)
        df_FSM = df[df['category'] != 'NNC'].copy()
        df_NNC = df[df['category'] == 'NNC'].copy()

        df_NNC['source'] = 'data'
        ref_anno['category'] = 'FSM'
        df_NNC = pd.concat([df_NNC, ref_anno], ignore_index=True).reset_index(drop=True)

        gene_clustering = GeneClustering(num_processes=self.num_processes)
        df_NNC = gene_clustering.cluster(df_NNC)
        
        # 过滤全是FSM的group
        df_NNC = df_NNC.groupby(['Chr', 'Strand', 'Group']).filter(lambda x: not (x['category'] == 'FSM').all())
        df_NNC['exonChain2'] = df_NNC['exonChain'].apply(lambda x: tuple(map(int, x.split('-'))))

        # 矫正和移除逻辑
        correct_dict = {}

        for _, df_group in df_NNC.groupby(['Chr', 'Strand', 'Group']):
            df_ref_exonChain = df_group[df_group['category'] == 'FSM']
            df_NNC_group = df_group[df_group['category'] != 'FSM']

            if len(df_ref_exonChain) > 0:
                ref_sites = set(site for exon in df_ref_exonChain['exonChain2'] for site in exon)

                for index, row in df_NNC_group.iterrows():
                    for ref_eC in df_ref_exonChain['exonChain2']:
                        # 假外显子矫正
                        if (row['Strand'] == '+' and row['exonChain2'][:-2] == ref_eC) or \
                           (row['Strand'] == '-' and row['exonChain2'][2:] == ref_eC):
                            if row['Strand'] == '+' and len(set(row['exonChain2'][-2:]) & ref_sites) == 0:
                                diff = np.mean(row['TrEnd']) - row['exonChain2'][-1]
                                if diff < self.nnc_fake_exon_max_diff:
                                    correct_dict[row['Chr'] + row['Strand'] + row['exonChain']] = '-'.join(map(str, ref_eC))
                                    break
                            if row['Strand'] == '-' and len(set(row['exonChain2'][:2]) & ref_sites) == 0:
                                diff = row['exonChain2'][0] - np.mean(row['TrStart'])
                                if diff < self.nnc_fake_exon_max_diff:
                                    correct_dict[row['Chr'] + row['Strand'] + row['exonChain']] = '-'.join(map(str, ref_eC))
                                    break

                        # 错配矫正
                        diff_sites_row = sorted(set(row['exonChain2']) - set(ref_eC))
                        diff_sites_ref = sorted(set(ref_eC) - set(row['exonChain2']))
                        if len(diff_sites_row) == len(diff_sites_ref) and len(diff_sites_row) <= 2:
                            diff_values = [abs(a - b) for a, b in zip(diff_sites_row, diff_sites_ref)]
                            average_diff = np.mean(diff_values)
                            potential_correct = '-'.join(map(str, ref_eC))
                            if average_diff < self.nnc_mismatch_max_diff and potential_correct in df_FSM['exonChain']:
                                correct_dict[row['Chr'] + row['Strand'] + row['exonChain']] = potential_correct
                                break

        # 应用矫正
        df_NNC['key'] = df_NNC['Chr'] + df_NNC['Strand'] + df_NNC['exonChain']
        df_NNC['exonChain'] = df_NNC.apply(lambda row: correct_dict.get(row['key'], row['exonChain']), axis=1)
        df_NNC['category'] = df_NNC.apply(lambda row: 'FSM' if row['key'] in correct_dict else row['category'], axis=1)
        df_NNC = df_NNC[df_NNC['source'] == 'data']

        # 合并FSM与NNC
        df = pd.concat([df_FSM[['Chr', 'Strand', 'exonChain', 'TrStart', 'TrEnd', 'category']],
                        df_NNC[['Chr', 'Strand', 'exonChain', 'TrStart', 'TrEnd', 'category']]],
                       ignore_index=True).reset_index(drop=True)

        # 数据整理
        df = df.explode(['TrStart', 'TrEnd'], ignore_index=True)
        df = df.groupby(['Chr', 'Strand', 'exonChain', 'category'], as_index=False).agg({
            'TrStart': tuple,
            'TrEnd': tuple
        }).reset_index(drop=True)

        df['frequency'] = df['TrStart'].apply(len)
        return df
    
    
    def nic_correct(self, df, ref_anno):
        """
        nic correct
        """
        # 初始化
        correct_dict = {}  # 存储需要矫正的exonChain {错误exonChain: 正确exonChain}
        df_FSM = df[df['category'] != 'NIC'].copy()
        df_NIC = df[df['category'] == 'NIC'].copy()

        # 标记数据来源
        df_NIC['source'] = 'data'
        ref_anno['source'] = 'ref'
        ref_anno['category'] = 'FSM'

        # 合并数据并进行聚类
        df_NIC = pd.concat([df_NIC, ref_anno], ignore_index=True).reset_index(drop=True)
        gene_clustering = GeneClustering(num_processes=self.num_processes)
        df_NIC = gene_clustering.cluster(df_NIC)

        # 过滤全是FSM的组
        df_NIC = df_NIC.groupby(['Chr', 'Strand', 'Group']).filter(lambda x: not (x['category'] == 'FSM').all())

        # 计算外显子元组和小外显子元组
        df_NIC['exonChain2'] = df_NIC['exonChain'].apply(lambda x: tuple(map(int, x.split('-'))))
        df_NIC['exon_tuple'] = df_NIC['exonChain2'].apply(
            lambda x: [(x[i], x[i + 1]) for i in range(1, len(x) - 2, 2)]
        )
        df_NIC['exon_tuple_filter'] = df_NIC['exon_tuple'].apply(
            lambda x: [exon for exon in x if abs(exon[1] - exon[0]) <= self.nic_little_exon_diff]
        )
        df_NIC['exon_tuple_filter'] = df_NIC['exon_tuple_filter'].apply(lambda x: x if x else np.nan)

        # 开始矫正
        for _, df_group in df_NIC.groupby(['Chr', 'Strand', 'Group']):
            df_ref_exonChain = df_group[df_group['category'] == 'FSM']
            df_NIC_group = df_group[df_group['category'] != 'FSM']

            # 获取参考小外显子
            if len(df_ref_exonChain[~df_ref_exonChain['exon_tuple_filter'].isna()]) > 0:
                little_exons = list(set(exon for l in df_ref_exonChain[~df_ref_exonChain['exon_tuple_filter'].isna()]['exon_tuple_filter'] for exon in l))

                # 对NIC组数据进行矫正
                for _, row in df_NIC_group.iterrows():
                    for exon in little_exons:
                        # 如果当前exonChain不包含该小外显子
                        if exon not in row['exon_tuple']:
                            corrected_exon_chain = '-'.join(map(str, sorted(list(row['exonChain2']) + list(exon))))
                            # 如果矫正后的exonChain在参考数据中
                            if corrected_exon_chain in df_ref_exonChain['exonChain'].tolist():
                                correct_dict[row['Chr'] + row['Strand'] + row['exonChain']] = corrected_exon_chain

        # 应用矫正
        df_NIC['key'] = df_NIC['Chr'] + df_NIC['Strand'] + df_NIC['exonChain']
        df_NIC['exonChain'] = df_NIC.apply(lambda row: correct_dict.get(row['key'], row['exonChain']), axis=1)
        df_NIC['category'] = df_NIC.apply(lambda row: 'FSM' if row['key'] in correct_dict else row['category'], axis=1)
        df_NIC = df_NIC[df_NIC['source'] == 'data']

        # 合并FSM和修正后的NIC数据
        df = pd.concat(
            [df_FSM[['Chr', 'Strand', 'exonChain', 'TrStart', 'TrEnd', 'category']],
             df_NIC[['Chr', 'Strand', 'exonChain', 'TrStart', 'TrEnd', 'category']]],
            ignore_index=True
        ).reset_index(drop=True)

        # 数据整理
        df = df.explode(['TrStart', 'TrEnd'], ignore_index=True)
        df = df.groupby(['Chr', 'Strand', 'exonChain', 'category'], as_index=False).agg({
            'TrStart': tuple,
            'TrEnd': tuple
        }).reset_index(drop=True)
        df['frequency'] = df['TrStart'].apply(len)

        return df
    
    
    def filter_groups(self, df, ref_anno):
        """
        碎片段过滤
        """
        gene_clustering = GeneClustering(num_processes=self.num_processes)
        df = gene_clustering.cluster(df)
        
        df['key'] = df['Chr'] + df['Strand'] + df['exonChain']
        df['frequency'] = df['TrStart'].apply(len)
        df['ref'] = df['key'].isin(ref_anno['key']).astype(int)
        
        # 初始化需要移除的组
        remove_group = []

        # 计算每个 Chr-Strand 的频率均值
        Chr_Strand_mean_Freq = pd.DataFrame(df.groupby(['Chr','Strand'])['frequency'].apply('mean')).reset_index()

        # 按 Chr-Strand-Group 分组
        for (chr_, strand_, group), df_group in df.groupby(['Chr', 'Strand', 'Group']):
            # 如果组内没有 ref == 1，则进行检查
            if 1 not in set(df_group['ref']):
                # 获取当前组的频率均值
                threshold_freq = Chr_Strand_mean_Freq[
                    (Chr_Strand_mean_Freq['Chr'] == chr_) & 
                    (Chr_Strand_mean_Freq['Strand'] == strand_)
                ]['frequency'].iloc[0]

                # 计算当前组的外显子链长度均值
                df_group['eC_len'] = df_group['exonChain'].apply(lambda x: len(x.split('-')))
                exon_chain_mean_length = df_group['eC_len'].mean()

                # 检查过滤条件
                if (
                    df_group['frequency'].sum() < threshold_freq * self.two_exon_FreqRatio and 
                    exon_chain_mean_length == 2 and 
                    len(df_group) == 1
                ):
                    remove_group.append(f"{chr_}{strand_}{group}")

        # 标记需要移除的组
        df['group_key'] = df.apply(lambda row: f"{row['Chr']}{row['Strand']}{row['Group']}", axis=1)

        # 过滤掉标记的组
        df = df[~df['group_key'].isin(remove_group)].copy()

        # 删除临时列
        df = df.drop(columns='group_key')

        return df
    

    def ism_filter(self, df):
        """
        ISM filter
        """
        truncationprocessor = TruncationProcessor(num_processes = self.num_processes)
        df = truncationprocessor.get_truncation(df)
        Chr_Strand_mean_Freq = pd.DataFrame(df.groupby(['Chr','Strand'])['frequency'].apply('mean')).reset_index()

        def remove_ism_for_row(row):
            threshold_freq = Chr_Strand_mean_Freq[
                (Chr_Strand_mean_Freq['Chr'] == row['Chr']) & 
                (Chr_Strand_mean_Freq['Strand'] == row['Strand'])
            ].iloc[0]['frequency']

            if row['truncation_source'] != 'n' and row['frequency'] <= 2:
                return 0
            if row['truncation_source'] == 'n' and row['frequency'] < threshold_freq * self.ism_freqRatio_notrun:
                return 0
            return 1

        df_ism = df[df['category'] == 'ISM'].copy()
        df_ism['remove'] = df_ism.apply(remove_ism_for_row, axis=1)
        remove_index = df_ism[df_ism['remove'] == 0].index
        df = df.drop(index=remove_index).copy()

        # 删除不再需要的列
        return df.drop(columns=['sourceIso_num', 'trun_source_freq', 'truncation_source'], errors='ignore')


    def nic_filter(self, df):
        """
        NIC filter
        """
        Chr_Strand_mean_Freq = pd.DataFrame(df.groupby(['Chr','Strand'])['frequency'].apply('mean')).reset_index()
        df['Group_freq'] = df.groupby(['Chr', 'Strand', 'Group'])['frequency'].transform('sum')

        def remove_nic_for_row(row):
            threshold_freq = Chr_Strand_mean_Freq[
                (Chr_Strand_mean_Freq['Chr'] == row['Chr']) & 
                (Chr_Strand_mean_Freq['Strand'] == row['Strand'])
            ].iloc[0]['frequency']

            if row['frequency'] < threshold_freq * self.nic_freqratio_mean or row['frequency'] < row['Group_freq'] * self.nic_freqratio_group:
                return 0
            return 1

        df_nic = df[df['category'] == 'NIC'].copy()
        df_nic['remove'] = df_nic.apply(remove_nic_for_row, axis=1)
        remove_index = df_nic[df_nic['remove'] == 0].index
        return df.drop(index=remove_index).copy()


    def nnc_filter(self, df):
        """
        NNC filter
        """
        Chr_Strand_mean_Freq = pd.DataFrame(df.groupby(['Chr','Strand'])['frequency'].apply('mean')).reset_index()
        df['Group_freq'] = df.groupby(['Chr','Strand','Group'])['frequency'].transform('sum')

        def remove_nnc_for_row(row):
            threshold_freq = Chr_Strand_mean_Freq[
                (Chr_Strand_mean_Freq['Chr'] == row['Chr']) & 
                (Chr_Strand_mean_Freq['Strand'] == row['Strand'])
            ].iloc[0]['frequency']

            if row['frequency'] < threshold_freq * self.nnc_freqratio_mean or row['frequency'] < row['Group_freq'] * self.nnc_freqratio_group:
                return 0
            return 1

        df_nnc = df[df['category'] == 'NNC'].copy()
        df_nnc['remove'] = df_nnc.apply(remove_nnc_for_row, axis=1)
        remove_index = df_nnc[df_nnc['remove'] == 0].index
        return df.drop(index=remove_index).copy()
    

    def anno_TS_prediction_correct(self, df, ref_anno):
        """
        处理终止位点的预测与参考位点的校正。
        """
        # 预处理 df
        df = df[['Chr', 'Strand', 'exonChain', 'TrStart', 'TrEnd', 'category', 'ref', 'key']].copy()
        df = df.explode(['TrStart', 'TrEnd'], ignore_index=True)
        df = df.groupby(['Chr', 'Strand', 'exonChain', 'category', 'key'], as_index=False).agg({
            'TrStart': tuple,
            'TrEnd': tuple
        }).reset_index(drop=True)

        # TS 预测
        terminalsitesprocessor = TerminalSitesProcessor(cluster_group_size = self.cluster_group_size,
                                                        eps = self.eps,
                                                        min_samples = self.min_samples,
                                                        num_processes = self.num_processes)
        df = terminalsitesprocessor.get_terminal_sites(df)

        # 获取参考位点字典
        ref_sites_dict = ref_anno.groupby('key').apply(
            lambda df: list(zip(df['TrStart'], df['TrEnd']))
        ).to_dict()

        # 添加参考终止位点
        df['terminal_sites'] = df.apply(lambda row: ref_sites_dict.get(row['key'], []), axis=1)

        # 校正函数
        def in_range(val, se, bp):
            return se - bp <= val <= se + bp

        def get_terminal_sites_ref(row):
            """
            获取预测和参考的终止位点，并返回匹配或校正的终止位点。

            参数:
            - row: 包含预测终止位点和参考终止位点信息的行。

            返回:
            - 返回匹配的 TSS 和 TES 位点，或者校正后的位点。
            """
            pred_TS = row['pred_terminal_sites']
            TSS, TES = pred_TS[0], pred_TS[1]
            ref_TS_list = row['terminal_sites']

            # 如果参考的终止位点为空，返回 0
            if pd.isna(ref_TS_list).any():
                return 0

            def get_corrected_site(site, ref_site, ref_list, bp=10):
                """
                校正位点。如果与参考位点范围内匹配，返回原位点。
                否则，根据条件选择最接近的参考位点或最小差值的位点。
                """
                if in_range(site, ref_site, bp=bp):
                    return site

                # 最小差值大于10 直接取ref
                min_diff = min(ref_list, key=lambda x: abs(x - ref_site))

                if min_diff > 10:
                    return ref_site

            TS_ref = []
            for ref in ref_TS_list:
                start, end = ref

                if in_range(TSS, start, bp=10) and in_range(TES, end, bp=10):
                    TS_ref.append((TSS, TES))
                    continue

                corrected_TSS = get_corrected_site(TSS, start, row['TrStart'], bp=10)
                corrected_TES = get_corrected_site(TES, end, row['TrEnd'], bp=10)
                TS_ref.append((corrected_TSS, corrected_TES))

            if len(ref_TS_list) == 0:
                return [(TSS, TES)]

            return TS_ref

        # 应用参考位点校正
        df['pred_terminal_sites_ref'] = df.apply(get_terminal_sites_ref, axis=1)
        df = df.explode(['pred_terminal_sites_ref'], ignore_index=True)
        df = df.drop(columns=['pred_terminal_sites'])
        df.rename(columns={'pred_terminal_sites_ref': 'pred_terminal_sites'}, inplace=True)

        # 清理重复数据
        df = df[['Chr', 'Strand', 'exonChain', 'TrStart', 'TrEnd', 'pred_terminal_sites']].drop_duplicates()

        return df
    
    
    def anno_process(self,df_raw,df,ref_anno):
        # 获取分类
        if 'category' in df.columns:
            df = df.drop(columns=['category'])

        isoformclassifier = IsoformClassifier(num_processes = self.num_processes)
        df = isoformclassifier.add_category(df,ref_anno)

        ################
        ## 1. rescue
        df = self.rescue_fsm(df_raw,df,ref_anno)
        df.to_parquet('SIRV_PB.parquet')
        
        #############
        ## 2. NNC correct
        df = self.nnc_correct(df,ref_anno)
        
        ##########                        
        ## 3. NIC correct
        df = self.nic_correct(df,ref_anno)
        
        #############
        ## 4. 碎片段
        df = self.filter_groups(df,ref_anno)
        
        #############
        ## 5. ISM filter
        df = self.ism_filter(df)
        
        #############
        ## 6. NIC filter
        df = self.nic_filter(df)
        
        #############
        ## 7. NNC filter
        df = self.nnc_filter(df)
        
        ##########################
        ## 8. TS prediction
        df = self.anno_TS_prediction_correct(df, ref_anno)
        
        return df
        
