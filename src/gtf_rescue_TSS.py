import pandas as pd
import numpy as np
from src.isoform_classify import IsoformClassifier
from src.gene_grouping import GeneClustering
from src.ISM_filter import TruncationProcessor
from src.get_terminal_sites import TerminalSitesProcessor
import multiprocessing as mp
from src.NNC_NIC_filter import NncNicGraphProcessor
import bisect
from collections import Counter


class GTFRescueTSS:
    def __init__(self, 
                 little_exon_bp=30,
                 mismatch_error_sites_bp=20,
                 mismatch_error_sites_groupfreq_multiple=0.25,
                 exon_excursion_diff_bp=20,
                 fake_exon_group_freq_multiple=0.1,
                 fake_exon_bp=50,
                 
                 two_exon_FreqRatio=1,
                 ism_freqRatio_notrun=4,
                 nic_freqratio_mean=0.1, nic_freqratio_group=0.25,
                 nnc_freqratio_mean=0.1, nnc_freqratio_group=0.25,
                 
                 cluster_group_size=1500, eps=10, min_samples=20,
                 num_processes = 10):
        """
        初始化参数
        """
        self.little_exon_bp = little_exon_bp
        self.mismatch_error_sites_bp=mismatch_error_sites_bp
        self.mismatch_error_sites_groupfreq_multiple=mismatch_error_sites_groupfreq_multiple
        self.exon_excursion_diff_bp=exon_excursion_diff_bp
        self.fake_exon_group_freq_multiple=fake_exon_group_freq_multiple
        self.fake_exon_bp=fake_exon_bp
        
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
    
    @staticmethod
    def rescue_fsm(df_raw, df, ref_anno):
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

        df = pd.concat([df[['Chr','Strand','exonChain','category','TrStart','TrEnd']],
                    df_rescue[['Chr','Strand','exonChain','category','TrStart','TrEnd']]],ignore_index = True).reset_index(drop = True)  
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
    

#     def ism_filter(self, df):
#         """
#         ISM filter
#         """   
#         truncationprocessor = TruncationProcessor(num_processes = self.num_processes)
#         df = truncationprocessor.get_truncation(df)
#         Chr_Strand_mean_Freq = pd.DataFrame(df.groupby(['Chr','Strand'])['frequency'].apply('mean')).reset_index()

#         def remove_ism_for_row(row):
#             threshold_freq = Chr_Strand_mean_Freq[
#                 (Chr_Strand_mean_Freq['Chr'] == row['Chr']) & 
#                 (Chr_Strand_mean_Freq['Strand'] == row['Strand'])
#             ].iloc[0]['frequency']

#             if row['truncation_source'] != 'n' and row['frequency'] <= 2:
#                 return 0
#             if row['truncation_source'] == 'n' and row['frequency'] < threshold_freq * self.ism_freqRatio_notrun:
#                 return 0
#             return 1

#         df_ism = df[df['category'] == 'ISM'].copy()
#         df_ism['remove'] = df_ism.apply(remove_ism_for_row, axis=1)
#         remove_index = df_ism[df_ism['remove'] == 0].index
#         df = df.drop(index=remove_index).copy()

#         # 删除不再需要的列
#         return df.drop(columns=['sourceIso_num', 'trun_source_freq', 'truncation_source'], errors='ignore')



#     def nic_filter(self, df):
#         """
#         NIC filter
#         """
#         Chr_Strand_mean_Freq = pd.DataFrame(df.groupby(['Chr','Strand'])['frequency'].apply('mean')).reset_index()
#         df['Group_freq'] = df.groupby(['Chr', 'Strand', 'Group'])['frequency'].transform('sum')

#         def remove_nic_for_row(row):
#             threshold_freq = Chr_Strand_mean_Freq.loc[
#                 (Chr_Strand_mean_Freq['Chr'] == row['Chr']) & 
#                 (Chr_Strand_mean_Freq['Strand'] == row['Strand']), 'frequency'
#             ].values[0]


#             if row['frequency'] < threshold_freq * self.nic_freqratio_mean or row['frequency'] < row['Group_freq'] * self.nic_freqratio_group:
#                 return 0
#             return 1

#         df_nic = df[df['category'] == 'NIC'].copy()
#         df_nic['remove'] = df_nic.apply(remove_nic_for_row, axis=1)
#         remove_index = df_nic[df_nic['remove'] == 0].index
#         return df.drop(index=remove_index).copy()


#     def nnc_filter(self, df):
#         """
#         NNC filter
#         """
#         Chr_Strand_mean_Freq = pd.DataFrame(df.groupby(['Chr','Strand'])['frequency'].apply('mean')).reset_index()
#         df['Group_freq'] = df.groupby(['Chr','Strand','Group'])['frequency'].transform('sum')

#         def remove_nnc_for_row(row):
#             threshold_freq = Chr_Strand_mean_Freq[
#                 (Chr_Strand_mean_Freq['Chr'] == row['Chr']) & 
#                 (Chr_Strand_mean_Freq['Strand'] == row['Strand'])
#             ].iloc[0]['frequency']

#             if row['frequency'] < threshold_freq * self.nnc_freqratio_mean or row['frequency'] < row['Group_freq'] * self.nnc_freqratio_group:
#                 return 0
#             return 1

#         df_nnc = df[df['category'] == 'NNC'].copy()
#         df_nnc['remove'] = df_nnc.apply(remove_nnc_for_row, axis=1)
#         remove_index = df_nnc[df_nnc['remove'] == 0].index
#         return df.drop(index=remove_index).copy()


    def ism_filter(self, df):
        """
        ISM filter
        """
        df = df[~((df['category'] =='ISM') & (df['Group_freq_ratio'] < 0.25))].copy()
        return df
    
    
    def nic_filter(self, df):
        """
        NIC filter
        """
        df = df[~((df['category'] =='NIC') & (df['Group_freq_ratio'] < 0.05))].copy()
        return df
    
    
    def nnc_filter(self, df):
        """
        NNC filter
        """
        df = df[~((df['category'] =='NNC') & (df['Group_freq_ratio'] < 0.05))].copy()
        return df
    

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
    
    
    def find_nearest_exon(self,ref_sites, query_sites):
        """使用二分查找找到每个 query 对应的最近匹配位点"""
        mapped_sites = []
        for q in query_sites:
            idx = bisect.bisect_left(ref_sites, q)

            # 处理边界情况
            if idx == 0:  # query 在最左边
                mapped_sites.append(ref_sites[0])
            elif idx == len(ref_sites):  # query 在最右边
                mapped_sites.append(ref_sites[-1])
            else:  # query 在中间，选择最近的点
                before = ref_sites[idx - 1]
                after = ref_sites[idx]
                mapped_sites.append(before if abs(q - before) <= abs(q - after) else after)

        return mapped_sites

    
    def is_within_littleExon_range(self,lst, lower, upper):
        """
        判断列表中的位点是否全部落在小exon区域内
        """
        return all(lower <= x <= upper for x in lst)
    
        
    def correct_exonChain(self, row, dict_correct, df1):
        """
        dict_correct ： {整条错误exonChain : 纠正exonChain}
        纠正exonChain： 满足以下即纠正：
            1. 纠正后为FSM
            2. 纠正后为FSM的截断
        """

        FSM_exonChain = df1[
            (df1['Chr'] == row['Chr']) & 
            (df1['Strand'] == row['Strand']) & 
            (df1['Group'] == row['Group']) & 
            (df1['category'] == 'FSM')
        ]['exonChain'].tolist()

        corrected_exonChain = dict_correct.get(row['key'])

        if corrected_exonChain:
            if corrected_exonChain in FSM_exonChain or any(corrected_exonChain in fsm for fsm in FSM_exonChain):
                return corrected_exonChain

        return row['exonChain']


    
    def correct_for_Chr(self, df1):
        """
        单条染色体进行 nnc/nic correct
        """
        df = df1.copy()
        
        # 过滤不需要矫正的group 过滤:不含FSM 或只含FSM和ISM
        df = df.groupby(["Chr", "Strand", "Group"]).filter(lambda g: any(g["category"] == 'FSM'))
        df = df.groupby(["Chr", "Strand", "Group"]).filter(lambda g: "NIC" in g["category"].values or "NNC" in g["category"].values)
    
        df['exonChain2'] = df['exonChain'].apply(lambda x: tuple(map(int, x.split('-'))))
        df['exonChain2_tuple'] = df['exonChain2'].apply(lambda x: [(x[i], x[i+1]) for i in range(0, len(x), 2)])
        
        correct_dict = {}
        for _,df_group in df.groupby(['Chr','Strand','Group']):
            df_ref_exonChain = df_group[df_group['category'] == 'FSM'].copy()
            df_group_NNC_NIC = df_group[df_group['category'].isin(['NNC','NIC'])].copy()

            # 获取小exon路径
            ref_exon = set([(exonC[i],exonC[i+1]) for exonC in df_ref_exonChain['exonChain2'].tolist() for i in range(1,len(exonC)-1,2)])
            ref_little_exon = [tup for tup in ref_exon if abs(tup[1] - tup[0]) < self.little_exon_bp]
            ref_little_exon_path = []
            for little_exon in ref_little_exon:
                for exonC2 in df_ref_exonChain['exonChain2']:
                    # ref的exonChain包含小exon
                    if len(set(little_exon) & set(exonC2)) == 2:
                        indx1 = exonC2.index(little_exon[0])
                        indx2 = exonC2.index(little_exon[1])
                        ref_little_exon_path.append(tuple(exonC2[indx1-1:indx2+2]))

            ref_little_exon_path = set(ref_little_exon_path)
            
            # 是否含有小exon
            correcting_NIC = False
            if len(ref_little_exon_path) > 0:
                correcting_NIC = True

            # 获取所有FSM的位点
            ref_sites = sorted(set(s for exonC in df_ref_exonChain['exonChain2'] for s in exonC))

            for _,row in df_group_NNC_NIC.iterrows():

                ###############
                # NIC校正  小exon跳跃
                if row['category'] == 'NIC' and correcting_NIC:
                    skip_row = False

                    # print('-'*10)
                    # print(ref_little_exon_path)
                    for little_exon_path in ref_little_exon_path:
                        # 小exon跳跃
                        if len(set(little_exon_path[1:-1]) - set(row['exonChain2'])) == 2 and \
                            len(set([little_exon_path[0],little_exon_path[-1]]) & set(row['exonChain2'])) == 2:

                            error_path = row['Chr'] + row['Strand'] + '_' + row['exonChain']
                            correct_path = row['exonChain'].replace('-'.join(list(map(str,[little_exon_path[0],little_exon_path[-1]]))),
                                                                    '-'.join(list(map(str,little_exon_path))))

                            correct_dict[error_path] = correct_path
                            skip_row = True

                            # print('小exon跳跃' ,error_path)
                            # print(correct_path)
                            break

                    if skip_row:
                        continue


                #################
                # NNC校正： 小exon错配，凹陷/凸起，exon偏移， 假exon
                elif row['category'] == 'NNC':
                    skip_row = False

                    # print('-'*10)
                    # print(ref_little_exon_path)
                    # 二分查找法匹配每个位点最近的匹配位点
                    NNCnic_exonChain = list(row['exonChain2'])
                    mapped = self.find_nearest_exon(ref_sites, NNCnic_exonChain)

                    # 和ref不同的sites
                    diff_sites = sorted(set(NNCnic_exonChain) - set(mapped))
                    # 对应的匹配位点
                    diff_sites_ref = [mapped[NNCnic_exonChain.index(diff)] for diff in diff_sites]

                    ##########################
                    # 1. 小exon错配
                    # ref存在小exon
                    if correcting_NIC:
                        for little_exon_path in ref_little_exon_path:
                            # 错误位点位于小exon区域
                            if self.is_within_littleExon_range(diff_sites, little_exon_path[0], little_exon_path[-1]):
                                error_path = row['Chr'] + row['Strand'] + '_' + row['exonChain']

                                # 单端错配
                                if len(diff_sites) == 1:
                                    # 左边错配
                                    if diff_sites_ref[0] == little_exon_path[0]:
                                        correct_path = row['exonChain'].replace('-'.join(list(map(str,[diff_sites[0],little_exon_path[-1]]))),
                                                                                '-'.join(list(map(str,little_exon_path))))
                                        correct_dict[error_path] = correct_path
                                        skip_row = True
                                        
                                        # print('左边错配',error_path)
                                        # print(correct_path)
                                        break


                                    # 右边错配
                                    elif diff_sites_ref[0] == little_exon_path[-1]:
                                        correct_path = row['exonChain'].replace('-'.join(list(map(str,[little_exon_path[0],diff_sites[0]]))),
                                                                                '-'.join(list(map(str,little_exon_path))))
                                        correct_dict[error_path] = correct_path
                                        skip_row = True
                                        
                                        # print('右边错配',error_path)
                                        # print(correct_path)
                                        break


                                # 双端错配
                                elif len(diff_sites) == 2:
                                    if diff_sites_ref[0] == little_exon_path[0] and diff_sites_ref[1] == little_exon_path[-1]:
                                        correct_path = row['exonChain'].replace('-'.join(list(map(str,diff_sites))),
                                                                                '-'.join(list(map(str,little_exon_path))))
                                        correct_dict[error_path] = correct_path
                                        skip_row = True

                                        # print('双端错配',error_path)
                                        # print(correct_path)
                                        break

                    if skip_row:
                        continue

                    ################
                    # 2. 凹陷/凸起  发生在非小exon附近
                    # 和匹配的ref相差 20bp； freq/group freq < 0.25
                    if len(diff_sites) == 1 and abs(diff_sites[0] - diff_sites_ref[0]) < self.mismatch_error_sites_bp:
                        group_freq = df_group[df_group['source'] == 'data']['frequency'].sum()

                        if row['frequency']/group_freq < self.mismatch_error_sites_groupfreq_multiple:
                            error_path = row['Chr'] + row['Strand'] + '_' + row['exonChain']
                            correct_path = '-'.join(list(map(str,mapped)))

                            correct_dict[error_path] = correct_path

                            # print('凹陷凸起错配',error_path)
                            # print(correct_path)
                            continue

                    ################
                    # 3. exon偏移
                    # 偏移量差值小于20bp
                    elif len(diff_sites) == 2 and '-'.join(list(map(str,diff_sites))) in row['exonChain']:
                        excursion = [a - b for a, b in zip(diff_sites, diff_sites_ref)]
                        if ((excursion[0] >0 and excursion[1] >0) or (excursion[0] <0 and excursion[1] <0)) and \
                            abs(excursion[0] - excursion[1]) < self.exon_excursion_diff_bp:

                            error_path = row['Chr'] + row['Strand'] + '_' + row['exonChain']
                            correct_path = row['exonChain'].replace('-'.join(list(map(str,diff_sites))),
                                                                    '-'.join(list(map(str,diff_sites_ref))))

                            correct_dict[error_path] = correct_path

                            # print('exon错配',error_path)
                            # print(correct_path)
                            continue

                    ################
                    # 4. 假exon
                    # 最后一个exon长度小于50, freq/group freq < 0.1
                    if len(diff_sites) >=2:
                        group_freq = df_group[df_group['source'] == 'data']['frequency'].sum()

                        counter = Counter(mapped)
                        duplicates = [k for k, v in counter.items() if v >= 2]

                        # 正链假exon
                        if len(duplicates) == 1 and row['Strand'] == '+':
                            index = mapped.index(duplicates[0]) # 第一个匹配的值
                            tail_sites = mapped[index:]

                            if len(set(tail_sites)) == 1 and tail_sites[0] == duplicates[0] and \
                                row['frequency'] / group_freq < self.fake_exon_group_freq_multiple and \
                                abs(row['exonChain2'][-1] - int(np.mean(row['TrEnd']))) < self.fake_exon_bp:

                                error_path = row['Chr'] + row['Strand'] + '_' + row['exonChain']
                                correct_path = '-'.join(list(map(str,mapped[:index+1])))
                                correct_dict[error_path] = correct_path

                                # print('假exon',error_path)
                                # print(correct_path)


                        # 负链假exon
                        if len(duplicates) == 1 and row['Strand'] == '-':
                            index = mapped.index(duplicates[0]) # 最后一个匹配的值
                            index = len(mapped) - 1 - mapped[::-1].index(duplicates[0])
                            tail_sites = mapped[:index]

                            if len(set(tail_sites)) == 1 and tail_sites[0] == duplicates[0] and \
                                row['frequency'] / group_freq < self.fake_exon_group_freq_multiple and \
                                abs(row['exonChain2'][0] - int(np.mean(row['TrStart']))) < self.fake_exon_bp:

                                error_path = row['Chr'] + row['Strand'] + '_' + row['exonChain']
                                correct_path = '-'.join(list(map(str,mapped[index:])))
                                correct_dict[error_path] = correct_path

                                # print('假exon',error_path)
                                # print(correct_path)
        
        
        # 校正
        df1['key'] = df1['Chr'] + df1['Strand'] + '_' + df1['exonChain']
        df1['exonChain'] = df1.apply(lambda row : self.correct_exonChain(row, correct_dict, df1),axis = 1)
        
        # 整理数据
        df1 = df1[df1['source'] == 'data'].copy()
        df1 = df1[['Chr','Strand','exonChain','TrStart','TrEnd','Group']].copy()
        
        df1 = df1.explode(['TrStart', 'TrEnd'], ignore_index=True)
        df1 = df1.groupby(['Chr', 'Strand', 'exonChain', 'Group'],
                        as_index=False).agg({
                            'TrStart': list,
                            'TrEnd': list
                        }).reset_index(drop=True)

        return df1
                    
            
            
            
    def correct(self, df, ref_anno):
        """
        nnc/nic correct
        """
        df = df[['Chr','Strand','exonChain','category','TrStart','TrEnd']].copy()
        df['key'] = df['Chr'] + df['Strand'] + df['exonChain']
        
        ref_anno = ref_anno.copy()
        ref_anno['key'] = ref_anno['Chr'] + ref_anno['Strand'] + ref_anno['exonChain']
        ref_anno['category'] = 'FSM'
        # 注释中未表达的exonChain添加到df中
        ref_anno = ref_anno[~ref_anno['key'].isin(df['key'])]
        
        df['source'] = 'data'
        ref_anno['source'] = 'ref'
        
        df['frequency'] = df['TrStart'].apply(len)
        ref_anno['frequency'] = 1
        
        df = pd.concat([df,ref_anno],ignore_index = True).reset_index(drop = True)  
        
        # 聚类
        gene_clustering = GeneClustering(num_processes=self.num_processes)
        df = gene_clustering.cluster(df)
        
        # 去除只含ref的group
        df = df.groupby(['Chr', 'Strand', 'Group']).filter(lambda g: not all(g['source'] == 'ref'))
        
        # 多线程
        chr_grouped = [chr_group for _, chr_group in df.groupby('Chr')]

        with mp.Pool(self.num_processes) as pool:
            results = pool.map(self.correct_for_Chr, chr_grouped)

        df = pd.concat(results, ignore_index=True)
        
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
        
        
        #############
        ## 2. NNC/NIC correct
        df = self.correct(df,ref_anno)
        
        isoformclassifier = IsoformClassifier(num_processes = self.num_processes)
        df = isoformclassifier.add_category(df,ref_anno)
        
        #############
        ## 3. 碎片段
        df = self.filter_groups(df,ref_anno)
        
        df['Group_freq'] = df.groupby(['Chr', 'Strand', 'Group'])['frequency'].transform('sum')
        df['Group_freq_ratio'] = df.apply(lambda row: row['frequency'] / row['Group_freq'] if row['Group_freq'] != 0 else 0, axis=1)
        
        #############
        ## 4. ISM filter
        df = self.ism_filter(df)
        
        #############
        ## 5. NIC filter
        df = self.nic_filter(df)
        
        #############
        ## 6. NNC filter
        df = self.nnc_filter(df)
        
        
        ##########################
        ## 7. TS prediction
        df = self.anno_TS_prediction_correct(df, ref_anno)
        
        return df
