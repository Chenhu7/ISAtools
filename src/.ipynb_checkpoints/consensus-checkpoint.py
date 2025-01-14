import pandas as pd
import multiprocessing as mp
from functools import partial


import pandas as pd
import multiprocessing as mp
from functools import partial


class ConsensusFilter:
    """
    处理Consensus的过滤逻辑。
    """

    def __init__(self, consensus_bp=10, consensus_multiple=0.1, num_processes=10):
        """
        初始化参数。
        :param consensus_bp: 相似位点的范围。
        :param consensus_multiple: 共识频率的倍数阈值。
        :param num_processes: 使用的进程数量，默认为 CPU 核心数。
        """
        self.consensus_bp = consensus_bp
        self.consensus_multiple = consensus_multiple
        self.num_processes = num_processes

    def _consensus_sites_for_chr(self, df):
        """
        单条 Chr 的 Consensus 筛选。
        :param df: 单个 Chr 的 DataFrame。
        """
        df = df.copy().reset_index(drop=True)
        df_grouped = df.groupby(['Chr', 'Strand', 'Group'])
        groups = [group for _, group in df_grouped]

        remove_index = []
        for df_group in groups:
            # 按频率和 exonChain 长度排序
            df_group['exonChain_list'] = df_group['exonChain'].str.split('-').apply(lambda x: list(map(int, x)))
            df_group['exonChain_length'] = df_group['exonChain_list'].apply(len)
            df_group = df_group.sort_values(by=['Chr', 'Strand', 'Group', 'frequency', 'exonChain_length'],
                                            ascending=[True, True, True, False, False])

            # 构建相似位点字典
            site_info_dict = {}
            for i, row in enumerate(df_group.itertuples()):
                for site in row.exonChain_list:
                    key_site = range(site - self.consensus_bp, site + self.consensus_bp + 1)
                    found = False
                    for key in site_info_dict:
                        if site in key:
                            site_info_dict[key].append((site, i, row.frequency))
                            found = True
                            break
                    if not found:
                        site_info_dict[key_site] = [(site, i, row.frequency)]

            # 获取 Consensus 位点
            l = []
            si = 1
            for key, similar_sites in site_info_dict.items():
                site_group_name = f"site{si}"
                l += [[site_group_name] + list(site_info) for site_info in similar_sites]
                si += 1

            df_similar_sites = pd.DataFrame(l, columns=['site_group', 'site', 'order', 'frequency'])
            df_similar_sites['group_count'] = df_similar_sites.groupby('site_group')['site_group'].transform('count')
            df_similar_sites['group_total_freq'] = df_similar_sites.groupby('site_group')['frequency'].transform('sum')
            df_similar_sites = df_similar_sites[
                ~(df_similar_sites['group_count'] == df_similar_sites['group_total_freq'])]

            grouped_df = df_similar_sites.groupby(['site_group', 'site']).agg(
                count=('site', 'size'),
                min_order=('order', 'min'),
                total_freq=('frequency', 'sum')
            ).reset_index()

            grouped_df = grouped_df.sort_values(['site_group', 'count', 'min_order'], ascending=[True, False, True])
            consensus_sites_list = []
            error_sites = []

            for _, group in grouped_df.groupby('site_group'):
                if len(group) > 1:
                    ref_freq = group.iloc[0]['total_freq']
                    for index, row in group.iloc[1:].iterrows():
                        if row['total_freq'] / ref_freq >= self.consensus_multiple:
                            consensus_sites_list.append(row['site'])
                        else:
                            error_sites.append(row['site'])
                else:
                    consensus_sites_list.append(group.iloc[0]['site'])

            for index, row in df_group.iterrows():
                if len(set(row['exonChain_list']) & set(error_sites)) >= 1:
                    remove_index.append(index)

        return df.drop(remove_index)

    
    def consensus(self, df):
        """
        对输入数据进行 Consensus 筛选。
        :param df: 输入 DataFrame。
        """
        # 按 Chr 分组
        df_grouped = df.groupby('Chr')
        df_chrs = [dfchr for _, dfchr in df_grouped]

        # 多进程处理
        partial_func = partial(self._consensus_sites_for_chr)
        with mp.Pool(self.num_processes) as pool:
            results = pool.map(partial_func, df_chrs)

        # 合并结果
        df_return = pd.concat(results, ignore_index=True)
        return df_return
    