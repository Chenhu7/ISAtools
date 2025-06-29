import pandas as pd
import multiprocessing as mp
from functools import partial

class ConsensusFilter:
    """
    consensus
    """

    def __init__(self, consensus_bp=10, consensus_multiple=0.1, num_processes=10):
        self.consensus_bp = consensus_bp
        self.consensus_multiple = consensus_multiple
        self.num_processes = num_processes

    def _consensus_sites_for_chr(self, df):
        df = df.copy().reset_index(drop=True)
        df_grouped = df.groupby(['Chr', 'Strand', 'Group'], observed=True)
        groups = [group for _, group in df_grouped]

        remove_index = []
        for df_group in groups:
            df_group['SSC_list'] = df_group['SSC'].str.split('-').apply(lambda x: list(map(int, x)))
            df_group['SSC_length'] = df_group['SSC_list'].apply(len)
            df_group = df_group.sort_values(by=['Chr', 'Strand', 'Group', 'frequency', 'SSC_length'],
                                            ascending=[True, True, True, False, False])

            site_info_dict = {}
            for i, row in enumerate(df_group.itertuples()):
                for site in row.SSC_list:
                    key_site = range(site - self.consensus_bp, site + self.consensus_bp + 1)
                    found = False
                    for key in site_info_dict:
                        if site in key:
                            site_info_dict[key].append((site, i, row.frequency))
                            found = True
                            break
                    if not found:
                        site_info_dict[key_site] = [(site, i, row.frequency)]

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

            grouped_df = grouped_df.groupby('site_group').filter(lambda x: len(x) > 1)
            grouped_df = grouped_df.sort_values(['site_group', 'total_freq','count'], ascending=[True, False,False])
            error_sites = []

            for _, group in grouped_df.groupby('site_group'):
                ref_freq = group.iloc[0]['total_freq']
                for index, row in group.iloc[1:].iterrows():
                    if row['total_freq'] / ref_freq < self.consensus_multiple:
                        error_sites.append(row['site'])

            for index, row in df_group.iterrows():
                if len(set(row['SSC_list']) & set(error_sites)) >= 1:
                    remove_index.append(index)

        return df.drop(remove_index)

    def consensus(self, df):
        df_grouped = df.groupby('Chr', observed=True)
        df_chrs = [dfchr for _, dfchr in df_grouped]

        partial_func = partial(self._consensus_sites_for_chr)
        with mp.Pool(self.num_processes) as pool:
            results = pool.map(partial_func, df_chrs)
            
        df_return = pd.concat(results, ignore_index=True)
        return df_return
