import pandas as pd
import numpy as np
from src.isoform_classify import IsoformClassifier
from src.gene_grouping import GeneClustering
from src.get_terminal_sites import TerminalSitesProcessor
import multiprocessing as mp
import bisect
from collections import Counter


class GTFRescueFiltering:
    def __init__(self, 
                 little_exon_bp=30,
                 mismatch_error_sites_bp=20,
                 mismatch_error_sites_groupfreq_multiple=0.25,
                 exon_excursion_diff_bp=20,
                 fake_exon_group_freq_multiple=0.1,
                 fake_exon_bp=50,
                 
                 ism_freqRatio_notrun=4,
                 
                 cluster_group_size=1500, eps=10, min_samples=20,
                 num_processes = 10):
        self.little_exon_bp = little_exon_bp
        self.mismatch_error_sites_bp=mismatch_error_sites_bp
        self.mismatch_error_sites_groupfreq_multiple=mismatch_error_sites_groupfreq_multiple
        self.exon_excursion_diff_bp=exon_excursion_diff_bp
        self.fake_exon_group_freq_multiple=fake_exon_group_freq_multiple
        self.fake_exon_bp=fake_exon_bp
        
        self.ism_freqRatio_notrun = ism_freqRatio_notrun
        
        self.cluster_group_size = cluster_group_size
        self.eps = eps
        self.min_samples = min_samples
        self.num_processes = num_processes
    
    @staticmethod
    def rescue_fsm(df_raw, df, ref_anno):
        ref_anno["key"] = ref_anno["Chr"].astype(str) + ref_anno["Strand"].astype(str) + ref_anno["SSC"].astype(str)
        df_raw["key"] = df_raw["Chr"].astype(str) + df_raw["Strand"].astype(str) + df_raw["SSC"].astype(str)
        df["key"] = df["Chr"].astype(str) + df["Strand"].astype(str) + df["SSC"].astype(str)

        df_rescue = df_raw[
            df_raw["key"].isin(ref_anno["key"]) & ~df_raw["key"].isin(df["key"])
        ].copy()
        df_rescue["category"] = "FSM"

        df = pd.concat([df[['Chr','Strand','SSC','category','TrStart_reads','TrEnd_reads']],
                    df_rescue[['Chr','Strand','SSC','category','TrStart_reads','TrEnd_reads']]],ignore_index = True).reset_index(drop = True)  
        return df
    
    def filter_groups(self, df, ref_anno):
        gene_clustering = GeneClustering(num_processes=self.num_processes)
        df = gene_clustering.cluster(df)

        df['frequency'] = df['TrStart_reads'].apply(len)

        df['key'] = df['Chr'].astype(str) + df['Strand'].astype(str) + df['SSC'].astype(str)
        ref_keys = set(ref_anno['key'])
        df['ref'] = df['key'].isin(ref_keys).astype(int)

        freq_threshold = df['frequency'].quantile(0.25)

        group_stats = df.groupby(['Chr', 'Strand', 'Group'], observed=True).agg(
            any_ref=('ref', 'any'),
            total_freq=('frequency', 'sum')
        )

        groups_to_keep = group_stats[
            (group_stats['any_ref']) | (group_stats['total_freq'] >= freq_threshold)
        ].index

        df = df.set_index(['Chr', 'Strand', 'Group'])
        df = df.loc[groups_to_keep].reset_index()

        return df.drop(columns=['key', 'ref'])
        
    def ism_filter(self, df):
        """
        ISM filter
        """
        df = df[~((df['category'] =='ISM') & (df['Group_freq_ratio'] < self.ism_freqRatio_notrun))]
        return df
    
    def find_nearest_exon(self,ref_sites, query_sites):
        mapped_sites = []
        for q in query_sites:
            idx = bisect.bisect_left(ref_sites, q)

            if idx == 0:
                mapped_sites.append(ref_sites[0])
            elif idx == len(ref_sites):
                mapped_sites.append(ref_sites[-1])
            else:
                before = ref_sites[idx - 1]
                after = ref_sites[idx]
                mapped_sites.append(before if abs(q - before) <= abs(q - after) else after)

        return mapped_sites

    def is_within_littleExon_range(self,lst, lower, upper):
        return all(lower <= x <= upper for x in lst)
    
    def correct_ssc(self, row, dict_correct, df1):

        FSM_SSC = df1[
            (df1['Chr'] == row['Chr']) & 
            (df1['Strand'] == row['Strand']) & 
            (df1['Group'] == row['Group']) & 
            (df1['category'] == 'FSM')
        ]['SSC'].tolist()

        corrected_SSC = dict_correct.get(row['key'])

        if corrected_SSC:
            if corrected_SSC in FSM_SSC or any(corrected_SSC in fsm for fsm in FSM_SSC):
                return corrected_SSC

        return row['SSC']

    def correct_for_chr(self, df: pd.DataFrame) -> pd.DataFrame:
        if df.empty:
            return df[['Chr', 'Strand', 'SSC', 'TrStart_reads', 'TrEnd_reads', 'Group']]

        df = df.copy()
        df['SSC2'] = df['SSC'].str.split('-').map(lambda x: tuple(map(int, x)))

        correct_dict = {}
        for (chr_, strand, group), df_group in df.groupby(['Chr', 'Strand', 'Group'], observed=True):
            df_ref_ssc = df_group[df_group['category'] == 'FSM']
            df_nnc_nic = df_group[df_group['category'].isin(['NNC', 'NIC'])]
            if df_ref_ssc.empty or df_nnc_nic.empty:
                continue

            ref_exons = {(exons[i], exons[i+1]) for exons in df_ref_ssc['SSC2'] for i in range(1, len(exons)-1, 2)}
            little_exons = [(s, e) for s, e in ref_exons if e - s < self.little_exon_bp]

            ref_little_exon_paths = set()
            for start, end in little_exons:
                mask = df_ref_ssc['SSC2'].map(lambda x: start in x and end in x)
                for ssc2 in df_ref_ssc.loc[mask, 'SSC2']:
                    try:
                        idx1, idx2 = ssc2.index(start), ssc2.index(end)
                        ref_little_exon_paths.add(ssc2[idx1-1:idx2+2])
                    except ValueError:
                        continue

            ref_sites = np.unique(np.concatenate(df_ref_ssc['SSC2'].values))

            for row in df_nnc_nic.itertuples():
                ssc2 = np.array(row.SSC2)
                mapped = np.array([ref_sites[np.abs(ref_sites - site).argmin()] for site in ssc2])
                diff_sites = np.setdiff1d(ssc2, mapped)
                diff_sites_ref = np.array([mapped[np.where(ssc2 == diff)[0][0]] for diff in diff_sites])

                error_key = f"{chr_}{strand}_{row.SSC}"
                group_freq = df_group.loc[df_group['source'] == 'data', 'frequency'].sum()

                if row.category == 'NIC' and ref_little_exon_paths:
                    for path in ref_little_exon_paths:
                        path_set = set(path[1:-1])
                        if (len(path_set - set(ssc2)) == 2 and set([path[0], path[-1]]).issubset(ssc2)):
                            correct_path = row.SSC.replace(
                                f"{path[0]}-{path[-1]}",
                                '-'.join(map(str, path))
                            )
                            correct_dict[error_key] = correct_path
                            break

                elif row.category == 'NNC':
                    # 1
                    if ref_little_exon_paths:
                        for path in ref_little_exon_paths:
                            if not self.is_within_littleExon_range(diff_sites, path[0], path[-1]):
                                continue
                            if len(diff_sites) == 1:
                                if diff_sites_ref[0] == path[0]:
                                    correct_path = row.SSC.replace(
                                        f"{diff_sites[0]}-{path[-1]}",
                                        '-'.join(map(str, path))
                                    )
                                    correct_dict[error_key] = correct_path
                                    break
                                elif diff_sites_ref[0] == path[-1]:
                                    correct_path = row.SSC.replace(
                                        f"{path[0]}-{diff_sites[0]}",
                                        '-'.join(map(str, path))
                                    )
                                    correct_dict[error_key] = correct_path
                                    break
                            elif len(diff_sites) == 2 and diff_sites_ref[0] == path[0] and diff_sites_ref[1] == path[-1]:
                                correct_path = row.SSC.replace(
                                    '-'.join(map(str, diff_sites)),
                                    '-'.join(map(str, path))
                                )
                                correct_dict[error_key] = correct_path
                                break

                    # 2
                    elif len(diff_sites) == 1 and abs(diff_sites[0] - diff_sites_ref[0]) < self.mismatch_error_sites_bp:
                        if row.frequency / group_freq < self.mismatch_error_sites_groupfreq_multiple:
                            correct_dict[error_key] = '-'.join(map(str, mapped))

                    # 3
                    elif len(diff_sites) == 2 and f"{diff_sites[0]}-{diff_sites[1]}" in row.SSC:
                        excursion = diff_sites - diff_sites_ref
                        if ((excursion[0] * excursion[1] > 0) and 
                            abs(excursion[0] - excursion[1]) < self.exon_excursion_diff_bp):
                            correct_dict[error_key] = row.SSC.replace(
                                f"{diff_sites[0]}-{diff_sites[1]}",
                                f"{diff_sites_ref[0]}-{diff_sites_ref[1]}"
                            )

                    # 4
                    elif len(diff_sites) >= 2:
                        counter = Counter(mapped)
                        duplicates = [k for k, v in counter.items() if v >= 2]
                        if len(duplicates) == 1:
                            index = np.where(mapped == duplicates[0])[0][0]
                            if strand == '+':
                                tail_sites = mapped[index:]
                                if (len(np.unique(tail_sites)) == 1 and 
                                    row.frequency / group_freq < self.fake_exon_group_freq_multiple and
                                    abs(ssc2[-1] - int(np.mean(row.TrEnd_reads))) < self.fake_exon_bp):
                                    correct_dict[error_key] = '-'.join(map(str, mapped[:index+1]))
                            else:
                                index = len(mapped) - 1 - np.where(mapped[::-1] == duplicates[0])[0][0]
                                tail_sites = mapped[:index]
                                if (len(np.unique(tail_sites)) == 1 and 
                                    row.frequency / group_freq < self.fake_exon_group_freq_multiple and
                                    abs(ssc2[0] - int(np.mean(row.TrStart_reads))) < self.fake_exon_bp):
                                    correct_dict[error_key] = '-'.join(map(str, mapped[index:]))

        df['key'] = df['Chr'].astype(str) + df['Strand'].astype(str) + '_' + df['SSC']
        df['SSC'] = df.apply(lambda row: self.correct_ssc(row, correct_dict, df), axis=1)

        df = df[df['source'] == 'data'][['Chr', 'Strand', 'SSC', 'TrStart_reads', 'TrEnd_reads', 'Group']]
        df[['Chr', 'Strand']] = df[['Chr', 'Strand']].astype(str)
        df = df.groupby(['Chr', 'Strand', 'SSC', 'Group'], as_index=False).agg({
            'TrStart_reads': lambda x: np.concatenate(x.values),
            'TrEnd_reads': lambda x: np.concatenate(x.values)
            })
        df['Chr'] = pd.Categorical(df['Chr'])
        df['Strand'] = pd.Categorical(df['Strand'])

        return df
                    
  
    def correct(self, df, ref_anno):
        df = df[['Chr', 'Strand', 'SSC', 'category', 'TrStart_reads', 'TrEnd_reads']].copy()
        df['key'] = df['Chr'].astype(str) + df['Strand'].astype(str) + '_' + df['SSC'].astype(str)

        ref_anno = ref_anno[['Chr', 'Strand', 'SSC']].copy()
        ref_anno['key'] = ref_anno['Chr'].astype(str) + ref_anno['Strand'].astype(str) + '_' + ref_anno['SSC'].astype(str)
        ref_anno['category'] = 'FSM'
        ref_anno['TrStart_reads'] = np.empty((len(ref_anno), 0)).tolist()
        ref_anno['TrEnd_reads'] = np.empty((len(ref_anno), 0)).tolist()

        ref_anno = ref_anno[~ref_anno['key'].isin(df['key'])]

        df['source'] = 'data'
        ref_anno['source'] = 'ref'
        df['frequency'] = df['TrStart_reads'].apply(len)
        ref_anno['frequency'] = 1

        df = pd.concat([df, ref_anno], ignore_index=True)
        gene_clustering = GeneClustering(num_processes=self.num_processes)
        df = gene_clustering.cluster(df)
        valid_groups = (
            df.groupby(['Chr', 'Strand', 'Group'], observed=True)['source']
            .agg(lambda x: not x.eq('ref').all())
            .loc[lambda x: x].index
            )
        df = df[df.set_index(['Chr', 'Strand', 'Group']).index.isin(valid_groups)]

        if df.empty:
            return df[['Chr', 'Strand', 'SSC', 'TrStart_reads', 'TrEnd_reads', 'Group']]

        chr_grouped = [chr_group for _, chr_group in df.groupby(['Chr', 'Strand'], observed=True)]
        with mp.Pool(self.num_processes) as pool:
            results = pool.map(self.correct_for_chr, chr_grouped)

        return pd.concat(results, ignore_index=True)
        
    
    def anno_process(self,df_raw,df,ref_anno):
        if 'category' in df.columns:
            df = df.drop(columns=['category'])

        isoformclassifier = IsoformClassifier(num_processes = self.num_processes)
        df = isoformclassifier.add_category(df,ref_anno)

        # 1. rescue
        df = self.rescue_fsm(df_raw,df,ref_anno)

        # 2. filter group
        df = self.filter_groups(df,ref_anno)

        # 3. filter ISM
        df['Group_freq'] = df.groupby(['Chr', 'Strand', 'Group'],observed=True)['frequency'].transform('sum')
        df['Group_freq_ratio'] = df.apply(lambda row: row['frequency'] / row['Group_freq'] if row['Group_freq'] != 0 else 0, axis=1)
        df = self.ism_filter(df)
        df = df.drop(columns = ['Group_freq','Group_freq_ratio'])

        # 4. NNC/NIC correct
        df = self.correct(df,ref_anno)
        df['frequency'] = df['TrStart_reads'].apply(len)

        # 5. TS prediction
        terminalsitesprocessor = TerminalSitesProcessor(cluster_group_size = self.cluster_group_size,
                                                        eps = self.eps,
                                                        min_samples = self.min_samples,
                                                        num_processes = self.num_processes)
        df = terminalsitesprocessor.get_terminal_sites(df)
        
        return df
