import os
import numpy as np
import pandas as pd
from multiprocessing import Pool
from typing import Optional, Dict, List
from src.gene_grouping import GeneClustering
import time

class IsoformAnnotator:
    def __init__(self, num_processes: int = 20):
        self.num_processes = num_processes

    def _annotate_novel_only(self, df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()
        df['Group'] = df['Chr'].astype(str) + '_' + df['Strand'].astype(str) + '_' + df['Group'].astype(str)

        Gi = 1  # Gene计数器
        gene_id_map = {}
        tr_id_map = {}
        gene_names_map = {}

        for _, group_df in df.groupby('Group',observed=True):
            Ti = 1  # Transcript计数器
            first_row = group_df.iloc[0]
            gene_id = f"{first_row['Chr']}_{first_row['Strand']}_Gene{Gi}"
            gene_name = gene_id
            Gi += 1

            for idx in group_df.index:
                gene_id_map[idx] = gene_id
                tr_id_map[idx] = f"{gene_id}_Tr{Ti}"
                gene_names_map[idx] = gene_name
                Ti += 1

        df['GeneID'] = pd.Series(gene_id_map)
        df['TrID'] = pd.Series(tr_id_map)
        df['GeneName'] = pd.Series(gene_names_map)

        return df.drop(columns='Group')

    def _vectorized_annotation_match(self, merged: pd.DataFrame) -> pd.DataFrame:
        ref_counts = merged['TrStart_ref'].str.len()
        max_ref = ref_counts.max()
        
        start_arr = np.stack(
            merged['TrStart_ref'].apply(
                lambda x: np.pad(x, (0, max_ref - len(x)), constant_values=np.inf)
            )
        )
        end_arr = np.stack(
            merged['TrEnd_ref'].apply(
                lambda x: np.pad(x, (0, max_ref - len(x)), constant_values=np.inf)
            )
        )
        diffs = (
            np.abs(merged['TrStart'].values[:, None] - start_arr) +
            np.abs(merged['TrEnd'].values[:, None] - end_arr)
        )
        
        with np.errstate(invalid='ignore'):
            min_indices = np.nanargmin(np.where(np.isfinite(diffs), diffs, np.nan), axis=1)
            valid_mask = ~np.all(np.isinf(diffs), axis=1)
        
        tr_ids = np.empty(len(merged), dtype=object)
        gene_ids = np.empty(len(merged), dtype=object)
        gene_names = np.empty(len(merged), dtype=object)
        
        for i, (valid, idx) in enumerate(zip(valid_mask, min_indices)):
            if valid and idx < len(merged['TrID'].iloc[i]):
                tr_ids[i] = merged['TrID'].iloc[i][idx]
                gene_ids[i] = merged['GeneID'].iloc[i][idx]
                gene_names[i] = merged['GeneName'].iloc[i][idx]
            else:
                tr_ids[i] = None
                gene_ids[i] = None
                gene_names[i] = None
        
        return pd.DataFrame({
            'TrID': tr_ids,
            'GeneID': gene_ids,
            'GeneName': gene_names
        })

    def _process_unannotated(self, merged: pd.DataFrame) -> pd.DataFrame:
        mask_none = merged['GeneID'].isna()
        if not mask_none.any():
            return merged

        merged['key'] = (
            merged['Chr'].astype(str) + 
            merged['Strand'].astype(str) + 
            merged['Group'].astype(str)
        )
        
        group_info = merged[~mask_none].groupby('key',observed=True).agg({
            'GeneID': 'first',
            'GeneName': 'first'
        })
        merged.loc[mask_none, 'GeneID'] = merged.loc[mask_none, 'key'].map(group_info['GeneID'])
        merged.loc[mask_none, 'GeneName'] = merged.loc[mask_none, 'key'].map(group_info['GeneName'])

        new_genes_mask = merged['GeneID'].isna()
        merged.loc[new_genes_mask, 'Gi'] = (
            merged[new_genes_mask]
            .groupby('key', observed=True)
            .ngroup() + 1
        )
        merged.loc[new_genes_mask, 'GeneID'] = (
            merged.loc[new_genes_mask, 'Chr'].astype(str) + "_" +
            merged.loc[new_genes_mask, 'Strand'].astype(str) + "_Gene" +
            merged.loc[new_genes_mask, 'Gi'].astype('int').astype(str)
        )
        merged.loc[new_genes_mask, 'GeneName'] = merged.loc[new_genes_mask, 'GeneID']
        merged.drop(columns='Gi', inplace=True)

        needs_tr_mask = merged['TrID'].isna()

        merged.loc[needs_tr_mask, 'Ti'] = (
            merged[needs_tr_mask]
            .groupby('GeneID', observed=True)
            .cumcount()
            .add(1)
            .astype('int')
        )

        merged.loc[needs_tr_mask, 'TrID'] = (
            merged.loc[needs_tr_mask, 'GeneID'] + 
            '_Tr' + 
            merged.loc[needs_tr_mask, 'Ti'].astype('int').astype(str)
            )
        merged.drop(columns='Ti', inplace=True)

        merged.drop(columns='key', inplace=True)
        return merged

    def annotate(self, df: pd.DataFrame, ref_anno: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        df['uniqueTr'] = 'Tr' + df.groupby(['Chr', 'Strand', 'SSC', 'TrStart', 'TrEnd'],observed=True).ngroup().astype(str)
        df_unique = df[['Chr', 'Strand', 'SSC', 'TrStart', 'TrEnd', 'uniqueTr']].drop_duplicates()
        
        gene_clustering = GeneClustering(num_processes=self.num_processes)
        df_unique = gene_clustering.cluster(df_unique)

        if ref_anno is None:
            df_unique = self._annotate_novel_only(df_unique)
        else:
            ref_anno = ref_anno[ref_anno['SSC'].isin(df_unique['SSC'].unique())].copy()
            
            merged = df_unique.merge(
                ref_anno,
                on=['Chr', 'Strand', 'SSC'],
                how='left',
                suffixes=('', '_ref')
            )
            
            dup_cols = ['Chr', 'Strand', 'SSC', 'TrStart', 'TrEnd', 'uniqueTr', 'Group']
            dup_rows = merged[merged.duplicated(subset=dup_cols, keep=False)]
            merged_nodup = merged[~merged.duplicated(subset=dup_cols, keep=False)]
            if not dup_rows.empty:
                dup_rows = dup_rows.groupby(dup_cols,observed=True).agg({
                    'TrID': list,
                    'GeneID': list,
                    'GeneName': list,
                    'TrStart_ref': list,
                    'TrEnd_ref': list
                }).reset_index()
                dup_rows[['TrID', 'GeneID', 'GeneName']] = self._vectorized_annotation_match(dup_rows)
            
            merged = pd.concat([merged_nodup, dup_rows]).pipe(self._process_unannotated)
            
            df_unique = df_unique.join(
                merged.set_index('uniqueTr')[['TrID', 'GeneID', 'GeneName']],
                on='uniqueTr'
            )
            
            group_cols = ['Chr', 'Strand', 'SSC', 'Group', 'TrID', 'GeneID', 'GeneName']
            group_sizes = df_unique.groupby(group_cols,observed=True)['TrID'].transform('size')
            group_idx = df_unique.groupby(group_cols,observed=True).cumcount() + 1
            df_unique['TrID'] = np.where(
                group_sizes > 1,
                df_unique['TrID'] + '_TssTes' + group_idx.astype(str),
                df_unique['TrID']
            )

        id_maps = df_unique.set_index('uniqueTr')[['TrID', 'GeneID', 'GeneName']].to_dict()
        for col in ['TrID', 'GeneID', 'GeneName']:
            df[col] = df['uniqueTr'].map(id_maps[col])
        
        return df

    def to_gtf(self, df: pd.DataFrame, output_path: str) -> None:
        df = df[['Chr', 'Strand', 'SSC','TrStart','TrEnd','TrID', 'GeneID', 'GeneName']].drop_duplicates()
        df['sites'] = df.apply(lambda r: sorted(
            list(map(int, r['SSC'].split('-'))) + [int(r['TrStart']), int(r['TrEnd'])]
        ), axis=1)
        
        with open(output_path, 'w') as f:
            for _, row in df.iterrows():
                f.write(
                    f"{row['Chr']}\tISAtools\ttranscript\t{row['sites'][0]}\t{row['sites'][-1]}\t.\t{row['Strand']}\t.\t"
                    f'gene_id "{row["GeneID"]}"; transcript_id "{row["TrID"]}"; gene_name "{row["GeneName"]}";\n'
                )
                for i in range(0, len(row['sites']), 2):
                    f.write(
                        f"{row['Chr']}\tISAtools\texon\t{row['sites'][i]}\t{row['sites'][i+1]}\t.\t{row['Strand']}\t.\t"
                        f'gene_id "{row["GeneID"]}"; transcript_id "{row["TrID"]}"; exon_number "{i//2 + 1}";\n'
                    )


    def quantify(self, df: pd.DataFrame) -> Dict[str, pd.DataFrame]:
        df['length'] = df['sites'].apply(
            lambda s: sum(s[i+1]-s[i] for i in range(0, len(s), 2))
        )

        df['rpk'] = df['quantification'] / (df['length'] / 1000)
        tpm_scale = 1e6 / df['rpk'].sum()
        df['TPM'] = (df['rpk'] * tpm_scale).round(2)
        
        return {
            "transcript_counts": df.pivot_table(
                index=['TrID', 'GeneID', 'GeneName'],
                columns='sample',
                values='quantification',
                fill_value=0
            ),
            "transcript_tpm": df.pivot_table(
                index=['TrID', 'GeneID', 'GeneName'],
                columns='sample',
                values='TPM',
                fill_value=0
            ),
            "gene_counts": df.groupby(['GeneID', 'GeneName', 'sample'],observed=True)['quantification'].sum().unstack(fill_value=0),
            "gene_tpm": df.groupby(['GeneID', 'GeneName', 'sample'],observed=True)['TPM'].sum().unstack(fill_value=0)
        }

    def save_results(self, df_results: pd.DataFrame, output_dir: str, ref_anno: Optional[pd.DataFrame] = None) -> None:
        if ref_anno is not None:
            ref_anno = ref_anno.fillna("-")
        total_start = time.time()

        annotated_df = self.annotate(df_results, ref_anno)
        self.to_gtf(annotated_df, os.path.join(output_dir, 'isatools.filtered.transcript.gtf'))
        df_results['sites'] = df_results.apply(lambda row: sorted(
            list(map(int, row['SSC'].split('-'))) + [int(row['TrStart']),int(row['TrEnd'])]
        ), axis=1)
        quant_tables = self.quantify(df_results)
        for name, table in quant_tables.items():
            table.to_csv(os.path.join(output_dir, f'isatools_{name}.tsv'), sep='\t')
            