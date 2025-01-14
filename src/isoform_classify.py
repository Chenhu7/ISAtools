import pandas as pd
import multiprocessing as mp


class IsoformClassifier:
    def __init__(self, num_processes=10):
        self.num_processes = num_processes

    @staticmethod
    def _get_isoform_category_for_row(row, df_ref):
        """
        获取单行 isoform 的分类
        """
        # 提取参考数据中的所有位点
        ref_sites = set(df_ref['exonChain'].str.split('-').explode().astype(int))
        row_sites = list(map(int, row['exonChain'].split('-')))

        # 直接检查 FSM
        if row['exonChain'] in df_ref['exonChain'].values:
            return 'FSM'

        # 检查 ISM：行的 exonChain 是否是某个参考链的前缀或后缀
        for ref_chain in df_ref['exonChain']:
            if ref_chain.startswith(row['exonChain']) or ref_chain.endswith(row['exonChain']):
                return 'ISM'

        # 检查 NIC：所有位点均在参考链中
        if set(row_sites).issubset(ref_sites):
            return 'NIC'

        # 如果有部分但不全位点匹配参考链，返回 NNC
        return 'NNC'

    @staticmethod
    def _get_isoform_category_for_chr(df_group):
        """
        对单个染色体组进行 isoform 分类
        """
        df_classification = []
        for _, group in df_group.groupby(['Chr', 'Strand', 'Group']):
            df_ref = group[group['source'] == 'ref']
            df_data = group[group['source'] == 'data']

            # 无参考链的情况，所有数据直接标记为 NNC
            if df_ref.empty:
                df_data['category'] = 'NNC'
                df_classification.append(df_data)
                continue

            # 跳过没有数据的情况
            if df_data.empty:
                continue

            # 分类
            df_data['category'] = df_data.apply(
                lambda row: IsoformClassifier._get_isoform_category_for_row(row, df_ref),
                axis=1
            )
            df_classification.append(df_data)

        # 合并所有分类结果
        return pd.concat(df_classification, ignore_index=True) if df_classification else pd.DataFrame()

    def get_isoform_category(self, df_combined):
        """
        多进程获取 isoform 分类
        """
        # 按染色体分组
        df_grouped = [group for _, group in df_combined.groupby('Chr')]

        # 使用多进程处理
        with mp.Pool(self.num_processes) as pool:
            results = pool.map(self._get_isoform_category_for_chr, df_grouped)

        # 合并所有结果
        return pd.concat(results, ignore_index=True).reset_index(drop=True)

    def add_category(self, df1, df_ref):
        """
        添加分类结果到原数据中
        """
        # 准备数据
        df1 = df1[['Chr', 'Strand', 'exonChain']].copy()
        df_ref = df_ref[['Chr', 'Strand', 'exonChain']].copy()
        df1['source'] = 'data'
        df_ref['source'] = 'ref'

        # 合并数据并进行聚类
        df_combined = pd.concat([df1, df_ref], axis=0, ignore_index=True)
        df_combined = self.cluster(df_combined)

        # 获取分类结果
        df_category = self.get_isoform_category(df_combined)

        # 合并分类结果到原数据
        return df1.merge(
            df_category[['Chr', 'Strand', 'exonChain', 'category']].drop_duplicates(),
            on=['Chr', 'Strand', 'exonChain'],
            how='left'
        )
