import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN
import multiprocessing as mp
from functools import partial


class TerminalSitesProcessor:
    def __init__(self, cluster_group_size=1500, eps=10, min_samples=20, num_processes = 10):
        """
        初始化参数
        :param cluster_group_size: 采样时的最大组大小
        :param eps: DBSCAN 的半径参数
        :param min_samples: DBSCAN 的最小样本数参数
        """
        self.cluster_group_size = cluster_group_size
        self.eps = eps
        self.min_samples = min_samples
        self.num_processes = num_processes

    def get_sites_DBSCAN(self, row, site_col=None):
        """
        使用 DBSCAN 获取单个目标 start/end 值
        :param row: 单行数据
        :param site_col: 列名 ('TrStart' 或 'TrEnd')
        :return: 聚类后的目标值
        """
        sites = pd.Series(row[site_col])
        sites_len = len(sites)

        # 如果 group 的数量大于阈值，进行采样
        if sites_len >= self.cluster_group_size:
            sites = sites.sample(n=self.cluster_group_size, replace=False, random_state=42)

        X = np.array(sites).reshape(-1, 1)

        # 动态调整 min_samples
        if len(X) >= self.min_samples:
            min_samples = int(len(X) * 0.1)
        else:
            min_samples = self.min_samples

        # 创建 DBSCAN 对象
        dbscan = DBSCAN(eps=self.eps, min_samples=min_samples)
        group = pd.DataFrame(X, columns=[site_col])
        group['cluster'] = dbscan.fit_predict(X)

        # 获取聚类簇
        clusters = [clu for clu in group['cluster'].unique() if clu != -1]

        # 聚类结果处理
        if len(clusters) > 0:
            cluster_avg = group[group['cluster'] != -1].groupby('cluster')[site_col].mean().reset_index()
            cluster_avg.columns = ['cluster', 'Average']

            if site_col == 'TrStart':
                target_cluster = cluster_avg.loc[cluster_avg['Average'].idxmin(), 'cluster']
                site_value = min(group[group['cluster'] == target_cluster][site_col].mode())
            elif site_col == 'TrEnd':
                target_cluster = cluster_avg.loc[cluster_avg['Average'].idxmax(), 'cluster']
                site_value = max(group[group['cluster'] == target_cluster][site_col].mode())
        else:
            # 无法聚类时取最值
            site_value = min(group[site_col]) if site_col == 'TrStart' else max(group[site_col])

        return site_value

    def get_terminal_sites_forChr(self, df):
        """
        处理单个 Chr 的数据，获取起始和结束位点
        :param df: 单个 Chr 的 DataFrame
        :return: 处理后的 DataFrame
        """
        df_result = df.copy()
        df_result['pred_Start'] = df_result.apply(
            lambda row: self.get_sites_DBSCAN(row, site_col='TrStart'), axis=1
        )
        df_result['pred_End'] = df_result.apply(
            lambda row: self.get_sites_DBSCAN(row, site_col='TrEnd'), axis=1
        )
        df_result['pred_terminal_sites'] = list(zip(df_result['pred_Start'], df_result['pred_End']))
        return df_result

    def get_terminal_sites(self, df):
        """
        按 Chr 分组并并行处理所有数据
        :param df: 输入 DataFrame
        :return: 处理后的 DataFrame
        """
        dfChr_list = [dfChr for _, dfChr in df.groupby('Chr')]
        partial_func = partial(self.get_terminal_sites_forChr)

        with mp.Pool(self.num_processes) as pool:
            results = pool.map(partial_func, dfChr_list)

        # 合并结果
        df_result = pd.concat(results, ignore_index=True).reset_index(drop=True)
        return df_result
