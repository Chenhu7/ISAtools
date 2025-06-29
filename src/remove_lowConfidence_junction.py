import numpy as np
import pandas as pd
import networkx as nx
from multiprocessing import Pool

class SpliceConsensusFilter:
    def __init__(self, threshold_lowWeight_edges, num_processes = 10):
        self.threshold_lowWeight_edges = threshold_lowWeight_edges
        self.num_processes = num_processes

    def get_error_paths(self, row, error_paths, error_sites_paths_freq):
        for error_path in error_paths:
            if error_path in row['SSC'] and row['frequency'] <= error_sites_paths_freq:
                return error_path
        return np.nan

    def get_edges_weight(self, df_group):
        dict_net = {}
        group_freq = df_group['frequency'].sum()

        for SSC, freq in zip(df_group['SSC2'], df_group['frequency']):
            freq_weight = freq / group_freq
            for i in range(len(SSC) - 1):
                edge = (SSC[i], SSC[i + 1])
                dict_net[edge] = dict_net.get(edge, [SSC[i], SSC[i + 1], 0])
                dict_net[edge][2] += freq_weight

        return dict_net

    def remove_lowWeight_edges_forgraph(self, subgraph, df_subgraph):
        edges_to_remove = [
            (u, v) for u, v, weight in subgraph.edges(data='weight') 
            if weight < self.threshold_lowWeight_edges and 
               (subgraph.out_degree(u) == 1 or subgraph.in_degree(u) == 1 or
                subgraph.out_degree(v) == 1 or subgraph.in_degree(v) == 1)
        ]
        remove_paths = [f"{u}-{v}" for u, v in edges_to_remove]
        df_subgraph['remove_edges'] = df_subgraph['SSC'].apply(
            lambda SSC: self.get_error_paths({'SSC': SSC, 'frequency': np.inf}, remove_paths, np.inf)
        )
        return df_subgraph

    def filter_edges_forchr(self, chr_group):
        _, df_chr = chr_group
        df_chr['SSC2'] = df_chr['SSC'].apply(lambda x: list(map(int, x.split('-'))))
        df_list = []

        for _, df_group in df_chr.groupby(['Strand', 'Group'], observed=True):
            dict_net = self.get_edges_weight(df_group)
            edges_list = [tuple(edge) for edge in dict_net.values()]
            G = nx.DiGraph()
            G.add_weighted_edges_from(edges_list)

            for component in nx.weakly_connected_components(G):
                subgraph = G.subgraph(component)
                subgraph_nodes = set(subgraph.nodes)
                subgraph_indices = df_group[df_group['SSC2'].apply(
                    lambda SSC: bool(set(SSC) & subgraph_nodes))].index
                df_subgraph = self.remove_lowWeight_edges_forgraph(subgraph, df_group.loc[subgraph_indices])
                df_list.append(df_subgraph)

        if not df_list:
            return pd.DataFrame(columns=df_chr.columns.tolist() + ['remove_edges'])

        return pd.concat(df_list, ignore_index=True)


    def remove_SJ(self, df):
        prev_length = -1
        while len(df) != prev_length:
            prev_length = len(df)
            with Pool(self.num_processes) as pool:
                chr_groups = list(df.groupby(['Chr','Strand'], observed=True))
                results = pool.map(self.filter_edges_forchr, chr_groups)

            results = [r for r in results if not r.empty]
            if not results:
                break

            df = pd.concat(results, ignore_index=True)
            df = df[(df['remove_edges'].isna()) | (df['frequency'] >= 5)]

        return df.drop(columns='remove_edges', errors='ignore')
