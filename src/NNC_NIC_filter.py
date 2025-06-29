import pandas as pd
import networkx as nx
from functools import partial
from multiprocessing import Pool


class NncNicGraphProcessor:
    def __init__(self, exon_excursion_diff_bp=10, error_sites_diff_bp=10, error_sites_multiple=0.01,little_exon_bp=30, 
                 little_exon_mismatch_diff_bp=10,Nolittle_exon_mismatch_diff_bp=20,little_exon_jump_multiple=0.1, num_processes=10):
        self.exon_excursion_diff_bp = exon_excursion_diff_bp
        self.error_sites_diff_bp = error_sites_diff_bp
        self.error_sites_multiple = error_sites_multiple
        self.little_exon_bp = little_exon_bp
        self.little_exon_mismatch_diff_bp = little_exon_mismatch_diff_bp
        self.Nolittle_exon_mismatch_diff_bp = Nolittle_exon_mismatch_diff_bp
        self.num_processes = num_processes
        
    def get_edges_weight_SSCfull(self, df_group):
        """
        获取df中所有边的权重  (frequency)/(group freq)   [node1,node2,weight]
        """
        dict_net = {}

        group_freq = df_group['frequency'].sum()
        for index, row in df_group.iterrows():
            SSC2 = row['SSC_full']
            freq_weight = row['frequency'] / group_freq

            for i in range(0, len(SSC2) - 1):
                node = (SSC2[i], SSC2[i + 1])
                if node not in dict_net:
                    dict_net[node] = [SSC2[i], SSC2[i + 1], freq_weight]
                else:
                    dict_net[node][2] += freq_weight
        return dict_net
    

    def get_edges_weight(self, df_group):
        """
        获取df中所有边的权重  (frequency)/(group freq)   [node1,node2,weight]
        """
        dict_net = {}

        group_freq = df_group['frequency'].sum()
        for index, row in df_group.iterrows():
            SSC2 = row['SSC2']
            freq_weight = row['frequency'] / group_freq

            for i in range(0, len(SSC2) - 1):
                node = (SSC2[i], SSC2[i + 1])
                if node not in dict_net:
                    dict_net[node] = [SSC2[i], SSC2[i + 1], freq_weight]
                else:
                    dict_net[node][2] += freq_weight
        return dict_net
    
    def correct_SSC(self, row, dict_correct, exon_chain_list):
        exon_chain = row['SSC']
        exon_chain_list = set(exon_chain_list)
        
        for error_path, correct_path in dict_correct.items():
            if error_path in exon_chain:
                corrected = exon_chain.replace(error_path, correct_path)
                if corrected in exon_chain_list:
                    exon_chain = corrected

        return exon_chain

    
    def corrected_df(self, df1, dict_correct):
        df = df1.copy()
        exon_chain_list = df['SSC'].tolist()
        df['SSC_corrected'] = df.apply(lambda row: self.correct_SSC(row, dict_correct, exon_chain_list), axis=1)

        df = df.explode(['TrStart_reads', 'TrEnd_reads'], ignore_index=True)

        df = df.groupby(['Chr', 'Strand', 'Group', 'SSC_corrected'],
                        as_index=False,observed=True).agg({
                            'TrStart_reads': list,
                            'TrEnd_reads': list
                        }).reset_index(drop=True)
        
        df.rename(columns={'SSC_corrected': 'SSC'}, inplace=True)
        df['SSC2'] = df['SSC'].apply(lambda x: list(map(int, x.split('-'))))
        df['frequency'] = df['TrStart_reads'].apply(len)
        
        df = df.merge(df1[['Chr', 'Strand', 'SSC', 'junction']],on = ['Chr', 'Strand', 'SSC'],how = 'left')
        
        return df
    
    def get_df_SSC_StartEnd(self,df1):
        df = df1.copy()
        df['SSC2'] = df['SSC'].apply(lambda x: list(map(int, x.split('-'))))

        extreme_values = (
            df.groupby(['Chr', 'Strand', 'Group'],observed=True)['SSC2']
            .agg(lambda x: (min([i for sublist in x for i in sublist]), 
                            max([i for sublist in x for i in sublist])))
            .apply(pd.Series)
            .rename(columns={0: 'min_value', 1: 'max_value'})
            .reset_index())

        df = df.merge(extreme_values, on=['Chr', 'Strand', 'Group'])

        df['SSC_full'] = df.apply(
            lambda row: [row['min_value'] - 100] + row['SSC2'] + [row['max_value'] + 100], axis=1
        )
        df = df.drop(columns = ['min_value','max_value'])

        df['SSC_full_str'] = df['SSC_full'].apply(lambda x: '-'.join(map(str, x)))

        return df
    

    def remove_misalignment(self, df_subgraph,subgraph):
        out_nodes = [node for node, out_degree in subgraph.out_degree() if out_degree > 1]
        in_nodes = [node2 for node2, in_degree in subgraph.in_degree() if in_degree > 1]

        forked_paths = []
        for out_node in out_nodes:
            for in_node in in_nodes:
                if nx.has_path(subgraph, out_node, in_node):
                    paths = nx.all_simple_paths(subgraph, source=out_node, target=in_node, cutoff=6)
                    forked_paths.extend(paths)

        SSCs = df_subgraph['SSC_full_str'].tolist()

        forked_paths_filtered = []
        for path in forked_paths:
            path_str = '-'.join(map(str, path))
            if any(path_str in chain for chain in SSCs):
                forked_paths_filtered.append(path)
                    
        dict_path = {}
        for path in forked_paths_filtered:
            start_end = (path[0],path[-1])
            if start_end not in dict_path:
                dict_path[start_end] = [path]
            else:
                dict_path[start_end] += [path]
        dict_path = {k: v for k, v in dict_path.items() if len(v) > 1}
        
        subgraph_total_freq = df_subgraph['frequency'].sum()
        dict_subgraph_SSC_freq = dict(zip(df_subgraph["SSC_full_str"], df_subgraph["frequency"]))

        for combine in dict_path:
            paths = dict_path[combine]
            for i in range(len(paths)):
                path = paths[i]
                path_weight = 0
                for SSC in dict_subgraph_SSC_freq:
                    if '-'.join(map(str,path)) in SSC:
                        path_weight += dict_subgraph_SSC_freq[SSC] / subgraph_total_freq

                dict_path[combine][i].append(path_weight)

        correct_paths = {}

        # 1
        dict_path_filtered = {key: [v for v in value if len(v) == 5] for key, value in dict_path.items()}
        dict_path_filtered = {k: v for k, v in dict_path_filtered.items() if len(v) > 1}
        
        if len(dict_path_filtered) > 0:
            for combine in dict_path_filtered:
                paths = dict_path_filtered[combine]
                paths = sorted(paths, key=lambda x: x[-1],reverse=True)

                ref_path = paths[0]
                ref_path_str = '-'.join(map(str,ref_path[:-1]))

                for path in paths[1:]:
                    path_str = '-'.join(map(str,path[:-1]))
              
                    diff1 = path[1] - ref_path[1]
                    diff2 = path[2] - ref_path[2]

                    if (diff1 > 0 and diff2 > 0) or (diff1 < 0 and diff2 < 0):
                        if abs(diff1 - diff2) <= self.exon_excursion_diff_bp:
                            error_exon = '-'.join(path_str.split('-')[1:-1])
                            ref_exon = '-'.join(ref_path_str.split('-')[1:-1])
                            correct_paths[error_exon] = ref_exon
        

        # 2
        sites_dict = {}
        for _,row in df_subgraph.iterrows():
            SSC2 = row['SSC2']
            freq = row['frequency']
            
            for site in SSC2:
                if site not in sites_dict:
                    sites_dict[site] = freq
                else:
                    sites_dict[site] += freq
        
        sites_dict_keys = sorted(list(sites_dict.keys()))
        sites_dict_keys = list(zip(sites_dict_keys, sites_dict_keys[1:]))
        
        potential_error_sites = [sites for sites in sites_dict_keys if abs(sites[1] - sites[0]) <= self.error_sites_diff_bp]
        
        error_sites = []
        for p_sites in potential_error_sites:
            site1_freq = sites_dict[p_sites[0]]
            site2_freq = sites_dict[p_sites[1]]
            
            if site1_freq < site2_freq and site1_freq/site2_freq <= self.error_sites_multiple:
                error_sites.append(p_sites[0])
            elif site1_freq > site2_freq and site2_freq/site1_freq <= self.error_sites_multiple:
                error_sites.append(p_sites[1])
        
        remove_index = []
        for index,row in df_subgraph.iterrows():
            if len(set(error_sites) & set(row['SSC2'])) > 0:
                remove_index.append(index)
        if len(remove_index) > 0:
            df_subgraph.drop(remove_index, inplace=True)
        
        # 3
        try:
            df_subgraph = self.get_df_SSC_StartEnd(df_subgraph)
        except:
            print(df_subgraph)
        
        node6_path = []
        node4_path = []
        
        for SSC in df_subgraph['SSC_full'].tolist():
            if len(SSC) >= 4:
                for i in range(0,len(SSC),2):
                    exon_bp = SSC[i+1] - SSC[i]

                    if exon_bp <= self.little_exon_bp:
                        node4_index_S = i - 1
                        node4_index_E = i + 3
                        node6_index_S = i - 2
                        node6_index_E = i + 4

                        try:
                            n4_path = SSC[node4_index_S:node4_index_E]
                            node4_path.append(tuple(n4_path))
                        except:
                            pass

                        try:
                            n6_path = SSC[node6_index_S:node6_index_E]
                            node6_path.append(tuple(n6_path))
                        except:
                            pass
        
        if len(node4_path) > 0:
            for Lexon in node4_path:
                error_path = '-'.join([str(s) for s in [Lexon[0],Lexon[-1]]])
                correct_path = '-'.join([str(s) for s in Lexon])
                correct_paths[error_path] = correct_path
            
        if len(node6_path) > 0:
            for N6_path in node6_path:
                for _,row in df_subgraph.iterrows():
                    row_SSCFull = row['SSC_full']
                    if N6_path[0] in row_SSCFull and N6_path[-1] in row_SSCFull:
                        row_path = row_SSCFull[row_SSCFull.index(N6_path[0]):row_SSCFull.index(N6_path[-1]) + 1]
                        if len(set(row_path) - set(N6_path)) > 0 and len(row_path) == 4:
                            N6_path_length = sum([N6_path[i+1] - N6_path[i] for i in range(0,len(N6_path),2)])
                            row_path_length = sum([row_path[i+1] - row_path[i] for i in range(0,len(row_path),2)])

                            if abs(N6_path_length - row_path_length) <= self.little_exon_mismatch_diff_bp:
                                error_path = '-'.join([str(row_path[1]),str(row_path[2])])
                                correct_path = '-'.join([str(s) for s in N6_path[1:-1]])
                                correct_paths[error_path] = correct_path
                            
        
        # 4
        no_little_exon_correct = {}
        forked_paths_4_6 = [path[:-1] for path in forked_paths_filtered if len(path) in [5,7]]

        dict_path46 = {}
        for path in forked_paths_4_6:
            start_end = (path[0],path[-1])
            if start_end not in dict_path46:
                dict_path46[start_end] = [path]
            else:
                dict_path46[start_end] += [path]
        dict_path46 = {k: v for k, v in dict_path46.items() if len(v) > 1}
        
        for combine in dict_path46:
            node6_paths = [path for path in dict_path46[combine] if len(path) == 6]
            node4_paths = [path for path in dict_path46[combine] if len(path) == 4]
            
            for node6 in node6_paths:
                for node4 in node4_paths:
                    node6_length = sum([node6[i+1] - node6[i] for i in range(0,len(node6),2)])
                    node4_length = sum([node4[i+1] - node4[i] for i in range(0,len(node4),2)])
                    
                    if abs(node4_length - node6_length) <= self.Nolittle_exon_mismatch_diff_bp:
                        error_path = '-'.join([str(node4[1]),str(node4[2])])
                        correct_path = '-'.join([str(s) for s in node6[1:-1]])
                        no_little_exon_correct[error_path] = correct_path
        
        for No_Lexon_error in no_little_exon_correct:
            for _,row in df_subgraph.iterrows():
                if No_Lexon_error in row['SSC']:
                    junctions = row['junction'].split(',')
                    if not all(x.upper() == 'GT-AG' for x in junctions):
                        correct_paths[No_Lexon_error] = no_little_exon_correct[No_Lexon_error]

        if len(correct_paths) > 0:
            remove_SSC = set(correct_paths.keys())
            df_subgraph = df_subgraph[~df_subgraph['SSC'].isin(remove_SSC)][['Chr','Strand','Group','SSC','TrStart_reads','TrEnd_reads','SSC2','frequency', 'junction']]
        else:
            df_subgraph = df_subgraph[['Chr','Strand','Group','SSC','TrStart_reads','TrEnd_reads','SSC2','frequency', 'junction']]

        return df_subgraph
    
    def nnc_nic_graph_forChr(self, df):
        df = df.copy()
        
        df = self.get_df_SSC_StartEnd(df)

        df_list = []
        i = 0
        for _, df_group in df.groupby(['Chr', 'Strand', 'Group'],observed=True):
            dict_net = self.get_edges_weight_SSCfull(df_group)

            edges_list = list(dict_net.values())
            edges_list = [tuple(l) for l in edges_list]
            G = nx.DiGraph()
            G.add_weighted_edges_from(edges_list)

            for component in nx.weakly_connected_components(G):
                subgraph = G.subgraph(component)

                subgraph_nodes = set(list(subgraph.nodes()))
                subgraph_index = []
                for index,row in df_group.iterrows():
                    SSC2 = set(row['SSC2'])
                    if len(subgraph_nodes & SSC2) >0:
                        subgraph_index.append(index)

                df_subgraph = df_group.loc[subgraph_index]

                df_subgraph['Group'] = i
                i += 1
                df_subgraph = self.remove_misalignment(df_subgraph,subgraph)

                df_list.append(df_subgraph)

        df = pd.concat(df_list,ignore_index=True)
        return df

    def nnc_nic_graph(self, df):
        dfChr_list = [dfChr for _,dfChr in df.groupby('Chr',observed=True)]

        partial_func = partial(self.nnc_nic_graph_forChr)

        with Pool(self.num_processes) as pool:
            results = pool.map(partial_func, dfChr_list)

        df_result = pd.concat(results,ignore_index= True).reset_index(drop=True)

        df_result = df_result.drop(columns = ['SSC2','junction'])
        return df_result
