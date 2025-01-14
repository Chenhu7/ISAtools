import pandas as pd
import networkx as nx
from functools import partial
from multiprocessing import Pool


class NncNicGraphProcessor:
    def __init__(self, exon_excursion_diff_bp=10, little_exon_bp=30, little_exon_jump_multiple=0.1, num_processes=10):
        self.exon_excursion_diff_bp = exon_excursion_diff_bp
        self.little_exon_bp = little_exon_bp
        self.little_exon_jump_multiple = little_exon_jump_multiple
        self.num_processes = num_processes
        

    def get_edges_weight_exonChainfull(self, df_group):
        """
        获取df中所有边的权重  (frequency)/(group freq)   [node1,node2,weight]
        """
        dict_net = {} # 收集所有边的总频率

        group_freq = df_group['frequency'].sum()
        # 为每条边分配权重
        for index, row in df_group.iterrows():
            exonChain2 = row['exonChain_full']
            freq_weight = row['frequency'] / group_freq
            # freq_weight = row['frequency']

            for i in range(0, len(exonChain2) - 1):
                node = (exonChain2[i], exonChain2[i + 1])
                if node not in dict_net:
                    dict_net[node] = [exonChain2[i], exonChain2[i + 1], freq_weight]
                else:
                    dict_net[node][2] += freq_weight
        return dict_net
    

    def get_edges_weight(self, df_group):
        """
        获取df中所有边的权重  (frequency)/(group freq)   [node1,node2,weight]
        """
        dict_net = {} # 收集所有边的总频率

        group_freq = df_group['frequency'].sum()
        # 为每条边分配权重
        for index, row in df_group.iterrows():
            exonChain2 = row['exonChain2']
            freq_weight = row['frequency'] / group_freq
            # freq_weight = row['frequency']

            for i in range(0, len(exonChain2) - 1):
                node = (exonChain2[i], exonChain2[i + 1])
                if node not in dict_net:
                    dict_net[node] = [exonChain2[i], exonChain2[i + 1], freq_weight]
                else:
                    dict_net[node][2] += freq_weight
        return dict_net
    

    def correct_exonChain(self, row, dict_correct, exon_chain_list):
        """纠正错误的外显子链"""
        exon_chain = row['exonChain']
        exon_chain_list = set(exon_chain_list)

        for error_path, correct_path in dict_correct.items():
            if error_path in exon_chain:
                corrected = exon_chain.replace(error_path, correct_path)
                if corrected in exon_chain_list:
                    exon_chain = corrected
        return exon_chain

    def corrected_df(self, df, dict_correct):
        """更新纠正后的DataFrame"""
        df = df.copy()
        exon_chain_list = df['exonChain'].tolist()
        df['exonChain_corrected'] = df.apply(lambda row: self.correct_exonChain(row, dict_correct, exon_chain_list), axis=1)

        df = df.explode(['TrStart', 'TrEnd'], ignore_index=True)
        # df = df.groupby(['Chr', 'Strand', 'Group', 'exonChain_corrected', 'junction', 'contains_non_conservative'],
        #                 as_index=False).agg({
        #                     'TrStart': list,
        #                     'TrEnd': list
        #                     # 'category': lambda x: sorted(x)[0]
        #                 }).reset_index(drop=True)
        
        df = df.groupby(['Chr', 'Strand', 'Group', 'exonChain_corrected', 'junction'],
                        as_index=False).agg({
                            'TrStart': list,
                            'TrEnd': list
                            # 'category': lambda x: sorted(x)[0]
                        }).reset_index(drop=True)
        
        df.rename(columns={'exonChain_corrected': 'exonChain'}, inplace=True)
        df['exonChain2'] = df['exonChain'].apply(lambda x: list(map(int, x.split('-'))))
        df['frequency'] = df['TrStart'].apply(len)
        return df    
    

    def remove_misalignment(self, df_subgraph,subgraph):
        ### 1. 获取中间段低频exon
        out_nodes = [node for node, out_degree in subgraph.out_degree() if out_degree > 1] # 起点 出度大于一的所有节点
        in_nodes = [node2 for node2, in_degree in subgraph.in_degree() if in_degree > 1] # 终点 入度大于一的所有节点

        # 存储所有起点-终点路径
        forked_paths = []

        # 获取所有起点到终点的所有路径
        for out_node in out_nodes:
            for in_node in in_nodes:
                all_path = list(nx.all_simple_paths(subgraph, source=out_node, target=in_node))
                forked_paths += all_path

        # 获取所有节点数为4的path
        forked_paths = [path for path in forked_paths if len(path) <= 6]


        # 路径有read支持
        exonChain_group = df_subgraph['exonChain_full_str'].tolist()
        forked_paths_filtered = []
        for path in forked_paths:
            for exonChain in exonChain_group:
                if '-'.join([str(p) for p in path]) in exonChain:
                    forked_paths_filtered.append(path)
                    break


        # 相同起止节点以及所有路径 相同起止位点至少2个路径
        dict_path = {}
        for path in forked_paths_filtered:
            start_end = (path[0],path[-1])
            if start_end not in dict_path:
                # dict_path[start_end] = ['-'.join([str(p) for p in path])]
                dict_path[start_end] = [path]
            else:
                # dict_path[start_end].append('-'.join([str(p) for p in path]))
                dict_path[start_end] += [path]
        dict_path = {k: v for k, v in dict_path.items() if len(v) > 1}

        ## 记录所有path的权重
        subgraph_total_freq = df_subgraph['frequency'].sum()
        dict_subgraph_exonChain_freq = dict(zip(df_subgraph["exonChain_full_str"], df_subgraph["frequency"]))

        for combine in dict_path:
            paths = dict_path[combine]
            # weights = []
            for i in range(len(paths)):
                path = paths[i]
                path_weight = 0 # path的总支持数
                for exonChain in dict_subgraph_exonChain_freq:
                    # if '-'.join([str(p) for p in path]) in exonChain:
                    # print('-'.join(map(str,path)))
                    if '-'.join(map(str,path)) in exonChain:
                        path_weight += dict_subgraph_exonChain_freq[exonChain] / subgraph_total_freq # path的权重为freq/group freq

                #直接将权重添加到path列表的末尾
                dict_path[combine][i].append(path_weight)


        ############
        ## exon错配
        ############
        correct_paths = {} # 纠正

        ### 1. 偏移错配: 偏移量小于10bp内 / 含有非保守剪接碱基
        dict_path_filtered = {key: [v for v in value if len(v) == 5] for key, value in dict_path.items()}
        dict_path_filtered = {k: v for k, v in dict_path_filtered.items() if len(v) > 1}

        if len(dict_path_filtered) > 0:
            for combine in dict_path_filtered:
                paths = dict_path_filtered[combine]
                paths = sorted(paths, key=lambda x: x[-1],reverse=True)
                # print(paths)

                ref_path = paths[0]
                ref_path_str = '-'.join(map(str,ref_path[:-1]))

                for path in paths[1:]:
                    path_str = '-'.join(map(str,path[:-1]))

                    # 计算偏移量                
                    diff1 = path[1] - ref_path[1] # 差异exon第一个位点差值
                    diff2 = path[2] - ref_path[2] # 差异exon第二个位点差值

                    # 有潜在偏移
                    if (diff1 > 0 and diff2 > 0) or (diff1 < 0 and diff2 < 0):
                        # 偏移量相差 bp
                        if abs(diff1 - diff2) <= self.exon_excursion_diff_bp:
                            correct_paths[path_str] = ref_path_str

        # 矫正
        if len(correct_paths) > 0:
            df_subgraph = self.corrected_df(df_subgraph,correct_paths)
        else:
            # 和已纠正的df保持一致
            # df_subgraph = df_subgraph[['Chr','Strand','Group','exonChain','TrStart','TrEnd','category','exonChain2','frequency','junction','contains_non_conservative']]
            df_subgraph = df_subgraph[['Chr','Strand','Group','exonChain','TrStart','TrEnd','exonChain2','frequency','junction']]


        ###############
        ##  小exon错配
        ##############
        correct_paths = {} # 纠正
        ## (1) 先找小exon  先筛选出含有小exon的path  4node path 和 6node path
        node6_path = [] # 含有小exon的6个位点path
        node4_path = [] # 含有小exon的4个位点path

        for exonChain in df_subgraph['exonChain2'].tolist():
            if len(exonChain) >= 4:
                for i in range(1,len(exonChain)-1,2):
                    exon_bp = exonChain[i+1] - exonChain[i]

                    if exon_bp <= self.little_exon_bp:
                        node4_index_S = i - 1
                        node4_index_E = i + 3
                        node6_index_S = i - 2
                        node6_index_E = i + 4

                        try:
                            n4_path = exonChain[node4_index_S:node4_index_E]
                            node4_path.append(tuple(n4_path))
                            # print('n4_path---',n4_path)
                        except:
                            pass

                        try:
                            n6_path = exonChain[node6_index_S:node6_index_E]
                            node6_path.append(tuple(n6_path))
                            # print('n6_path---',n6_path)
                        except:
                            pass

        ## (2) 含有小exon的path和相应的情况分别处理
        node4_path = list(set(node4_path))
        node6_path = list(set(node6_path))

        # 筛选过滤dict_path
        node4_combine = [(path[0], path[-1]) for path in node4_path if len(path) > 0]
        node6_combine = [(path2[0], path2[-1]) for path2 in node6_path if len(path2) > 0]

        combines = list(set(node4_combine + node6_combine))

        if len(combines) > 0: # 存在小exon
            # print('========',combines)

            ####### 重新构建图
            dict_subgraph = self.get_edges_weight(df_subgraph)

            edges_list = list(dict_subgraph.values())
            edges_list = [tuple(l) for l in edges_list]
            subgraph = nx.DiGraph()  # 有向图
            subgraph.add_weighted_edges_from(edges_list)

            ### 1. 获取中间段低频exon
            out_nodes = [node for node, out_degree in subgraph.out_degree() if out_degree > 1] # 起点 出度大于一的所有节点
            in_nodes = [node2 for node2, in_degree in subgraph.in_degree() if in_degree > 1] # 终点 入度大于一的所有节点

            # 存储所有起点-终点路径
            forked_paths = []

            # 获取所有起点到终点的所有路径
            for out_node in out_nodes:
                for in_node in in_nodes:
                    all_path = list(nx.all_simple_paths(subgraph, source=out_node, target=in_node))
                    forked_paths += all_path

            forked_paths = [path for path in forked_paths if len(path) <= 6]


            # 路径有read支持
            exonChain_group = df_subgraph['exonChain'].tolist()
            forked_paths_filtered = []
            for path in forked_paths:
                for exonChain in exonChain_group:
                    if '-'.join([str(p) for p in path]) in exonChain:
                        forked_paths_filtered.append(path)
                        break

            # 相同起止节点以及所有路径 相同起止位点至少2个路径
            dict_path = {}
            for path in forked_paths_filtered:
                start_end = (path[0],path[-1])
                if start_end not in dict_path:
                    dict_path[start_end] = [path]
                else:
                    dict_path[start_end] += [path]
            dict_path = {k: v for k, v in dict_path.items() if len(v) > 1}


            ## 记录所有path的权重
            subgraph_total_freq = df_subgraph['frequency'].sum()
            dict_subgraph_exonChain_freq = dict(zip(df_subgraph["exonChain"], df_subgraph["frequency"]))

            for combine in dict_path:
                paths = dict_path[combine]
                # weights = []
                for i in range(len(paths)):
                    path = paths[i]
                    path_weight = 0 # path的总支持数
                    for exonChain in dict_subgraph_exonChain_freq:
                        # if '-'.join([str(p) for p in path]) in exonChain:
                        # print('-'.join(map(str,path)))
                        if '-'.join(map(str,path)) in exonChain:
                            path_weight += dict_subgraph_exonChain_freq[exonChain] / subgraph_total_freq # path的权重为freq/group freq

                    #直接将权重添加到path列表的末尾
                    dict_path[combine][i].append(path_weight)


            dict_path_filtered = {k: v for k, v in dict_path.items() if k in combines}
            dict_path_filtered = {k: v for k, v in dict_path_filtered.items() if len(v) > 1}

            dict_node4 = dict(zip([(n[0], n[-1]) for n in node4_path if len(n) > 0], node4_path))

            if len(dict_path_filtered) > 0:
                correct_paths = {} # 纠正
                for combine in dict_path_filtered:
                    # 1. 小exon跳跃: 2node path : 4node path
                    if combine in node4_combine:
                        paths = dict_path_filtered[combine]
                        # print('ppppaths  ',paths)
                        paths2 = [p for p in paths if len(p) == 3]
                        # paths3 = [p for p in paths if len(p) == 4]
                        paths4 = [p for p in paths if len(p) == 5]

                        node4_sites = dict_node4[combine]
                        node4_sites_str = '-'.join(map(str,list(node4_sites)))
                        node4_weight = 0
                        for p in paths4:
                            if len(set(node4_sites) & set(p)) == 4:
                                node4_weight = p[-1]

                        # 小跳跃
                        for p2 in paths2:
                            if p2[-1] < self.little_exon_jump_multiple:
                                correct_paths['-'.join(map(str,list(p2[:-1])))] = node4_sites_str


                    # 3. 小exon错配 4node path : 6node path
                    if combine in node6_combine:
                        paths = dict_path_filtered[combine]
                        # print('ppppaths6666  ',paths)
                        paths4 = [p for p in paths if len(p) == 5]
                        paths6 = [p for p in paths if len(p) == 7]

                        for p4 in paths4:
                            p4_exonLength = p4[1] - p4[0] + p4[3] - p4[2]

                            for p6 in paths6:
                                p6_exonLength = p6[1] - p6[0] + p6[3] - p6[2] + p6[5] - p6[4]

                                if abs(p4_exonLength - p6_exonLength) <= 5:
                                    correct_paths['-'.join(map(str,list(p4[:-1])))] = '-'.join(map(str,list(p6[:-1])))

        # 矫正
        # print('==============2',correct_paths)
        if len(correct_paths) > 0:
            df_subgraph = self.corrected_df(df_subgraph,correct_paths)
        else:
            # 和已纠正的df保持一致
            # df_subgraph = df_subgraph[['Chr','Strand','Group','exonChain','TrStart','TrEnd','category','exonChain2','frequency','junction','contains_non_conservative']]
            df_subgraph = df_subgraph[['Chr','Strand','Group','exonChain','TrStart','TrEnd','exonChain2','frequency','junction']]

        return df_subgraph
    
    
    def nnc_nic_graph_forChr(self, df):
        """
        图论：过滤由于错配引起的NNC
        1. exon错配偏移
        2. 小exon跳跃
        3. 小exon错配：凸起  错配到两个exon

        原理：
        1. 添加初始起止位点： exonChain极值±100
        2. 给每个group构建子图，并获取所有边的权重(freq/group freq)
        3. 获取所有相同起止node的path:
            (1) 4 node的paths： 处理exon错配偏移
            (2) 2/4 node的paths: 处理小exon跳跃
            (3) 3/5 node的paths: 处理小exon 凸起
            (4) 4/6 node的paths: 处理小exon 错配到两个exon

        """
        df = df.copy()
        df['exonChain2'] = df['exonChain'].apply(lambda x: list(map(int, x.split('-'))))

        ## 添加初始Tr起止位点
        extreme_values = (
            df.groupby(['Chr', 'Strand', 'Group'])['exonChain2']
            .agg(lambda x: (min([i for sublist in x for i in sublist]), 
                            max([i for sublist in x for i in sublist])))
            .apply(pd.Series)
            .rename(columns={0: 'min_value', 1: 'max_value'})
            .reset_index())

        df = df.merge(extreme_values, on=['Chr', 'Strand', 'Group'])

        df['exonChain_full'] = df.apply(
            lambda row: [row['min_value'] - 100] + row['exonChain2'] + [row['max_value'] + 100], axis=1
        )
        df = df.drop(columns = ['min_value','max_value'])

        df['exonChain_full_str'] = df['exonChain_full'].apply(lambda x: '-'.join(map(str, x)))


        df_list = [] # 存储df
        i = 0 # 记录每个子图的Group
        for _, df_group in df.groupby(['Chr', 'Strand', 'Group']):
            # 获取图边和权重
            dict_net = self.get_edges_weight_exonChainfull(df_group) # 收集所有边的总频率

            # print('====权重freq:',dict_net)

            # 构建图
            edges_list = list(dict_net.values())
            edges_list = [tuple(l) for l in edges_list]
            G = nx.DiGraph()  # 有向图
            G.add_weighted_edges_from(edges_list)


            i = 1
            ############## 4. 处理每一个子图
            for component in nx.weakly_connected_components(G):
                subgraph = G.subgraph(component).copy()

                # 获取子图df
                subgraph_nodes = set(list(subgraph.nodes()))
                subgraph_index = []  # 子图的行索引
                for index,row in df_group.iterrows():
                    exonChain2 = set(row['exonChain2'])
                    if len(subgraph_nodes & exonChain2) >0:
                        subgraph_index.append(index)
                # 子图df
                df_subgraph = df_group.loc[subgraph_index]

                ################
                ## 重新分配group
                ################
                df_subgraph['Group'] = i
                i += 1

                ###############
                ## NNC图过滤错配
                ###############
                df_subgraph = self.remove_misalignment(df_subgraph,subgraph)


                ##########
                # 汇总子图
                df_list.append(df_subgraph)

        df = pd.concat(df_list,ignore_index=True)
        return df


    def nnc_nic_graph(self, df):
        # 按Chr对df进行分组
        dfChr_list = [dfChr for _,dfChr in df.groupby('Chr')]

        partial_func = partial(self.nnc_nic_graph_forChr)

        with Pool(self.num_processes) as pool:
            results = pool.map(partial_func, dfChr_list)

        # 合并结果
        df_result = pd.concat(results,ignore_index= True).reset_index(drop=True)

        return df_result