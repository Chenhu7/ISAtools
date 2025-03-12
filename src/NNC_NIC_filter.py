import pandas as pd
import networkx as nx
from functools import partial
from multiprocessing import Pool


class NncNicGraphProcessor:
    def __init__(self, exon_excursion_diff_bp=10, error_sites_diff_bp=10, error_sites_multiple=0.01,little_exon_bp=30, 
                 little_exon_mismatch_diff_bp=10,Nolittle_exon_mismatch_diff_bp=20,little_exon_jump_multiple=0.1, num_processes=10):
        self.exon_excursion_diff_bp = exon_excursion_diff_bp  # exon偏移的差值
        self.error_sites_diff_bp = error_sites_diff_bp # 错误位点的差值
        self.error_sites_multiple = error_sites_multiple # 错误位点的比值
        self.little_exon_bp = little_exon_bp # 小exon的长度
        # self.little_exon_jump_multiple = little_exon_jump_multiple
        self.little_exon_mismatch_diff_bp = little_exon_mismatch_diff_bp # 小exon错配差值
        self.Nolittle_exon_mismatch_diff_bp = Nolittle_exon_mismatch_diff_bp # 非小exon错配差值
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
        """纠正错误的外显子链, 纠正后的exonChain保证在df的exonChain中"""
        exon_chain = row['exonChain']
        exon_chain_list = set(exon_chain_list)
        
        for error_path, correct_path in dict_correct.items():
            if error_path in exon_chain:
                corrected = exon_chain.replace(error_path, correct_path)
                if corrected in exon_chain_list:
                    exon_chain = corrected

        return exon_chain

    
    def corrected_df(self, df1, dict_correct):
        """更新纠正后的DataFrame"""
        df = df1.copy()
        exon_chain_list = df['exonChain'].tolist()
        df['exonChain_corrected'] = df.apply(lambda row: self.correct_exonChain(row, dict_correct, exon_chain_list), axis=1)

        df = df.explode(['TrStart', 'TrEnd'], ignore_index=True)
        # df = df.groupby(['Chr', 'Strand', 'Group', 'exonChain_corrected', 'junction', 'contains_non_conservative'],
        #                 as_index=False).agg({
        #                     'TrStart': list,
        #                     'TrEnd': list
        #                     # 'category': lambda x: sorted(x)[0]
        #                 }).reset_index(drop=True)
        
        df = df.groupby(['Chr', 'Strand', 'Group', 'exonChain_corrected'],
                        as_index=False).agg({
                            'TrStart': list,
                            'TrEnd': list
                            # 'category': lambda x: sorted(x)[0]
                        }).reset_index(drop=True)
        
        df.rename(columns={'exonChain_corrected': 'exonChain'}, inplace=True)
        df['exonChain2'] = df['exonChain'].apply(lambda x: list(map(int, x.split('-'))))
        df['frequency'] = df['TrStart'].apply(len)
        
        # 添加junction
        df = df.merge(df1[['Chr', 'Strand', 'exonChain', 'junction']],on = ['Chr', 'Strand', 'exonChain'],how = 'left')
        
        return df
    
    
    def get_df_exonChain_StartEnd(self,df):
        """
        对分组内的exonChain添加虚拟起止位点  极值±100bp
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

        return df
    

    def remove_misalignment(self, df_subgraph,subgraph):
        ### 1. 获取中间段低频exon
        out_nodes = [node for node, out_degree in subgraph.out_degree() if out_degree > 1] # 起点 出度大于一的所有节点
        in_nodes = [node2 for node2, in_degree in subgraph.in_degree() if in_degree > 1] # 终点 入度大于一的所有节点

        print_w = 1 # 调试输出开关 0为输出
        
        # 存储所有起点-终点路径
        forked_paths = []

        # 获取所有起点到终点的所有路径
        for out_node in out_nodes:
            for in_node in in_nodes:
                all_path = list(nx.all_simple_paths(subgraph, source=out_node, target=in_node))
                forked_paths += all_path

        # 所有候选路径
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
                    if '-'.join(map(str,path)) in exonChain:
                        path_weight += dict_subgraph_exonChain_freq[exonChain] / subgraph_total_freq # path的权重为freq/group freq

                #直接将权重添加到path列表的末尾
                dict_path[combine][i].append(path_weight)
        
        if print_w == 0:
            print('============')
            print(dict_path)
            print('============')

        ############
        ## exon错配
        ############
        correct_paths = {} # 纠正

        ### 1. 偏移错配: 偏移量差值小于10bp内 / 含有非保守剪接碱基  只针对四个位点的path
        dict_path_filtered = {key: [v for v in value if len(v) == 5] for key, value in dict_path.items()}
        dict_path_filtered = {k: v for k, v in dict_path_filtered.items() if len(v) > 1}
        
        if print_w == 0:
            print('======偏移错配======')
            print(dict_path_filtered)
            # print('============')

        if len(dict_path_filtered) > 0:
            for combine in dict_path_filtered:
                paths = dict_path_filtered[combine]
                paths = sorted(paths, key=lambda x: x[-1],reverse=True)

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
                            error_exon = '-'.join(path_str.split('-')[1:-1])
                            ref_exon = '-'.join(ref_path_str.split('-')[1:-1])
                            correct_paths[error_exon] = ref_exon
                            # print(f'偏移纠正,{path_str} --> {ref_path_str}')
        
        if print_w == 0:
            print(correct_paths)

        
        ###############
        ##  错误位点
        ##############
        # 直接过滤掉 "相差10bp内"  "低频"的单个错误位点exonChain
        if print_w == 0:
            print('======单个错误位点======')
        # display(df_subgraph)
        # print('============')
        
        # 获取所有位点的权重
        sites_dict = {}
        for _,row in df_subgraph.iterrows():
            exonChain2 = row['exonChain2']
            freq = row['frequency']
            
            for site in exonChain2:
                if site not in sites_dict:
                    sites_dict[site] = freq
                else:
                    sites_dict[site] += freq
        
        # 排序并筛选相差10bp内的位点
        sites_dict_keys = sorted(list(sites_dict.keys()))
        sites_dict_keys = list(zip(sites_dict_keys, sites_dict_keys[1:]))
        
        # 筛选候选错误位点元组
        potential_error_sites = [sites for sites in sites_dict_keys if abs(sites[1] - sites[0]) <= self.error_sites_diff_bp]
        
        error_sites = []
        for p_sites in potential_error_sites:
            site1_freq = sites_dict[p_sites[0]]
            site2_freq = sites_dict[p_sites[1]]
            
            if site1_freq < site2_freq and site1_freq/site2_freq <= self.error_sites_multiple:
                error_sites.append(p_sites[0])
                print('p_sites[0]',p_sites[0])
            elif site1_freq > site2_freq and site2_freq/site1_freq <= self.error_sites_multiple:
                error_sites.append(p_sites[1])
                print('p_sites[0]',p_sites[1])
        
        # 删除exonChain
        remove_index = []
        for index,row in df_subgraph.iterrows():
            if len(set(error_sites) & set(row['exonChain2'])) > 0:
                remove_index.append(index)
        if len(remove_index) > 0:
            df_subgraph.drop(remove_index, inplace=True)
        
        if print_w == 0:
            print(error_sites)
        
        
        ###############
        ##  小exon错配
        ##############
        # 添加起止位点
        df_subgraph = self.get_df_exonChain_StartEnd(df_subgraph)
        
        # correct_paths = {} # 纠正
        ## (1) 先找小exon  先筛选出含有小exon的path  4node path 和 6node path
        node6_path = [] # 含有小exon的6个位点path
        node4_path = [] # 含有小exon的4个位点path
        
        for exonChain in df_subgraph['exonChain_full'].tolist():
            if len(exonChain) >= 4:
                for i in range(0,len(exonChain),2):
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
        
        ## 小exon跳跃,补全小exon, 如果不补全后的exonChain存在于df中,则纠正,否则保留原有小exon跳跃的exonChain
        if len(node4_path) > 0:
            for Lexon in node4_path:
                error_path = '-'.join([str(s) for s in [Lexon[0],Lexon[-1]]])
                correct_path = '-'.join([str(s) for s in Lexon])
                correct_paths[error_path] = correct_path
            
            if print_w == 0:
                print('======小exon跳跃======')
                print(correct_paths)
        
        
        ## 小exon错配, 错配差值小于<10bp, 且纠正后的exonChain在df中
        if len(node6_path) > 0:
            for N6_path in node6_path:
                for _,row in df_subgraph.iterrows():
                    row_exonChainFull = row['exonChain_full']
                    if N6_path[0] in row_exonChainFull and N6_path[-1] in row_exonChainFull:
                        row_path = row_exonChainFull[row_exonChainFull.index(N6_path[0]):row_exonChainFull.index(N6_path[-1]) + 1]
                        if len(set(row_path) - set(N6_path)) > 0 and len(row_path) == 4:
                            N6_path_length = sum([N6_path[i+1] - N6_path[i] for i in range(0,len(N6_path),2)])
                            row_path_length = sum([row_path[i+1] - row_path[i] for i in range(0,len(row_path),2)])

                            # 长度差值小于10 即纠正
                            if abs(N6_path_length - row_path_length) <= self.little_exon_mismatch_diff_bp:
                                error_path = '-'.join([str(row_path[1]),str(row_path[2])])
                                correct_path = '-'.join([str(s) for s in N6_path[1:-1]])
                                # print(N6_path,'--->',row_path)
                                # print(error_path,'--->',correct_path)
                                correct_paths[error_path] = correct_path
                            
        if print_w == 0:
            print('===小exon错配===')
            print(correct_paths)
        
        ###############
        ##  非小exon错配
        ##############
        # 长度差值<20 & 不包含GT-AG & 纠正后在df中
        no_little_exon_correct = {} # 暂时的错配纠正

        ## 1.记录只有4个和6个位点的path
        forked_paths_4_6 = [path[:-1] for path in forked_paths_filtered if len(path) in [5,7]]
        
        # 相同起止节点以及所有路径 相同起止位点至少2个路径
        dict_path46 = {}
        for path in forked_paths_4_6:
            start_end = (path[0],path[-1])
            if start_end not in dict_path46:
                # dict_path46[start_end] = ['-'.join([str(p) for p in path])]
                dict_path46[start_end] = [path]
            else:
                # dict_path46[start_end].append('-'.join([str(p) for p in path]))
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
        
        # 判断是否含有非GT-AG
        for No_Lexon_error in no_little_exon_correct:
            for _,row in df_subgraph.iterrows():
                if No_Lexon_error in row['exonChain']:
                    junctions = row['junction'].split(',')
                    if not all(x.upper() == 'GT-AG' for x in junctions):
                        correct_paths[No_Lexon_error] = no_little_exon_correct[No_Lexon_error]
        if print_w == 0:
            print('===非小exon错配===')
            print(correct_paths)
                        
        # 矫正
        # print('==============2',correct_paths)
        if len(correct_paths) > 0:
            df_subgraph = self.corrected_df(df_subgraph,correct_paths)
        else:
            # 和已纠正的df保持一致
            # df_subgraph = df_subgraph[['Chr','Strand','Group','exonChain','TrStart','TrEnd','category','exonChain2','frequency','junction','contains_non_conservative']]
            df_subgraph = df_subgraph[['Chr','Strand','Group','exonChain','TrStart','TrEnd','exonChain2','frequency', 'junction']]

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
        
        df = self.get_df_exonChain_StartEnd(df)

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