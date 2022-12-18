import random
import pandas as pd
import numpy as np
import networkx as nx
from sklearn.model_selection import KFold
from tqdm import trange
from joblib import Parallel, delayed

# 根据 SIR 模型，更新节点的状态; beta: 被感染的概率, gamma:恢复的概率
def updateNodeState(G, node, beta, gamma, day):
    if G.nodes[node]["state"] == "I": #感染者
        #p = random.random() # 生成一个0到1的随机数
        if day > gamma:   # gamma的概率恢复
            G.nodes[node]["state"] = "R" #将节点状态设置成“R”
    elif G.nodes[node]["state"] == "S": #易感者
        p = random.random() # 生成一个0到1的随机数
        k = 0  # 计算邻居中的感染者数量
        for neibor in G.adj[node]: # 查看所有邻居状态，遍历邻居用 G.adj[node]
            if G.nodes[neibor]["state"] == "I": #如果这个邻居是感染者，则k加1
                k = k + 1
        if p < 1 - (1 - beta)**k:  # 易感者被感染
            G.nodes[node]["state"] = "I"

#网络状态进行更新
def updateNetworkState(G, beta, gamma,days):
    for node in G: #遍历图中节点，每一个节点状态进行更新
        updateNodeState(G,node, beta, gamma, days[node])

# 计算三类人群的数量；S: 易感者；I: 感染者; R: 恢复者（第三个输出参数）
def countSIR(G):
    S = 0;I = 0
    for node in G:
        if G.nodes[node]["state"] == "S":
            S = S + 1
        elif G.nodes[node]["state"] == "I":
            I = I + 1
    return S,I, len(G.nodes) - S - I

def find_I(G):
    nodes_I = []
    nodes_R = []
    for node in G:
        if G.nodes[node]["state"] == "I":
            nodes_I.append(node)
        elif G.nodes[node]["state"] == "R":
            nodes_R.append(node)
    return nodes_I, nodes_R

#设置不同人群的显示颜色，易感者为橘色，感染者为红色，恢复者为绿色
color_dict = {"S":"orange","I":"red","R":"green"}
def get_node_color(G): #返回每一个节点的颜色组成的列表
    color_list = []
    for node in G:
        #使用我们前面创建的状态到颜色的映射字典 color_dict
        color_list.append(color_dict[G.nodes[node]["state"]])
    return color_list

#返回网络节点的度
def degreeNetwork(G):
    neighbors = []
    for node in G.nodes(): #遍历图中节点，每一个节点状态进行更新
        neighbors.append(G.degree(node))
    return neighbors

#引入感染者，增加密接网络.kI:感染者数目, kCC：平均密接数目
def addClose_Contacts(tmpG,kI,kCC):
# import_inf：输入型感染者；CC： 密切接触者。
    from scipy import stats
    # num_CC  = list(np.rint(stats.norm.rvs(kCC,0.5,kI)).astype(int))
    CC = []
    import_inf = []
    for i in range(kI):
        # CC.append(random.sample(list(tmpG.nodes()),num_CC[i]))
        CC.append(random.sample(list(tmpG.nodes()),kCC))
    for i in range(kI):
        import_inf.append(max(list(tmpG.nodes()))+1)
        tmpG.add_node(max(list(tmpG.nodes()))+1)
        tmpG.nodes[max(list(tmpG.nodes()))]["state"] = "I"
        tmpG.add_edges_from(list(zip([max(list(tmpG.nodes()))] * len(CC[i]),CC[i])))
    return import_inf, CC


# 简单核酸检测，仅返回阳性节点
def simple_nucleic_acid_testing(G,Fn,test_node):
    #Fn:假阴性率;p_nat:参加核酸检测的人员比率,
    #infected_node:被检测出来的阳性节点
    infected_node = []
    for node in test_node:
        if G.nodes[node]["state"] == "I":
            infected_node.append(node)
    if len(infected_node)>0:
        randnum = np.random.rand(len(infected_node))
        randnum = np.random.rand(len(infected_node))
        del_ind = list(np.where(randnum<Fn)[0])
        #print(del_ind)
        infected_node = [n for i,n in enumerate(infected_node) if i not in del_ind]
    return infected_node


#比较概率抽样检测与常态化7天1检
#固定网络
#当gamma大于指定天数，感染者恢复
def compare_top_multiply_p_sample_simple_test(nodenum=50000,avg_degree=4,beta=0.1,
                            gamma=7,FN=0.3,s_p = 1/7,infect_num=1):
    ba = nx.barabasi_albert_graph(nodenum,int(avg_degree/2))
    for mynode in ba.nodes():
        ba.nodes[mynode]["state"] = "S"
    days_infect = [0]*nodenum
    #引入感染者
    import_inf = []
    CC = []
    import_inf, CC = addClose_Contacts(ba,infect_num,avg_degree)
    days_infect = days_infect+[1]*infect_num
    days_infect = np.array(days_infect)
    updateNetworkState(ba,beta,gamma,list(days_infect))
    mynodes_I,mynodes_R = find_I(ba)
    days_infect[mynodes_I] = days_infect[mynodes_I]+1
#    mynodes_I,mynodes_R = find_I(ba)
    #SIR_list = []
    #SIR_list.append(list(countSIR(ba)))
    loops = 3 #3轮循环，每轮循环7天
    infected_node1 = []
    count_simple = 0
    count_simple_days = 0
    num_inf_simple = []
    start_simple = True
    start_hub = True
    infected_node2 = []
    count_hub = 0
    count_hub_days = 0
    num_inf_hub = []
    n_d = dict(nx.degree(ba))
    degrees = np.array(list(n_d.values()))
    degrees = degrees/degrees.sum()#按度 设置权重
    nodes = list(n_d.keys())
    #sample_num = int(len(nodes)/7)#暂定人数和常态化一样
    sample_num = int(len(nodes)*s_p)
    #a = np.arange(0, 7)
    #tmp_p = np.cumsum(a/np.sum(a))
    selected_nodes = []
    #last_selected_nodes = pd.DataFrame(columns = ['nodes','test_oder'])
    test_days = np.zeros(len(nodes),dtype = int)
    weights = np.zeros(len(nodes))
    for t in range(0,loops):
        random.shuffle(nodes)
        kf = KFold(n_splits=7)# 常态化7天一检
        for tmptr, tmpte in kf.split(nodes):
            if start_simple:
                test_node = list(np.array(nodes)[tmpte])
                count_simple += len(test_node)
                count_simple_days += 1
                infected_node1 = simple_nucleic_acid_testing(ba,FN,test_node)
                if len(infected_node1) > 0:
                    start_simple = False
                    #df = pd.DataFrame(SIR_list,columns=["S","I","R"])
                    #num_inf_simple = list(df['I'])[-1] + list(df['R'])[-1]
                    num_inf_simple = len(mynodes_I)+len(mynodes_R)
            if start_hub:
                count_hub_days += 1
                if count_hub_days == 1:
                    selected_nodes = list(np.argsort(degrees)[int(-1*sample_num):])
                    test_days[selected_nodes] = 1
                else:
                    tmp_w = test_days - 1
                    weights += (1/21)*tmp_w
                    weights = np.where(tmp_w==0, 0, weights)
                    weights = np.where(tmp_w<0, 1, weights)
                    Integrated_w = (weights/weights.sum())*degrees
                    Integrated_w = Integrated_w/Integrated_w.sum()
                    selected_nodes = list(np.argsort(Integrated_w)[int(-1*sample_num):])
                    test_days = np.where(test_days>=1, test_days+1, test_days)
                    test_days = np.where(test_days>=7, 0, test_days)
                    test_days[selected_nodes] = 1
                #tmp_snode = pd.DataFrame(selected_node,columns=['nodes'])
                #tmp['test_oder'] = np.zeros(len(selected_node),dtype=int)+2
#                last_selected_nodes = last_selected_nodes.append(pd.DataFrame(zip(selected_nodes,
#                                                                                  list(np.zeros(len(selected_nodes),dtype=int)+count_hub_days)),columns=['nodes','test_oder']))
                count_hub += len(selected_nodes)
                infected_node2 = simple_nucleic_acid_testing(ba,FN,selected_nodes)
                if len(infected_node2) > 0:
                    start_hub = False
                    #df = pd.DataFrame(SIR_list,columns=["S","I","R"])
                    #num_inf_hub = list(df['I'])[-1] + list(df['R'])[-1]
                    num_inf_hub = len(mynodes_I)+len(mynodes_R)
            updateNetworkState(ba,beta,gamma,list(days_infect))
            mynodes_I,mynodes_R = find_I(ba)#list
            days_infect[mynodes_I] = days_infect[mynodes_I]+1
            #updateNetworkState(ba,beta,gamma)
            #SIR_list.append(list(countSIR(ba)))
            if not(start_simple) and not(start_hub):
                break
        if not(start_simple) and not(start_hub):
            break
    tmprlt = [count_simple,count_simple_days,num_inf_simple,count_hub,count_hub_days,num_inf_hub]
    # rlt = []
    # rlt.append(tmprlt)
    # rlt = pd.DataFrame(rlt)
    # rlt.columns=['norm_total','norm_days','norm_I','hub_total','hub_days','hub_I']
    return tmprlt

if __name__ == "__main__":
    # nodenum=50000
    avg_degree=4
    beta=0.1
    gamma = 7
    s_p = 1/14
    for nodenum in [100000,500000,1000000,5000000,10000000]:
        #tested_nodes = pd.DataFrame(columns=['nodes','test_oder','repeat_id'])
        tmp = Parallel(n_jobs=10)(
            delayed( compare_top_multiply_p_sample_simple_test)(
                nodenum=nodenum,
                avg_degree=avg_degree,
                beta=beta,
                gamma=gamma,
                FN=0.3,s_p = s_p
            )  for _ in trange(1000)
        )
        rlt = pd.DataFrame(tmp,columns=['norm_total','norm_days','norm_I','hub_total','hub_days','hub_I'])
        rlt.to_csv(f'Number_of_people_{nodenum}_{avg_degree}_{beta}_{gamma}_0.3_{int(1/s_p)}.csv')