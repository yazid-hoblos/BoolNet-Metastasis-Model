import pickle
import networkx as nx
with open('model.pkl', 'rb') as f:
    model = pickle.load(f)
print('model loaded')
scc = list(map(frozenset, nx.strongly_connected_components(model)))

comps=[]
for s in scc:
    if len(s) > 1:
        print(s, len(s))
        comps.append(s)
    # print(s, len(s))
    
print(len(scc), 'strongly connected components')
print(len(comps), 'strongly connected components with more than one node')
