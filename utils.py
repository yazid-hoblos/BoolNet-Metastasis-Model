from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
import numpy as np
from matplotlib import pyplot as plt
import networkx as nx

def plot_stable_states(df, name, show=False):
    df_sorted = df.copy()
    df_sorted = df_sorted.astype(int)
    # active_gene_counts = df_sorted.sum(axis=0)
    # col_order = active_gene_counts.sort_values().index    
    row_linkage = linkage(pdist(df_sorted), method='average')
    row_order = dendrogram(row_linkage, no_plot=True)['leaves']
    df_sorted = df_sorted.iloc[row_order, :]  
    # df_sorted = df_sorted.loc[:, col_order]   
    
    plt.figure(figsize=(14, 12))
    cmap = plt.cm.RdBu_r  
    plt.imshow(df_sorted, cmap=cmap, interpolation='none', aspect='auto', vmin=0, vmax=1)
    
    col_labels = [f"{df.columns[i]} ({df_sorted.sum(axis=0)[col]})" for i,col in enumerate(df_sorted.columns)]
    
    # plt.colorbar(label='Node State (0=inactive, 1=active)')
    plt.colorbar(label='Node State (0=inactive, 1=active)', ticks=[0, 1], orientation='vertical')
    plt.xticks(range(len(df_sorted.columns)), col_labels, rotation=90)
    plt.yticks(range(len(df_sorted.index)), df_sorted.index)
    plt.title('Heatmap of Boolean Network States (Sorted by Active Variables Count)')
    
    ax = plt.gca()
    ax.set_xticks(np.arange(-.5, len(df_sorted.columns), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(df_sorted.index), 1), minor=True)
    ax.grid(which='minor', color='black', linestyle='-', linewidth=0.5)
    
    plt.tight_layout()
    plt.savefig(f'plots/{name}_stable_states_heatmap.png', dpi=300)
    
    if show:
        plt.show()
    plt.close()
    
    plt.figure(figsize=(10, 8))
    activation_freq = df_sorted.mean(axis=1).sort_values(ascending=False)
    activation_freq.plot(kind='bar', color='darkred')
    plt.title('Node Activation Frequency Across States')
    plt.xlabel('Node')
    plt.ylabel('Fraction of States Active')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(f'plots/{name}_node_activation_frequency.png', dpi=300)
    if show:
        plt.show()
    plt.close()

def rearrange_columns(df):
    column_activity = {}
    for col in df.columns:
        column_activity[col] = df[col].sum()  
    sorted_columns = sorted(column_activity.keys(), key=lambda col: column_activity[col])
    df = df.reindex(columns=sorted_columns)
    return df

def handle_input_variables(stable, vars):  
    for i,state in enumerate(stable):
        for node in vars:
            if node not in state:
                stable[i][node] = False
                new_state = state.copy()
                new_state[node] = True
                stable.append(new_state)
                
def identify_active_nodes(df):
    print("\n--- Stable states sorted by number of active nodes ---")
    for col in df.columns:
        active_nodes = df[col].sum()
        active_node_names = [str(node) for node in df.index[df[col] == True]]
        print(f"State {col}: {active_nodes} active nodes - {', '.join(active_node_names)}")
        

def draw_interaction_graph(graph, name, show=False):    
    plt.figure(figsize=(21, 15))  
    # pos = nx.spring_layout(graph, seed=42)  
    pos = nx.kamada_kawai_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=1000, node_color="lightblue", edgecolors="black")

    activation_edges = []
    inhibition_edges = []
    
    for u, v, data in graph.edges(data=True):
        if 'sign' in data:
            if data['sign'] > 0:
                activation_edges.append((u, v)) 
            elif data['sign'] < 0:
                inhibition_edges.append((u, v)) 

    nx.draw_networkx_edges(graph, pos, edgelist=activation_edges, edge_color="green", arrows=True, width=2)
    nx.draw_networkx_edges(graph, pos, edgelist=inhibition_edges, edge_color="red", arrows=True, width=2, style="dashed")
    nx.draw_networkx_labels(graph, pos, font_size=10, font_weight="bold")

    plt.title("Influence Network")
    plt.savefig(f'plots/{name}_interactions_network.png', dpi=300)
    if show:
        plt.show()
    plt.close()

def draw_act_inh_seperately(G, name, show=False):
    G_activation = nx.DiGraph() 
    G_inhibition = nx.DiGraph()  

    G_activation.add_nodes_from(G.nodes())
    G_inhibition.add_nodes_from(G.nodes())

    for u, v, data in G.edges(data=True):
        if data['sign'] > 0:
            G_activation.add_edge(u, v, color="green")
        elif data['sign'] < 0:
            G_inhibition.add_edge(u, v, color="red")

    G_activation.remove_nodes_from(list(nx.isolates(G_activation)))
    G_inhibition.remove_nodes_from(list(nx.isolates(G_inhibition)))

    fig, ax = plt.subplots(1, 2, figsize=(21, 10))
    # pos = nx.spring_layout(G)  
    pos = nx.kamada_kawai_layout(G)

    nx.draw(G_activation, pos, ax=ax[0], with_labels=True, edge_color="green", node_color="lightblue")
    ax[0].set_title("Activation Network")

    nx.draw(G_inhibition, pos, ax=ax[1], with_labels=True, edge_color="red", node_color="lightblue")
    ax[1].set_title("Inhibition Network")

    plt.savefig(f'plots/{name}_split_interactions_network.png', dpi=300)
    if show:
        plt.show()
    plt.close()


def draw_network_interactive(G, name='network_visualization'):
    import os
    from pyvis.network import Network

    # os.makedirs('plots', exist_ok=True)
    
    net = Network(notebook=False, cdn_resources='remote', height="750px", width="100%", 
                  bgcolor="#e8f1ff", font_color="#000c1f", 
                  select_menu=True, filter_menu=True)
    
    net.toggle_hide_edges_on_drag(False)
    net.force_atlas_2based(spring_length=100, spring_strength=0.15, 
                           damping=0.9, gravity=-50)
    
    for node in G.nodes():
        phenotypes = ['EMT', 'Metastasis', 'Apoptosis', 'CellCycleArrest']
        special_nodes = ['p53', 'CDH1', 'DNAdamage', 'ECM', 'ECMicroenv']
        
        if str(node) in phenotypes:
            color = "#ff5500"  # Orange for phenotypes
            size = 35
        elif str(node) in special_nodes:
            color = "#0055ff"  # Blue for special nodes
            size = 30
        else:
            color = "#888888"  # Gray for other nodes
            size = 20 + 5 * G.degree[node]
            
        net.add_node(str(node), label=str(node), title=str(node), 
                     color=color, size=size, value=G.degree[node])
    
    for u, v, data in G.edges(data=True):
        if 'sign' in data:
            if data['sign'] > 0:
                edge_color = "#00aa00"  # Green
                title = f"{u} activates {v}"
                arrows = "to"
                group = "activation"
            else:
                edge_color = "#cc0000"  # Red
                title = f"{u} inhibits {v}"
                arrows = "to;dash"
                group = "inhibition"
                
            edge_width = 1 + 2 * np.sqrt(G.degree[u]) / 5
                
            net.add_edge(str(u), str(v), color=edge_color, title=title, 
                        width=edge_width, arrows=arrows, physics=True, group=group)

    output_path = f'plots/{name}_interactive_interactions_network.html'
    net.save_graph(output_path)
    print(f"Interactive network saved to {output_path}")
    
    # Optional: Display network in notebook if running in jupyter
    from IPython import get_ipython
    ipython = get_ipython()
    notebook = True if ipython is not None and 'IPKernelApp' in ipython.config else False
    if notebook:
        from IPython.display import IFrame, display, HTML
        display(HTML(f'<a href="{output_path}" target="_blank">Open Network Visualization</a>'))
        display(IFrame(output_path, width="100%", height=750))

