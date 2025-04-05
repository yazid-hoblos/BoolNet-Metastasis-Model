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
        phenotypes = ['Metastasis', 'Apoptosis', 'CellCycleArrest', 'CCA']
        special_nodes = ['Invasion', 'Migration', 'Proliferation', 'EMT']
        input_nodes = ['DNAdamage', 'ECM', 'ECMicroenv']
        
        if str(node) in phenotypes:
            color = "#ff0000"  # Red for phenotypes
            size = 35
        elif str(node) in special_nodes:
            color = "#0055ff"  # Blue for special nodes
            size = 30
        elif str(node) in input_nodes:
            color = "#ffff00" # Yellow for input nodes
            size = 25
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
    
    
def read_controllability_file(file):
    import re, ast
    with open(f'data_files/{file}_controllability_analysis.txt', 'r') as file:
        data = file.read()

    control_nodes = []
    states_affected = {}

    current_node = None
    for line in data.split('\n'):
        line = line.strip()
        if '--------' in line and 'Controllability Analysis' not in line:
            current_node = re.sub(r'-+', '', line).strip()
            control_nodes.append(current_node)
            states_affected[current_node] = {'OFF': [], 'ON': []}
        elif line.startswith('OFF:'):
            states = line.replace('OFF:', '').strip()
            states = ast.literal_eval(states) if states else []
            states_affected[current_node]['OFF'] = states
        elif line.startswith('ON:'):
            states = line.replace('ON:', '').strip()
            states = ast.literal_eval(states) if states else []
            states_affected[current_node]['ON'] = states

    all_states = sorted({state for v in states_affected.values() for mode in v.values() for state in mode})
    return control_nodes, states_affected, all_states

def plot_controllability_results(file):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import matplotlib.patches as mpatches
    from matplotlib.colors import ListedColormap

    control_nodes, states_affected, all_states = read_controllability_file(file)
    on_matrix = np.zeros((len(control_nodes), len(all_states)))
    off_matrix = np.zeros((len(control_nodes), len(all_states)))

    for i, node in enumerate(control_nodes):
        for j, state in enumerate(all_states):
            if state in states_affected[node]['ON']:
                on_matrix[i, j] = 1
            if state in states_affected[node]['OFF']:
                off_matrix[i, j] = 2

    control_matrix = np.zeros_like(on_matrix)
    for i in range(len(control_nodes)):
        for j in range(len(all_states)):
            has_on = on_matrix[i, j] == 1
            has_off = off_matrix[i, j] == 2
            control_matrix[i, j] = 3 if has_on and has_off else (1 if has_on else (2 if has_off else 0))

    df_combined = pd.DataFrame(control_matrix, index=control_nodes, columns=all_states)

    colors = ['white', 'green', 'red', 'purple']  # None, ON, OFF, Both
    cmap = ListedColormap(colors)

    plt.figure(figsize=(12, 12))
    ax = sns.heatmap(df_combined, cmap=cmap, linewidths=0.5, linecolor='gray',
                    cbar=False, square=True)

    legend_patches = [
        mpatches.Patch(color='white', label='None'),
        mpatches.Patch(color='green', label='ON'),
        mpatches.Patch(color='red', label='OFF'),
        mpatches.Patch(color='purple', label='ON/OFF')
    ]
    plt.legend(handles=legend_patches, title="Effect Type",
            bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0.)

    plt.title('States Resulting from Nodes Control', fontsize=16)
    plt.xlabel('Phenotypic States', fontsize=14)
    plt.ylabel('Controlled Node', fontsize=14)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f'plots/{file}_controllability_analysis_heatmap.png', dpi=300)
    
def plot_controllability_bar_plot(filename):
    control_nodes, states_affected, all_states = read_controllability_file(filename)
    
    on_matrix = np.zeros((len(control_nodes), len(all_states)))
    off_matrix = np.zeros((len(control_nodes), len(all_states)))

    for i, node in enumerate(control_nodes):
        for j, state in enumerate(all_states):
            if state in states_affected[node]['ON']:
                on_matrix[i, j] = 1
            if state in states_affected[node]['OFF']:
                off_matrix[i, j] = 1
    
    on_counts = [sum(on_matrix[i]) for i in range(len(control_nodes))]
    off_counts = [sum(off_matrix[i]) for i in range(len(control_nodes))]

    x = np.arange(len(control_nodes))
    width = 0.35

    fig, ax = plt.subplots(figsize=(14, 8))
    rects1 = ax.bar(x - width/2, on_counts, width, label='ON effects', color='green', alpha=0.7)
    rects2 = ax.bar(x + width/2, off_counts, width, label='OFF effects', color='red', alpha=0.7)

    ax.set_ylabel('Number of Reached States')
    ax.set_title('Number of States Reached by Controlling Nodes (one at a time)')
    ax.set_xticks(x)
    ax.set_xticklabels(control_nodes, rotation=45, ha='right')
    ax.legend()

    def autolabel(rects):
        for rect in rects:
            height = rect.get_height()
            ax.annotate(f'{int(height)}',xy=(rect.get_x() + rect.get_width() / 2, height),xytext=(0, 3), textcoords="offset points",
                        ha='center', va='bottom')
    autolabel(rects1)
    autolabel(rects2)
    fig.tight_layout()
    plt.savefig(f'plots/{filename}_controllability_analysis_bar_plot.png', dpi=300)


