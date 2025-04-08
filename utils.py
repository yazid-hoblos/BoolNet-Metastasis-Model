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

    phenotypes = ['Metastasis', 'Apoptosis', 'CellCycleArrest', 'CCA', 'Invasion', 'Migration', 'EMT']
    input_nodes = ['DNAdamage', 'ECM', 'ECMicroenv', 'GF']
    for node in graph.nodes():
        if str(node) in phenotypes:
            nx.draw_networkx_nodes(graph, pos, nodelist=[node], node_color="grey", node_size=1000, edgecolors="black")
        elif str(node) in input_nodes:
            nx.draw_networkx_nodes(graph, pos, nodelist=[node], node_color="yellow", node_size=1000, edgecolors="black")
        else:
            nx.draw_networkx_nodes(graph, pos, nodelist=[node], node_color="lightblue", node_size=1000, edgecolors="black")
    
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

    phenotypes = ['Metastasis', 'Apoptosis', 'CellCycleArrest', 'CCA', 'Invasion', 'Migration', 'EMT']
    input_nodes = ['DNAdamage', 'ECM', 'ECMicroenv', 'GF']
    
    # G_activation.remove_nodes_from(list(nx.isolates(G_activation)))
    # G_inhibition.remove_nodes_from(list(nx.isolates(G_inhibition)))

    # pos = nx.spring_layout(G)  
    pos = nx.kamada_kawai_layout(G)

    def get_node_colors(graph):
        colors = []
        for node in graph.nodes():
            if str(node) in phenotypes:
                colors.append("grey")  # Phenotypes in red
            elif str(node) in input_nodes:
                colors.append("yellow")  # Input nodes in yellow
            else:
                colors.append("lightblue")  # Other nodes in light blue
        return colors

    node_colors = get_node_colors(G) 
        
    fig, ax = plt.subplots(1, 2, figsize=(21, 10))

    nx.draw(G_activation, pos, ax=ax[0], with_labels=True, edge_color="green", node_color=node_colors)
    ax[0].set_title("Activation Network")

    nx.draw(G_inhibition, pos, ax=ax[1], with_labels=True, edge_color="red", node_color=node_colors)
    ax[1].set_title("Inhibition Network")

    plt.savefig(f'plots/{name}_split_interactions_network.png', dpi=300)
    if show:
        plt.show()
    plt.close()


def draw_network_interactive(G, name='network_visualization'):
    print("Drawing interactive network...")
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
        phenotypes = ['Metastasis', 'Apoptosis', 'CellCycleArrest', 'CCA', 'Invasion', 'Migration', 'EMT']
        input_nodes = ['DNAdamage', 'ECM', 'ECMicroenv', 'GF']
        
        if str(node) in phenotypes:
            color = "#ff0000"  # Red for phenotypes
            size = 35
        elif str(node) in input_nodes:
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
    # net.save_graph(output_path)
    # print(f"Interactive network saved to {output_path}")
    net.show(os.path.join(os.getcwd(), output_path), notebook=False)
    # Optional: Display network in notebook if running in jupyter
    # from IPython import get_ipython
    # ipython = get_ipython()
    # notebook = True #if ipython is not None and 'IPKernelApp' in ipython.config else False
    # if notebook:
    #     from IPython.display import IFrame, display, HTML
    #     display(HTML(f'<a href="{output_path}" target="_blank">Open Network Visualization</a>'))
    #     display(IFrame(output_path, width="100%", height=750))
    # display(IFrame(output_path, width=800, height=600))
    
    
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

    blue = plt.cm.RdBu_r(0.01)
    red = plt.cm.RdBu_r(0.9)   

    colors = ['white', 'green', red, blue]  # None, ON, OFF, Both
    cmap = ListedColormap(colors)

    plt.figure(figsize=(12, 12))
    ax = sns.heatmap(df_combined, cmap=cmap, linewidths=0.5, linecolor='gray',
                    cbar=False, square=True)

    legend_patches = [
        mpatches.Patch(color='white', label='None'),
        mpatches.Patch(color='green', label='ON'),
        mpatches.Patch(color=red, label='OFF'),
        mpatches.Patch(color=blue, label='ON/OFF')
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


def evaluate_rule(rule_expr, input_values):
    import re
    """
    Evaluates a boolean rule expression like "WNT_pthw & (Notch_pthw | TGFb_pthw | GF | EMTreg) & ~miRNA & ~p53 & ~Ecadh"
    
    Parameters:
    -----------
    rule_expr : str
        Boolean rule expression
    input_values : dict
        Dictionary mapping node names to their boolean values
    
    Returns:
    --------
    bool : Result of rule evaluation
    """
    # Replace operators with Python equivalents
    rule_py = rule_expr.replace('&', ' and ').replace('|', ' or ').replace('~', ' not ')
    
    # Replace node names with their values
    for node, value in input_values.items():
        node_str = str(node)
        # Replace whole words only to avoid partial matches
        rule_py = re.sub(r'\b' + re.escape(node_str) + r'\b', str(value), rule_py)
    
    # Replace any remaining node names (not in input_values) with False
    # This handles nodes that might be in the rule but not in the predecessors
    remaining_nodes = re.findall(r'\b[A-Za-z0-9_]+\b', rule_py)
    for node in remaining_nodes:
        if node not in ('True', 'False', 'and', 'or', 'not'):
            rule_py = re.sub(r'\b' + node + r'\b', 'False', rule_py)
    
    try:
        # Evaluate the transformed expression
        return bool(eval(rule_py))
    except Exception as e:
        print(f"Error evaluating rule: {rule_expr} -> {rule_py}")
        print(f"Exception: {e}")
        return False
    
    
def simulate_boolean_network_for_gephi(G, model, initial_state=None, steps=100, output_dir='gephi_data'):
    
    import os
    import random
    import csv
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    s=model.symbols
    # Initialize state randomly if not provided
    if initial_state is None:
        initial_state = {}
        for node in G.nodes():
            initial_state[node] = bool(random.getrandbits(1))
    
    state_history = [initial_state.copy()]
    current_state = initial_state.copy()
    stable_state_reached = False
    stable_at_step = None
    
    print("Running simulation...")
    for step in range(steps):
        next_state = {}
        
        for node in G.nodes():
            node_str = s[str(node)]
            if node_str in model.model.desc:
                rule = str(model.model.desc[node_str])
                input_values = {str(n): current_state.get(n, False) for n in G.nodes()}
                
                try:
                    next_state[node] = evaluate_rule(rule, input_values)
                except Exception as e:
                    print(f"Error evaluating rule for {node}: {e}")
                    next_state[node] = current_state.get(node, False)
            else:
                next_state[node] = current_state.get(node, False)
        
        if next_state == current_state:
            stable_state_reached = True
            stable_at_step = step
            print(f"Stable state reached at step {step}")
            # print all nodes that are active in the stable state
            active_nodes = [str(node) for node, state in next_state.items() if state]
            print(len(active_nodes), "active nodes in stable state")
            print(f"Active nodes in stable state: {', '.join(active_nodes)}")
            break
        
        cycle_detected = False
        for i, prev_state in enumerate(state_history):
            if next_state == prev_state:
                cycle_start = i
                cycle_length = len(state_history) - i
                print(f"Cycle detected: States repeat from step {cycle_start} with period {cycle_length}")
                cycle_detected = True
                break
        
        if cycle_detected:
            break
        
        current_state = next_state.copy()
        state_history.append(current_state.copy())
    
    print(f"Simulation completed with {len(state_history)} states")
    
    # Create node attributes file for Gephi
    nodes_file = os.path.join(output_dir, 'nodes.csv')
    with open(nodes_file, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Write header
        header = ['Id', 'Label', 'Type']
        # Add columns for each time step
        for t in range(len(state_history)):
            header.append(f'State_t{t}')
        writer.writerow(header)
        
        for node in G.nodes():
            node_str = str(node)
            
            phenotypes = ['Metastasis', 'Apoptosis', 'CellCycleArrest', 'CCA', 'Invasion', 'Migration', 'EMT']
            input_nodes = ['DNAdamage', 'ECM', 'ECMicroenv', 'GF']
            
            if node_str in phenotypes:
                node_type = 'Phenotype'
            elif node_str in input_nodes:
                node_type = 'Input'
            else:
                node_type = 'Regular'
                
            # Create row with node data
            row = [node_str, node_str, node_type]
            
            # Add state at each time step
            for state in state_history:
                row.append(1 if state.get(node, False) else 0)
                
            writer.writerow(row)
    
    # Create edges file for Gephi
    edges_file = os.path.join(output_dir, 'edges.csv')
    with open(edges_file, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Write header
        writer.writerow(['Source', 'Target', 'Type', 'Sign', 'Weight'])
        
        # Write data for each edge
        for u, v, data in G.edges(data=True):
            edge_type = 'Directed'
            sign = data.get('sign', 0)
            weight = 1.0
            
            writer.writerow([str(u), str(v), edge_type, sign, weight])
    
    # Create summary file with simulation info
    summary_file = os.path.join(output_dir, 'simulation_info.txt')
    with open(summary_file, 'w') as f:
        f.write(f"Boolean Network Simulation Summary\n")
        f.write(f"================================\n\n")
        f.write(f"Total time steps: {len(state_history)}\n")
        
        if stable_state_reached:
            f.write(f"Stable state reached at step: {stable_at_step}\n")
        elif cycle_detected:
            f.write(f"Cycle detected starting at step {cycle_start} with period {cycle_length}\n")
        else:
            f.write(f"No stable state or cycle detected within {steps} steps\n")
            
        f.write("\nInitial state:\n")
        for node, state in initial_state.items():
            f.write(f"  {node}: {1 if state else 0}\n")
            
        if stable_state_reached:
            f.write("\nFinal stable state:\n")
            for node, state in state_history[-1].items():
                f.write(f"  {node}: {1 if state else 0}\n")
    
    timeseries_file = os.path.join(output_dir, 'timeseries.csv')
    with open(timeseries_file, 'w', newline='') as f:
        writer = csv.writer(f)
        
        header = ['Timestamp']
        header.extend([str(node) for node in G.nodes()])
        writer.writerow(header)
        
        for t, state in enumerate(state_history):
            row = [t]
            for node in G.nodes():
                row.append(1 if state.get(node, False) else 0)
            writer.writerow(row)
    
    # print(f"Files exported to {output_dir} directory:")
    # print(f"  - nodes.csv: Node data with states at each time step")
    # print(f"  - edges.csv: Edge relationships")
    # print(f"  - simulation_info.txt: Summary of the simulation")
    # print(f"  - timeseries.csv: Time series data for dynamic visualization")
    # print("\nImport instructions for Gephi:")
    # print("1. Import nodes.csv as nodes table")
    # print("2. Import edges.csv as edges table")
    # print("3. Use the 'State_t#' columns for dynamic coloring of nodes")
    # print("4. For dynamic visualization, use the timeseries.csv with the Dynamic Filter plugin")

# Helper function to evaluate boolean rules
def evaluate_rule(rule_expr, input_values):
    """
    Evaluates a boolean rule expression like "WNT_pthw & (Notch_pthw | TGFb_pthw | GF | EMTreg) & ~miRNA & ~p53 & ~Ecadh"
    
    Parameters:
    -----------
    rule_expr : str
        Boolean rule expression
    input_values : dict
        Dictionary mapping node names to their boolean values
    
    Returns:
    --------
    bool : Result of rule evaluation
    """
    import re
    
    # Replace operators with Python equivalents
    rule_py = rule_expr.replace('&', ' and ').replace('|', ' or ').replace('~', ' not ')
    
    # Replace node names with their values
    for node, value in input_values.items():
        node_str = str(node)
        # Replace whole words only to avoid partial matches
        rule_py = re.sub(r'\b' + re.escape(node_str) + r'\b', str(value), rule_py)
    
    # Replace any remaining node names (not in input_values) with False
    # This handles nodes that might be in the rule but not in the predecessors
    remaining_nodes = re.findall(r'\b[A-Za-z0-9_]+\b', rule_py)
    for node in remaining_nodes:
        if node not in ('True', 'False', 'and', 'or', 'not'):
            rule_py = re.sub(r'\b' + node + r'\b', 'False', rule_py)
    
    try:
        # Evaluate the transformed expression
        return bool(eval(rule_py))
    except Exception as e:
        print(f"Error evaluating rule: {rule_expr} -> {rule_py}")
        print(f"Exception: {e}")
        return False
    




from pyvis.network import Network
import pandas as pd
import networkx as nx

def simulate_evolution(model):
    simulate_boolean_network_for_gephi(model.model.interaction_graph, model, steps=1000, output_dir='gephi_data')
    # Load data
    df = pd.read_csv('gephi_data/timeseries.csv')
    timestamps = df['Timestamp']
    df_states = df.drop(columns='Timestamp')

    # Network
    G = model.model.interaction_graph
    net = Network(height='750px', width='100%', directed=True, notebook=False)
    
    import networkx as nx
    pos = nx.kamada_kawai_layout(G)
    # for node in G.nodes():
    #     init_state = int(df_states.iloc[0][str(node)])
    #     color = 'limegreen' if init_state == 1 else 'lightgray'
    #     x, y = pos[node]
    #     net.add_node(str(node), label=str(node), color=color, x=x*1000, y=y*1000, fixed=True)

    for node in G.nodes():
        x,y= pos[node]
        net.add_node(str(node), label=str(node),x=x*500, y=y*500, fixed=True)

    for source, target,data in G.edges(data=True):
        if 'sign' in data:
            if data['sign'] > 0:
                net.add_edge(str(source), str(target), color='green') 
            else:
                net.add_edge(str(source), str(target), color='red')

    # Generate states JS data
    js_data = []
    for i, row in df_states.iterrows():
        state = {str(k): int(v) for k, v in row.items()}
        js_data.append(state)

    # Add JavaScript to control state changes
    net.set_options("""
    var options = {
      "nodes": {
        "shape": "dot",
        "size": 25,
        "font": { "size": 14, "color": "#343434" }
      },
      "edges": {
        "arrows": { "to": { "enabled": true } },
        "color": { "inherit": true },
        "smooth": false
      },
      "physics": {
        "enabled": true,
        "stabilization": { "iterations": 100 }
      }
    }
    """)

    # Save and manually inject custom JS
    net.show("network_simulation.html", notebook=False)

    # Inject JS for slider and animation
    with open("network_simulation.html", "r") as f:
        html = f.read()

    import json
    js_json = json.dumps(js_data)

    js_script = f"""
    <script type="text/javascript">
    let states = {js_json};
    let curr = 0;
    let interval = null;

    function updateNodes(stateIndex) {{
        let state = states[stateIndex];
        for (let nodeId in state) {{
            let newColor = state[nodeId] === 1 ? 'limegreen' : 'lightgray';
            nodes.update({{id: nodeId, color: {{background: newColor}} }});
        }}
        document.getElementById("stateSlider").value = stateIndex;
        document.getElementById("stepLabel").innerText = "Time Step: " + stateIndex;
    }}

    function play() {{
        if (interval) return;
        interval = setInterval(() => {{
            curr = (curr + 1) % states.length;
            updateNodes(curr);
        }}, 1000);
    }}

    function pause() {{
        clearInterval(interval);
        interval = null;
    }}

    function setupControls() {{
        // Slider
        let slider = document.createElement("input");
        slider.type = "range";
        slider.min = 0;
        slider.max = states.length - 1;
        slider.value = 0;
        slider.id = "stateSlider";
        slider.oninput = function() {{
            curr = parseInt(this.value);
            updateNodes(curr);
        }};

        // Label
        let label = document.createElement("div");
        label.id = "stepLabel";
        label.style = "margin-top:10px; font-weight:bold;";
        label.innerText = "Time Step: 0";

        // Play/Pause buttons
        let playBtn = document.createElement("button");
        playBtn.innerText = "▶ Play";
        playBtn.onclick = play;
        playBtn.style = "margin-right: 10px;";

        let pauseBtn = document.createElement("button");
        pauseBtn.innerText = "⏸ Pause";
        pauseBtn.onclick = pause;

        // Container
        let container = document.createElement("div");
        container.style = "padding: 10px;";
        container.appendChild(slider);
        container.appendChild(label);
        container.appendChild(playBtn);
        container.appendChild(pauseBtn);

        document.body.insertBefore(container, document.body.firstChild);
        updateNodes(0);
    }}

    window.addEventListener('load', setupControls);
    </script>
    """
    # Inject JS before </body>
    html = html.replace("</body>", js_script + "\n</body>")

    # Save updated HTML
    with open("network_simulation.html", "w") as f:
        f.write(html)

    print("✅ Interactive simulation saved to network_simulation.html")

def git_simulation(model):
    import pandas as pd
    import networkx as nx
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from metastasisModel import MetastasisModel

    df = pd.read_csv('gephi_data/timeseries.csv')
    timestamps = df['Timestamp']
    df_states = df.drop(columns='Timestamp')

    G = model.model.interaction_graph

    pos = nx.kamada_kawai_layout(G) 
    fig, ax = plt.subplots(figsize=(12, 8))

    for node, edge, data in G.edges(data=True):
        if data['sign'] > 0:
            G.edges[node, edge]['sign'] = 1
        else:
            G.edges[node, edge]['color'] = -1
    edge_colors = ['green' if G.edges[node, edge]['sign'] > 0 else 'red' for node, edge in G.edges()]


    def update(frame):
        ax.clear()
        state = df_states.iloc[frame]
        colors = ['limegreen' if state[str(n)] == 1 else 'lightgray' for n in G.nodes]
        
        nx.draw(
            G, pos, ax=ax, with_labels=True, node_color=colors,
            node_size=800, font_size=10, edge_color=edge_colors
        )
        ax.set_title(f"Network State at Timestamp {timestamps[frame]}", fontsize=14)


    ani = animation.FuncAnimation(fig, update, frames=len(df_states), interval=1000, repeat=True)
    plt.tight_layout()
    # save the animation as a gif
    ani.save('network_animation.gif', writer='pillow', fps=1)
    # plt.show()
