from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
import numpy as np
from matplotlib import pyplot as plt

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
    plt.savefig(f'plots/{name}_heatmap.png', dpi=300)
    
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