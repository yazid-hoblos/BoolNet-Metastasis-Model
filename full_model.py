from tabulate import tabulate

from utils import *
import pandas as pd 
from metastasisModel import MetastasisModel

full_model = MetastasisModel(modular=False) 

draw_network_interactive(full_model.model.interaction_graph, filename='full_model_interactive')
draw_seperately(full_model.model.interaction_graph, 'full_model')
draw_interaction_graph(full_model.model.interaction_graph, 'full_model')

# print(met_model)
stable = full_model.model.stable_states
# print(tabulate(stable, headers='keys', tablefmt='dpsl'))
handle_input_variables(stable, full_model.variables)
            
df = pd.DataFrame(stable)
df_T = df.T
# print(tabulate(df_T, headers='keys', tablefmt='dpsl'))

df = full_model.identify_stable_states(df_T)
df = rearrange_columns(df)

identify_active_nodes(df)

# df.astype(int).to_csv('data_files/stable_states.csv', index=True, header=True)
# plot_stable_states(df, 'model', show=True)

# boon.control(frozenfalse={DNAdamage},frozentrue={ECM})