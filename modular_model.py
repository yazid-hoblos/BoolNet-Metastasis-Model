from BooN import *
from utils import *
from sympy import symbols
from tabulate import tabulate
import pandas as pd

from metastasisModel import MetastasisModel

reduced_model = MetastasisModel(modular=True)

# draw_seperately(reduced_model.model.interaction_graph, 'reduced_model', show=False)
# draw_network_interactive(reduced_model.model.interaction_graph, filename='reduced_model_interactive')
# draw_interaction_graph(reduced_model.model.interaction_graph, 'reduced_model', show=False)

stable = reduced_model.model.stable_states

print(f"Number of Modules: {len(reduced_model.variables)}")

handle_input_variables(stable, reduced_model.variables)

df = pd.DataFrame(stable)
df_T = df.T
# print(tabulate(df_T, headers='keys', tablefmt='dpsl'))

df = reduced_model.identify_stable_states(df_T)

df = rearrange_columns(df)

identify_active_nodes(df)

# df.astype(int).to_csv('data_files/reduced_model_stable_states.csv', index=True, header=True)
# plot_stable_states(df, 'reduced_model')
