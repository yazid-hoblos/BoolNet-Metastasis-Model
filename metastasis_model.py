from tabulate import tabulate

from utils import *
import pandas as pd 
from model import Model

met_model = Model(modular=False) 
     
# print(met_model)
stable = met_model.model.stable_states
# print(tabulate(stable, headers='keys', tablefmt='dpsl'))
handle_input_variables(stable, met_model.variables)
            
df = pd.DataFrame(stable)
df_T = df.T
# print(tabulate(df_T, headers='keys', tablefmt='dpsl'))

df = met_model.identify_stable_states(df_T)
df = rearrange_columns(df)

identify_active_nodes(df)

# df.astype(int).to_csv('data_files/stable_states.csv', index=True, header=True)
# plot_stable_states(df, 'model', show=True)

# boon.control(frozenfalse={DNAdamage},frozentrue={ECM})