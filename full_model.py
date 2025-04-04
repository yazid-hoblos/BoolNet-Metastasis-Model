from utils import *
from metastasisModel import MetastasisModel

full_model = MetastasisModel(modular=False) 
print(f"Number of Modules: {len(full_model.variables)}")

full_model.draw_interaction_graph('full_model', split=True, interactive=True, show=True)
full_model.get_stable_states_df(display=True)
full_model.identify_active_nodes()

full_model.write_stable_states('full_model')

full_model.plot_stable_states('full_model', show=True)

# boon.control(frozenfalse={DNAdamage},frozentrue={ECM})