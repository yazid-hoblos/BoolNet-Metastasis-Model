from utils import *
from metastasisModel import MetastasisModel

full_model = MetastasisModel(modular=False) 
# print(full_model)

# full_model.draw_interaction_graph('full_model', split=True, interactive=True, show=True)
full_model.get_stable_states_df(display=True)
# full_model.identify_active_nodes()

# full_model.write_stable_states('full_model')

# full_model.plot_stable_states('full_model', show=True)

# mutated_model = full_model.control(frozentrue={'CTNNB1'})
# mutated_model.get_stable_states_df(display=True)
# mutated_model.write_stable_states('mutated_model')

# mutated_model.plot_stable_states('mutated_model', show=False)