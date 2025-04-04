from metastasisModel import MetastasisModel

reduced_model = MetastasisModel(modular=True)
# print(f"Number of Modules: {len(reduced_model.variables)}")

# reduced_model.draw_interaction_graph('reduced_model', split=True, interactive=True, show=True)
# reduced_model.get_stable_states_df(display=True)
# reduced_model.identify_active_nodes()

# reduced_model.write_stable_states('reduced_model')

# reduced_model.plot_stable_states('reduced_model', show=True)

mutated_model = reduced_model.control(frozenfalse={'p53','Ecadh'})
mutated_model.get_stable_states_df(display=True)
mutated_model.write_stable_states('second_mutated_model')

mutated_model.plot_stable_states('second_mutated_model')