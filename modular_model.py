from metastasisModel import MetastasisModel

reduced_model = MetastasisModel(modular=True)
print(f"Number of Modules: {len(reduced_model.variables)}")

reduced_model.draw_interaction_graph('reduced_model', split=True, interactive=True, show=True)
reduced_model.get_stable_states_df(display=True)
reduced_model.identify_active_nodes()

reduced_model.write_stable_states('reduced_model')

# reduced_model.plot_stable_states('reduced_model', show=True)