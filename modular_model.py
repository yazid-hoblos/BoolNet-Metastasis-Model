from metastasisModel import MetastasisModel

reduced_model = MetastasisModel(modular=True)
# mutated_model =reduced_model.control(frozentrue={'Metastasis','Invasion','Migration','EMT','CCA'}, frozenfalse={'Apoptosis'})
# mutated_model.get_stable_states_df(display=True)
# reduced_model.draw_interaction_graph('reduced_model', split=True, interactive=False, show=True)
# reduced_model.necessary(trueset={'Metastasis','Invasion','Migration','EMT','CCA'}, falseset={'Apoptosis'},max_cnf=20000000, trace=True) 

# print(f"Number of Modules: {len(reduced_model.variables)}")

# reduced_model.draw_interaction_graph('reduced_model', split=True, interactive=True, show=True)
# reduced_model.get_stable_states_df(display=True)
# reduced_model.identify_active_nodes()

# reduced_model.write_stable_states('reduced_model')

# reduced_model.plot_stable_states('reduced_model', show=True)

# mutated_model = reduced_model.control(frozenfalse={'p53','Ecadh'})
# mutated_model.get_stable_states_df(display=True)
# mutated_model.write_stable_states('second_mutated_model')

# mutated_model.plot_stable_states('second_mutated_model')

# reduced_model.controllability_analysis('reduced_model', prevent_duplicates=True, plot=False)

# import pickle
# with open('model.pkl', 'rb') as f:
#     model = pickle.load(f)

# def synchronous(variables: list | set) -> frozenset:
#     return frozenset({frozenset({*variables})})
# eqs=reduced_model.model.equilibria(model=model,trace=True)[1]

# for eq in eqs:
#     print(eq, len(eq))

# with open('eqs.txt', 'w') as f:
    # for eq in eqs:
        # f.write(f"{eq}\n")