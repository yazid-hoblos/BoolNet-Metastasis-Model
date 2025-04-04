from BooN import *
from utils import *
from sympy import symbols
from tabulate import tabulate
import pandas as pd

AKT1, AKT2, Ecadh, WNT_pthw, ERK_pthw, GF, miRNA, Notch_pthw, p53, p63_73, TGFb_pthw, EMTreg, CCA, Apoptosis, EMT, Invasion, Migration, Metastasis, DNAdamage, ECMicroenv = symbols(
    'AKT1 AKT2 Ecadh WNT_pthw ERK_pthw GF miRNA Notch_pthw p53 p63_73 TGFb_pthw EMTreg CCA Apoptosis EMT Invasion Migration Metastasis DNAdamage ECMicroenv')

reduced_model = BooN(descriptor={
    AKT1: WNT_pthw & (Notch_pthw | TGFb_pthw | GF | EMTreg) & ~miRNA & ~p53 & ~Ecadh,
    AKT2: (TGFb_pthw | GF | Notch_pthw | EMTreg) & EMTreg & ~miRNA & ~p53,
    Ecadh: ~AKT2 & ~EMTreg,
    WNT_pthw: ~Notch_pthw & ~EMTreg & ~miRNA & ~p53 & ~p63_73 & ~AKT1 & ~Ecadh & ~WNT_pthw,
    ERK_pthw: (TGFb_pthw | Notch_pthw | GF | EMTreg) & ~AKT1,
    GF: (GF | EMTreg) & ~Ecadh,
    miRNA: (p53 | p63_73) & ~AKT2 & ~EMTreg & ~AKT1,
    Notch_pthw: ECMicroenv & ~p53 & ~p63_73 & ~miRNA,
    p53: (Notch_pthw | DNAdamage | WNT_pthw) & ~AKT1 & ~AKT2 & ~p63_73 & ~EMTreg,
    p63_73: ~Notch_pthw & ~p53 & DNAdamage & ~AKT2 & ~AKT1 & ~EMTreg,
    TGFb_pthw: (Notch_pthw | ECMicroenv) & ~WNT_pthw & ~miRNA,
    EMTreg: (Notch_pthw | WNT_pthw | EMTreg) & ~miRNA & ~p53,
    CCA: (((p53 | p63_73 | (TGFb_pthw & Notch_pthw) | AKT2) & ~ERK_pthw) | miRNA | EMTreg) & ~AKT1,
    Apoptosis: ~ERK_pthw & ~AKT1 & ~EMTreg & (miRNA | p63_73 | p53),
    EMT: ~Ecadh & EMTreg,
    Invasion: (TGFb_pthw & EMTreg) | WNT_pthw,
    Migration: EMT & ERK_pthw & AKT2 & Invasion & ~AKT1 & ~miRNA & ~p63_73,
    Metastasis: Migration,
    DNAdamage: DNAdamage,
    ECMicroenv: ECMicroenv
})

stable = reduced_model.stable_states

print(len(reduced_model.variables))

handle_input_variables(stable, reduced_model.variables)

df = pd.DataFrame(stable)
df_T = df.T
# print(tabulate(df_T, headers='keys', tablefmt='dpsl'))

state_names = []
for col in df_T.columns:
    state = df_T[col]
    if state[Ecadh] and state.sum() == 1:
        name = "HS"
    elif state[Apoptosis]:
        if state[p53]:
            name = "Apop1" if state[ECMicroenv] == False else "Apop3"
        else:
            name = "Apop2" if state[ECMicroenv] == False else "Apop4"
    elif state[EMT] and not state[Metastasis]:
        name = "EMT1" if state[DNAdamage] == True else "EMT2"
    else:
        name = "M1" if state[DNAdamage] == True else "M2"
    state_names.append(name)
df_T.columns = state_names
df_T.index.name = 'Variables'

df_T = rearrange_columns(df_T)

df_T.astype(int).to_csv('data_files/reduced_model_stable_states.csv', index=True, header=True)

# plot_stable_states(df_T, 'reduced_model')
