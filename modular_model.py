from BooN import *
from utils import *
from sympy import symbols
from tabulate import tabulate
import pandas as pd

AKT1, AKT2, Ecadh, WNT_pthw, ERK_pthw, GF, miRNA, Notch_pthw, p53, p63_73, TGFb_pthw, EMTreg, CCA, Apoptosis, EMT, Invasion, Migration, Metastasis, DNAdamage, ECMicroenv = symbols(
    'AKT1 AKT2 Ecadh WNT_pthw ERK_pthw GF miRNA Notch_pthw p53 p63_73 TGFb_pthw EMTreg CCA Apoptosis EMT Invasion Migration Metastasis DNAdamage ECMicroenv')

boon = BooN(descriptor={
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

stable = boon.stable_states

print(len(boon.variables))

all_nodes = set(boon.variables)
for i,state in enumerate(stable):
    for node in all_nodes:
        if node not in state:
            stable[i][node] = False
            new_state = state.copy()
            new_state[node] = True
            stable.append(new_state)

# print(tabulate(stable, headers='keys', tablefmt='dpsl'))
df = pd.DataFrame(stable)
df_T = df.T
print(tabulate(df_T, headers='keys', tablefmt='dpsl'))

stable_state_names=['HS','Apop2','Apop1','EMT2', 'EMT1','Apop4','Apop3','M2','M1']

plot_stable_states(df_T, stable_state_names, 'reduced_model')

