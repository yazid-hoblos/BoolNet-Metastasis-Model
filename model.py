import matplotlib.pyplot as plt
import networkx as nx
from sympy.abc import w, x, y, z, v
from tabulate import tabulate
from BooN import *

from sympy import symbols  

AKT1, AKT2, CDH1, CDH2, CTNNB1, DKK1, ERK, GF, miR200, miR203, miR34, NICD, p21, p53, p63, p73, SMAD, SNAI1, SNAI2, TGFbeta, TWIST1, VIM, ZEB1, ZEB2, CellCycleArrest, Apoptosis, EMT, Invasion, Migration, Metastasis, DNAdamage, ECM = symbols(
    'AKT1 AKT2 CDH1 CDH2 CTNNB1 DKK1 ERK GF miR200 miR203 miR34 NICD p21 p53 p63 p73 SMAD SNAI1 SNAI2 TGFbeta TWIST1 VIM ZEB1 ZEB2 CellCycleArrest Apoptosis EMT Invasion Migration Metastasis DNAdamage ECM')

boon = BooN({
    AKT1: CTNNB1 & (NICD | TGFbeta | GF | CDH2) & ~p53 & ~miR34 & ~CDH1,
    AKT2: TWIST1 & (TGFbeta | GF | CDH2) & ~(miR203 | miR34 | p53),
    CDH1: ~TWIST1 & ~SNAI2 & ~ZEB1 & ~ZEB2 & ~SNAI1 & ~AKT2,
    CDH2: TWIST1,
    CTNNB1: ~DKK1 & ~p53 & ~AKT1 & ~miR34 & ~miR200 & ~CDH1 & ~CDH2 & ~p63,
    DKK1: CTNNB1 | NICD,
    ERK: (SMAD | CDH2 | GF | NICD) & ~AKT1,
    GF: ~CDH1 & (GF | CDH2),
    miR200: (p63 | p53 | p73) & ~(AKT2 | SNAI1 | SNAI2 | ZEB1 | ZEB2),
    miR203: p53 & ~(SNAI1 | ZEB1 | ZEB2),
    miR34: ~(SNAI1 | ZEB1 | ZEB2) & (p53 | p73) & AKT2 & ~p63 & ~AKT1,
    NICD: ~p53 & ~p63 & ~p73 & ~miR200 & ~miR34 & ECM,
    p21: ((SMAD & NICD) | p63 | p53 | p73 | AKT2) & ~(AKT1 | ERK),
    p53: (DNAdamage | CTNNB1 | NICD | miR34) & ~SNAI2 & ~p73 & ~AKT1 & ~AKT2,
    p63: DNAdamage & ~NICD & ~AKT1 & ~AKT2 & ~p53 & ~miR203,
    p73: DNAdamage & ~p53 & ~ZEB1 & ~AKT1 & ~AKT2,
    SMAD: TGFbeta & ~miR200 & ~miR203,
    SNAI1: (NICD | TWIST1) & ~miR203 & ~miR34 & ~p53 & ~CTNNB1,
    SNAI2: (TWIST1 | CTNNB1 | NICD) & ~miR200 & ~p53 & ~miR203,
    TGFbeta: (ECM | NICD) & ~CTNNB1,
    TWIST1: CTNNB1 | NICD | SNAI1,
    VIM: CTNNB1 | ZEB2,
    ZEB1: ((TWIST1 & SNAI1) | CTNNB1 | SNAI2 | NICD) & ~miR200,
    ZEB2: (SNAI1 | (SNAI2 & TWIST1) | NICD) & ~miR200 & ~miR203,
    CellCycleArrest: (miR203 | miR200 | miR34 | ZEB2 | p21) & ~AKT1,
    Apoptosis: (p53 | p63 | p73 | miR200 | miR34) & ~ZEB2 & ~AKT1 & ~ERK,
    EMT: CDH2 & ~CDH1,
    Invasion: (SMAD & CDH2) | CTNNB1,
    Migration: VIM & AKT2 & ERK & ~miR200 & ~AKT1 & EMT & Invasion & ~p63,
    Metastasis: Migration,
    DNAdamage: False,
    ECM: ECM
})

import pandas as pd
# compute the stable states
print("- STABLE STATES -")
stable = boon.stable_states
print(len(stable))
for s in stable:
    for k, v in s.items():
        if v == True:
            s[k] = 1
        elif v == False:
            s[k] = 0
# print(tabulate(stable, headers='keys', tablefmt='dpsl'))

boon = boon.cnf()
df = pd.DataFrame(stable)
df_T = df.T
print(tabulate(df_T, headers='keys', tablefmt='dpsl'))
