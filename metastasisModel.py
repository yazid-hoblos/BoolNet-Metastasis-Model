from sympy import symbols
from BooN import BooN

from utils import *

class MetastasisModel:
    def __init__(self, modular=False):
        
        AKT1, AKT2, CDH1, CDH2, CTNNB1, DKK1, ERK, GF, miR200, miR203, miR34, NICD, p21, p53, p63, p73, SMAD, SNAI1, SNAI2, TGFbeta, TWIST1, VIM, ZEB1, ZEB2, CellCycleArrest, Apoptosis, EMT, Invasion, Migration, Metastasis, DNAdamage, ECM = symbols(
        'AKT1 AKT2 CDH1 CDH2 CTNNB1 DKK1 ERK GF miR200 miR203 miR34 NICD p21 p53 p63 p73 SMAD SNAI1 SNAI2 TGFbeta TWIST1 VIM ZEB1 ZEB2 CellCycleArrest Apoptosis EMT Invasion Migration Metastasis DNAdamage ECM')
        
        self.symbols = {'AKT1': AKT1, 'AKT2': AKT2, 'CDH1': CDH1, 'CDH2': CDH2, 'CTNNB1': CTNNB1,
            'DKK1': DKK1, 'ERK': ERK, 'GF': GF, 'miR200': miR200, 'miR203': miR203,
            'miR34': miR34, 'NICD': NICD, 'p21': p21, 'p53': p53, 'p63': p63, 'p73': p73,
            'SMAD': SMAD, 'SNAI1': SNAI1, 'SNAI2': SNAI2, 'TGFbeta': TGFbeta, 'TWIST1': TWIST1,
            'VIM': VIM, 'ZEB1': ZEB1, 'ZEB2': ZEB2, 'CellCycleArrest': CellCycleArrest,
            'Apoptosis': Apoptosis, 'EMT': EMT, 'Invasion': Invasion, 'Migration': Migration,
            'Metastasis': Metastasis, 'DNAdamage': DNAdamage, 'ECM': ECM,}
        
        if modular:
            Ecadh, WNT_pthw, ERK_pthw, miRNA, Notch_pthw, p63_73, TGFb_pthw, EMTreg, CCA, ECMicroenv = symbols('Ecadh WNT_pthw ERK_pthw miRNA Notch_pthw p63_73 TGFb_pthw EMTreg CCA ECMicroenv')
            
            self.symbols.update({'Ecadh': Ecadh, 'WNT_pthw': WNT_pthw, 'ERK_pthw': ERK_pthw, 'miRNA': miRNA,
                'Notch_pthw': Notch_pthw, 'p63_73': p63_73, 'TGFb_pthw': TGFb_pthw,
                'EMTreg': EMTreg, 'CCA': CCA, 'ECMicroenv': ECMicroenv})

            desc={AKT1: WNT_pthw & (Notch_pthw | TGFb_pthw | GF | EMTreg) & ~miRNA & ~p53 & ~Ecadh,
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
            ECMicroenv: ECMicroenv}
            
        else:
            desc={AKT1: CTNNB1 & (NICD | TGFbeta | GF | CDH2) & ~p53 & ~miR34 & ~CDH1,
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
            DNAdamage: DNAdamage,
            ECM: ECM}
        
        self.modular = modular
        self.model = BooN(descriptor=desc)
        self.variables = self.model.variables
        self.stable_states_df = None
    
    def __str__(self):
        output = f"Number of Variables: {len(self.variables)}\n"
        output += f"Variables: {self.variables}\n"
        output += self.model.__str__() 
        return output
     
    def identify_stable_states(self):
        state_names = []
        s=self.symbols
        df = self.stable_states_df
        CDH = s['CDH1'] if s['CDH1'] in df.index.tolist() else s['Ecadh']
        ECM = s['ECM'] if s['ECM'] in df.index.tolist() else s['ECMicroenv']
        CCA = s['CellCycleArrest'] if s['CellCycleArrest'] in df.index.tolist() else s['CCA']
        for col in df.columns:
            state = df[col]
            if state[CDH] and state.sum() == 1:
                name = "HS"
            elif state[s['Apoptosis']] and state[CCA]:
                if state[s['p53']]:
                    name = "Apop1" if not state[ECM] else "Apop3"
                else:
                    name = "Apop2" if not state[ECM] else "Apop4"
            elif state[s['EMT']] and not state[s['Metastasis']] and state[CCA]:
                name = "EMT1" if state[s['DNAdamage']] else "EMT2"
            elif state[s['Metastasis']] and state[s['Invasion']] and state[s['Migration']] and state[s['EMT']]:
                name = "M1" if state[s['DNAdamage']] else "M2"
            else:
                if state[s['Invasion']]:
                    name = "Inv"
                elif state[CCA]:
                    name = "CCA"
                else:
                    name = "Novel"
            state_names.append(name)
        
        name_counts = {}
        unique_names = []
        for name in state_names:
            if name in name_counts:
                name_counts[name] += 1
                unique_name = f"{name}_{name_counts[name]}"
            else:
                name_counts[name] = 0
                unique_name = name
            unique_names.append(unique_name)
            
        df.columns = unique_names
        self.stable_states_df = df
        return df

    def get_stable_states_df(self, display=False):
        import pandas as pd
        from tabulate import tabulate
        stable = self.model.stable_states
        handle_input_variables(stable, self.variables)    
        df = pd.DataFrame(stable)
        df = df.T
        df = rearrange_columns(df)
        self.stable_states_df = df
        df = self.identify_stable_states()
        df.index.name = 'Variables'
        if display:
            print(tabulate(df, headers='keys', tablefmt='dpsl'))
        return df
    
    def identify_active_nodes(self):
        if self.stable_states_df is not None:
            identify_active_nodes(self.stable_states_df)
        else:
            print("Please call get_stable_states_df() first")
            
    def write_stable_states(self, name):
        if self.stable_states_df is not None:
            self.stable_states_df.astype(int).to_csv(f"data_files/{name}_stable_states.csv", index=True, header=True)
        else:
            print("Please call get_stable_states_df() first")
            
    def plot_stable_states(self, name, show=False):
        if self.stable_states_df is not None:
            plot_stable_states(self.stable_states_df, name, show)
        else:
            print("Please call get_stable_states_df() first")
            
    def draw_interaction_graph(self, name, split=False, interactive=False, show=False):
        draw_interaction_graph(self.model.interaction_graph, name, show)
        if split:
            draw_act_inh_seperately(self.model.interaction_graph, name, show)
        if interactive:
            draw_network_interactive(self.model.interaction_graph, name)
            
    def control(self, frozenfalse=None, frozentrue=None):
        s=self.symbols
        ff=set([s[var] for var in frozenfalse]) if frozenfalse else set()
        ft=set([s[var] for var in frozentrue]) if frozentrue else set()
        model_replicate = MetastasisModel(modular=self.modular)
        for var in ff:
            model_replicate.model.desc[var]=False
        for var in ft:
            model_replicate.model.desc[var]=True
        return model_replicate
    
    def controllability_analysis(self, name):
        vars = sorted(self.variables, key=lambda x: str(x))
        with open(f'data_files/{name}_controllability_analysis.txt', 'w') as f:
            f.write("------------ Controllability Analysis ------------\n")
            for var in vars:
                LoF=list(self.control(frozenfalse={str(var)}).get_stable_states_df().columns)
                GoF=list(self.control(frozentrue={str(var)}).get_stable_states_df().columns)
                f.write(f"\t\t --------{var}---------\n")
                f.write(f"OFF: {LoF}\n")
                f.write(f"ON: {GoF}\n\n")
            
        


            
            