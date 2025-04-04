from sympy import symbols
from BooN import *

class Model:
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
        
        self.model = BooN(descriptor=desc)
        self.variables = self.model.variables
    
    def __str__(self):
        return self.model.__str__()
     
    def identify_stable_states(self, df):
        state_names = []
        s=self.symbols
        CDH = s['CDH1'] if s['CDH1'] in df.index.tolist() else s['Ecadh']
        ECM_ = s['ECM'] if s['ECM'] in df.index.tolist() else s['ECMicroenv']
        for col in df.columns:
            state = df[col]
            if state[CDH] and state.sum() == 1:
                name = "HS"
            elif state[s['Apoptosis']]:
                if state[s['p53']]:
                    name = "Apop1" if not state[ECM_] else "Apop3"
                else:
                    name = "Apop2" if not state[ECM_] else "Apop4"
            elif state[s['EMT']] and not state[s['Metastasis']]:
                name = "EMT1" if state[s['DNAdamage']] else "EMT2"
            else:
                name = "M1" if state[s['DNAdamage']] else "M2"
            state_names.append(name)
        df.columns = state_names
        df.index.name = 'Variables'
        return df

            
            