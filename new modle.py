#%%

import numpy as np
import scipy
import matplotlib.pyplot as plt
import pandas as pd
import math
import os
from sympy import symbols, Eq, solve, nsolve, re,S, Abs
from matplotlib import colors
colors.XKCD_COLORS
from matplotlib import patheffects
from scipy.optimize import curve_fit
from matplotlib import cm, ticker
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm
from matplotlib import font_manager
from matplotlib import rcParams
from matplotlib.patches import FancyBboxPatch
from scipy.optimize import minimize
from prettytable import PrettyTable


#%%
DIC = 2.5e-3       # mol/L
Ca_tot = 10.52e-3    
Mg_tot = 54.13e-3
k_co2 = 10**(-6.4)
k_co3 = 10**(-10.3)

#%%
#original ocean concentration:
i = 8.1
H = 10 ** (-i)
OH = 10 ** (-14) / H
CO2_0 = DIC / (1 + k_co2 / H + k_co2*k_co3 / (H) ** 2)
HCO3_0 = DIC / (1 + H / k_co2 + k_co3 / H)
CO3_0 = DIC / (1 + H / k_co3 + (H) ** 2 / k_co2*k_co3)

print(CO2_0)
print(HCO3_0)
print(CO3_0)

#%%

#Dolomite(ordered) CaMg(CO3)2
ksp_camgco3_2 = 10**(-17.09)

#Brucite Mg(OH)2
ksp_mgoh2 = 10**(-11.14)

#Artinite Mg2(CO3)(OH)2
ksp_mg_2co3oh_2 = 10**(-18.4)
 
#Portlandite Ca(OH)2
ksp_caoh2 = 10**(-5.18)

#Hydromagnesite Mg5(CO3)4(OH)2•4H2O
ksp_mg5co34oh2 = 10**(-36.8)

#Hunitite Mg3Ca(CO3)4
ksp_mg3caco34 = 10**(-29.967)

#Calcite CaCO3
ksp_caco3 = 10**(-8.48)

#Magnesite MgCO3
ksp_mgco3 = 10**(-8.03)

F = 96485.3321 #s A / mol

pH = []
CO2_conc= []
HCO3_conc= []
CO3_conc = []
Mg_conc = []
Ca_conc= []

CaMgCO32_conc= []
CaCO3_conc = []
MgCO3_conc = []
MgOH2_conc= []
Mg2CO3OH2_conc= []
CaOH2_conc = []
Mg5CO34OH2_conc = []

Q_list= []
Total_Pre_conc = []
CO3_Pre_conc = []
OH_Pre_conc = []
Carbon_release = []



#%%
for i in np.arange (8.8,9.5,0.1):
    H = 10 ** (-i)
    OH = 10 ** (-14) / H
    CO2_i = DIC / (1 + k_co2 / H + k_co2*k_co3 / (H) ** 2)
    HCO3_i = DIC / (1 + H / k_co2 + k_co3 / H)
    CO3_i = DIC / (1 + H / k_co3 + (H) ** 2 / k_co2*k_co3)
    
    n0_camgco3_2 = Ca_tot * Mg_tot * (CO3_i) ** 2 - ksp_camgco3_2
    n0_mgoh2 = Mg_tot* (OH**2) -ksp_mgoh2
    n0_mg_2co3oh_2 = Mg_tot**2 *CO3_i* (OH**2) - ksp_mg_2co3oh_2
    n0_caoh2 = Ca_tot * (OH**2) -ksp_caoh2
    n0_mg5co34oh2 = Mg_tot**5 *CO3_i**4 *OH**2 - ksp_mg5co34oh2
    n0_mg3caco34 = Mg_tot**3*Ca_tot * CO3_i**4 -ksp_mg3caco34
    n0_caco3 = Ca_tot*CO3_i - ksp_caco3
    n0_mgco3 = Mg_tot *CO3_i - ksp_mgco3
    values = {
    'n0_camgco3_2': n0_camgco3_2,
    'n0_mgoh2': n0_mgoh2,
    'n0_mg_2co3oh_2': n0_mg_2co3oh_2,
    'n0_caoh2': n0_caoh2,
    'n0_mg5co34oh2': n0_mg5co34oh2,
    'n0_mg3caco34': n0_mg3caco34,
    'n0_caco3': n0_caco3,
    'n0_mgco3': n0_mgco3
}
    print(i)
    table = PrettyTable(['i', 'name', 'value'])

    for i, (name, value) in enumerate(values.items()):
        if value >= 0:
            table.add_row([i, name, value])
    print(table)


# %%


df = pd.DataFrame(columns=['pH', 'CO2', 'HCO3', 'CO3', 'Ca', 'Mg', 'CaMgCO32', 'MgOH2', 'CaCO3', 'CaOH2', 'Mg5CO34OH2', 'Mg3CaCO34', 'MgCO3'])

for i in np.arange (0,14.01,0.001):
    H = 10 ** (-i)
    OH = 10 ** (-14) / H
    CO2_i = DIC / (1 + k_co2 / H + k_co2*k_co3 / (H) ** 2)
    HCO3_i = DIC / (1 + H / k_co2 + k_co3 / H)
    CO3_i = DIC / (1 + H / k_co3 + (H) ** 2 / k_co2*k_co3)

    n0_camgco3_2 = Ca_tot * Mg_tot * (CO3_i) ** 2 - ksp_camgco3_2
    n0_mgoh2 = Mg_tot* (OH**2) -ksp_mgoh2
    n0_mg_2co3oh_2 = Mg_tot**2 *CO3_i* (OH**2) - ksp_mg_2co3oh_2
    n0_caoh2 = Ca_tot * (OH**2) -ksp_caoh2
    n0_mg5co34oh2 = Mg_tot**5 *CO3_i**4 *OH**2 - ksp_mg5co34oh2
    n0_mg3caco34 = Mg_tot**3*Ca_tot * CO3_i**4 -ksp_mg3caco34
    n0_caco3 = Ca_tot*CO3_i - ksp_caco3
    n0_mgco3 = Mg_tot *CO3_i - ksp_mgco3



    if n0_camgco3_2 >0 :
        print(i,"CaMgCO32 appears")
        CO2, HCO3, CO3, Ca, Mg, CaMgCO32= symbols('CO2 HCO3 CO3 Ca Mg CaMgCO32' )
        eq1= Eq(k_co2 * CO2 , H* HCO3)
        eq2 = Eq(k_co3 * HCO3, H* CO3)
        eq3 = Eq(ksp_camgco3_2, Ca* Mg * CO3**2)
        eq4 = Eq(Ca_tot, Ca + CaMgCO32)
        eq5 = Eq(Mg_tot, Mg + CaMgCO32)
        eq6 = Eq(DIC, CO2 + HCO3 + CO3 + 2*CaMgCO32)
        sol = solve((eq1 , eq2, eq3, eq4, eq5, eq6), (CO2, HCO3, CO3, Ca, Mg, CaMgCO32), dict = True)

        if sol:  
            print(i,"single CaMgCO32 has solution")
            real_solutions = [{k: re(v) for k, v in solution.items()} for solution in sol]
            positive_solutions = [solution for solution in real_solutions if all(v > 0 for v in solution.values())]
            if positive_solutions :
                print(i, "single CaMgCO32 has positive solution")
                solution = positive_solutions[0]
                print(i,solution)

                Data = {
                    "pH": i,
                    "CO2": solution[CO2],
                    "HCO3": solution[HCO3],
                    "CO3": solution[CO3],
                    "Ca": solution[Ca],
                    "Mg": solution[Mg],
                    "CaMgCO32": solution[CaMgCO32],
                    "MgOH2": 0,
                    "CaCO3": 0,
                    "CaOH2": 0,
                    "Mg5CO34OH2": 0,
                    "Mg3CaCO34": 0,
                    "MgCO3": 0,
                }
                df.loc[len(df)] = Data
                # print(i,solution)
                n1_mg2co3oh2 = solution[Mg]**2 * solution[CO3] * (OH**2) - ksp_mg_2co3oh_2
                n1_caoh2 = solution[Ca] * (OH**2) - ksp_caoh2
                n1_mgoh2 = solution[Mg] * (OH**2) - ksp_mgoh2
                n1_mg5co34oh2 = solution[Mg]**5 * solution[CO3]**4 * OH**2 - ksp_mg5co34oh2
                n1_mg3caco34 = solution[Mg]**3 * solution[Ca] * solution[CO3]**4 - ksp_mg3caco34
                n1_caco3 = solution[Ca] * solution[CO3] - ksp_caco3
                n1_mgco3 = solution[Mg] * solution[CO3] - ksp_mgco3
                values = {
                    "n1_mg2co3oh2": n1_mg2co3oh2,
                    "n1_caoh2": n1_caoh2,
                    "n1_mgoh2": n1_mgoh2,
                    "n1_mg5co34oh2": n1_mg5co34oh2,
                    "n1_mg3caco34": n1_mg3caco34,
                    "n1_caco3": n1_caco3,
                    "n1_mgco3": n1_mgco3
                        }
                for name, value in values.items():
                    if value > 0:
                        print(i, name, value)
                
                # result: mgoh2  is the next precipitation
                if n1_mgoh2 >0:
                    print(i,"CaMgCO32 +MgOH2")
                    df = df.drop(df.index[-1]) 
                    CO2, HCO3, CO3, Ca, Mg, CaMgCO32, MgOH2= symbols('CO2 HCO3 CO3 Ca Mg CaMgCO32 MgOH2' )
                    eq1= Eq(k_co2 * CO2 , H* HCO3)
                    eq2 = Eq(k_co3 * HCO3, H* CO3)
                    eq3 = Eq(ksp_camgco3_2, Ca* Mg * CO3**2)
                    eq4 = Eq(ksp_mgoh2, Mg* (OH**2))
                    eq5 = Eq(Ca_tot, Ca + CaMgCO32)
                    eq6 = Eq(Mg_tot, Mg + CaMgCO32 + MgOH2)
                    eq7 = Eq(DIC, CO2 + HCO3 + CO3 + 2*CaMgCO32)
                    sol = solve((eq1 , eq2, eq3, eq4, eq5, eq6, eq7), (CO2, HCO3, CO3, Ca, Mg, CaMgCO32, MgOH2 ), dict = True)
                    if sol:
                        real_solutions = [{k: re(v) for k, v in solution.items()} for solution in sol]
                        positive_solutions = [solution for solution in real_solutions if all(v > 0 for v in solution.values())]
                        if positive_solutions :
                            solution = positive_solutions[0]
                            print(i,solution)
                            Data = {
                                "pH": i,
                                "CO2": solution[CO2],
                                "HCO3": solution[HCO3],
                                "CO3": solution[CO3],
                                "Ca": solution[Ca],
                                "Mg": solution[Mg],
                                "CaMgCO32": solution[CaMgCO32],
                                "MgOH2": solution[MgOH2],
                                "CaCO3": 0,
                                "CaOH2": 0,
                                "Mg5CO34OH2": 0,
                                "Mg3CaCO34": 0,
                                "MgCO3": 0,
                            }
                            df.loc[len(df)] = Data 
                            n2_mg2co3oh2 = solution[Mg]**2 * solution[CO3] * (OH**2) - ksp_mg_2co3oh_2
                            n2_caoh2 = solution[Ca] * (OH**2) - ksp_caoh2
                            n2_mg5co34oh2 = solution[Mg]**5 * solution[CO3]**4 * OH**2 - ksp_mg5co34oh2
                            n2_mg3caco34 = solution[Mg]**3 * solution[Ca] * solution[CO3]**4 - ksp_mg3caco34
                            n2_caco3 = solution[Ca] * solution[CO3] - ksp_caco3
                            n2_mgco3 = solution[Mg] * solution[CO3] - ksp_mgco3
                            values = {
                                "n2_mg2co3oh2": n2_mg2co3oh2,
                                "n2_caoh2": n2_caoh2,
                                "n2_mg5co34oh2": n2_mg5co34oh2,
                                "n2_mg3caco34": n2_mg3caco34,
                                "n2_caco3": n2_caco3,
                                "n2_mgco3": n2_mgco3
                                    }
                            for name, value in values.items():
                                if value > 0:
                                    print(i, name, value)

                            if n2_caco3 >=0:
                                print(i,"CaMgCO32 +MgOH2 + CaCO3")
                                df = df.drop(df.index[-1]) 
                                CO2, HCO3, CO3, Ca, Mg, CaMgCO32, MgOH2, CaCO3= symbols('CO2 HCO3 CO3 Ca Mg CaMgCO32 MgOH2 CaCO3' )
                                eq1= Eq(k_co2 * CO2 , H* HCO3)
                                eq2 = Eq(k_co3 * HCO3, H* CO3)
                                eq3 = Eq(ksp_camgco3_2, Ca* Mg * CO3**2)
                                eq4= Eq(ksp_mgoh2, Mg* (OH**2))
                                eq5 = Eq(ksp_caco3, Ca* CO3)
                                eq6 = Eq(Ca_tot, Ca + CaCO3 + CaMgCO32)
                                eq7 = Eq(Mg_tot, Mg + MgOH2+ CaMgCO32)
                                eq8 = Eq(DIC, CO2 + HCO3 + CO3 + CaCO3 + 2*CaMgCO32)
                                sol = solve((eq1 , eq2, eq3, eq4, eq5, eq6, eq7,eq8), (CO2, HCO3, CO3, Ca, Mg, CaMgCO32, MgOH2,CaCO3), dict = True)
                                if sol:
                                    real_solutions = [{k: re(v) for k, v in solution.items()} for solution in sol]
                                    positive_solutions = [solution for solution in real_solutions if all(v > 0 for v in solution.values())]
                                    if positive_solutions :
                                        solution = positive_solutions[0]
                                        print(i, solution)
                                    else:
                                        print(i,"CaMgCO32 +MgOH2 + CaCO3 not have positive solution")

                                        # ##if MgOH2 is dissolved
                                        # print(i,"444")
                                        # CO2, HCO3, CO3, Ca, Mg, CaMgCO32, CaCO3= symbols('CO2 HCO3 CO3 Ca Mg CaMgCO32 CaCO3' )
                                        # eq1= Eq(k_co2 * CO2 , H* HCO3)
                                        # eq2 = Eq(k_co3 * HCO3, H* CO3)
                                        # eq3 = Eq(ksp_camgco3_2, Ca* Mg * CO3**2)
                                        # eq4 = Eq(ksp_caco3, Ca* CO3)
                                        # eq5 = Eq(Ca_tot, Ca + CaCO3 + CaMgCO32)
                                        # eq6 = Eq(Mg_tot, Mg + CaMgCO32)
                                        # eq7 = Eq(DIC, CO2 + HCO3 + CO3 + CaCO3 + 2*CaMgCO32)
                                        # sol = solve((eq1 , eq2, eq3, eq4, eq5, eq6, eq7), (CO2, HCO3, CO3, Ca, Mg, CaMgCO32, CaCO3), dict = True)
                                        # if sol:
                                        #     print(sol)


                                        #if CaMgCO32 is dissolved

                                        CO2, HCO3, CO3, Ca, Mg, MgOH2, CaCO3= symbols('CO2 HCO3 CO3 Ca Mg MgOH2 CaCO3' )
                                        eq1= Eq(k_co2 * CO2 , H* HCO3)
                                        eq2 = Eq(k_co3 * HCO3, H* CO3)
                                        eq3 = Eq(ksp_mgoh2, Mg* (OH**2))
                                        eq4 = Eq(ksp_caco3, Ca* CO3)
                                        eq5 = Eq(Ca_tot, Ca + CaCO3)
                                        eq6 = Eq(Mg_tot, Mg + MgOH2)
                                        eq7 = Eq(DIC, CO2 + HCO3 + CO3 + CaCO3)
                                        sol = solve((eq1 , eq2, eq3, eq4, eq5, eq6, eq7), (CO2, HCO3, CO3, Ca, Mg, MgOH2,CaCO3), dict = True)
                                        if sol:
                                            real_solutions = [{k: re(v) for k, v in solution.items()} for solution in sol]
                                            positive_solutions = [solution for solution in real_solutions if all(v > 0 for v in solution.values())]
                                            if positive_solutions :
                                                solution = positive_solutions[0]
                                                print(i, solution)
                                                Data = {
                                                    "pH": i,
                                                    "CO2": solution[CO2],
                                                    "HCO3": solution[HCO3],
                                                    "CO3": solution[CO3],
                                                    "Ca": solution[Ca],
                                                    "Mg": solution[Mg],
                                                    "CaMgCO32": 0,
                                                    "MgOH2": solution[MgOH2],
                                                    "CaCO3": solution[CaCO3],
                                                    "CaOH2": 0,
                                                    "Mg5CO34OH2": 0,
                                                    "Mg3CaCO34": 0,
                                                    "MgCO3": 0,
                                                }
                                                df.loc[len(df)] = Data 

                                                n3_mg2co3oh2 = solution[Mg]**2 * solution[CO3] * (OH**2) - ksp_mg_2co3oh_2
                                                n3_caoh2 = solution[Ca] * (OH**2) - ksp_caoh2
                                                n3_mg5co34oh2 = solution[Mg]**5 * solution[CO3]**4 * OH**2 - ksp_mg5co34oh2
                                                n3_mg3caco34 = solution[Mg]**3 * solution[Ca] * solution[CO3]**4 - ksp_mg3caco34
                                                n3_mgco3 = solution[Mg] * solution[CO3] - ksp_mgco3
                                                values = {
                                                    "n3_mg2co3oh2": n3_mg2co3oh2,
                                                    "n3_caoh2": n3_caoh2,
                                                    "n3_mg5co34oh2": n3_mg5co34oh2,
                                                    "n3_mg3caco34": n3_mg3caco34,
                                                    "n3_mgco3": n3_mgco3
                                                        }
                                                for name, value in values.items():
                                                    if value > 0:
                                                        print(i, name, value)

                                                if n3_caoh2 >=0:
                                                    print(i,"MgOH2 + CaCO3 + CaOH2")
                                                    df = df.drop(df.index[-1]) 
                                                    CO2, HCO3, CO3, Ca, Mg, MgOH2, CaCO3, CaOH2= symbols('CO2 HCO3 CO3 Ca Mg MgOH2 CaCO3 CaOH2' )
                                                    eq1= Eq(k_co2 * CO2 , H* HCO3)
                                                    eq2 = Eq(k_co3 * HCO3, H* CO3)
                                                    eq3 = Eq(ksp_mgoh2, Mg* (OH**2))
                                                    eq4 = Eq(ksp_caco3, Ca* CO3)
                                                    eq5 = Eq(ksp_caoh2, Ca * (OH**2))
                                                    eq6 = Eq(Ca_tot, Ca + CaCO3 + CaOH2)
                                                    eq7 = Eq(Mg_tot, Mg + MgOH2)
                                                    eq8 = Eq(DIC, CO2 + HCO3 + CO3 + CaCO3)
                                                    sol = solve((eq1 , eq2, eq3, eq4, eq5, eq6, eq7, eq8), (CO2, HCO3, CO3, Ca, Mg, MgOH2,CaCO3, CaOH2), dict = True)
                                                    if sol:
                                                        real_solutions = [{k: re(v) for k, v in solution.items()} for solution in sol]
                                                        positive_solutions = [solution for solution in real_solutions if all(v > 0 for v in solution.values())]
                                                    if positive_solutions :
                                                        solution = positive_solutions[0]
                                                        print(i, solution)

                                                        Data = {
                                                            "pH": i,
                                                            "CO2": solution[CO2],
                                                            "HCO3": solution[HCO3],
                                                            "CO3": solution[CO3],
                                                            "Ca": solution[Ca],
                                                            "Mg": solution[Mg],
                                                            "CaMgCO32": 0,
                                                            "MgOH2": solution[MgOH2],
                                                            "CaCO3": solution[CaCO3],
                                                            "CaOH2": solution[CaOH2],
                                                            "Mg5CO34OH2": 0,
                                                            "Mg3CaCO34": 0,
                                                            "MgCO3": 0,
                                                        }
                                                        df.loc[len(df)] = Data
                                                        n4_mg2co3oh2 = solution[Mg]**2 * solution[CO3] * (OH**2) - ksp_mg_2co3oh_2
                                                        n4_mg5co34oh2 = solution[Mg]**5 * solution[CO3]**4 * OH**2 - ksp_mg5co34oh2
                                                        n4_mg3caco34 = solution[Mg]**3 * solution[Ca] * solution[CO3]**4 - ksp_mg3caco34
                                                        n4_mgco3 = solution[Mg] * solution[CO3] - ksp_mgco3
                                                        values = {
                                                            "n4_mg2co3oh2": n4_mg2co3oh2,
                                                            "n4_mg5co34oh2": n4_mg5co34oh2,
                                                            "n4_mg3caco34": n4_mg3caco34,
                                                            "n4_mgco3": n4_mgco3
                                                                }
                                                        # for name, value in values.items():
                                                        #     if value > 0:
                                                        #         print(i, name, value)
                        else:
                            ##CaMgCO32 is dissolved
                            print(i,"CaMgCO32 +MgOH2 no positibe solution")
                            CO2, HCO3, CO3, Mg, MgOH2= symbols('CO2 HCO3 CO3 Mg MgOH2' )
                            eq1= Eq(k_co2 * CO2 , H* HCO3)
                            eq2 = Eq(k_co3 * HCO3, H* CO3)
                            eq3 = Eq(ksp_mgoh2, Mg* (OH**2))
                            eq4 = Eq(Mg_tot, Mg + MgOH2)
                            eq5 = Eq(DIC, CO2 + HCO3 + CO3)
                            sol = solve((eq1 , eq2, eq3, eq4, eq5), (CO2, HCO3, CO3, Mg, MgOH2), dict = True)
                            if sol:
                                real_solutions = [{k: re(v) for k, v in solution.items()} for solution in sol]
                                positive_solutions = [solution for solution in real_solutions if all(v > 0 for v in solution.values())]
                                if positive_solutions :
                                    solution = positive_solutions[0]
                                    print(i,solution)
                                Data = {
                                    "pH": i,
                                    "CO2": solution[CO2],
                                    "HCO3": solution[HCO3],
                                    "CO3": solution[CO3],
                                    "Ca": Ca_tot,
                                    "Mg": solution[Mg],
                                    "CaMgCO32":0,
                                    "MgOH2": solution[MgOH2],
                                    "CaCO3": 0,
                                    "CaOH2": 0,
                                    "Mg5CO34OH2": 0,
                                    "Mg3CaCO34": 0,
                                    "MgCO3": 0,
                                }
                                df.loc[len(df)] = Data
                                # print(i,solution)
                                n5_camgco3_2 = Ca_tot * solution[Mg] * (CO3) ** 2 - ksp_camgco3_2
                                n5_mg2co3oh2 = solution[Mg]**2 * solution[CO3] * (OH**2) - ksp_mg_2co3oh_2
                                n5_caoh2 = Ca_tot * (OH**2) - ksp_caoh2
                                n5_mg5co34oh2 = solution[Mg]**5 * solution[CO3]**4 * OH**2 - ksp_mg5co34oh2
                                n5_mg3caco34 = solution[Mg]**3 * Ca_tot * solution[CO3]**4 - ksp_mg3caco34
                                n5_caco3 = Ca_tot* solution[CO3] - ksp_caco3
                                n5_mgco3 = solution[Mg] * solution[CO3] - ksp_mgco3
                                # if n5_camgco3_2 >=0:
                                #     print(i,n5_camgco3_2)
                                # if n5_mg2co3oh2 >=0:
                                #     print(i,n5_mg2co3oh2)
                                # if n5_caoh2 >=0:
                                #     print(i,n5_caoh2)
                                # if n5_mg5co34oh2 >=0:
                                #     print(i,n5_mg5co34oh2)
                                # if n5_mg3caco34 >=0:
                                #     print(i,n5_mg3caco34)
                                # if n5_caco3 >=0:
                                #     print(i,n5_caco3)
                                # if n5_mgco3 >=0:
                                #     print(i,n5_mgco3)
                                ##result: CaCO3 + CaOH2
                                if n5_caco3 >=0 and n5_caoh2 >=0:
                                    print(i,"MgOH2+CaOH2 + CaCO3")
                                    df = df.drop(df.index[-1]) 
                                    CO2, HCO3, CO3, Ca, Mg, MgOH2, CaCO3, CaOH2= symbols('CO2 HCO3 CO3 Ca Mg MgOH2 CaCO3 CaOH2' )
                                    eq1= Eq(k_co2 * CO2 , H* HCO3)
                                    eq2 = Eq(k_co3 * HCO3, H* CO3)
                                    eq3 = Eq(ksp_mgoh2, Mg* (OH**2))
                                    eq4 = Eq(ksp_caco3, Ca* CO3)
                                    eq5 = Eq(ksp_caoh2, Ca * (OH**2))
                                    eq6 = Eq(Ca_tot, Ca + CaCO3 + CaOH2)
                                    eq7 = Eq(Mg_tot, Mg + MgOH2)
                                    eq8 = Eq(DIC, CO2 + HCO3 + CO3 + CaCO3)
                                    sol = solve((eq1 , eq2, eq3, eq4, eq5, eq6, eq7, eq8), (CO2, HCO3, CO3, Ca, Mg, MgOH2,CaCO3, CaOH2), dict = True)
                                    if sol:
                                        real_solutions = [{k: re(v) for k, v in solution.items()} for solution in sol]
                                        positive_solutions = [solution for solution in real_solutions if all(v > 0 for v in solution.values())]
                                        if positive_solutions :
                                            solution = positive_solutions[0]
                                            print(i, solution)
                                            Data = {
                                                            "pH": i,
                                                            "CO2": solution[CO2],
                                                            "HCO3": solution[HCO3],
                                                            "CO3": solution[CO3],
                                                            "Ca": solution[Ca],
                                                            "Mg": solution[Mg],
                                                            "CaMgCO32": 0,
                                                            "MgOH2": solution[MgOH2],
                                                            "CaCO3": solution[CaCO3],
                                                            "CaOH2": solution[CaOH2],
                                                            "Mg5CO34OH2": 0,
                                                            "Mg3CaCO34": 0,
                                                            "MgCO3": 0,
                                                        }
                                            df.loc[len(df)] = Data
                                            n6_mg2co3oh2 = solution[Mg]**2 * solution[CO3] * (OH**2) - ksp_mg_2co3oh_2
                                            n6_mg5co34oh2 = solution[Mg]**5 * solution[CO3]**4 * OH**2 - ksp_mg5co34oh2
                                            n6_mg3caco34 = solution[Mg]**3 * solution[Ca] * solution[CO3]**4 - ksp_mg3caco34
                                            n6_mgco3 = solution[Mg] * solution[CO3] - ksp_mgco3
                                            values = {
                                                            "n6_mg2co3oh2": n6_mg2co3oh2,
                                                            "n6_mg5co34oh2": n6_mg5co34oh2,
                                                            "n6_mg3caco34": n6_mg3caco34,
                                                            "n6_mgco3": n6_mgco3
                                                                }
                                            for name, value in values.items():
                                                if value > 0:
                                                    print(i, name, value)
            else:                
                Data = {
                    "pH": i,
                    "CO2": CO2_i,
                    "HCO3": HCO3_i,
                    "CO3": CO3_i,
                    "Ca": Ca_tot,
                    "Mg": Mg_tot,
                    "CaMgCO32": 0,
                    "MgOH2": 0,
                    "CaCO3": 0,
                    "CaOH2": 0,
                    "Mg5CO34OH2": 0,
                    "Mg3CaCO34": 0,
                    "MgCO3": 0,
                }
                df.loc[len(df)] = Data
    else:
        Data = {
                    "pH": i,
                    "CO2": CO2_i,
                    "HCO3": HCO3_i,
                    "CO3": CO3_i,
                    "Ca": Ca_tot,
                    "Mg": Mg_tot,
                    "CaMgCO32": 0,
                    "MgOH2": 0,
                    "CaCO3": 0,
                    "CaOH2": 0,
                    "Mg5CO34OH2": 0,
                    "Mg3CaCO34": 0,
                    "MgCO3": 0,
                }
        df.loc[len(df)] = Data







df.to_excel("base & acid_mol.xlsx", index=False)


# %%
