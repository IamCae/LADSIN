"""
Created on Fri Mar 28 00:17:54 2025

@author: cae
"""
import numpy as np
import pandas as pd

# Constantes
k = 0.4  # Constante de von Kármán
g = 9.81  # Gravidade (m/s²)
mi = 1200  # Coeficiente para rugosidade de Carmo et al. (2022)
zo_tab_maior = 0.0002  # Rugosidade tabelada (ajustável)
L = 100  # Comprimento de Monin-Obukhov (ajustável)

# Entradas manuais
Hs = float(input("Digite Hs (altura significativa da onda em metros): "))
Tp = float(input("Digite Tp (período de pico da onda em segundos): "))
vento10 = float(input("Digite a velocidade do vento a 10 m (m/s): "))

# Funções para cálculo
def beta(h, zo):
    return np.log(h / zo)

def ufric(h, uh, zo):
    return (k * uh) / beta(h, zo)

# Alturas para cálculo
niveis = np.arange(10, 130, 10)

# Criando DataFrame para armazenar resultados
df = pd.DataFrame({"Altura": niveis})

# METODOLOGIA 2 - Lei Logarítmica com Rugosidade Tabelada (DNV)
df["MET2"] = (ufric(10, vento10, zo_tab_maior) / k) * beta(df["Altura"], zo_tab_maior)

# METODOLOGIA 5 - Carmo et al. (2022) - Perfil Neutro
Lp = (Tp**2 * g) / (2 * np.pi)
zo_car22 = mi * Hs * (Hs / Lp) ** 4.5
ufric_car22 = (k * vento10) / np.log(10 / zo_car22)
df["MET5"] = (ufric_car22 / k) * np.log(df["Altura"] / zo_car22)

# METODOLOGIA 6 - Carmo et al. (2022) - Perfil Estável
df["MET6"] = (ufric_car22 / k) * (np.log(df["Altura"] / zo_car22) + (4.7 * df["Altura"]) / L)

# Salvando CSV
df.to_csv("/media/ladsin/Seagate Expansion Drive/SIEME/Perfis/Recife/PERFIL4_Recife_boia.csv", index=False)
print("Resultados salvos")
