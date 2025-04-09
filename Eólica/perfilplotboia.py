"""
Created on Fri Mar 28 00:56:53 2025

@author: cae
"""

import pandas as pd
import matplotlib.pyplot as plt

# Carregar o arquivo CSV
file_path = "/media/ladsin/Seagate Expansion Drive/SIEME/Perfis/Recife/PERFIL4_Recife_boia.csv"  # Substitua pelo caminho correto do arquivo
df = pd.read_csv(file_path)

# Definir as alturas de 10 em 10 até 120
alturas = list(range(10, 130, 10))  # Lista de alturas de 10 a 120

# Filtrar o DataFrame para garantir que contenha apenas essas alturas
df = df[df["Altura"].isin(alturas)]

# Criar o gráfico
plt.figure(figsize=(6, 8))
plt.plot(df["MET2"], df["Altura"], label="DNV", marker="o", linestyle="-")
plt.plot(df["MET5"], df["Altura"], label="Carmo 22 (Neutro)", marker="o", linestyle="-")
plt.plot(df["MET6"], df["Altura"], label="Carmo 22 (Estável)", marker="o", linestyle="-")

# Configurações do gráfico
plt.xlabel("Velocidade do Vento (m/s)")
plt.ylabel("Altura (m)")
plt.title("Perfil do Vento (boia) Recife 29/12/2013 09:00:00")
plt.yticks(alturas)  # Define os intervalos de 10 em 10 no eixo Y
plt.legend()
plt.grid(True)

# Mostrar o gráfico
plt.show()
