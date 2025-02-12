#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""
Created on Wed Jan 22 11:53:49 2025

@author: Caetano Pereira
""

import numpy as np
import pandas as pd
import xarray as xr
from tqdm import tqdm
import matplotlib.pyplot as plt
#from matplotlib.windrose import WindroseAxes
from windrose import WindroseAxes

# Diretorios (caminhos) principais
pathname='/home/ladsin/'
pathname2 = '/media/operacao/Seagate Expansion Drive/ERA5VENTO(NOVO)'
#pathname3=glob.glob('/home/ladsin/Caetano/boias/*.csv')
pathname4= pathname2+'/Caetano/TSM_MUR.2/data'
pathname5='/media/ladsin/Seagate Expansion Drive/Caetano/Energia/INTERPOLADOS'

# Constantes
k = 0.4  # Constante de Von Karman
g = 9.81  # Aceleração gravitacional (m/s²)
alfa_oc = 0.14  # Coeficiente da Lei da Potência (ajustado de 0.12 para 0.14)
zo_tab_maior = 0.001  # Rugosidade para Lei Logarítmica 1
epsilon= 1200 ##const de ajuste da estabilidade pela T-TSM
Cp = 1.67  # Constante (P) Termodinâmica
mi = 1200  # Coeficiente para Carmo et al. 2022
L = 90  # Comprimento de Monin-Obukhov (mudar conforme a estabilidade, valor escolhido por causa do artigo de Schalldemose )
niveis=[10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170]

# Abrindo dados .nc

ds = xr.open_mfdataset(pathname2+'/ERA5_WIND_INTERPOLADO.nc')
ds1 = xr.open_mfdataset(pathname2+'/ERA5_WAVE_INTERPOLADO.nc')
ds2 = xr.open_mfdataset(pathname2+'/ERA5_T_INTERPOLADO_CELSIUS.nc')
ds3 = xr.open_mfdataset(pathname5+'/WAVERYS_INTERPOLADO_6H.nc')
ds4 = xr.open_mfdataset(pathname5+'/OISST_INTERPOLADO_6H.nc')
ds5 = xr.open_mfdataset(pathname5+'/TSM_MUR_INTERPOLADO_6H.nc')

# Transformando lat e lon em matriz
lons, lats = np.meshgrid(ds.lon, ds.lat)
lons1, lats1 = np.meshgrid(ds1.lon, ds1.lat)
lons2, lats2 = np.meshgrid(ds2.lon, ds2.lat)
lons3, lats3 = np.meshgrid(ds3.lon, ds3.lat)
lons4, lats4 = np.meshgrid(ds4.lon, ds4.lat)
lons5, lats5 = np.meshgrid(ds5.lon, ds5.lat)

# Coordenadas e intervalo de tempo
def filtrar_dados(ds, lon, lat, start_date, end_date):
   # Seleciona o ponto espacial mais próximo
    ds_filtrado = ds.sel(lon=lon, lat=lat, method='nearest')
    # Aplica filtro de temporal
    ds_filtrado = ds_filtrado.sel(time=slice(start_date, end_date))
    return ds_filtrado

lon, lat = -40.04083, -22.43278
start_date, end_date = '2003-01-01 00:00', '2018-12-31 18:00'

ERA5_wind = filtrar_dados(ds, lon, lat, start_date, end_date)
ERA5_wave = filtrar_dados(ds1, lon, lat, start_date, end_date)
ERA5_T = filtrar_dados(ds2, lon, lat, start_date, end_date)
WAVERYS = filtrar_dados(ds3, lon, lat, start_date, end_date)
OISST = filtrar_dados(ds4, lon, lat, start_date, end_date)
MUR = filtrar_dados(ds5, lon, lat, start_date, end_date)

# Variáveis de vento, onda e temperatura (ERA5)
ERA5_wind['mag10_ERA5'] = ((ERA5_wind['u10'] ** 2) + (ERA5_wind['v10'] ** 2)) ** 0.5
ERA5_wind['dir10_ERA5'] = np.degrees(np.arctan2(-ERA5_wind['u10'], -ERA5_wind['v10']))
ERA5_wind.loc[ERA5_wind['dir10_ERA5'] < 0, 'dir10_ERA5'] += 360
ERA5_wind['Hs_ERA5'] = ERA5_wave['swh']
ERA5_wind['Tp_ERA5'] = ERA5_wave['pp1d']
ERA5_wind['Permed_ERA5'] = ERA5_wave['mwp']
ERA5_wind['T2m_ERA5'] = ERA5_T['t2m']
ERA5_wind['TSM_ERA5'] = ERA5_T['sst']

# Variáveis de onda (WAVERYS)
WAVERYS['Permed_WAV'] = WAVERYS['VTM01_WW'] 
WAVERYS['Tp_WAV'] = WAVERYS['VTPK']
WAVERYS['Hs_WAV'] = WAVERYS['VHM0']

# Variáveis de TSM (NOAA)
OISST['TSM_OISST'] = OISST['sst']

# Variáveis de TSM (NASA)
MUR['TSM_MUR'] = MUR['analysed_sst']

# Validação inicial de dados
def validar_dados(df):
    """Verifica se há valores nulos e retorna informações básicas."""
    if df.isnull().sum().any():
        print("Aviso: Existem valores nulos no DataFrame. Verifique antes de continuar.")
    print(df.info())
    return df

# Funções auxiliares
def beta(h, zo):
    """Calcula o fator beta para a velocidade de fricção."""
    return np.log(h / zo)

def ufric(h, uh, zo):
    """Calcula a velocidade de fricção."""
    return (k * uh) / beta(h, zo)

# Metodologias
def metodologia1(df, niveis):
    """Lei da Potência"""
    for i in niveis:
        df[f'mag{i}_ERA5_MET1'] = df['mag10_ERA5'] * ((i / 10) ** alfa_oc)
    return df

def metodologia2(df, niveis):
    """DNV (Rugosidade com valor tabelado)"""
    for i in niveis:
        df[f'mag{i}_ERA5_MET2'] = (ufric(10, df['mag10_ERA5'], zo_tab_maior) / k) * beta(i, zo_tab_maior)
    return df

def metodologia3(df, niveis):
    """Perfil Neutro (PSI=0) + Rugosidade De Donelan 90"""
    df['zo_Do_90'] = 0.033 * (df['Hs_ERA5'] / 4)
    for i in niveis:
        df[f'mag{i}_ERA5_MET3'] = (ufric(10, df['mag10_ERA5'], df['zo_Do_90']) / k) * beta(i, df['zo_Do_90'])
    return df

def metodologia4(df, niveis):
    """Perfil Neutro (PSI=0) + Rugosidade De Donelan et al. 93"""
    df['zo_Do_93'] = (df['Hs_ERA5'] / 4) * (6.7e-4) * ((df['mag10_ERA5'] / Cp) ** 2.6)
    for i in niveis:
        df[f'mag{i}_ERA5_MET4'] = (ufric(10, df['mag10_ERA5'], df['zo_Do_93']) / k) * beta(i, df['zo_Do_93'])
    return df

def metodologia5(df, niveis):
    """Perfil Neutro (PSI=0) + Carmo et al. 2022 Rugosidade de Taylor & Yelland Corrigida"""
    df['Lp'] = (df['Tp_ERA5'] ** 2 * g) / (2 * np.pi)
    df['zo_car22'] = mi * df['Hs_ERA5'] * ((df['Hs_ERA5'] / df['Lp']) ** 4.5)
    df['ufric_car22'] = (k * df['mag10_ERA5']) / np.log(10 / df['zo_car22'])
    for i in niveis:
        df[f'mag{i}_ERA5_MET5'] = (df['ufric_car22'] / k) * np.log(i / df['zo_car22'])
    return df

def metodologia6(df, niveis):
    """Perfil Estável (PSI=4.7Z/L) + Carmo et al. 2022 Rugosidade de Taylor & Yelland Corrigida"""
    for i in niveis:
        df[f'mag{i}_ERA5_MET6'] = (df['ufric_car22'] / k) * (np.log(i / df['zo_car22']) + (4.7 * i / L))
    return df

# Processamento principal
def processar_metodologias(df, niveis):
    """Aplica todas as metodologias e retorna o DataFrame processado."""
    with tqdm(total=6, desc="Aplicando metodologias") as pbar:
        df = metodologia1(df, niveis)
        pbar.update(1)
        df = metodologia2(df, niveis)
        pbar.update(1)
        df = metodologia3(df, niveis)
        pbar.update(1)
        df = metodologia4(df, niveis)
        pbar.update(1)
        df = metodologia5(df, niveis)
        pbar.update(1)
        df = metodologia6(df, niveis)
        pbar.update(1)
    return df

# Geração de relatório
def gerar_relatorio(df, niveis, metodologias, output_path):
    """Gera um relatório com o resumo dos resultados."""
    resumo = {}
    for met in metodologias:
        resumo[met] = [df[f'mag{nivel}_ERA5_{met}'].mean() for nivel in niveis]
    relatorio = pd.DataFrame(resumo, index=niveis)
    relatorio.to_csv(output_path)
    print(f"Relatório gerado em {output_path}")

# Exemplo de uso
pathname2 = '/caminho/para/salvar/dados'
ERA5_wind = pd.read_csv('/caminho/para/arquivo/ERA5.csv')  # Ajustar caminho
niveis = [20, 50, 100]  # Exemplos de níveis (ajustar conforme necessário)

# Processar
ERA5_wind = processar_metodologias(ERA5_wind, niveis)

# Salvar resultados
ERA5_wind.to_csv(f'{pathname2}/ERA5_Celsiusv2.csv', index=False)


# CRIANDO DF VAZIO PARA CALCULAR OS PERFIS
df_perfis = pd.DataFrame(np.random.randn(17, 1), index='10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170'.split(),
                         columns=['z', 'Met1_ERA5', 'Met2_ERA5', 'Met3_ERA5', 'Met4_ERA5', 'Met5_ERA5', 'Met6_ERA5'])
dados = ['ERA5']

# Loop para preencher os valores das colunas 'Met1_ERA5' até 'Met6_ERA5'
for i in range(len(niveis)):
    index_value = df_perfis.index[i]  # Acessando o índice apenas uma vez
    df_perfis.loc[i, 'Met1_ERA5'] = ERA5_wind[f'mag{index_value}_ERA5_MET1'].mean()
    df_perfis.loc[i, 'Met2_ERA5'] = ERA5_wind[f'mag{index_value}_ERA5_MET2'].mean()
    df_perfis.loc[i, 'Met3_ERA5'] = ERA5_wind[f'mag{index_value}_ERA5_MET3'].mean()
    df_perfis.loc[i, 'Met4_ERA5'] = ERA5_wind[f'mag{index_value}_ERA5_MET4'].mean()
    df_perfis.loc[i, 'Met5_ERA5'] = ERA5_wind[f'mag{index_value}_ERA5_MET5'].mean()
    df_perfis.loc[i, 'Met6_ERA5'] = ERA5_wind[f'mag{index_value}_ERA5_MET6'].mean()

# Função para plotar gráficos de perfis de vento
def plot_perfil_vento(df, ax, title, xlabel, ylabel, filename):
    ax.set_xlim(6, 16)
    ax.legend()
    # Adicionando as diferentes metodologias de vento no gráfico
    for met in range(1, 7):  # De Met1 a Met6
        ax.plot(df[f'Met{met}_ERA5'], df.index, label=f'Met{met}', linestyle='-', color=plt.cm.Set1(met/7))  # Usando diferentes cores
    for k in range(len(df.index)):
        ax.quiver(0, df.index[k], df.iloc[k, 0], 0, units="dots", width=2.5, headwidth=2, headlength=4, 
                  headaxislength=5, color='black', scale=0.0225)
    ax.set_title(title, fontsize=20)
    ax.set_ylabel(ylabel, fontsize=20)
    ax.set_xlabel(xlabel, fontsize=20)
    ax.grid(True)
    ax.legend()
    plt.savefig(filename)
    plt.clf()

# Função para plotar a rosa dos ventos
def plot_rosa_dos_ventos(df_wind, ax, title, filename):
    ax.bar(df_wind['dir10_ERA5'], df_wind['mag10_ERA5'], normed=True, opening=0.8, edgecolor='black', 
           cmap=plt.cm.winter, lw=1, bins=np.arange(0.01, 16, 2))
    ax.set_title(title, fontsize=12)
    ax.set_legend(loc=2, decimal_places=1)
    plt.savefig(filename)

# Função para plotar os gráficos de potência eólica
def plot_potencial_eolico(df, ax, title, xlabel, ylabel, filename):
    ax.set_xlim(0, 3)
    ax.legend()
    ax.set_title(title, fontsize=20)
    ax.set_ylabel(ylabel, fontsize=20)
    ax.set_xlabel(xlabel, fontsize=20)
    plt.savefig(filename)
    plt.clf()

# Gráfico de Perfis de Vento
fig, ax = plt.subplots(figsize=(8, 15))
plot_perfil_vento(df_perfis, ax, 'PERFIL DE VENTO CEARÁ - 2003 a 2018 (ERA5)', 'Velocidade (m/s)', 'Altura (m)', 
                  '/media/ladsin/Seagate Expansion Drive/Caetano/Energia/Plots/perfil_vento2_CE_cpam.png')

# Gráfico Rosa dos Ventos
ax = WindroseAxes.from_ax()
plot_rosa_dos_ventos(ERA5_wind, ax, 'Rosa dos ventos para RS - 2003 a 2018 (ERA5)', 
                     '/media/ladsin/Seagate Expansion Drive/Caetano/Energia/Plots/P1_rosa_dos_ventos_RS_Cpam.png')

# Cálculo do Potencial Eólico
Drotor_3000kW = 82
cp_3000kW = 0.45
dens_ar = 1.22  # kgm-³
A = 5281  # área do rotor

# Calculando a potência
for met in ['Met2_ERA5', 'Met5_ERA5', 'Met6_ERA5']:
    df_perfis[met] = (0.5 * dens_ar * A * (df_perfis[met] ** 3) * cp_3000kW) / 1000000

# Gráfico de Potencial Eólico
fig, ax = plt.subplots(figsize=(8, 15))
plot_potencial_eolico(df_perfis, ax, 'Potencial Eólico Rio Grande do Sul - 2003 a 2018', 'Potência (MW)', 'Altura (m)', 
                      '/media/ladsin/Seagate Expansion Drive/Caetano/Energia/Plots/P1_Pot_Eol_RS.png')

# Preenchendo as colunas de Met2_ERA5, Met5_ERA5, Met6_ERA5 com base nos dados
for local in dados:
    for met in ['Met2_ERA5', 'Met5_ERA5', 'Met6_ERA5']:
        df_perfis[f'{met}_{local}'] = np.nan

# Atualizando as colunas de acordo com as médias dos dados
for i in range(len(df_perfis.index)):
    for met in ['Met2_ERA5', 'Met5_ERA5', 'Met6_ERA5']:
        df_perfis.loc[df_perfis.index[i], f'Met{met}_ERA5'] = ERA5_wind[f'mag{i+1}_ERA5_{met}'].mean()
