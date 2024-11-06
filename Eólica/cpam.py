#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 16:52:04 2024

@author: ladsin
"""

# ------------------------------------------------------------------------------------------------------------
import os
import glob
import pandas as pd
import xarray as xr
import numpy as np
from windrose import WindroseAxes
import colormap as cm
import matplotlib.pyplot as plt
randn = np.random.randn
from PIL import Image
from numpy import asarray, max, arange, round, insert, radians
from numpy import ceil, ma, cumsum, array, argmin
from numpy import linspace, meshgrid, histogram2d, flipud, size, sum
from numpy import nanmax, nanmean, nansum
import xlsxwriter
from bokeh.plotting import figure, output_file, show
#from scipy.stats import pearsonr
#import skill_metrics as sm
#import seaborn as sns
import cartopy.crs as ccrs
import cartopy.feature as cfeature
#from cartopy.io.shapereader import Reader
#from cartopy.feature import ShapelyFeature
#from cartopy.mpl.gridliner import lon_FORMATTER, lat_FORMATTER
#import cartopy.io.shapereader as shpreader # Import shapefiles

# ------------------------------------------------------------------------------------------------------------
#                        CPAM 2024
#usar DNV 0.14power law 
#comparar ERA5, WAVERYS, TSM-OISST, TSM-MUR e PNBOIA (interpolados temporalmente de 6 em 6 horas)
#dado observado, ponto da plataforma (?)
#fazer perfil de vento
#calcular potecial eólico
#aplicar estatísticas comparativas
#comparar SODAR (omega) superficie a 300m <------- pedir dado
# ------------------------------------------------------------------------------------------------------------
#                                                  OBJETIVO
#PERFIL 1 -> DNV vs CARMO 22 (NEUTRO)                | ERA5 EM AMBOS
#PERFIL 2 -> CARMO 22 (NEUTRO) vs CARMO 22 (ESTÁVEL) | ERA5 EM AMBOS 
#PERFIL 3 -> CARMO 22 (NEUTRO) vs CARMO 22 (NEUTRO)  | ERA5 vs WAVERYS



# ------------------------------------------------------------------------------------------------------------

#diretorios (caminhos) principais
pathname='/home/ladsin/'
pathname2='/media/ladsin/Seagate Expansion Drive/'
#pathname3=glob.glob('/home/ladsin/Caetano/boias/*.csv')
pathname4= pathname2+'/Caetano/TSM_MUR.2/data'
pathname5='/media/ladsin/Seagate Expansion Drive/Caetano/Energia/INTERPOLADOS'

# ------------------------------------------------------------------------------------------------------------

#                                          CONSTANTES
dens_ag=1025 #kgm-³
dens_ar=1.22 #kgm-³
alfa_oc=0.14#Termo de ajuste da lei da potencia pelo metodo da DNV ----ajustado de 0.12 para 0.14-----
alfa_Powerlaw=0.14 #alfa do power law 1/7
beta=0.7
zo_tab_menor=0.0001  #tabelaDNV - open sea without waves
zo_tab_maior=0.001  #tabelaDNV - Open sea with waves
epsilon=1200 ##const de ajuste da estabilidade pela T-TSM
E=0.033
mi=1200 #constante de ajuste de CARMO et al. 22 - ajustado do metodo de Taylor & YElland 2001
a=2.7
b=0.142
c=13.09
k=0.4 ##constsnte de von-karman
g=9.80665 #gravidade
Cp=1.67 #constante (P) - termodinamica
fc=0.5
L=90 #AQUI EH O COMPRIMENTO DE MONIN OBUKHOV. AQUI VOCE MUDAR O VALOR CONFORME A ESTABILIDADE. FOI ESCOLHIDO ESSE VALOR
#PARA IR DE ENCONTRO COM O ARTIGO DE SCHALDEMOSE
zo_tab_menor=0.0001  #tabelaDNV - open sea without waves
zo_tab_maior=0.001  #tabelaDNV - Open sea with waves

niveis=[10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170]

#                                 ABRINDO DADOS .nc

## ERA5 (1993-2018) 
ds = xr.open_dataset(pathname2+'/ERA5VENTO(NOVO)/ERA5_WIND_INTERPOLADO.nc') #abrindo os arquivos nc

ds1 = xr.open_dataset(pathname2+'/ERA5VENTO(NOVO)/ERA5_WAVE_INTERPOLADO.nc') #abrindo os arquivos nc

## WAVERYS -> ECWMF(1993 - 2018)
#ds3 = xr.open_mfdataset(pathname5+'/WAVERYS_INTERPOLADO_6H.nc')

# TESTANDO NC ALTERADO
ds3= xr.open_dataset('/media/ladsin/Seagate Expansion Drive/Caetano/Energia/waverys6h_v2.nc')
# ^     FUNCIONOU!     ^

# utilizar xr.open_mfdataset para abrir multiplos nc's

# ------------------------------------------------------------------------------------------------------------
#                                            TRATANDO ÁREAS DE INTERESSE (COMENTAR E DESCOMENTAR LINHAS PARA TIRAR PERFIL DE ESTADO)

#Transformando lat e lon em matriz
lons,lats = np.meshgrid(ds.lon,ds.lat)
lons1,lats1 = np.meshgrid(ds1.lon,ds1.lat)
lons3,lats3 = np.meshgrid(ds3.lon,ds3.lat)

#Fortaleza - Ceará (BOIA - PNBOIA)
#ERA5_wind=ds.sel(lon=-38.4325,lat=-3.2136, method='nearest').to_dataframe()
#ERA5_wave=ds1.sel(lon=-38.4325,lat=-3.2136, method='nearest').to_dataframe()
#WAVERYS=ds3.sel(lon=-38.4325,lat=-3.2136, method='nearest').to_dataframe()


#Santos - São Paulo (BOIA - PNBOIA)
# ERA5_wind = ds.sel(lon=-45.0361,lat=-25.433944, method='nearest').to_dataframe()
# ERA5_wave = ds1.sel(lon=-45.0361,lat=-25.433944, method='nearest').to_dataframe()
# WAVERYS=ds3.sel(lon=-49.51,lat=-31.32, method='nearest').to_dataframe()

#Rio Grande do Sul (BOIA - PNBOIA (boia Rio Grande))
#ERA5_wind = ds.sel(lon=-49.51,lat=-31.32, method='nearest').to_dataframe()
#ERA5_wave = ds1.sel(lon=-49.51,lat=-31.32, method='nearest').to_dataframe()
#WAVERYS   = ds3.sel(lon=-49.51,lat=-31.32, method='nearest').to_dataframe()

#latlon -> maranhão <- alto potencial eolico
ERA5_wind=ds.sel(lon=-42.63,lat=-2.31,method='nearest').to_dataframe()#mudei a lon de 40 para 45
ERA5_wave=ds1.sel(lon=-42.63,lat=-2.31,method='nearest').to_dataframe()
#ERA5_T=ds2.sel(lon=-42.63,lat=-2.31,method='nearest').to_dataframe()
WAVERYS=ds3.sel(lon=-42.63,lat=-2.31,method='nearest').to_dataframe()


# ------------------------------------------------------------------------------------------------------------
#                                    VARIÁVEIS GERAIS

#separa em intervalo de tempo desejado
ERA5_wind=ERA5_wind[(ERA5_wind.index>='2018-01-01 06:00') & (ERA5_wind.index<='2018-12-31 18:00')] 
ERA5_wave=ERA5_wave[(ERA5_wave.index>='2018-01-01 06:00') & (ERA5_wave.index<='2018-12-31 18:00')]
WAVERYS=WAVERYS[(WAVERYS.index>='2018-01-01 06:00') & (WAVERYS.index<='2018-12-31 18:00')]

# variaveis de vento, onda e temperatura (ERA5)
ERA5_wind['mag10_ERA5']=np.nan
ERA5_wind['mag10_ERA5']=((ERA5_wind['u10']**2)+(ERA5_wind['v10']**2))**0.5#1.7
ERA5_wind['dir10_ERA5']=np.nan
ERA5_wind['dir10_ERA5']=np.degrees(np.arctan2(-ERA5_wind['u10'],-ERA5_wind['v10']))
#ERA5_wind['dir10_ERA5'][ERA5_wind['dir10_ERA5']<0]=ERA5_wind['dir10_ERA5']+360
ERA5_wind.loc[ERA5_wind['dir10_ERA5'] < 0, 'dir10_ERA5'] = ERA5_wind['dir10_ERA5'] + 360
ERA5_wind['Hs_ERA5']=np.nan
ERA5_wind['Hs_ERA5']=ERA5_wave['swh'] # altura significativa de onda
ERA5_wind['Tp_ERA5']=np.nan
ERA5_wind['Tp_ERA5']=ERA5_wave['pp1d'] #periodo de pico da onda
ERA5_wind['Permed_ERA5']=ERA5_wave['mwp'] #periodo medio de vento

# variaveis de onda (WAVERYS)
WAVERYS['Permed_WAV']=np.nan
WAVERYS['Permed_WAV']= WAVERYS['VTM01_WW'] #periodo medio das ondas de vento
WAVERYS['Tp_WAV']=np.nan
WAVERYS['Tp_WAV']= WAVERYS['VTPK'] #periodo de pico da onda
#WAVERYS['Hs_WAV']=np.nan
#WAVERYS['Hs_WAV']= WAVERYS['VHM0'] # altura significativa de onda

#_________________________________________________________________________________________________________________________________

#                                    METODOLOGIAS

#METODOLOGIA 2 - LEI LOGARITMA 1 (DNV (RUGOSIDADE COM VALOR TABELADO))
#METODOLOGIA 5 - LEI LOGARITMA 4 (PERFIL NEUTRO (PSI=0) + CARMO ET AL. 2022 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA))
#METODOLOGIA 6 - LEI LOGARITMA 5 (PERFIL ESTÁVEL (PSI=4.7Z/L) + CARMO ET AL. 2022 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA))

##Funcoes para calcular o para beta para poder calcular a velocidade de friccao 
def beta(h,zo):
    B=np.log(h/zo)
    return B

##nessa funcao calcula-se a velocidade de friccao. eh necessario fornecer os seguintes parametros:
##h - altura do valor do vento "observado"
##uh - eh o proprio vento
##zo - eh o metodo de rugosidade que voce vai usar
def ufric(h,uh,zo):
    ufric=(k*uh)/beta(h,zo)
    return ufric

#METODOLOGIA 1 - LEI DA POTENCIA    
for i in niveis:
    ERA5_wind['mag10_ERA5_MET1']=ERA5_wind['mag10_ERA5']
    ERA5_wind['mag{}_ERA5_MET1'.format(i)]=ERA5_wind['mag10_ERA5']*((i/10)**(alfa_oc))  
    ERA5_wind['mag10_ERA5_MET1']


#METODOLOGIA 2 - LEI LOGARITMA 1 (DNV (RUGOSIDADE COM VALOR TABELADO))    
for i in niveis:
    ERA5_wind['mag10_ERA5_MET2']=ERA5_wind['mag10_ERA5']
    ERA5_wind['mag{}_ERA5_MET2'.format(i)]=(ufric(10,ERA5_wind['mag10_ERA5'],zo_tab_maior)/k)*(beta(i,zo_tab_maior))


#METODOLOGIA 3 - LEI LOGARITMA 2 (PERFIL NEUTRO (PSI=0) + RUGOSIDADE DE DONELAN 90)    
#CALCULO DA RUGOSIDADE PELO Hs - Donelan 90
ERA5_wind['zo_Do_90']=0.033*(ERA5_wind['Hs_ERA5']/4)

for i in niveis:
     ERA5_wind['mag10_ERA5_MET3']=ERA5_wind['mag10_ERA5']
     ERA5_wind['mag{}_ERA5_MET3'.format(i)]=(ufric(10,ERA5_wind['mag10_ERA5'],ERA5_wind['zo_Do_90'])/k)*(beta(i,ERA5_wind['zo_Do_90']))


#METODOLOGIA 4 - LEI LOGARITMA 3 (PERFIL NEUTRO (PSI=0) + DONELAN ET AL. 93)    
#CALCULO DA RUGOSIDADE PELO Hs - DOnelan 93
ERA5_wind['zo_Do_93']=(ERA5_wind['Hs_ERA5']/4)*(6.7*(10**-4))*((ERA5_wind['mag10_ERA5']/Cp)**2.6)

for i in niveis:
    ERA5_wind['mag10_ERA5_MET4']=ERA5_wind['mag10_ERA5']
    ERA5_wind['mag{}_ERA5_MET4'.format(i)]=(ufric(10,ERA5_wind['mag10_ERA5'],ERA5_wind['zo_Do_93'])/k)*(beta(i,ERA5_wind['zo_Do_93']))

#METODOLOGIA 5 - LEI LOGARITMA 4 (PERFIL NEUTRO (PSI=0) + CARMO ET AL. 2022 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA))
#calculando Zo pelo metodo de CARMO 22 - adaptado de taylor e yelland
ERA5_wind['Lp']=(ERA5_wind['Tp_ERA5']**2*g)/(2*np.pi)  
ERA5_wind['ufric_car22']=(k*ERA5_wind['mag10_ERA5'])/np.log(10/(mi*ERA5_wind['Hs_ERA5']*((ERA5_wind['Hs_ERA5']/ERA5_wind['Lp'])**4.5)))
ERA5_wind['zo_car22']=(mi*ERA5_wind['Hs_ERA5']*((ERA5_wind['Hs_ERA5']/ERA5_wind['Lp'])**4.5))

for i in niveis:
    ERA5_wind['mag10_ERA5_MET5']=ERA5_wind['mag10_ERA5']
    ERA5_wind['mag{}_ERA5_MET5'.format(i)]=(ERA5_wind['ufric_car22']/k)*((np.log((i)/ERA5_wind['zo_car22'])))


#METODOLOGIA 6 - LEI LOGARITMA 5 (PERFIL ESTÁVEL (PSI=4.7Z/L) + CARMO ET AL. 2022 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA))
for i in niveis:
    ERA5_wind['mag10_ERA5_MET6']=ERA5_wind['mag10_ERA5']
    ERA5_wind['mag{}_ERA5_MET6'.format(i)]=(ERA5_wind['ufric_car22']/k)*((np.log((i)/ERA5_wind['zo_car22']))+((4.7*(i))/L))
     
ERA5_wind.to_csv(pathname2+'/Caetano/Energia/PROCESSADOS/SAZ/PERFIL_ERA5_MA.csv')    


# ------------------------------------------------------------------------------------------------------------
#                                CRIANDO DF VAZIO PARA CALCULAR OS PERFIS

df_perfis = pd.DataFrame(randn(17,1),index='10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170'.split(),
                columns='z'.split())
dados =['ERA5'] 



for local in dados:
#     df_perfis['Met1_{}'.format(local)]=np.nan  #zo_tab1 = 0.001
     df_perfis['Met2_{}'.format(local)]=np.nan   #Lei da potencia alfa
#     df_perfis['Met3_{}'.format(local)]=np.nan   #zo_Don1990
#     df_perfis['Met4_{}'.format(local)]=np.nan   #zo_Don1993
     df_perfis['Met5_{}'.format(local)]=np.nan   #zo_Charnock1955 - solução Oost 2002
     df_perfis['Met6_{}'.format(local)]=np.nan   #zo_He2019
#     df_perfis['Met7_{}'.format(local)]=np.nan   #zo_Don90 - Est
    

i=0
while i < len(niveis):
#     df_perfis['Met1_ERA5'][i]=ERA5_wind['mag{}_ERA5_MET1'.format(df_perfis.index[i])].mean()
#     df_perfis['Met1_ERAW'][i]=dfERAW_wind['mag{}_dfERAW_MET1'.format(df_perfis.index[i])].mean()
#     df_perfis['Met1_ERAWO'][i]=dfERAWO_wind['mag{}_dfERAWO_MET1'.format(df_perfis.index[i])].mean()
#     df_perfis['Met1_ERAWM'][i]=dfERAWM_wind['mag{}_dfERAWM_MET1'.format(df_perfis.index[i])].mean()
    
     df_perfis['Met2_ERA5'][i]=ERA5_wind['mag{}_ERA5_MET2'.format(df_perfis.index[i])].mean()
#     df_perfis['Met2_ERAW'][i]=dfERAW_wind['mag{}_dfERAW_MET2'.format(df_perfis.index[i])].mean()
#     df_perfis['Met2_ERAWO'][i]=dfERAWO_wind['mag{}_dfERAWO_MET2'.format(df_perfis.index[i])].mean()
#     df_perfis['Met2_ERAWM'][i]=dfERAWM_wind['mag{}_dfERAWM_MET2'.format(df_perfis.index[i])].mean()
    
#     df_perfis['Met3_ERA5'][i]=ERA5_wind['mag{}_ERA5_MET3'.format(df_perfis.index[i])].mean()
#     df_perfis['Met3_ERAW'][i]=dfERAW_wind['mag{}_dfERAW_MET3'.format(df_perfis.index[i])].mean()
#     df_perfis['Met3_ERAWO'][i]=dfERAWO_wind['mag{}_dfERAWO_MET3'.format(df_perfis.index[i])].mean()
#     df_perfis['Met3_ERAWM'][i]=dfERAWM_wind['mag{}_dfERAWM_MET3'.format(df_perfis.index[i])].mean()
    
#     df_perfis['Met4_ERA5'][i]=ERA5_wind['mag{}_ERA5_MET4'.format(df_perfis.index[i])].mean()
#     df_perfis['Met4_ERAW'][i]=dfERAW_wind['mag{}_dfERAW_MET4'.format(df_perfis.index[i])].mean()
#     df_perfis['Met4_ERAWO'][i]=dfERAWO_wind['mag{}_dfERAWO_MET4'.format(df_perfis.index[i])].mean()
#     df_perfis['Met4_ERAWM'][i]=dfERAWM_wind['mag{}_dfERAWM_MET4'.format(df_perfis.index[i])].mean()
    
     df_perfis['Met5_ERA5'][i]=ERA5_wind['mag{}_ERA5_MET5'.format(df_perfis.index[i])].mean()
#     df_perfis['Met5_ERAW'][i]=dfERAW_wind['mag{}_dfERAW_MET5'.format(df_perfis.index[i])].mean()
#     df_perfis['Met5_ERAWO'][i]=dfERAWO_wind['mag{}_dfERAWO_MET5'.format(df_perfis.index[i])].mean()
#     df_perfis['Met5_ERAWM'][i]=dfERAWM_wind['mag{}_dfERAWM_MET5'.format(df_perfis.index[i])].mean()
    
     df_perfis['Met6_ERA5'][i]=ERA5_wind['mag{}_ERA5_MET6'.format(df_perfis.index[i])].mean()
#     df_perfis['Met6_ERAW'][i]=dfERAW_wind['mag{}_dfERAW_MET6'.format(df_perfis.index[i])].mean()
#     df_perfis['Met6_ERAWO'][i]=dfERAWO_wind['mag{}_dfERAWO_MET6'.format(df_perfis.index[i])].mean()
#     df_perfis['Met6_ERAWM'][i]=dfERAWM_wind['mag{}_dfERAWM_MET6'.format(df_perfis.index[i])].mean()
    
     i=i+1

# ------------------------------------------------------------------------------------------------------------

#                                                 GRAFICOS DIREÇÃO 

#    ERA5
ax = WindroseAxes.from_ax()
ax.bar(ERA5_wind['dir10_ERA5'], ERA5_wind['mag10_ERA5'], normed=True, opening=0.8, edgecolor='black', cmap=plt.cm.winter, lw=1, bins=np.arange(0.01,16,2))
ax.set_title('Rosa dos ventos para RS  - 2003 to 2018 (ERA5)', fontsize=12)
ax.set_legend(loc=2,decimal_places=1)
#logo = plt.imread('ladsin_2.png')
#ax.figure.figimage(logo, xo = 60, yo = 70)
plt.savefig('/media/ladsin/Seagate Expansion Drive/Caetano/Energia/Plots/P1_rosa_dos_ventos_RS_Cpam.png')


# #   ERA5 + WAVERYS 
# ax = WindroseAxes.from_ax()
# ax.bar(dfERAW_wind['dir10_dfERAW'], dfERAW_wind['mag10_dfERAW'], normed=True, opening=0.8, edgecolor='black', cmap=plt.cm.winter, lw=1, bins=np.arange(0.01,16,2))
# ax.set_title('Rosa dos ventos para P18 - 1997 a 2018 (ERA5 + WAVERYS)', fontsize=12)
# ax.set_legend(loc=2,decimal_places=1)
# #logo = plt.imread('ladsin_2.png')
# #ax.figure.figimage(logo, xo = 60, yo = 70)
# plt.savefig('/media/ladsin/Seagate Expansion Drive/Caetano/energia/Plots/rosa_dos_ventos_ERAW.png', fontsize=15)

# #   ERA5 + WAVERYS + OISST
# ax = WindroseAxes.from_ax()
# ax.bar(dfERAWO_wind['dir10_dfERAWO'], dfERAWO_wind['mag10_dfERAWO'], normed=True, opening=0.8, edgecolor='black', cmap=plt.cm.winter, lw=1, bins=np.arange(0.01,16,2))
# ax.set_title('Rosa dos ventos para P18 - 1997 a 2018 (ERA5 + WAVERYS + OISST)', fontsize=12)
# ax.set_legend(loc=2,decimal_places=1)
# #logo = plt.imread('ladsin_2.png')
# #ax.figure.figimage(logo, xo = 60, yo = 70)
# plt.savefig('/media/ladsin/Seagate Expansion Drive/Caetano/energia/Plots/rosa_dos_ventos_ERAWO.png', fontsize=15)

# #   ERA5 + WAVERYS + MUR
# ax = WindroseAxes.from_ax()
# ax.bar(dfERAW_wind['dir10_dfERAWM'], dfERAW_wind['mag10_dfERAWM'], normed=True, opening=0.8, edgecolor='black', cmap=plt.cm.winter, lw=1, bins=np.arange(0.01,16,2))
# ax.set_title('Rosa dos ventos para P18 - 1997 a 2018 (ERA5 + WAVERYS + MUR)', fontsize=12)
# ax.set_legend(loc=2,decimal_places=1)
# #logo = plt.imread('ladsin_2.png')
# #ax.figure.figimage(logo, xo = 60, yo = 70)
# plt.savefig('/media/ladsin/Seagate Expansion Drive/Caetano/energia/Plots/rosa_dos_ventos_ERAWM.png', fontsize=15)


# # ------------------------------------------------------------------------------------------------------------

# #                                                PERFIS DE VENTO

# #      PERFIS DE VENTO - ERA5
fig,ax=plt.subplots(figsize=(8,15))
# line1= ax.plot(df_perfis['Met1_ERA5'],df_perfis.index,color='red',label='1/7 PowerLaw', linestyle='dotted')
#line2= ax.plot(df_perfis['Met2_ERA5'],df_perfis.index,color='orange',label='DNV', linestyle='-.')
# line3= ax.plot(df_perfis['Met3_ERA5'],df_perfis.index,color='green',label='Zo Donelan 90', linestyle=':')
# line4= ax.plot(df_perfis['Met4_ERA5'],df_perfis.index, color='purple',label='Zo Donelan 93', linestyle='-')
line5= ax.plot(df_perfis['Met5_ERA5'],df_perfis.index, color='gray',label='Carmo 22 - Neu', linestyle='dashed')
line6= ax.plot(df_perfis['Met6_ERA5'],df_perfis.index, color='blue',label='Carmo 22 - Est', linestyle='-')
# line7= ax.plot(df_perfis['Met7_ERA5'],df_perfis.index, color='yellow',label='Carmo 22 - Inst', linestyle='-') #estourando muito no plot
# #line7= ax.plot(df_perfis['Met9_ERA5'],df_perfis.index, color='red',label='Zo TayYel01 - Stable', linestyle='-.')
# #line8=plt.scatter(df_perfis['SODAR'],df_perfis.index, marker="*", color='black', label='SODAR')

ax.set_xlim(6,16)
ax.legend()
#criando o for para fazer o vetor - ERA5
for k in range (0,np.size(df_perfis.index)):
    qax = plt.quiver(0,df_perfis.index[k],df_perfis['Met6_ERA5'][k],0, units="dots",
              width=2.5,
              headwidth = 2,
              headlength = 4,
              headaxislength = 5,
              color='black',
              scale=0.0225)
ax.set_title('PERFIL DE VENTO CEARÁ - 2003 a 2018 (ERA5)', fontsize=20)
ax.set_ylabel('Altura (m)', fontsize=20)
ax.set_xlabel('Velocidade (m/s)', fontsize=20)
plt.savefig('/media/ladsin/Seagate Expansion Drive/Caetano/Energia/Plots/perfil_vento2_CE_cpam.png') #coloca 2 digitos
plt.clf() #fecha o arquivo para nao plotar por cima

# ------------------------------------------------------------------------------------------------------------

#                            Calculo do PotenciaL Eólico ERA5
Drotor_3000kW=82
cp_3000kW=0.45 
dens_ar=1.22 #kgm-³
A=5281 #ou (np.pi*(Drotor_3000kW**2))/4
#df_perfis['Pot_MET1_ERA5']=(0.5*dens_ar*A*(df_perfis['Met1_ERA5']**3)*cp_3000kW)/1000000
df_perfis['Pot_MET2_ERA5']=(0.5*dens_ar*A*(df_perfis['Met2_ERA5']**3)*cp_3000kW)/1000000
#df_perfis['Pot_MET3_ERA5']=(0.5*dens_ar*A*(df_perfis['Met3_ERA5']**3)*cp_3000kW)/1000000
#df_perfis['Pot_MET4_ERA5']=(0.5*dens_ar*A*(df_perfis['Met4_ERA5']**3)*cp_3000kW)/1000000
df_perfis['Pot_MET5_ERA5']=(0.5*dens_ar*A*(df_perfis['Met5_ERA5']**3)*cp_3000kW)/1000000
df_perfis['Pot_MET6_ERA5']=(0.5*dens_ar*A*(df_perfis['Met6_ERA5']**3)*cp_3000kW)/1000000

#                            PERFIS DE VENTO DO POTENCIAL EÓLICO - ERA5
fig,ax=plt.subplots(figsize=(8,15))
#line1= ax.plot(df_perfis['Pot_MET1_ERA5'],df_perfis.index,color='red',label='1/7 PowerLaw', linestyle='dotted')
line2= ax.plot(df_perfis['Pot_MET2_ERA5'],df_perfis.index,color='orange',label='DNV', linestyle='-.')
#line3= ax.plot(df_perfis['Pot_MET3_ERA5'],df_perfis.index,color='green',label='Zo Donelan 90', linestyle=':')
#line4= ax.plot(df_perfis['Pot_MET4_ERA5'],df_perfis.index, color='purple',label='Zo Donelan 93', linestyle='-')
line5= ax.plot(df_perfis['Pot_MET5_ERA5'],df_perfis.index, color='gray',label='Carmo 22 - Neu', linestyle='dashed')
line6= ax.plot(df_perfis['Pot_MET6_ERA5'],df_perfis.index, color='blue',label='Carmo 22 - Est', linestyle='-')
#line7= ax.plot(df_perfis['Met9_ERA5'],df_perfis.index, color='red',label='Zo TayYel01 - Stable', linestyle='-.')
#line8=plt.scatter(df_perfis['SODAR'],df_perfis.index, marker="*", color='black', label='SODAR')


ax.set_xlim(0,3)
ax.legend()
#criando o for para fazer o vetor - ERA5
ax.set_title('Potencial Eólico Rio Grande do Sul - 2003 a 2018 ', fontsize=20)
ax.set_ylabel('Altura (m)', fontsize=20)
ax.set_xlabel('Potência (MW)', fontsize=20)
plt.savefig('/media/ladsin/Seagate Expansion Drive/Caetano/Energia/Plots/P1_Pot_Eol_RS.png') #coloca 2 digitos
plt.clf() #fecha o arquivo para nao plotar por cima

  