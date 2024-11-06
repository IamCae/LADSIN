#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 11:53:12 2023

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
import json
from bokeh.plotting import figure, output_file, show
from scipy.stats import pearsonr
#import skill_metrics as sm
import seaborn as sns
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
#from cartopy.mpl.gridliner import lon_FORMATTER, lat_FORMATTER
import cartopy.io.shapereader as shpreader # Import shapefiles

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
#niveis2=[5,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160]
#metodologias=['MET1.1_EST','MET1.1_NEU','MET1.2_EST','MET1.2_NEU',
#         'MET2','MET2','MET3_EST','MET3_NEU','MET4_EST',
#         'MET4_NEU','MET5_EST','MET5_NEU','MET3_NEU_P18',
#         'MET3_EST_P18','MET4_NEU_P18','MET4_EST_P18',
#         'MET5_NEU_P18','MET5_EST_P18','MET6_NEU_P18',
#         'MET6A_NEU_P18']

# -------------------------------------------------- ----------------------------------------------------------

#                                 ABRINDO DADOS DE BOIA

#pnb1 = '/home/ladsin/Caetano/boias/ondas_fortaleza_01.csv'
#pnb2 = '/home/ladsin/Caetano/boias/ondas_rio_grande_01.csv'
#pnb3 = '/home/ladsin/Caetano/boias/ondas_santos_01.csv'

#Fortaleza = pd.read_csv(pnb1, sep=';', decimal='.')
#RioGrande = pd.read_csv(pnb2, sep=';', decimal='.')
#Santos = pd.read_csv(pnb3, sep=';', decimal='.')


#                                 ABRINDO DADOS .nc

## ERA5 (1993-2018) 
ds = xr.open_mfdataset(pathname2+'/ERA5VENTO(NOVO)/ERA5_WIND_INTERPOLADO.nc') #abrindo os arquivos nc

ds1 = xr.open_mfdataset(pathname2+'/ERA5VENTO(NOVO)/ERA5_WAVE_INTERPOLADO.nc') #abrindo os arquivos nc

ds2 = xr.open_mfdataset(pathname2+'/ERA5VENTO(NOVO)/ERA5_T_INTERPOLADO_CELSIUS.nc') #abrindo os arquivos nc, convertido para celcius

## WAVERYS -> ECWMF(1993 - 2018)
ds3 = xr.open_mfdataset(pathname5+'/WAVERYS_INTERPOLADO_6H.nc')
#ds2 = xr.open_mfdataset(pathname2+'/ERA5VENTO(NOVO)/Waverys_interp_2013.nc')
# TESTANDO NC ALTERADO
ds3= xr.open_mfdataset('/media/ladsin/Seagate Expansion Drive/hs_waverys/HS_interp_2003.nc')
# ^     FUNCIONOU!     ^

## OISST -> NOAA (1993 - 2023)
ds4 = xr.open_mfdataset(pathname5+'/OISST_INTERPOLADO_6H.nc')

## TSM MUR -> NASA (09/2002 - 08/2023)
ds5= xr.open_mfdataset(pathname5+'/TSM_MUR_INTERPOLADO_6H.nc')

# ------------------------------------------------------------------------------------------------------------


#Transformando lat e lon em matriz
lons,lats = np.meshgrid(ds.lon,ds.lat)
lons1,lats1 = np.meshgrid(ds1.lon,ds1.lat)
lons2,lats2 = np.meshgrid(ds2.lon,ds2.lat)
lons3,lats3 = np.meshgrid(ds3.lon,ds3.lat)
lons4,lats4 = np.meshgrid(ds4.lon,ds4.lat)
lons5,lats5 = np.meshgrid(ds5.lon,ds5.lat)



##PONTO PRA-1 - ou outro ponto a definir (boia, plataforma ou LIDAR)
#ERA5=ds.sel(lat=-22.134251,lon=-38.412515,method='nearest').to_dataframe()
#ERA5_wave=ds1.sel(lat=-22.134251,lon=-38.412515,method='nearest').to_dataframe()
#ERA5_T=ds2.sel(lat=-22.134251,lon=-38.412515,method='nearest').to_dataframe()
#### descobrir qual latlon acima

#latlon -> PlataformaP18 - Bacia de Campos
ERA5_wind=ds.sel(lon=-40.04083,lat=-22.43278,method='nearest').to_dataframe()
ERA5_wave=ds1.sel(lon=-40.04083,lat=-22.43278,method='nearest').to_dataframe()
ERA5_T=ds2.sel(lon=-40.04083,lat=-22.43278,method='nearest').to_dataframe()
WAVERYS=ds3.sel(lon=-40.04083,lat=-22.43278,method='nearest').to_dataframe()
OISST=ds4.sel(lon=-40.04083,lat=-22.43278,method='nearest').to_dataframe()
MUR=ds5.sel(lon=-40.04083,lat=-22.43278,method='nearest').to_dataframe()



#Fortaleza - Ceará (BOIA - PNBOIA)
#dfERAW_wind_=ds.sel(lon=-38.4325,lat=-3.2136, method='nearest').to_dataframe()
#dfERAW_wave=ds1.sel(lon=-38.4325,lat=-3.2136, method='nearest').to_dataframe()
#dfERAW_T=ds2.sel(lon=-38.4325,lat=-3.2136, method='nearest').to_dataframe()

#Santos - São Paulo (BOIA - PNBOIA)
#dfERAWO_wind = ds.sel(lon=-45.0361,lat=-25.433944, method='nearest').to_dataframe()
#dfERAWO_wave = ds.sel(lon=-45.0361,lat=-25.433944, method='nearest').to_dataframe()
#dfERAWO_T=ds2.sel(lon=-45.0361,lat=-25.433944, method='nearest').to_dataframe()

#Rio Grande do Sul (BOIA - PNBOIA (boia Rio Grande))
#dfERAWM_wind = ds.sel(lon=-49.51,lat=-31.32, method='nearest').to_dataframe()
#dfERAWM_wave = ds.sel(lon=-49.51,lat=-31.32, method='nearest').to_dataframe()
#dfERAWM_T=ds2.sel(lon=-49.51,lat=-31.32, method='nearest').to_dataframe()

#ERA5_wind['T2m_ERA5']=np.nan
#ERA5_wind['T2m_ERA5']=ERA5_T['t2m'] # temperatura do ar a 2 metros
#ERA5_wind['TSM_ERA5']=np.nan
#ERA5_wind['TSM_ERA5']=ERA5_T['sst'] #TSM
#ERA5=ERA5.drop(columns=['lon'])
#ERA5=ERA5.drop(columns=['lat'])

#salvando o arquivo "ERA5" em csv, para nao precisar fazer sempre esse procedimento
#ERA5_wind.to_csv(pathname2+'/Caetano/energia/PROCESSADOS/proc1.csv')
################
#ABRINDO O CSV SALVO
#ERA5 = pd.read_csv(pathname2+'/Caetano/energia/PROCESSADOS/proc1.csv')

# ------------------------------------------------------------------------------------------------------------
#                                    VARIÁVEIS GERAIS


#df_RS = ds_wind.sel(lon=-49.834522,lat=-31.551188, method='nearest').to_dataframe()
#df_RS_wave = ds_wave.sel(lon=-49.834522,lat=-31.551188, method='nearest').to_dataframe()

#separa em intervalo de tempo desejado
ERA5_wind=ERA5_wind[(ERA5_wind.index>='2003-01-01 00:00') & (ERA5_wind.index<='2018-12-31 18:00')] 
ERA5_wave=ERA5_wave[(ERA5_wave.index>='2003-01-01 00:00') & (ERA5_wave.index<='2018-12-31 18:00')]
ERA5_T=ERA5_T[(ERA5_T.index>='2003-01-01 00:00') & (ERA5_T.index<='2018-12-31 18:00')]
WAVERYS=WAVERYS[(WAVERYS.index>='2003-01-01 00:00') & (WAVERYS.index<='2018-12-31 18:00')]
OISST=OISST[(OISST.index>='2003-01-01 00:00') & (OISST.index<='2018-12-31 18:00')]
MUR=MUR[(MUR.index>='2003-01-01 00:00') & (MUR.index<='2018-12-31 18:00')]




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
ERA5_wind['T2m_ERA5']=np.nan
ERA5_wind['T2m_ERA5']=ERA5_T['t2m'] # temperatura do ar a 2 metros 
ERA5_wind['TSM_ERA5']=np.nan
ERA5_wind['TSM_ERA5']=ERA5_T['sst'] #TSM

# variaveis de onda (WAVERYS)
WAVERYS['Permed_WAV']=np.nan
WAVERYS['Permed_WAV']= WAVERYS['VTM01_WW'] #periodo medio da onda
WAVERYS['Tp_WAV']=np.nan
WAVERYS['Tp_WAV']= WAVERYS['VTPK'] #periodo de pico da onda
WAVERYS['Hs_WAV']=np.nan
WAVERYS['Hs_WAV']= WAVERYS['VHM0'] # altura significativa da onda

# variaveis de tsm (NOAA)
OISST['TSM_OISST']=np.nan
OISST['TSM_OISST']=OISST['sst'] #TSM

# variaveis de tsm (NASA)
MUR['TSM_MUR']=np.nan
MUR['TSM_MUR']=MUR['analysed_sst'] #TSM

#_________________________________________________________________________________________________________________________________

#                                     ERA5

#METODOLOGIA 1 - LEI DA POTENCIA
#METODOLOGIA 2 - LEI LOGARITMA 1 (DNV (RUGOSIDADE COM VALOR TABELADO))
#METODOLOGIA 3 - LEI LOGARITMA 2 (PERFIL NEUTRO (PSI=0) + RUGOSIDADE DE DONELAN 90)
#METODOLOGIA 4 - LEI LOGARITMA 3 (PERFIL NEUTRO (PSI=0) + DONELAN ET AL. 93)
#METODOLOGIA 5 - LEI LOGARITMA 4 (PERFIL NEUTRO (PSI=0) + CARMO ET AL. 2022 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA))
#METODOLOGIA 6 - LEI LOGARITMA 5 (PERFIL ESTÁVEL (PSI=4.7Z/L) + CARMO ET AL. 2022 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA))
#METODOLOGIA 7 - LEI LOGARITMA 6 (CARMO ET AL. 2023 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA + CALCULO DE PSI PELA T e TSM))

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
    ERA5_wind['mag10_ER A5_MET1']


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
     
   # nao usae met7 no cpam

#            SOLUCIONAR O PROBLEMA DO CALCULO DA TSM
#METODOLOGIA 7 - LEI LOGARITMA 6 (CARMO ET AL. 2023 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA + CALCULO DE PSI PELA T e TSM))
#ERA5_wind['EST_CORRIGIDO']=((np.abs(ERA5_wind['TSM_ERA5']*ERA5_wind['T2m_ERA5']*0.7))/epsilon)  ## TERMO DE AJUSTE DA ESTABILIDADE PELA TSM E T
#ERA5_wind['zo_car23']=(mi*ERA5_wind['Hs_ERA5']*((ERA5_wind['Hs_ERA5']/ERA5_wind['Lp'])**4.5)) #DE REPENTE, PODE FAZER UMA CORRECAO NO ZO

#for i in niveis:
 #   ERA5_wind['mag10_ERA5_MET7']=ERA5_wind['mag10_ERA5']
  #  ERA5_wind['mag{}_ERA5_MET7'.format(i)]=(ERA5_wind['ufric_car22']/k)*((np.log((i)/ERA5_wind['zo_car23']))+(((i)/ERA5_wind['EST_CORRIGIDO'])))

ERA5_wind.to_csv(pathname2+'/Caetano/Energia/PROCESSADOS/CPAM/ERA5_Celsiusv2.csv') #V2 pela por alterar latlon, comparar

# ------------------------------------------------------------------------------------------------------------

#                                ERA5+WAVERYS 

dfERAW_wind=ds.sel(lon=-40.04083,lat=-22.43278, method='nearest').to_dataframe()     #ERA5
dfERAW_wave=ds3.sel(lon=-40.04083,lat=-22.43278, method='nearest').to_dataframe()    #WAVEYRS
dfERAW_T=ds2.sel(lon=-40.04083,lat=-22.43278, method='nearest').to_dataframe()       #ERA5

#separa em intervalo de tempo desejado
dfERAW_wind=dfERAW_wind[(dfERAW_wind.index>='2003-06-01 06:00') & (dfERAW_wind.index<='2018-12-31 18:00')]  #TESTE TEMPO
dfERAW_wave=dfERAW_wave[(dfERAW_wave.index>='2003-06-01 06:00') & (dfERAW_wave.index<='2018-12-31 18:00')]
dfERAW_T=dfERAW_T[(dfERAW_T.index>='2003-06-01 06:00') & (dfERAW_T.index<='2018-12-31 18:00')] 

# variaveis de vento, onda e temperatura
dfERAW_wind['mag10_dfERAW']=np.nan
dfERAW_wind['mag10_dfERAW']=((ERA5_wind['u10']**2)+(ERA5_wind['v10']**2))**0.5#1.7
dfERAW_wind['dir10_dfERAW']=np.nan
dfERAW_wind['dir10_dfERAW']=np.degrees(np.arctan2(-ERA5_wind['u10'],-ERA5_wind['v10']))
#dfERAW_wind['dir10_dfERAW'][dfERAW_wind['dir10_dfERAW']<0]=dfERAW_wind['dir10_dfERAW']+360
dfERAW_wind.loc[dfERAW_wind['dir10_dfERAW'] < 0, 'dir10_dfERAW'] =ERA5_wind['dir10_ERA5'] + 360
dfERAW_wind['Hs_dfERAW']=np.nan
dfERAW_wind['Hs_dfERAW']=ERA5_wave['swh']                       #WAVERYS['VTPK'] # altura significativa de onda ||| ver se estou chamando certo, se não, testar =WAVERYS['VTPK']
dfERAW_wind['Tp_dfERAW']=np.nan
dfERAW_wind['Tp_dfERAW']=WAVERYS['VTM01_WW'] #periodo de pico da onda
dfERAW_wind['Permed_dfERAW']=np.nan
dfERAW_wind['Permed_dfERAW']=WAVERYS['VMDR'] #periodo medio de vento
dfERAW_wind['T2m_ERA5']=np.nan
dfERAW_wind['T2m_ERA5']=ERA5_T['t2m'] # temperatura do ar a 2 metros 
dfERAW_wind['TSM_ERA5']=np.nan
dfERAW_wind['TSM_ERA5']=ERA5_T['sst'] #TSM


#METODOLOGIA 1 - LEI DA POTENCIA
#METODOLOGIA 2 - LEI LOGARITMA 1 (DNV (RUGOSIDADE COM VALOR TABELADO))
#METODOLOGIA 3 - LEI LOGARITMA 2 (PERFIL NEUTRO (PSI=0) + RUGOSIDADE DE DONELAN 90)
#METODOLOGIA 4 - LEI LOGARITMA 3 (PERFIL NEUTRO (PSI=0) + DONELAN ET AL. 93)
#METODOLOGIA 5 - LEI LOGARITMA 4 (PERFIL NEUTRO (PSI=0) + CARMO ET AL. 2022 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA))
#METODOLOGIA 6 - LEI LOGARITMA 5 (PERFIL ESTÁVEL (PSI=4.7Z/L) + CARMO ET AL. 2022 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA))
#METODOLOGIA 7 - LEI LOGARITMA 6 (CARMO ET AL. 2023 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA + CALCULO DE PSI PELA T e TSM))

# ##Funcoes para calcular o para beta para poder calcular a velocidade de friccao 
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
     dfERAW_wind['mag10_dfERAW_MET1']=dfERAW_wind['mag10_dfERAW']
     dfERAW_wind['mag{}_dfERAW_MET1'.format(i)]=dfERAW_wind['mag10_dfERAW']*((i/10)**(alfa_oc))


#METODOLOGIA 2 - LEI LOGARITMA 1 (DNV (RUGOSIDADE COM VALOR TABELADO))    
for i in niveis:
     dfERAW_wind['mag10_dfERAW_MET2']=dfERAW_wind['mag10_dfERAW']
     dfERAW_wind['mag{}_dfERAW_MET2'.format(i)]=(ufric(10,dfERAW_wind['mag10_dfERAW'],zo_tab_maior)/k)*(beta(i,zo_tab_maior))


#METODOLOGIA 3 - LEI LOGARITMA 2 (PERFIL NEUTRO (PSI=0) + RUGOSIDADE DE DONELAN 90)    
#CALCULO DA RUGOSIDADE PELO Hs - Donelan 90
dfERAW_wind['zo_Do_90']=0.033*(dfERAW_wind['Hs_dfERAW']/4)

for i in niveis:
     dfERAW_wind['mag10_dfERAW_MET3']=dfERAW_wind['mag10_dfERAW']
     dfERAW_wind['mag{}_dfERAW_MET3'.format(i)]=(ufric(10,dfERAW_wind['mag10_dfERAW'],dfERAW_wind['zo_Do_90'])/k)*(beta(i,dfERAW_wind['zo_Do_90']))

#METODOLOGIA 4 - LEI LOGARITMA 3 (PERFIL NEUTRO (PSI=0) + DONELAN ET AL. 93)    
#CALCULO DA RUGOSIDADE PELO Hs - DOnelan 93
dfERAW_wind['zo_Do_93']=(dfERAW_wind['Hs_dfERAW']/4)*(6.7*(10**-4))*((dfERAW_wind['mag10_dfERAW']/Cp)**2.6)

for i in niveis:
     dfERAW_wind['mag10_dfERAW_MET4']=dfERAW_wind['mag10_dfERAW']
     dfERAW_wind['mag{}_dfERAW_MET4'.format(i)]=(ufric(10,dfERAW_wind['mag10_dfERAW'],dfERAW_wind['zo_Do_93'])/k)*(beta(i,dfERAW_wind['zo_Do_93']))


#METODOLOGIA 5 - LEI LOGARITMA 4 (PERFIL NEUTRO (PSI=0) + CARMO ET AL. 2022 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA))
#calculando Zo pelo metodo de CARMO 22 - adaptado de taylor e yelland
dfERAW_wind['Lp']=(dfERAW_wind['Tp_dfERAW']**2*g)/(2*np.pi)   #descobrir oq é 'lp' e se tem q mudar algo nas variáveis
dfERAW_wind['ufric_car22']=(k*dfERAW_wind['mag10_dfERAW'])/np.log(10/(mi*dfERAW_wind['Hs_dfERAW']*((dfERAW_wind['Hs_dfERAW']/dfERAW_wind['Lp'])**4.5))) 
dfERAW_wind['zo_car22']=(mi*dfERAW_wind['Hs_dfERAW']*((dfERAW_wind['Hs_dfERAW']/dfERAW_wind['Lp'])**4.5))
                                                      # Evitar divisão por zero na zo_car22
#dfERAW_wind['zo_car22'] = np.where(dfERAW_wind['zo_car22'] > 0, dfERAW_wind['zo_car22'], np.nan)


for i in niveis:
     dfERAW_wind['mag10_dfERAW_MET5']=dfERAW_wind['mag10_dfERAW']
     dfERAW_wind['mag{}_dfERAW_MET5'.format(i)]=(dfERAW_wind['ufric_car22']/k)*((np.log((i)/dfERAW_wind['zo_car22'])))


#METODOLOGIA 6 - LEI LOGARITMA 5 (PERFIL ESTÁVEL (PSI=4.7Z/L) + CARMO ET AL. 2022 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA))
for i in niveis:
     dfERAW_wind['mag10_dfERAW_MET6']=dfERAW_wind['mag10_dfERAW']
     dfERAW_wind['mag{}_dfERAW_MET6'.format(i)]=(dfERAW_wind['ufric_car22']/k)*((np.log((i)/dfERAW_wind['zo_car22']))+((4.7*(i))/L))

#            SOLUCIONAR O PROBLEMA DO CALCULO DA TSM
#METODOLOGIA 7 - LEI LOGARITMA 6 (CARMO ET AL. 2023 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA + CALCULO DE PSI PELA T e TSM))
#dfERAW_wind['EST_CORRIGIDO']=((np.abs(dfERAW_wind['TSM_ERA5']*dfERAW_wind['T2m_ERA5']*0.7))/epsilon) ## TERMO DE AJUSTE DA ESTABILIDADE PELA TSM E T
#dfERAW_wind['zo_car23']=(mi*dfERAW_wind['Hs_dfERAW']*((dfERAW_wind['Hs_dfERAW']/dfERAW_wind['Lp'])**4.5)) #DE REPENTE, PODE FAZER UMA CORRECAO NO ZO

#for i in niveis:
 #    dfERAW_wind['mag10_dfERAW_MET7']=dfERAW_wind['mag10_dfERAW']
  #   dfERAW_wind['mag{}_dfERAW_MET7'.format(i)]=(dfERAW_wind['ufric_car22']/k)*((np.log((i)/dfERAW_wind['zo_car23']))+(((i)/dfERAW_wind['EST_CORRIGIDO'])))

# Salvando o DataFrame em CSV
dfERAW_wind.to_csv(pathname2+'/Caetano/Energia/PROCESSADOS/CPAM/ERAW_hsERA5.csv')#CONSERTAR DEPOIS

# ------------------------------------------------------------------------------------------------------------

#                            ERA5+WAVERYS+OISST 

dfERAWO_wind = ds.sel(lon=-40.04083,lat=-22.43278, method='nearest').to_dataframe()    #ERA5
dfERAWO_wave = ds3.sel(lon=-40.04083,lat=--22.43278, method='nearest').to_dataframe()  #WAVERYS
dfERAWO_T=ds2.sel(lon=-40.04083,lat=-22.43278, method='nearest').to_dataframe()        #ERA5
dfERAWO_TSM=ds4.sel(lon=-40.04083,lat=-22.43278, method='nearest').to_dataframe()      #OISST  

#separa em intervalo de tempo desejado
dfERAWO_wind=dfERAWO_wind[(dfERAWO_wind.index>='2003-01-01 00:00') & (dfERAWO_wind.index<='2018-12-31 18:00')] 
dfERAWO_wave=dfERAWO_wave[(dfERAWO_wave.index>='2003-01-01 00:00') & (dfERAWO_wave.index<='2018-12-31 18:00')]
dfERAWO_T=dfERAWO_T[(dfERAWO_T.index>='2003-01-01 00:00') & (dfERAWO_T.index<='2018-12-31 18:00')] 
dfERAWO_TSM=dfERAWO_TSM[(dfERAWO_TSM.index>='2003-01-01 00:00') & (dfERAWO_TSM.index<='2018-12-31 18:00')] 


# variaveis de vento, onda e temperatura
dfERAWO_wind['mag10_dfERAWO']=np.nan
dfERAWO_wind['mag10_dfERAWO']=((ERA5_wind['u10']**2)+(ERA5_wind['v10']**2))**0.5#1.7
dfERAWO_wind['dir10_dfERAWO']=np.nan
dfERAWO_wind['dir10_dfERAWO']=np.degrees(np.arctan2(-ERA5_wind['u10'],-ERA5_wind['v10']))
#dfERAWO_wind['dir10_dfERAWO'][dfERAWO_wind['dir10_dfERAWO']<0]=dfERAWO_wind['dir10_dfERAWO']+360
dfERAWO_wind.loc[dfERAWO_wind['dir10_dfERAWO'] < 0, 'dir10_dfERAWO'] = ERA5_wind['dir10_ERA5'] + 360
dfERAWO_wind['T2m_ERA5']=np.nan
dfERAWO_wind['T2m_ERA5']=ERA5_T['t2m'] # temperatura do ar a 2 metros 
dfERAWO_wind['TSM_OISST']=np.nan
dfERAWO_wind['TSM_OISST']=OISST['sst'] #TSM
dfERAWO_wind['Hs_dfERAWO']=np.nan
dfERAWO_wind['Hs_dfERAWO']=WAVERYS['VTPK'] # altura significativa de onda
dfERAWO_wind['Tp_dfERAWO']=np.nan
dfERAWO_wind['Tp_dfERAWO']=WAVERYS['VTM01_WW'] #periodo de pico da onda
dfERAWO_wind['Permed_dfERAWO']=np.nan
dfERAWO_wind['Permed_dfERAWO']=WAVERYS['VMDR'] #periodo medio de vento


#METODOLOGIA 1 - LEI DA POTENCIA
#METODOLOGIA 2 - LEI LOGARITMA 1 (DNV (RUGOSIDADE COM VALOR TABELADO))
#METODOLOGIA 3 - LEI LOGARITMA 2 (PERFIL NEUTRO (PSI=0) + RUGOSIDADE DE DONELAN 90)
#METODOLOGIA 4 - LEI LOGARITMA 3 (PERFIL NEUTRO (PSI=0) + DONELAN ET AL. 93)
#METODOLOGIA 5 - LEI LOGARITMA 4 (PERFIL NEUTRO (PSI=0) + CARMO ET AL. 2022 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA))
#METODOLOGIA 6 - LEI LOGARITMA 5 (PERFIL ESTÁVEL (PSI=4.7Z/L) + CARMO ET AL. 2022 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA))
#METODOLOGIA 7 - LEI LOGARITMA 6 (CARMO ET AL. 2023 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA + CALCULO DE PSI PELA T e TSM))

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
    dfERAWO_wind['mag10_dfERAWO_MET1']=dfERAWO_wind['mag10_dfERAWO']
    dfERAWO_wind['mag{}_dfERAWO_MET1'.format(i)]=dfERAWO_wind['mag10_dfERAWO']*((i/10)**(alfa_oc))


#METODOLOGIA 2 - LEI LOGARITMA 1 (DNV (RUGOSIDADE COM VALOR TABELADO))    
for i in niveis:
    dfERAWO_wind['mag10_dfERAWO_MET2']=dfERAWO_wind['mag10_dfERAWO']
    dfERAWO_wind['mag{}_dfERAWO_MET2'.format(i)]=(ufric(10,dfERAWO_wind['mag10_dfERAWO'],zo_tab_maior)/k)*(beta(i,zo_tab_maior))


#METODOLOGIA 3 - LEI LOGARITMA 2 (PERFIL NEUTRO (PSI=0) + RUGOSIDADE DE DONELAN 90)    
#CALCULO DA RUGOSIDADE PELO Hs - Donelan 90
dfERAWO_wind['zo_Do_90']=0.033*(dfERAWO_wind['Hs_dfERAWO']/4)

for i in niveis:
    dfERAWO_wind['mag10_dfERAWO_MET3']=dfERAWO_wind['mag10_dfERAWO']
    dfERAWO_wind['mag{}_dfERAWO_MET3'.format(i)]=(ufric(10,dfERAWO_wind['mag10_dfERAWO'],dfERAWO_wind['zo_Do_90'])/k)*(beta(i,dfERAWO_wind['zo_Do_90']))


#METODOLOGIA 4 - LEI LOGARITMA 3 (PERFIL NEUTRO (PSI=0) + DONELAN ET AL. 93)    
#CALCULO DA RUGOSIDADE PELO Hs - DOnelan 93
dfERAWO_wind['zo_Do_93']=(dfERAWO_wind['Hs_dfERAWO']/4)*(6.7*(10**-4))*((dfERAWO_wind['mag10_dfERAWO']/Cp)**2.6)

for i in niveis:
    dfERAWO_wind['mag10_dfERAWO_MET4']=dfERAWO_wind['mag10_dfERAWO']
    dfERAWO_wind['mag{}_dfERAWO_MET4'.format(i)]=(ufric(10,dfERAWO_wind['mag10_dfERAWO'],dfERAWO_wind['zo_Do_93'])/k)*(beta(i,dfERAWO_wind['zo_Do_93']))


#METODOLOGIA 5 - LEI LOGARITMA 4 (PERFIL NEUTRO (PSI=0) + CARMO ET AL. 2022 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA))
#calculando Zo pelo metodo de CARMO 22 - adaptado de taylor e yelland
dfERAWO_wind['Lp']=(dfERAWO_wind['Tp_dfERAWO']**2*g)/(2*np.pi)   #descobrir oq é 'lp' e se tem q mudar algo nas variáveis
dfERAWO_wind['ufric_car22']=(k*dfERAWO_wind['mag10_dfERAWO'])/np.log(10/(mi*dfERAWO_wind['Hs_dfERAWO']*((dfERAWO_wind['Hs_dfERAWO']/dfERAWO_wind['Lp'])**4.5)))
dfERAWO_wind['zo_car22']=(mi*dfERAWO_wind['Hs_dfERAWO']*((dfERAWO_wind['Hs_dfERAWO']/dfERAWO_wind['Lp'])**4.5))

for i in niveis:
    dfERAWO_wind['mag10_dfERAWO_MET5']=dfERAWO_wind['mag10_dfERAWO']
    dfERAWO_wind['mag{}_dfERAWO_MET5'.format(i)]=(dfERAWO_wind['ufric_car22']/k)*((np.log((i)/dfERAWO_wind['zo_car22'])))


#METODOLOGIA 6 - LEI LOGARITMA 5 (PERFIL ESTÁVEL (PSI=4.7Z/L) + CARMO ET AL. 2022 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA))
for i in niveis:
    dfERAWO_wind['mag10_dfERAWO_MET6']=dfERAWO_wind['mag10_dfERAWO']
    dfERAWO_wind['mag{}_dfERAWO_MET6'.format(i)]=(dfERAWO_wind['ufric_car22']/k)*((np.log((i)/dfERAWO_wind['zo_car22']))+((4.7*(i))/L))

#            SOLUCIONAR O PROBLEMA DO CALCULO DA TSM
#METODOLOGIA 7 - LEI LOGARITMA 6 (CARMO ET AL. 2023 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA + CALCULO DE PSI PELA T e TSM))
dfERAWO_wind['EST_CORRIGIDO']=((np.abs(dfERAWO_wind['TSM_OISST']*dfERAWO_wind['T2m_ERA5']*0.7))/epsilon) ## TERMO DE AJUSTE DA ESTABILIDADE PELA TSM E T
dfERAWO_wind['zo_car23']=(mi*dfERAWO_wind['Hs_dfERAWO']*((dfERAWO_wind['Hs_dfERAWO']/dfERAWO_wind['Lp'])**4.5)) #DE REPENTE, PODE FAZER UMA CORRECAO NO ZO

for i in niveis:
    dfERAWO_wind['mag10_dfERAWO_MET7']=dfERAWO_wind['mag10_dfERAWO']
    dfERAWO_wind['mag{}_dfERAWO_MET7'.format(i)]=(dfERAWO_wind['ufric_car22']/k)*((np.log((i)/dfERAWO_wind['zo_car23']))+(((i)/dfERAWO_wind['EST_CORRIGIDO'])))

dfERAWO_wind.to_csv(pathname2+'/Caetano/Energia/PROCESSADOS/CPAM/ERAWO.csv')#COSNERTAR DEPOIS

# ------------------------------------------------------------------------------------------------------------

#                            ERA5+WAVERYS+MUR 

dfERAWM_wind = ds.sel(lon=-40.04083,lat=-22.43278, method='nearest').to_dataframe()    #ERA5
dfERAWM_wave = ds3.sel(lon=-40.04083,lat=-22.43278, method='nearest').to_dataframe()   #WAVERYS
dfERAWM_T    = ds2.sel(lon=-40.04083,lat=-22.43278, method='nearest').to_dataframe()   #ERA5
dfERAWM_TSM  = ds5.sel(lon=-40.04083,lat=-22.43278, method='nearest').to_dataframe()   #MUR  


#separa em intervalo de tempo desejado
dfERAWM_wind=dfERAWM_wind[(dfERAWM_wind.index>='2003-01-01 00:00') & (dfERAWM_wind.index<='2018-12-31 18:00')] 
dfERAWM_wave=dfERAWM_wave[(dfERAWM_wave.index>='2003-01-01 00:00') & (dfERAWM_wave.index<='2018-12-31 18:00')]
dfERAWM_T=dfERAWM_T[(dfERAWM_T.index>='2003-01-01 00:00') & (dfERAWM_T.index<='2018-12-31 18:00')] 
dfERAWM_TSM=dfERAWM_TSM[(dfERAWM_TSM.index>='2003-01-01 00:00') & (dfERAWM_TSM.index<='2018-12-31 18:00')]

# variaveis de vento, onda e temperatura
dfERAWM_wind['mag10_dfERAWM']=np.nan
dfERAWM_wind['mag10_dfERAWM']=((ERA5_wind['u10']**2)+(ERA5_wind['v10']**2))**0.5#1.7
dfERAWM_wind['dir10_dfERAWM']=np.nan
dfERAWM_wind['dir10_dfERAWM']=np.degrees(np.arctan2(-ERA5_wind['u10'],-ERA5_wind['v10']))
#dfERAWM_wind['dir10_dfERAWM'][dfERAWM_wind['dir10_dfERAWM']<0]=ERA5_wind['dir10_ERA5'] + 360
dfERAWM_wind.loc[dfERAWM_wind['dir10_dfERAWM'] < 0, 'dir10_dfERAWM'] = ERA5_wind['dir10_ERA5'] + 360
dfERAWM_wind['T2m_ERA5']=np.nan
dfERAWM_wind['T2m_ERA5']=ERA5_T['t2m'] # temperatura do ar a 2 metros 
dfERAWM_wind['TSM_MUR']=np.nan
dfERAWM_wind['TSM_MUR']=MUR['analysed_sst'] #TSM
dfERAWM_wind['Hs_dfERAWM']=np.nan
dfERAWM_wind['Hs_dfERAWM']=WAVERYS['VTPK'] # altura significativa de onda
dfERAWM_wind['Tp_dfERAWM']=np.nan
dfERAWM_wind['Tp_dfERAWM']=WAVERYS['VTM01_WW'] #periodo de pico da onda
dfERAWM_wind['Permed_dfERAWM']=np.nan
dfERAWM_wind['Permed_dfERAWM']=WAVERYS['VMDR'] #periodo medio de vento
#METODOLOGIA 1 - LEI DA POTENCIA
#METODOLOGIA 2 - LEI LOGARITMA 1 (DNV (RUGOSIDADE COM VALOR TABELADO))
#METODOLOGIA 3 - LEI LOGARITMA 2 (PERFIL NEUTRO (PSI=0) + RUGOSIDADE DE DONELAN 90)
#METODOLOGIA 4 - LEI LOGARITMA 3 (PERFIL NEUTRO (PSI=0) + DONELAN ET AL. 93)
#METODOLOGIA 5 - LEI LOGARITMA 4 (PERFIL NEUTRO (PSI=0) + CARMO ET AL. 2022 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA))
#METODOLOGIA 6 - LEI LOGARITMA 5 (PERFIL ESTÁVEL (PSI=4.7Z/L) + CARMO ET AL. 2022 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA))
#METODOLOGIA 7 - LEI LOGARITMA 6 (CARMO ET AL. 2023 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA + CALCULO DE PSI PELA T e TSM))

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
    dfERAWM_wind['mag10_dfERAWM_MET1']=dfERAWM_wind['mag10_dfERAWM']
    dfERAWM_wind['mag{}_dfERAWM_MET1'.format(i)]=dfERAWM_wind['mag10_dfERAWM']*((i/10)**(alfa_oc))


#METODOLOGIA 2 - LEI LOGARITMA 1 (DNV (RUGOSIDADE COM VALOR TABELADO))    
for i in niveis:
    dfERAWM_wind['mag10_dfERAWM_MET2']=dfERAWM_wind['mag10_dfERAWM']
    dfERAWM_wind['mag{}_dfERAWM_MET2'.format(i)]=(ufric(10,dfERAWM_wind['mag10_dfERAWM'],zo_tab_maior)/k)*(beta(i,zo_tab_maior))


#METODOLOGIA 3 - LEI LOGARITMA 2 (PERFIL NEUTRO (PSI=0) + RUGOSIDADE DE DONELAN 90)    
#CALCULO DA RUGOSIDADE PELO Hs - Donelan 90
dfERAWM_wind['zo_Do_90']=0.033*(dfERAWM_wind['Hs_dfERAWM']/4)

for i in niveis:
    dfERAWM_wind['mag10_dfERAWM_MET3']=dfERAWM_wind['mag10_dfERAWM']
    dfERAWM_wind['mag{}_dfERAWM_MET3'.format(i)]=(ufric(10,dfERAWM_wind['mag10_dfERAWM'],dfERAWM_wind['zo_Do_90'])/k)*(beta(i,dfERAWM_wind['zo_Do_90']))


#METODOLOGIA 4 - LEI LOGARITMA 3 (PERFIL NEUTRO (PSI=0) + DONELAN ET AL. 93)    
#CALCULO DA RUGOSIDADE PELO Hs - DOnelan 93
dfERAWM_wind['zo_Do_93']=(dfERAWM_wind['Hs_dfERAWM']/4)*(6.7*(10**-4))*((dfERAWM_wind['mag10_dfERAWM']/Cp)**2.6)

for i in niveis:
    dfERAWM_wind['mag10_dfERAWM_MET4']=dfERAWM_wind['mag10_dfERAWM']
    dfERAWM_wind['mag{}_dfERAWM_MET4'.format(i)]=(ufric(10,dfERAWM_wind['mag10_dfERAWM'],dfERAWM_wind['zo_Do_93'])/k)*(beta(i,dfERAWM_wind['zo_Do_93']))


#METODOLOGIA 5 - LEI LOGARITMA 4 (PERFIL NEUTRO (PSI=0) + CARMO ET AL. 2022 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA))
#calculando Zo pelo metodo de CARMO 22 - adaptado de taylor e yelland
dfERAWM_wind['Lp']=(dfERAWM_wind['Tp_dfERAWM']**2*g)/(2*np.pi)   #descobrir oq é 'lp' e se tem q mudar algo nas variáveis
dfERAWM_wind['ufric_car22']=(k*dfERAWM_wind['mag10_dfERAWM'])/np.log(10/(mi*dfERAWM_wind['Hs_dfERAWM']*((dfERAWM_wind['Hs_dfERAWM']/dfERAWM_wind['Lp'])**4.5)))
dfERAWM_wind['zo_car22']=(mi*dfERAWM_wind['Hs_dfERAWM']*((dfERAWM_wind['Hs_dfERAWM']/dfERAWM_wind['Lp'])**4.5))

for i in niveis:
    dfERAWM_wind['mag10_dfERAWM_MET5']=dfERAWM_wind['mag10_dfERAWM']
    dfERAWM_wind['mag{}_dfERAWM_MET5'.format(i)]=(dfERAWM_wind['ufric_car22']/k)*((np.log((i)/dfERAWM_wind['zo_car22'])))


#METODOLOGIA 6 - LEI LOGARITMA 5 (PERFIL ESTÁVEL (PSI=4.7Z/L) + CARMO ET AL. 2022 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA))
for i in niveis:
    dfERAWM_wind['mag10_dfERAWM_MET6']=dfERAWM_wind['mag10_dfERAWM']
    dfERAWM_wind['mag{}_dfERAWM_MET6'.format(i)]=(dfERAWM_wind['ufric_car22']/k)*((np.log((i)/dfERAWM_wind['zo_car22']))+((4.7*(i))/L))

#            SOLUCIONAR O PROBLEMA DO CALCULO DA TSM
#METODOLOGIA 7 - LEI LOGARITMA 6 (CARMO ET AL. 2023 (RUGOSIDADE DE TAYLOR & YELLAND CORRIGIDA + CALCULO DE PSI PELA T e TSM))
dfERAWM_wind['EST_CORRIGIDO']=((np.abs(dfERAWM_wind['TSM_MUR']*dfERAWM_wind['T2m_ERA5']*0.7))/epsilon) ## TERMO DE AJUSTE DA ESTABILIDADE PELA TSM E T
dfERAWM_wind['zo_car23']=(mi*dfERAWM_wind['Hs_dfERAWM']*((dfERAWM_wind['Hs_dfERAWM']/dfERAWM_wind['Lp'])**4.5)) #DE REPENTE, PODE FAZER UMA CORRECAO NO ZO

for i in niveis:
    dfERAWM_wind['mag10_dfERAWM_MET7']=dfERAWM_wind['mag10_dfERAWM']
    dfERAWM_wind['mag{}_dfERAWM_MET7'.format(i)]=(dfERAWM_wind['ufric_car22']/k)*((np.log((i)/dfERAWM_wind['zo_car23']))+(((i)/dfERAWM_wind['EST_CORRIGIDO'])))

dfERAWM_wind.to_csv(pathname2+'/Caetano/Energia/PROCESSADOS/CPAM/ERAWM.csv') #feito até aqui

# ------------------------------------------------------------------------------------------------------------
#                                CRIANDO DF VAZIO PARA CALCULAR OS PERFIS

df_perfis = pd.DataFrame(randn(17,1),index='10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170'.split(),
                  columns='z'.split())
dados =['ERA5','ERAW','ERAWO','ERAWM'] 



for local in dados:
    df_perfis['Met1_{}'.format(local)]=np.nan  #zo_tab1 = 0.001
    df_perfis['Met2_{}'.format(local)]=np.nan   #Lei da potencia alfa
    df_perfis['Met3_{}'.format(local)]=np.nan   #zo_Don1990
    df_perfis['Met4_{}'.format(local)]=np.nan   #zo_Don1993
    df_perfis['Met5_{}'.format(local)]=np.nan   #zo_Charnock1955 - solução Oost 2002
    df_perfis['Met6_{}'.format(local)]=np.nan   #zo_He2019
    df_perfis['Met7_{}'.format(local)]=np.nan   #zo_Don90 - Est
    

i=0
while i < len(niveis):
    df_perfis['Met1_ERA5'][i]=ERA5_wind['mag{}_ERA5_MET1'.format(df_perfis.index[i])].mean()
    df_perfis['Met1_ERAW'][i]=dfERAW_wind['mag{}_dfERAW_MET1'.format(df_perfis.index[i])].mean()
    df_perfis['Met1_ERAWO'][i]=dfERAWO_wind['mag{}_dfERAWO_MET1'.format(df_perfis.index[i])].mean()
    df_perfis['Met1_ERAWM'][i]=dfERAWM_wind['mag{}_dfERAWM_MET1'.format(df_perfis.index[i])].mean()
    
    df_perfis['Met2_ERA5'][i]=ERA5_wind['mag{}_ERA5_MET2'.format(df_perfis.index[i])].mean()
    df_perfis['Met2_ERAW'][i]=dfERAW_wind['mag{}_dfERAW_MET2'.format(df_perfis.index[i])].mean()
    df_perfis['Met2_ERAWO'][i]=dfERAWO_wind['mag{}_dfERAWO_MET2'.format(df_perfis.index[i])].mean()
    df_perfis['Met2_ERAWM'][i]=dfERAWM_wind['mag{}_dfERAWM_MET2'.format(df_perfis.index[i])].mean()
    
    df_perfis['Met3_ERA5'][i]=ERA5_wind['mag{}_ERA5_MET3'.format(df_perfis.index[i])].mean()
    df_perfis['Met3_ERAW'][i]=dfERAW_wind['mag{}_dfERAW_MET3'.format(df_perfis.index[i])].mean()
    df_perfis['Met3_ERAWO'][i]=dfERAWO_wind['mag{}_dfERAWO_MET3'.format(df_perfis.index[i])].mean()
    df_perfis['Met3_ERAWM'][i]=dfERAWM_wind['mag{}_dfERAWM_MET3'.format(df_perfis.index[i])].mean()
    
    df_perfis['Met4_ERA5'][i]=ERA5_wind['mag{}_ERA5_MET4'.format(df_perfis.index[i])].mean()
    df_perfis['Met4_ERAW'][i]=dfERAW_wind['mag{}_dfERAW_MET4'.format(df_perfis.index[i])].mean()
    df_perfis['Met4_ERAWO'][i]=dfERAWO_wind['mag{}_dfERAWO_MET4'.format(df_perfis.index[i])].mean()
    df_perfis['Met4_ERAWM'][i]=dfERAWM_wind['mag{}_dfERAWM_MET4'.format(df_perfis.index[i])].mean()
    
    df_perfis['Met5_ERA5'][i]=ERA5_wind['mag{}_ERA5_MET5'.format(df_perfis.index[i])].mean()
    df_perfis['Met5_ERAW'][i]=dfERAW_wind['mag{}_dfERAW_MET5'.format(df_perfis.index[i])].mean()
    df_perfis['Met5_ERAWO'][i]=dfERAWO_wind['mag{}_dfERAWO_MET5'.format(df_perfis.index[i])].mean()
    df_perfis['Met5_ERAWM'][i]=dfERAWM_wind['mag{}_dfERAWM_MET5'.format(df_perfis.index[i])].mean()
    
    df_perfis['Met6_ERA5'][i]=ERA5_wind['mag{}_ERA5_MET6'.format(df_perfis.index[i])].mean()
    df_perfis['Met6_ERAW'][i]=dfERAW_wind['mag{}_dfERAW_MET6'.format(df_perfis.index[i])].mean()
    df_perfis['Met6_ERAWO'][i]=dfERAWO_wind['mag{}_dfERAWO_MET6'.format(df_perfis.index[i])].mean()
    df_perfis['Met6_ERAWM'][i]=dfERAWM_wind['mag{}_dfERAWM_MET6'.format(df_perfis.index[i])].mean()
    
    #df_perfis['Met7_ERA5'][i]=ERA5_wind['mag{}_ERA5_MET7'.format(df_perfis.index[i])].mean()
   # df_perfis['Met7_ERAW'][i]=dfERAW_wind['mag{}_dfERAW_MET7'.format(df_perfis.index[i])].mean()
   # df_perfis['Met7_ERAWO'][i]=dfERAWO_wind['mag{}_dfERAWO_MET7'.format(df_perfis.index[i])].mean()
    #df_perfis['Met7_ERAWM'][i]=dfERAWM_wind['mag{}_dfERAWM_MET7'.format(df_perfis.index[i])].mean()
    
    
    i=i+1
  
# ------------------------------------------------------------------------------------------------------------

#                                                 GRAFICOS DIREÇÃO 

#    ERA5
ax = WindroseAxes.from_ax()
ax.bar(ERA5_wind['dir10_ERA5'], ERA5_wind['mag10_ERA5'], normed=True, opening=0.8, edgecolor='black', cmap=plt.cm.winter, lw=1, bins=np.arange(0.01,16,2))
ax.set_title('Rosa dos ventos para P18  - 2003 to 2018 (ERA5)', fontsize=12)
ax.set_legend(loc=2,decimal_places=1)
#logo = plt.imread('ladsin_2.png')
#ax.figure.figimage(logo, xo = 60, yo = 70)
plt.savefig('/media/ladsin/Seagate Expansion Drive/Caetano/Energia/Plots/rosa_dos_ventos_ERA5_Cpam.png', fontsize=15)

#   ERA5 + WAVERYS 
ax = WindroseAxes.from_ax()
ax.bar(dfERAW_wind['dir10_dfERAW'], dfERAW_wind['mag10_dfERAW'], normed=True, opening=0.8, edgecolor='black', cmap=plt.cm.winter, lw=1, bins=np.arange(0.01,16,2))
ax.set_title('Rosa dos ventos para P18 - 1997 a 2018 (ERA5 + WAVERYS)', fontsize=12)
ax.set_legend(loc=2,decimal_places=1)
#logo = plt.imread('ladsin_2.png')
#ax.figure.figimage(logo, xo = 60, yo = 70)
plt.savefig('/media/ladsin/Seagate Expansion Drive/Caetano/energia/Plots/rosa_dos_ventos_ERAW.png', fontsize=15)

#   ERA5 + WAVERYS + OISST
ax = WindroseAxes.from_ax()
ax.bar(dfERAWO_wind['dir10_dfERAWO'], dfERAWO_wind['mag10_dfERAWO'], normed=True, opening=0.8, edgecolor='black', cmap=plt.cm.winter, lw=1, bins=np.arange(0.01,16,2))
ax.set_title('Rosa dos ventos para P18 - 1997 a 2018 (ERA5 + WAVERYS + OISST)', fontsize=12)
ax.set_legend(loc=2,decimal_places=1)
#logo = plt.imread('ladsin_2.png')
#ax.figure.figimage(logo, xo = 60, yo = 70)
plt.savefig('/media/ladsin/Seagate Expansion Drive/Caetano/energia/Plots/rosa_dos_ventos_ERAWO.png', fontsize=15)

#   ERA5 + WAVERYS + MUR
ax = WindroseAxes.from_ax()
ax.bar(dfERAW_wind['dir10_dfERAWM'], dfERAW_wind['mag10_dfERAWM'], normed=True, opening=0.8, edgecolor='black', cmap=plt.cm.winter, lw=1, bins=np.arange(0.01,16,2))
ax.set_title('Rosa dos ventos para P18 - 1997 a 2018 (ERA5 + WAVERYS + MUR)', fontsize=12)
ax.set_legend(loc=2,decimal_places=1)
#logo = plt.imread('ladsin_2.png')
#ax.figure.figimage(logo, xo = 60, yo = 70)
plt.savefig('/media/ladsin/Seagate Expansion Drive/Caetano/energia/Plots/rosa_dos_ventos_ERAWM.png', fontsize=15)


# ------------------------------------------------------------------------------------------------------------

#                                                PERFIS DE VENTO

#      PERFIS DE VENTO - ERA5
fig,ax=plt.subplots(figsize=(8,15))
line1= ax.plot(df_perfis['Met1_ERA5'],df_perfis.index,color='red',label='1/7 PowerLaw', linestyle='dotted')
line2= ax.plot(df_perfis['Met2_ERA5'],df_perfis.index,color='orange',label='DNV', linestyle='-.')
line3= ax.plot(df_perfis['Met3_ERA5'],df_perfis.index,color='green',label='Zo Donelan 90', linestyle=':')
line4= ax.plot(df_perfis['Met4_ERA5'],df_perfis.index, color='purple',label='Zo Donelan 93', linestyle='-')
line5= ax.plot(df_perfis['Met5_ERA5'],df_perfis.index, color='gray',label='Carmo 22 - Neu', linestyle='dashed')
line6= ax.plot(df_perfis['Met6_ERA5'],df_perfis.index, color='blue',label='Carmo 22 - Est', linestyle='-')
line7= ax.plot(df_perfis['Met7_ERA5'],df_perfis.index, color='yellow',label='Carmo 22 - Inst', linestyle='-') #estourando muito no plot
#line7= ax.plot(df_perfis['Met9_ERA5'],df_perfis.index, color='red',label='Zo TayYel01 - Stable', linestyle='-.')
#line8=plt.scatter(df_perfis['SODAR'],df_perfis.index, marker="*", color='black', label='SODAR')

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
ax.set_title('P18 - 2003 a 2018 (ERA5)', fontsize=20)
ax.set_ylabel('Altura (m)', fontsize=20)
ax.set_xlabel('Velocidade (m/s)', fontsize=20)
plt.savefig('/media/ladsin/Seagate Expansion Drive/Caetano/Energia/Plots/perfil_vento_ERA5_cpam2.png') #coloca 2 digitos
plt.clf() #fecha o arquivo para nao plotar por cima


#      PERFIS DE VENTO - ERA5 + WAVERYS
fig,ax=plt.subplots(figsize=(8,15))
line1= ax.plot(df_perfis['Met1_ERAW'],df_perfis.index,color='red',label='1/7 PowerLaw', linestyle='dotted')
line2= ax.plot(df_perfis['Met2_ERAW'],df_perfis.index,color='orange',label='DNV', linestyle='-.')
line3= ax.plot(df_perfis['Met3_ERAW'],df_perfis.index,color='green',label='Zo Donelan 90', linestyle=':')
line4= ax.plot(df_perfis['Met4_ERAW'],df_perfis.index, color='purple',label='Zo Donelan 93', linestyle='-')
line5= ax.plot(df_perfis['Met5_ERAW'],df_perfis.index, color='gray',label='Carmo 22 - Neu', linestyle='dashed')
line6= ax.plot(df_perfis['Met6_ERAW'],df_perfis.index, color='blue',label='Carmo 22 - Est', linestyle='-')
#line7= ax.plot(df_perfis['Met9_ERAW'],df_perfis.index, color='red',label='Zo TayYel01 - Stable', linestyle='-.')

ax.set_xlim(6,16)
ax.legend()
#criando o for para fazer o vetor - ERA5 + WAVERYS
for k in range (0,np.size(df_perfis.index)):
    qax = plt.quiver(0,df_perfis.index[k],df_perfis['Met6_ERAW'][k],0, units="dots",
             width=2.5,
             headwidth = 2,
             headlength = 4,
             headaxislength = 5,
             color='black',
             scale=0.0225)
ax.set_title('P18 - 2003 a 2018 (ERA5 + WAVERYS)', fontsize=20)
ax.set_ylabel('Altura (m)', fontsize=20)
ax.set_xlabel('Velocidade (m/s)', fontsize=20)
plt.savefig('/media/ladsin/Seagate Expansion Drive/Caetano/Energia/Plots/perfil_vento_ERAW_cpam_HsERA5.png') #coloca 2 digitos
plt.clf() #fecha o arquivo para nao plotar por cima



#      PERFIS DE VENTO - ERA5 + WAVERYS + OISST
fig,ax=plt.subplots(figsize=(8,15))
line1= ax.plot(df_perfis['Met1_ERAWO'],df_perfis.index,color='red',label='1/7 PowerLaw', linestyle='dotted')
line2= ax.plot(df_perfis['Met2_ERAWO'],df_perfis.index,color='orange',label='DNV', linestyle='-.')
line3= ax.plot(df_perfis['Met3_ERAWO'],df_perfis.index,color='green',label='Zo Donelan 90', linestyle=':')
line4= ax.plot(df_perfis['Met4_ERAWO'],df_perfis.index, color='purple',label='Zo Donelan 93', linestyle='-')
line5= ax.plot(df_perfis['Met5_ERAWO'],df_perfis.index, color='gray',label='Carmo 22 - Neu', linestyle='dashed')
line6= ax.plot(df_perfis['Met6_ERAWO'],df_perfis.index, color='blue',label='Carmo 22 - Est', linestyle='-')
#line7= ax.plot(df_perfis['Met9_ERAWO'],df_perfis.index, color='red',label='Zo TayYel01 - Stable', linestyle='-.')

ax.set_xlim(6,16)
ax.legend()
#criando o for para fazer o vetor - ERA5 + WAVERYS + OISST
for k in range (0,np.size(df_perfis.index)):
    qax = plt.quiver(0,df_perfis.index[k],df_perfis['Met6_ERAWO'][k],0, units="dots",
             width=2.5,
             headwidth = 2,
             headlength = 4,
             headaxislength = 5,
             color='gray',
             scale=0.0225)
ax.set_title('P18 - 2003 a 2018 (ERA5 + WAVERYS + OISST)', fontsize=20)
ax.set_ylabel('Altura (m)', fontsize=20)
ax.set_xlabel('Velocidade (m/s)', fontsize=20)
plt.savefig('/media/ladsin/Seagate Expansion Drive/Caetano/Energia/Plots/perfil_vento_ERAWO_cpam.png') #coloca 2 digitos
plt.clf() #fecha o arquivo para nao plotar por cima

#      PERFIS DE VENTO - ERA5 + WAVERYS + MUR
fig,ax=plt.subplots(figsize=(8,15))
line1= ax.plot(df_perfis['Met1_ERAWM'],df_perfis.index,color='red',label='1/7 PowerLaw', linestyle='dotted')
line2= ax.plot(df_perfis['Met2_ERAWM'],df_perfis.index,color='orange',label='DNV ', linestyle='-.')
line3= ax.plot(df_perfis['Met3_ERAWM'],df_perfis.index,color='green',label='Zo Donelan 90', linestyle=':')
line4= ax.plot(df_perfis['Met4_ERAWM'],df_perfis.index, color='purple',label='Zo Donelan 93', linestyle='-')
line5= ax.plot(df_perfis['Met5_ERAWM'],df_perfis.index, color='gray',label='Carmo 22 - Neu', linestyle='dashed')
line6= ax.plot(df_perfis['Met6_ERAWM'],df_perfis.index, color='blue',label='Carmo 22 - Est', linestyle='-')
#line7= ax.plot(df_perfis['Met9_ERAWM'],df_perfis.index, color='red',label='Zo TayYel01 - Stable', linestyle='-.')

ax.set_xlim(6,16)
ax.legend()
#criando o for para fazer o vetor - ERA5 + WAVERYS + MUR
for k in range (0,np.size(df_perfis.index)):
    qax = plt.quiver(0,df_perfis.index[k],df_perfis['Met6_ERAWM'][k],0, units="dots",
             width=2.5,
             headwidth = 2,
             headlength = 4,
             headaxislength = 5,
             color='black',
             scale=0.0225)
ax.set_title('P18 - 2003 a 2018 (ERA5 + WAVERYS + MUR)', fontsize=20)
ax.set_ylabel('Altura (m)', fontsize=20)
ax.set_xlabel('Velocidade (m/s)', fontsize=20)
plt.savefig('/media/ladsin/Seagate Expansion Drive/Caetano/Energia/Plots/perfil_vento_ERAWM_cpam.png') #coloca 2 digitos
plt.clf() #fecha o arquivo para nao plotar por cima

# ------------------------------------------------------------------------------------------------------------

#                            Calculo do PotenciaL Eólico ERA5
Drotor_3000kW=82
cp_3000kW=0.45 
dens_ar=1.22 #kgm-³
A=5281 #ou (np.pi*(Drotor_3000kW**2))/4
df_perfis['Pot_MET1_ERA5']=(0.5*dens_ar*A*(df_perfis['Met1_ERA5']**3)*cp_3000kW)/1000000
df_perfis['Pot_MET2_ERA5']=(0.5*dens_ar*A*(df_perfis['Met2_ERA5']**3)*cp_3000kW)/1000000
df_perfis['Pot_MET3_ERA5']=(0.5*dens_ar*A*(df_perfis['Met3_ERA5']**3)*cp_3000kW)/1000000
df_perfis['Pot_MET4_ERA5']=(0.5*dens_ar*A*(df_perfis['Met4_ERA5']**3)*cp_3000kW)/1000000
df_perfis['Pot_MET5_ERA5']=(0.5*dens_ar*A*(df_perfis['Met5_ERA5']**3)*cp_3000kW)/1000000
df_perfis['Pot_MET6_ERA5']=(0.5*dens_ar*A*(df_perfis['Met6_ERA5']**3)*cp_3000kW)/1000000

#                            PERFIS DE VENTO DO POTENCIAL EÓLICO - ERA5
fig,ax=plt.subplots(figsize=(8,15))
line1= ax.plot(df_perfis['Pot_MET1_ERA5'],df_perfis.index,color='red',label='1/7 PowerLaw', linestyle='dotted')
line2= ax.plot(df_perfis['Pot_MET2_ERA5'],df_perfis.index,color='orange',label='DNV', linestyle='-.')
line3= ax.plot(df_perfis['Pot_MET3_ERA5'],df_perfis.index,color='green',label='Zo Donelan 90', linestyle=':')
line4= ax.plot(df_perfis['Pot_MET4_ERA5'],df_perfis.index, color='purple',label='Zo Donelan 93', linestyle='-')
line5= ax.plot(df_perfis['Pot_MET5_ERA5'],df_perfis.index, color='gray',label='Carmo 22 - Neu', linestyle='dashed')
line6= ax.plot(df_perfis['Pot_MET6_ERA5'],df_perfis.index, color='blue',label='Carmo 22 - Est', linestyle='-')
#line7= ax.plot(df_perfis['Met9_ERA5'],df_perfis.index, color='red',label='Zo TayYel01 - Stable', linestyle='-.')
#line8=plt.scatter(df_perfis['SODAR'],df_perfis.index, marker="*", color='black', label='SODAR')


ax.set_xlim(0,3)
ax.legend()
#criando o for para fazer o vetor - ERA5
ax.set_title('Potencial Eólico P18 - 1997 a 2014 (ERA5)', fontsize=20)
ax.set_ylabel('Altura (m)', fontsize=20)
ax.set_xlabel('Potência (MW)', fontsize=20)
plt.savefig('/media/ladsin/Seagate Expansion Drive/Caetano/energia/Plots/Pot_Eol_ERA5.png') #coloca 2 digitos
plt.clf() #fecha o arquivo para nao plotar por cima

# ------------------------------------------------------------------------------------------------------------

#                            Calculo do PotenciaL Eólico ERA5 + WAVERYS
Drotor_3000kW=82
cp_3000kW=0.45 
dens_ar=1.22 #kgm-³
A=5281 #ou (np.pi*(Drotor_3000kW**2))/4
df_perfis['Pot_MET1_ERAW']=(0.5*dens_ar*A*(df_perfis['Met1_ERAW']**3)*cp_3000kW)/1000000
df_perfis['Pot_MET2_ERAW']=(0.5*dens_ar*A*(df_perfis['Met2_ERAW']**3)*cp_3000kW)/1000000
df_perfis['Pot_MET3_ERAW']=(0.5*dens_ar*A*(df_perfis['Met3_ERAW']**3)*cp_3000kW)/1000000
df_perfis['Pot_MET4_ERAW']=(0.5*dens_ar*A*(df_perfis['Met4_ERAW']**3)*cp_3000kW)/1000000
df_perfis['Pot_MET5_ERAW']=(0.5*dens_ar*A*(df_perfis['Met5_ERAW']**3)*cp_3000kW)/1000000
df_perfis['Pot_MET6_ERAW']=(0.5*dens_ar*A*(df_perfis['Met6_ERAW']**3)*cp_3000kW)/1000000

#                            PERFIS DE VENTO - ERA5 + WAVERYS
fig,ax=plt.subplots(figsize=(8,15))
line1= ax.plot(df_perfis['Pot_MET1_ERAW'],df_perfis.index,color='red',label='1/7 PowerLaw', linestyle='dotted')
line2= ax.plot(df_perfis['Pot_MET2_ERAW'],df_perfis.index,color='orange',label='DNV', linestyle='-.')
line3= ax.plot(df_perfis['Pot_MET3_ERAW'],df_perfis.index,color='green',label='Zo Donelan 90', linestyle=':')
line4= ax.plot(df_perfis['Pot_MET4_ERAW'],df_perfis.index, color='purple',label='Zo Donelan 93', linestyle='-')
line5= ax.plot(df_perfis['Pot_MET5_ERAW'],df_perfis.index, color='gray',label='Carmo 22 - Neu', linestyle='dashed')
#line6= ax.plot(df_perfis['Pot_MET6_ERAW'],df_perfis.index, color='blue',label='Carmo 22 - Est', linestyle='-')
#line7= ax.plot(df_perfis['Met9_CE'],df_perfis.index, color='red',label='Zo TayYel01 - Stable', linestyle='-.')

ax.set_xlim(0,3)
ax.legend()
#criando o for para fazer o vetor - ERA5 + WAVERYS
ax.set_title('Pontencial Eólico P18 - 1997 a 2018 (ERAW)', fontsize=20)
ax.set_ylabel('Altura (m)', fontsize=20)
ax.set_xlabel('Potência (MW)', fontsize=20)
plt.savefig('/media/ladsin/Seagate Expansion Drive/Caetano/energia/Plots/Pot_Eol_ERAW.png') #coloca 2 digitos
plt.clf() #fecha o arquivo para nao plotar por cima

# ------------------------------------------------------------------------------------------------------------

#                            Calculo do PotenciaL Eólico ERA5 + WAVERYS + OISST
Drotor_3000kW=82
cp_3000kW=0.45 
dens_ar=1.22 #kgm-³
A=5281 #ou (np.pi*(Drotor_3000kW**2))/4
df_perfis['Pot_MET1_ERAWO']=(0.5*dens_ar*A*(df_perfis['Met1_ERAWO']**3)*cp_3000kW)/1000000
df_perfis['Pot_MET2_ERAWO']=(0.5*dens_ar*A*(df_perfis['Met2_ERAWO']**3)*cp_3000kW)/1000000
df_perfis['Pot_MET3_ERAWO']=(0.5*dens_ar*A*(df_perfis['Met3_ERAWO']**3)*cp_3000kW)/1000000
df_perfis['Pot_MET4_ERAWO']=(0.5*dens_ar*A*(df_perfis['Met4_ERAWO']**3)*cp_3000kW)/1000000
df_perfis['Pot_MET5_ERAWO']=(0.5*dens_ar*A*(df_perfis['Met5_ERAWO']**3)*cp_3000kW)/1000000
df_perfis['Pot_MET6_ERAWO']=(0.5*dens_ar*A*(df_perfis['Met6_ERAWO']**3)*cp_3000kW)/1000000

#                            PERFIS DE VENTO - ERA5 + WAVERYS + OISST
fig,ax=plt.subplots(figsize=(8,15))
line1= ax.plot(df_perfis['Pot_MET1_ERAWO'],df_perfis.index,color='red',label='1/7 PowerLaw', linestyle='dotted')
line2= ax.plot(df_perfis['Pot_MET2_ERAWO'],df_perfis.index,color='orange',label='DNV', linestyle='-.')
line3= ax.plot(df_perfis['Pot_MET3_ERAWO'],df_perfis.index,color='green',label='Zo Donelan 90', linestyle=':')
line4= ax.plot(df_perfis['Pot_MET4_ERAWO'],df_perfis.index, color='purple',label='Zo Donelan 93', linestyle='-')
line5= ax.plot(df_perfis['Pot_MET5_ERAWO'],df_perfis.index, color='gray',label='Carmo 22 - Neu', linestyle='dashed')
line6= ax.plot(df_perfis['Pot_MET6_ERAWO'],df_perfis.index, color='blue',label='Carmo 22 - Est', linestyle='-')
#line7= ax.plot(df_perfis['Met9_ERAWO'],df_perfis.index, color='red',label='Zo TayYel01 - Stable', linestyle='-.')

ax.set_xlim(0,3)
ax.legend()
#criando o for para fazer o vetor - ERA5 + WAVERYS + OISST
ax.set_title('Potencial Eólico Santos - SP - 1997 a 2018 (PNBOIA)', fontsize=20)
ax.set_ylabel('Altura (m)', fontsize=20)
ax.set_xlabel('Potência (MW)', fontsize=20)
plt.savefig('/media/ladsin/Seagate Expansion Drive/Caetano/energia/Plots/pOT_Eol_ERAWO.png') #coloca 2 digitos
plt.clf() #fecha o arquivo para nao plotar por cima

# ------------------------------------------------------------------------------------------------------------

#                           Calculo do PotenciaL Eólico ERA5 + WAVERYS + MUR
Drotor_3000kW=82
cp_3000kW=0.45 
dens_ar=1.22 #kgm-³
A=5281 #ou (np.pi*(Drotor_3000kW**2))/4
df_perfis['Pot_MET1_ERAWM']=(0.5*dens_ar*A*(df_perfis['Met1_ERAWM']**3)*cp_3000kW)/1000000
df_perfis['Pot_MET2_ERAWM']=(0.5*dens_ar*A*(df_perfis['Met2_ERAWM']**3)*cp_3000kW)/1000000
df_perfis['Pot_MET3_ERAWM']=(0.5*dens_ar*A*(df_perfis['Met3_ERAWM']**3)*cp_3000kW)/1000000
df_perfis['Pot_MET4_ERAWM']=(0.5*dens_ar*A*(df_perfis['Met4_ERAWM']**3)*cp_3000kW)/1000000
df_perfis['Pot_MET5_ERAWM']=(0.5*dens_ar*A*(df_perfis['Met5_ERAWM']**3)*cp_3000kW)/1000000
df_perfis['Pot_MET6_ERAWM']=(0.5*dens_ar*A*(df_perfis['Met6_ERAWM']**3)*cp_3000kW)/1000000

#                           PERFIS DE VENTO - ERA5 + WAVERYS + MUR
fig,ax=plt.subplots(figsize=(8,15))
line1= ax.plot(df_perfis['Pot_MET1_ERAWM'],df_perfis.index,color='red',label='1/7 PowerLaw', linestyle='dotted')
line2= ax.plot(df_perfis['Pot_MET2_ERAWM'],df_perfis.index,color='orange',label='DNV ', linestyle='-.')
line3= ax.plot(df_perfis['Pot_MET3_ERAWM'],df_perfis.index,color='green',label='Zo Donelan 90', linestyle=':')
line4= ax.plot(df_perfis['Pot_MET4_ERAWM'],df_perfis.index, color='purple',label='Zo Donelan 93', linestyle='-')
line5= ax.plot(df_perfis['Pot_MET5_ERAWM'],df_perfis.index, color='gray',label='Carmo 22 - Neu', linestyle='dashed')
line6= ax.plot(df_perfis['Pot_MET6_ERAWM'],df_perfis.index, color='blue',label='Carmo 22 - Est', linestyle='-')
#line7= ax.plot(df_perfis['Met9_ERAWM'],df_perfis.index, color='red',label='Zo TayYel01 - Stable', linestyle='-.')

ax.set_xlim(0,3)
ax.legend()
#criando o for para fazer o vetor - ERA5 + WAVERYS + MUR
ax.set_title('Potencial Eólico P18 - 1990 a 2019 (ERA5 + WAVERYS + MUR)', fontsize=20)
ax.set_ylabel('Altura (m)', fontsize=20)
ax.set_xlabel('Potência (MW)', fontsize=20)
plt.savefig('/media/ladsin/Seagate Expansion Drive/Caetano/energia/Plots/Pot_Eol_ERAWM.png') #coloca 2 digitos
plt.clf() #fecha o arquivo para nao plotar por cima

