#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 08:41:48 2023

@author: Luiz do Carmo
"""

#### procesando dados concomitantes

import pandas as pd
import xarray as xr
import numpy as np
from windrose import WindroseAxes
import matplotlib.pyplot as plt
randn = np.random.randn
from PIL import Image
from numpy import asarray, max, arange, round, insert, radians
from numpy import ceil, ma, cumsum, array, argmin
from numpy import linspace, meshgrid, histogram2d, flipud, size, sum
from numpy import nanmax, nanmean, nansum
import xlsxwriter
import glob
import json
from bokeh.plotting import figure, output_file, show
from scipy.stats.stats import pearsonr
import skill_metrics as sm
import seaborn as sns
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.io.shapereader as shpreader # Import shapefiles


#diretorio principal
pathname='/home/ladsin/Documents/ladsin/Petrobras/LIDAR/'


#CONSTANTES
dens_ag=1025 #kgm-³
dens_ar=1.22 #kgm-³
alfa_oc=0.12
beta=0.7
epsilon=1200
E=0.033
mi=1200
a=2.7
b=0.142
c=13.09
k=0.4
g=9.80665
Cp=1.67
fc=0.5
#niveissodar=['40','50','60','80','100','120','140','160']
niveis=[10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170]
#niveis2=[5,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160]
#metodologias=['MET1.1_EST','MET1.1_NEU','MET1.2_EST','MET1.2_NEU',
#         'MET2','MET2','MET3_EST','MET3_NEU','MET4_EST',
#         'MET4_NEU','MET5_EST','MET5_NEU','MET3_NEU_P18',
#         'MET3_EST_P18','MET4_NEU_P18','MET4_EST_P18',
#         'MET5_NEU_P18','MET5_EST_P18','MET6_NEU_P18',
#         'MET6A_NEU_P18']


##ARQUIVO GERAL - CONTÉM TODOS AS DATAS (EXCETO JANEIRO)

###DADOS LIDAR PRA-1

LIDAR_10M = pd.read_csv(pathname+'DADOS_TODOS/LIDAR/lidar_PROC.csv',sep='\t')
#colocando a data no indice
LIDAR_10M.index=LIDAR_10M['Date']+' '+LIDAR_10M['hour']
#Transformando o indice de string pra datetime
LIDAR_10M.index=pd.to_datetime(LIDAR_10M.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME
LIDAR_1H=LIDAR_10M.resample('1H').last() #reamostrando para 1H


#ARQUIVO FEVEREIRO - CASO DO DIA 02 A 05 DE FEVEREIRO DE 2020
#DADOS LIDAR PRA-1

LIDAR_FEV_10M = pd.read_csv(pathname+'DADOS_TODOS/LIDAR/LiDAR_PRA01_02-05_fev_2020.csv',sep='\t')
#colocando a data no indice
LIDAR_FEV_10M.index=LIDAR_FEV_10M['Date']+' '+LIDAR_FEV_10M['hour']
#Transformando o indice de string pra datetime
LIDAR_FEV_10M.index=pd.to_datetime(LIDAR_FEV_10M.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME
LIDAR_FEV_1H=LIDAR_FEV_10M.resample('1H').last() #reamostrando para 1H

#ARQUIVO MARÇO - CASO DO DIA 01 A 05 DE MARÇO DE 2020
#DADOS LIDAR PRA-1

LIDAR_MAR_10M = pd.read_csv(pathname+'DADOS_TODOS/LIDAR/LiDAR_PRA01_01-05_mar_2020.csv',sep='\t')
#colocando a data no indice
LIDAR_MAR_10M.index=LIDAR_MAR_10M['Date']+' '+LIDAR_MAR_10M['hour']
#Transformando o indice de string pra datetime
LIDAR_MAR_10M.index=pd.to_datetime(LIDAR_MAR_10M.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME
LIDAR_MAR_1H=LIDAR_MAR_10M.resample('1H').last() #reamostrando para 1H

#ARQUIVO MAIO - CASO DO DIA 20 A 21 DE MAIO DE 2020
#DADOS LIDAR PRA-1

LIDAR_MAI_10M = pd.read_csv(pathname+'DADOS_TODOS/LIDAR/LiDAR_PRA01_20-21_mai_2020.csv',sep='\t')
#colocando a data no indice
LIDAR_MAI_10M.index=LIDAR_MAI_10M['Date']+' '+LIDAR_MAI_10M['hour']
#Transformando o indice de string pra datetime
LIDAR_MAI_10M.index=pd.to_datetime(LIDAR_MAI_10M.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME
LIDAR_MAI_1H=LIDAR_MAI_10M.resample('1H').last() #reamostrando para 1H

#ARQUIVO JUNHO - CASO DO DIA 16 A 20 DE JUNHO DE 2020
#DADOS LIDAR PRA-1

LIDAR_JUN_10M = pd.read_csv(pathname+'DADOS_TODOS/LIDAR/LiDAR_PRA01_16-20_jun_2020.csv',sep='\t')
#colocando a data no indice
LIDAR_JUN_10M.index=LIDAR_JUN_10M['Date']+' '+LIDAR_JUN_10M['hour']
#Transformando o indice de string pra datetime
LIDAR_JUN_10M.index=pd.to_datetime(LIDAR_JUN_10M.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME
LIDAR_JUN_1H=LIDAR_JUN_10M.resample('1H').last() #reamostrando para 1H

#ARQUIVO JULHO - CASO DO DIA 13 A 20 DE JULHO DE 2020
#DADOS LIDAR PRA-1

LIDAR_JUL_10M = pd.read_csv(pathname+'DADOS_TODOS/LIDAR/LiDAR_PRA01_13-20_jul_2020.csv',sep='\t')
#colocando a data no indice
LIDAR_JUL_10M.index=LIDAR_JUL_10M['Date']+' '+LIDAR_JUL_10M['hour']
#Transformando o indice de string pra datetime
LIDAR_JUL_10M.index=pd.to_datetime(LIDAR_JUL_10M.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME
LIDAR_JUL_1H=LIDAR_JUL_10M.resample('1H').last() #reamostrando para 1H

#ARQUIVO SETEMBRO - CASO DO DIA 02 A 05 DE SETEMBRO DE 2020
#DADOS LIDAR PRA-1

LIDAR_SET_10M = pd.read_csv(pathname+'DADOS_TODOS/LIDAR/LiDAR_PRA01_02-05_set_2020.csv',sep='\t')
#colocando a data no indice
LIDAR_SET_10M.index=LIDAR_SET_10M['Date']+' '+LIDAR_SET_10M['hour']
#Transformando o indice de string pra datetime
LIDAR_SET_10M.index=pd.to_datetime(LIDAR_SET_10M.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME
LIDAR_SET_1H=LIDAR_SET_10M.resample('1H').last() #reamostrando para 1H

#ARQUIVO OUTUBRO - CASO DO DIA 20 23 DE OUTUBRO DE 2020
#DADOS LIDAR PRA-1

LIDAR_OUT_10M = pd.read_csv(pathname+'DADOS_TODOS/LIDAR/LiDAR_PRA01_20-23_out_2020.csv',sep='\t')
#colocando a data no indice
LIDAR_OUT_10M.index=LIDAR_OUT_10M['Date']+' '+LIDAR_OUT_10M['hour']
#Transformando o indice de string pra datetime
LIDAR_OUT_10M.index=pd.to_datetime(LIDAR_OUT_10M.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME
LIDAR_OUT_1H=LIDAR_OUT_10M.resample('1H').last() #reamostrando para 1H


##DADOS DO ERA5
###Abrindo dados do ERA5
ds = xr.open_mfdataset(pathname+'DADOS_TODOS/ERA5/ERA5_wind*.nc') #abrindo os 2 arquivos nc
ds1 = xr.open_mfdataset(pathname+'DADOS_TODOS/ERA5/ERA5_wave*.nc') #abrindo os 2 arquivos nc
ds2 = xr.open_mfdataset(pathname+'DADOS_TODOS/ERA5/ERA5_TSM*.nc') #abrindo os 2 arquivos nc


#Transformando lat e lon em matriz
lons,lats = np.meshgrid(ds.longitude,ds.latitude)
lons1,lats1 = np.meshgrid(ds1.longitude,ds1.latitude)
lons2,lats2 = np.meshgrid(ds2.longitude,ds2.latitude)

##PONTO PRA-1
ERA5_wind=ds.sel(latitude=-22.134251,longitude=-40.412515,method='nearest').to_dataframe()
ERA5_wave=ds1.sel(latitude=-22.134251,longitude=-40.412515,method='nearest').to_dataframe()
ERA5_T=ds2.sel(latitude=-22.134251,longitude=-40.412515,method='nearest').to_dataframe()


ERA5_wind=ERA5_wind[(ERA5_wind.index>='2020-01-01 00:00') & (ERA5_wind.index<='2020-12-31 23:00')] 
ERA5_wave=ERA5_wave[(ERA5_wave.index>='2020-01-01 00:00') & (ERA5_wave.index<='2020-12-31 23:00')]
ERA5_T=ERA5_T[(ERA5_T.index>='2020-01-01 00:00') & (ERA5_T.index<='2020-12-31 23:00')] 
 
DADOS_1H=ERA5_wind
DADOS_1H['mag10_ERA5']=np.nan
DADOS_1H['mag10_ERA5']=((DADOS_1H['u10']**2)+(DADOS_1H['v10']**2))**0.5
DADOS_1H['dir10_ERA5']=np.nan
DADOS_1H['dir10_ERA5']=np.degrees(np.arctan2(-DADOS_1H['u10'],-DADOS_1H['v10']))
DADOS_1H['dir10_ERA5'][DADOS_1H['dir10_ERA5']<0]=DADOS_1H['dir10_ERA5']+360
DADOS_1H['Hs_ERA5']=np.nan
DADOS_1H['Hs_ERA5']=ERA5_wave['swh']
DADOS_1H['Tp_ERA5']=np.nan
DADOS_1H['Tp_ERA5']=ERA5_wave['pp1d']
DADOS_1H['T2m_ERA5']=np.nan
DADOS_1H['T2m_ERA5']=ERA5_T['t2m']
DADOS_1H['TSM_ERA5']=np.nan
DADOS_1H['TSM_ERA5']=ERA5_T['sst']
DADOS_1H=DADOS_1H.drop(columns=['longitude'])
DADOS_1H=DADOS_1H.drop(columns=['latitude'])


####DADOS DE BOIAS - P25 E PRA-1

##CASO - 2 a 5 FEVEREIRO
#PRA01
PRA01_FEV = pd.read_csv(pathname+'/DADOS_TODOS/boias_P25_PRA01/02_02a05_final_PRA01.csv',sep=',')
#colocando a data no indice
PRA01_FEV.index=PRA01_FEV['Data']+' '+PRA01_FEV['Hora']
#Transformando o indice de string pra datetime
PRA01_FEV.index=pd.to_datetime(PRA01_FEV.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME
PRA01_FEV_10M=PRA01_FEV.resample('10min').last() #reamostrando para 1H
PRA01_FEV_1H=PRA01_FEV.resample('1H').last() #reamostrando para 1H

#P25
P25_FEV = pd.read_csv(pathname+'/DADOS_TODOS/boias_P25_PRA01/02_02a05_TEMP_final_P25.csv',sep=',')
#colocando a data no indice
P25_FEV.index=P25_FEV['Data']+' '+P25_FEV['Hora']
#Transformando o indice de string pra datetime
P25_FEV.index=pd.to_datetime(P25_FEV.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME
P25_FEV_10M=P25_FEV.resample('10min').last() #reamostrando para 1H
P25_FEV_1H=P25_FEV.resample('1H').last() #reamostrando para 1H


##CASO - 1 a 5 MARCO
#PRA01
PRA01_MAR = pd.read_csv(pathname+'/DADOS_TODOS/boias_P25_PRA01/03_01a05_final_PRA01.csv',sep=',')
#colocando a data no indice
PRA01_MAR.index=PRA01_MAR['Data']+' '+PRA01_MAR['Hora']
#Transformando o indice de string pra datetime
PRA01_MAR.index=pd.to_datetime(PRA01_MAR.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME
PRA01_MAR_10M=PRA01_MAR.resample('10min').last() #reamostrando para 1H
PRA01_MAR_1H=PRA01_MAR.resample('1H').last() #reamostrando para 1H

#P25
P25_MAR = pd.read_csv(pathname+'/DADOS_TODOS/boias_P25_PRA01/03_01a05_TEMP_final_P25.csv',sep=',')
#colocando a data no indice
P25_MAR.index=P25_MAR['Data']+' '+P25_MAR['Hora']
#Transformando o indice de string pra datetime
P25_MAR.index=pd.to_datetime(P25_MAR.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME
P25_MAR_10M=P25_MAR.resample('10min').last() #reamostrando para 1H
P25_MAR_1H=P25_MAR.resample('1H').last() #reamostrando para 1H



##CASO - 20 a 21 MAIO
#PRA01
PRA01_MAI = pd.read_csv(pathname+'/DADOS_TODOS/boias_P25_PRA01/05_20a21_final_PRA01.csv',sep=',')
#colocando a data no indice
PRA01_MAI.index=PRA01_MAI['Data']+' '+PRA01_MAI['Hora']
#Transformando o indice de string pra datetime
PRA01_MAI.index=pd.to_datetime(PRA01_MAI.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME
PRA01_MAI_10M=PRA01_MAI.resample('10min').last() #reamostrando para 1H
PRA01_MAI_1H=PRA01_MAI.resample('1H').last() #reamostrando para 1H


#P25
P25_MAI = pd.read_csv(pathname+'/DADOS_TODOS/boias_P25_PRA01/05_20a21_TEMP_final_P25.csv',sep=',')
#colocando a data no indice
P25_MAI.index=P25_MAI['Data']+' '+P25_MAI['Hora']
#Transformando o indice de string pra datetime
P25_MAI.index=pd.to_datetime(P25_MAI.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME
P25_MAI_10M=P25_MAI.resample('10min').last() #reamostrando para 1H
P25_MAI_1H=P25_MAI.resample('1H').last() #reamostrando para 1H



##CASO - 16 a 20 JUNHO
#PRA01
PRA01_JUN = pd.read_csv(pathname+'/DADOS_TODOS/boias_P25_PRA01/06_16a20_final_PRA01.csv',sep=',')
#colocando a data no indice
PRA01_JUN.index=PRA01_JUN['Data']+' '+PRA01_JUN['Hora']
#Transformando o indice de string pra datetime
PRA01_JUN.index=pd.to_datetime(PRA01_JUN.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME
PRA01_JUN_10M=PRA01_JUN.resample('10min').last() #reamostrando para 1H
PRA01_JUN_1H=PRA01_JUN.resample('1H').last() #reamostrando para 1H


#P25
P25_JUN = pd.read_csv(pathname+'/DADOS_TODOS/boias_P25_PRA01/06_16a20_TEMP_final_P25.csv',sep=',')
#colocando a data no indice
P25_JUN.index=P25_JUN['Data']+' '+P25_JUN['Hora']
#Transformando o indice de string pra datetime
P25_JUN.index=pd.to_datetime(P25_JUN.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME
P25_JUN_10M=P25_JUN.resample('10min').last() #reamostrando para 1H
P25_JUN_1H=P25_JUN.resample('1H').last() #reamostrando para 1H



##CASO - 13 a 20 JULHO
#PRA01
PRA01_JUL = pd.read_csv(pathname+'/DADOS_TODOS/boias_P25_PRA01/07_13a20_final_PRA01.csv',sep=',')
#colocando a data no indice
PRA01_JUL.index=PRA01_JUL['Data']+' '+PRA01_JUL['Hora']
#Transformando o indice de string pra datetime
PRA01_JUL.index=pd.to_datetime(PRA01_JUL.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME
PRA01_JUL_10M=PRA01_JUL.resample('10min').last() #reamostrando para 1H
PRA01_JUL_1H=PRA01_JUL.resample('1H').last() #reamostrando para 1H

#P25
P25_JUL = pd.read_csv(pathname+'/DADOS_TODOS/boias_P25_PRA01/07_13a20_TEMP_final_P25.csv',sep=',')
#colocando a data no indice
P25_JUL.index=P25_JUL['Data']+' '+P25_JUL['Hora']
#Transformando o indice de string pra datetime
P25_JUL.index=pd.to_datetime(P25_JUL.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME
P25_JUL_10M=P25_JUL.resample('10min').last() #reamostrando para 1H
P25_JUL_1H=P25_JUL.resample('1H').last() #reamostrando para 1H



##CASO - 02 a 05 SETEMBRO
#PRA01
PRA01_SET = pd.read_csv(pathname+'/DADOS_TODOS/boias_P25_PRA01/09_02a05_final_PRA01.csv',sep=',')
#colocando a data no indice
PRA01_SET.index=PRA01_SET['Data']+' '+PRA01_SET['Hora']
#Transformando o indice de string pra datetime
PRA01_SET.index=pd.to_datetime(PRA01_SET.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME
PRA01_SET_10M=PRA01_SET.resample('10min').last() #reamostrando para 1H
PRA01_SET_1H=PRA01_SET.resample('1H').last() #reamostrando para 1H


#P25
P25_SET = pd.read_csv(pathname+'/DADOS_TODOS/boias_P25_PRA01/09_02a05_TEMP_final_P25.csv',sep=',')
#colocando a data no indice
P25_SET.index=P25_SET['Data']+' '+P25_SET['Hora']
#Transformando o indice de string pra datetime
P25_SET.index=pd.to_datetime(P25_SET.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME
P25_SET_10M=P25_SET.resample('10min').last() #reamostrando para 1H
P25_SET_1H=P25_SET.resample('1H').last() #reamostrando para 1H



##CASO - 20 a 23 OUTUBRO
#PRA01
PRA01_OUT = pd.read_csv(pathname+'/DADOS_TODOS/boias_P25_PRA01/10_20a23_final_PRA01.csv',sep=',')
#colocando a data no indice
PRA01_OUT.index=PRA01_OUT['Data']+' '+PRA01_OUT['Hora']
#Transformando o indice de string pra datetime
PRA01_OUT.index=pd.to_datetime(PRA01_OUT.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME
PRA01_OUT_10M=PRA01_OUT.resample('10min').last() #reamostrando para 1H
PRA01_OUT_1H=PRA01_OUT.resample('1H').last() #reamostrando para 1H


#P25
P25_OUT = pd.read_csv(pathname+'/DADOS_TODOS/boias_P25_PRA01/10_20a23_TEMP_final_P25.csv',sep=',')
#colocando a data no indice
P25_OUT.index=P25_OUT['Data']+' '+P25_OUT['Hora']
#Transformando o indice de string pra datetime
P25_OUT.index=pd.to_datetime(P25_OUT.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME
P25_OUT_10M=P25_OUT.resample('10min').last() #reamostrando para 1H
P25_OUT_1H=P25_OUT.resample('1H').last() #reamostrando para 1H



###ABRINDO P98
P98_HS = pd.read_csv(pathname+'/DADOS_TODOS/P98_DADOS/P98_HS.csv',sep=',')
P98_TP = pd.read_csv(pathname+'/DADOS_TODOS/P98_DADOS/P98_TP.csv',sep=',')
P98_HS.index=P98_HS['time']
P98_TP.index=P98_TP['time']
P98_HS.index=pd.to_datetime(P98_HS.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME
P98_TP.index=pd.to_datetime(P98_TP.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME

###WAVERYS

##PONTO PRA-1
ds_waverys = xr.open_mfdataset(pathname+'DADOS_TODOS/waverys/waverys_PRA-1.nc') #abrindo os 2 arquivos nc

#Transformando lat e lon em matriz
lons3,lats3 = np.meshgrid(ds_waverys.longitude,ds_waverys.latitude)
WAVERYS=ds_waverys.sel(latitude=-22.134251,longitude=-40.412515,method='nearest').to_dataframe()
WAVERYS['DATA']=WAVERYS.index

#COLOCANDO OS DADOS JUNTOS
#1H - Waverys
#DADOS_1H['Hs_WAVERYS']=-9999.00
#DADOS_1H['Tp_WAVERYS']=-9999.00

#PERIODO DE PICO
#i=0
#j=0
#while j<len(WAVERYS):
#    while i< len(DADOS_1H):
#        if WAVERYS.index[j]==DADOS_1H.index[i]:
#            DADOS_1H['Tp_WAVERYS'][i]=WAVERYS['wav_tp'][j]
#        else:
#            print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#        i=i+1
#    i=0
#    j=j+1

#ALTURA SIGNIFICATIVA DE ONDA
#i=0
#j=0
#while j<len(WAVERYS):
#    while i< len(DADOS_1H):
#        if WAVERYS.index[j]==DADOS_1H.index[i]:
#            DADOS_1H['Hs_WAVERYS'][i]=WAVERYS['wav_hs'][j]
#        else:
#            print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#        i=i+1
#    i=0
#    j=j+1

#SALVANDO EM CSV - TROCAR VERSAO
#DADOS_1H.to_csv(pathname+'DADOS_1H_v1.csv')

#continuar daqui

#MAGNITUDE E DIRECAO DO VENTO
#i=0
#j=0
#DADOS_1H['mag10m_WAVERYS']=-9999.00

#while j<len(WAVERYS):
#    while i< len(DADOS_1H):
#        if WAVERYS.index[j]==DADOS_1H.index[i]:
#            DADOS_1H['mag10m_WAVERYS'][i]=WAVERYS['atm_wnd_spd_10m'][j]
#        else:
#            print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
#        
#        i=i+1
#    i=0
#    j=j+1

#SALVANDO EM CSV - TROCAR VERSAO
#DADOS_1H.to_csv(pathname+'DADOS_1H_v2.csv')


#i=0
#j=0
#DADOS_1H['dir10m_WAVERYS']=-9999.00
#
#while j<len(WAVERYS):
#    while i< len(DADOS_1H):
#        if WAVERYS.index[j]==DADOS_1H.index[i]:
#            DADOS_1H['dir10m_WAVERYS'][i]=WAVERYS['atm_wnd_dir_10m'][j]
#        else:
#            print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
#        
#        i=i+1
#    i=0
#    j=j+1



#SALVANDO EM CSV - TROCAR VERSAO
#DADOS_1H.to_csv(pathname+'DADOS_1H_v3.csv')

#1H - p98 - variaveis: mag10, dir10, Hs, Tp, T, TSM
#DADOS_1H['mag10_P98']=-9999.00
#DADOS_1H['dir10_P98']=-9999.00

#MAGNITUDE P98
#i=0
#j=0
#while j<len(P98_HS):
#    while i< len(DADOS_1H):
#        if P98_HS.index[j]==DADOS_1H.index[i]:
#            DADOS_1H['mag10_P98'][i]=P98_HS['mag10'][j]
#        else:
#            print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
#        
#        i=i+1
#    i=0
#    j=j+1

#SALVANDO EM CSV - TROCAR VERSAO
#DADOS_1H.to_csv(pathname+'DADOS_1H_v4.csv')


#DIRECAO P98
#i=0
#j=0
#while j<len(P98_HS):
#    while i< len(DADOS_1H):
#        if P98_HS.index[j]==DADOS_1H.index[i]:
#            DADOS_1H['dir10_P98'][i]=P98_HS['dir10'][j]
#        else:
#            DADOS_1H['dir10_P98'][i]=-9999.00
#            print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
#        
#        i=i+1
#    i=0
#    j=j+1

#SALVANDO EM CSV - TROCAR VERSAO
#DADOS_1H.to_csv(pathname+'DADOS_1H_v5.csv')

#HSP98
#DADOS_1H['Hs_P98']=-9999.00
#i=0
#j=0
#while j<len(P98_HS):
#    while i< len(DADOS_1H):
#        if P98_HS.index[j]==DADOS_1H.index[i]:
#            DADOS_1H['Hs_P98'][i]=P98_HS['Hs'][j]
#        else:
#            print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
#        
#        i=i+1
#    i=0
#    j=j+1

#SALVANDO EM CSV - TROCAR VERSAO
#DADOS_1H.to_csv(pathname+'DADOS_1H_v6.csv')


#TPP98
#DADOS_1H['Tp_P98']=-9999.00
#i=0
#j=0
#while j<len(P98_HS):
#    while i< len(DADOS_1H):
#        if P98_HS.index[j]==DADOS_1H.index[i]:
#            DADOS_1H['Tp_P98'][i]=P98_HS['Tp'][j]
#        else:
#            print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#        i=i+1
#    i=0
#    j=j+1

#SALVANDO EM CSV - TROCAR VERSAO
#DADOS_1H.to_csv(pathname+'DADOS_1H_v7.csv')


#HSP98
#DADOS_1H['T_P98']=-9999.00
#i=0
#j=0
#while j<len(P98_HS):
#    while i< len(DADOS_1H):
#        if P98_HS.index[j]==DADOS_1H.index[i]:
#            DADOS_1H['T_P98'][i]=P98_HS['T'][j]
#        else:
#            print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
#        
#        i=i+1
#    i=0
#    j=j+1

#SALVANDO EM CSV - TROCAR VERSAO
#DADOS_1H.to_csv(pathname+'DADOS_1H_v8.csv')

#TSMP98
#DADOS_1H['TSM_P98']=-9999.00
#i=0
#j=0
#while j<len(P98_HS):
#    while i< len(DADOS_1H):
#        if P98_HS.index[j]==DADOS_1H.index[i]:
#            DADOS_1H['TSM_P98'][i]=P98_HS['TSM'][j]
#        else:
#            print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
#        
#        i=i+1
#    i=0
#    j=j+1

#SALVANDO EM CSV - TROCAR VERSAO
#DADOS_1H.to_csv(pathname+'DADOS_1H_v9.csv')



#TLIDAR
#DADOS_1H['T_LIDAR']=-9999.00
#i=0
#j=0
#while j<len(LIDAR_1H):
#    while i< len(DADOS_1H):
#        if LIDAR_1H.index[j]==DADOS_1H.index[i]:
#            DADOS_1H['T_LIDAR'][i]=LIDAR_1H['T'][j]
#        else:
#            print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
#        
#        i=i+1
#    i=0
#    j=j+1

#SALVANDO EM CSV - TROCAR VERSAO
#DADOS_1H.to_csv(pathname+'DADOS_1H_v10.csv')


#TLIDAR
#DADOS_1H['P_LIDAR']=-9999.00
#i=0
#j=0
#while j<len(LIDAR_1H):
#    while i< len(DADOS_1H):
#        if LIDAR_1H.index[j]==DADOS_1H.index[i]:
#            DADOS_1H['P_LIDAR'][i]=LIDAR_1H['P'][j]
#        else:
#            print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
#        
#        i=i+1
#    i=0
#    j=j+1

#SALVANDO EM CSV - TROCAR VERSAO
#DADOS_1H.to_csv(pathname+'DADOS_1H_v11.csv')


#TLIDAR
#DADOS_1H['UR_LIDAR']=-9999.00
#i=0
#j=0
#while j<len(LIDAR_1H):
#    while i< len(DADOS_1H):
#        if LIDAR_1H.index[j]==DADOS_1H.index[i]:
#            DADOS_1H['UR_LIDAR'][i]=LIDAR_1H['UR'][j]
#        else:
#            print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
#        
#        i=i+1
#    i=0
#    j=j+1

#SALVANDO EM CSV - TROCAR VERSAO
#DADOS_1H.to_csv(pathname+'DADOS_1H_v12.csv')


#TLIDAR
#DADOS_1H['mag90_LIDAR']=-9999.00
#i=0
#j=0
#while j<len(LIDAR_1H):
#    while i< len(DADOS_1H):
#        if LIDAR_1H.index[j]==DADOS_1H.index[i]:
#            DADOS_1H['mag90_LIDAR'][i]=LIDAR_1H['mag90'][j]
#        else:
#            print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
#        
#        i=i+1
#    i=0
#    j=j+1

#SALVANDO EM CSV - TROCAR VERSAO
#DADOS_1H.to_csv(pathname+'DADOS_1H_v13.csv')


#TLIDAR
#DADOS_1H['dir90_LIDAR']=-9999.00
#i=0
#j=0
#while j<len(LIDAR_1H):
#    while i< len(DADOS_1H):
#        if LIDAR_1H.index[j]==DADOS_1H.index[i]:
#            DADOS_1H['dir90_LIDAR'][i]=LIDAR_1H['dir90'][j]
#        else:
#            print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
#        
#        i=i+1
#    i=0
#    j=j+1

#SALVANDO EM CSV - TROCAR VERSAO
#DADOS_1H.to_csv(pathname+'DADOS_1H_v14.csv')


#TLIDAR
#DADOS_1H['mag100_LIDAR']=-9999.00
#j=0
#while j<len(LIDAR_1H):
#    while i< len(DADOS_1H):
#        if LIDAR_1H.index[j]==DADOS_1H.index[i]:
#            DADOS_1H['mag100_LIDAR'][i]=LIDAR_1H['mag100'][j]
#        else:
#            print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
#        
#        i=i+1
#    i=0
#    j=j+1

#SALVANDO EM CSV - TROCAR VERSAO
#DADOS_1H.to_csv(pathname+'DADOS_1H_v15.csv')


#TLIDAR
#DADOS_1H['dir100_LIDAR']=-9999.00
#i=0
#j=0
#while j<len(LIDAR_1H):
#    while i< len(DADOS_1H):
#        if LIDAR_1H.index[j]==DADOS_1H.index[i]:
#            DADOS_1H['dir100_LIDAR'][i]=LIDAR_1H['dir100'][j]
#        else:
#            print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
#        
#        i=i+1
#    i=0
#    j=j+1

#SALVANDO EM CSV - TROCAR VERSAO
#DADOS_1H.to_csv(pathname+'DADOS_1H_v16.csv')


#TLIDAR
#DADOS_1H['mag110_LIDAR']=-9999.00
#i=0
#j=0
#while j<len(LIDAR_1H):
#    while i< len(DADOS_1H):
#        if LIDAR_1H.index[j]==DADOS_1H.index[i]:
#            DADOS_1H['mag110_LIDAR'][i]=LIDAR_1H['mag110'][j]
#        else:
#            print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
#        
#        i=i+1
#    i=0
#    j=j+1

#SALVANDO EM CSV - TROCAR VERSAO
#DADOS_1H.to_csv(pathname+'DADOS_1H_v17.csv')



#TLIDAR
#DADOS_1H['dir110_LIDAR']=-9999.00
#i=0
#j=0
#while j<len(LIDAR_1H):
#    while i< len(DADOS_1H):
#        if LIDAR_1H.index[j]==DADOS_1H.index[i]:
#            DADOS_1H['dir110_LIDAR'][i]=LIDAR_1H['dir110'][j]
#        else:
#            print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
#        
#        i=i+1
#    i=0
#    j=j+1

#SALVANDO EM CSV - TROCAR VERSAO
#DADOS_1H.to_csv(pathname+'DADOS_1H_v18.csv')




#TLIDAR
#DADOS_1H['mag120_LIDAR']=-9999.00
#i=0
#j=0
#while j<len(LIDAR_1H):
#    while i< len(DADOS_1H):
#        if LIDAR_1H.index[j]==DADOS_1H.index[i]:
#            DADOS_1H['mag120_LIDAR'][i]=LIDAR_1H['mag120'][j]
#        else:
#            print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
#        
#        i=i+1
#    i=0
#    j=j+1

#SALVANDO EM CSV - TROCAR VERSAO
#DADOS_1H.to_csv(pathname+'DADOS_1H_v19.csv')


#TLIDAR
#DADOS_1H['dir120_LIDAR']=-9999.00
#i=0
#j=0
#while j<len(LIDAR_1H):
#    while i< len(DADOS_1H):
#        if LIDAR_1H.index[j]==DADOS_1H.index[i]:
#            DADOS_1H['dir120_LIDAR'][i]=LIDAR_1H['dir120'][j]
#        else:
#            print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
#        
#        i=i+1
#    i=0
#    j=j+1

#SALVANDO EM CSV - TROCAR VERSAO
#DADOS_1H.to_csv(pathname+'DADOS_1H_v20.csv')



#TLIDAR
#DADOS_1H['mag150_LIDAR']=-9999.00
#i=0
#j=0
#while j<len(LIDAR_1H):
#    while i< len(DADOS_1H):
#        if LIDAR_1H.index[j]==DADOS_1H.index[i]:
#            DADOS_1H['mag150_LIDAR'][i]=LIDAR_1H['mag150'][j]
#        else:
#            print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
#        
#        i=i+1
#    i=0
#    j=j+1

#SALVANDO EM CSV - TROCAR VERSAO
#DADOS_1H.to_csv(pathname+'DADOS_1H_v21.csv')


# #TLIDAR
# DADOS_1H['dir150_LIDAR']=-9999.00
# i=0
# j=0
# while j<len(LIDAR_1H):
#     while i< len(DADOS_1H):
#         if LIDAR_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['dir150_LIDAR'][i]=LIDAR_1H['dir150'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v22.csv')



# #TLIDAR
# DADOS_1H['mag170_LIDAR']=-9999.00
# i=0
# j=0
# while j<len(LIDAR_1H):
#     while i< len(DADOS_1H):
#         if LIDAR_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['mag170_LIDAR'][i]=LIDAR_1H['mag170'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v23.csv')


# #TLIDAR
# DADOS_1H['dir170_LIDAR']=-9999.00
# i=0
# j=0
# while j<len(LIDAR_1H):
#     while i< len(DADOS_1H):
#         if LIDAR_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['dir170_LIDAR'][i]=LIDAR_1H['dir170'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v24.csv')



# #TLIDAR
# DADOS_1H['mag200_LIDAR']=-9999.00
# i=0
# j=0
# while j<len(LIDAR_1H):
#     while i< len(DADOS_1H):
#         if LIDAR_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['mag200_LIDAR'][i]=LIDAR_1H['mag200'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v25.csv')


# #TLIDAR
# DADOS_1H['dir200_LIDAR']=-9999.00
# i=0
# j=0
# while j<len(LIDAR_1H):
#     while i< len(DADOS_1H):
#         if LIDAR_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['dir200_LIDAR'][i]=LIDAR_1H['dir200'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v26.csv')

# #PRA-FEV
# DADOS_1H['mag10_PRA01']=-9999.00
# i=0
# j=0
# while j<len(PRA01_FEV_1H):
#     while i< len(DADOS_1H):
#         if PRA01_FEV_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['mag10_PRA01'][i]=PRA01_FEV_1H['mag'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v27.csv')

# #PRA-FEV
# DADOS_1H['dir10_PRA01']=-9999.00
# i=0
# j=0
# while j<len(PRA01_FEV_1H):
#     while i< len(DADOS_1H):
#         if PRA01_FEV_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['dir10_PRA01'][i]=PRA01_FEV_1H['dir'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v28.csv')


# #PRA-FEV
# DADOS_1H['T_PRA01']=-9999.00
# i=0
# j=0
# while j<len(PRA01_FEV_1H):
#     while i< len(DADOS_1H):
#         if PRA01_FEV_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['T_PRA01'][i]=PRA01_FEV_1H['T'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v29.csv')


# #PRA-FEV
# DADOS_1H['Td_PRA01']=-9999.00
# i=0
# j=0
# while j<len(PRA01_FEV_1H):
#     while i< len(DADOS_1H):
#         if PRA01_FEV_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['Td_PRA01'][i]=PRA01_FEV_1H['Td'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v30.csv')


# #PRA-FEV
# DADOS_1H['UR_PRA01']=-9999.00
# i=0
# j=0
# while j<len(PRA01_FEV_1H):
#     while i< len(DADOS_1H):
#         if PRA01_FEV_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['UR_PRA01'][i]=PRA01_FEV_1H['UR'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v31.csv')


# #PRA-MAR
# i=0
# j=0
# while j<len(PRA01_MAR_1H):
#     while i< len(DADOS_1H):
#         if PRA01_MAR_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['mag10_PRA01'][i]=PRA01_MAR_1H['mag'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v32.csv')

# #PRA-MAR
# i=0
# j=0
# while j<len(PRA01_MAR_1H):
#     while i< len(DADOS_1H):
#         if PRA01_MAR_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['dir10_PRA01'][i]=PRA01_MAR_1H['dir'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v33.csv')


# #PRA-MAR
# i=0
# j=0
# while j<len(PRA01_MAR_1H):
#     while i< len(DADOS_1H):
#         if PRA01_MAR_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['T_PRA01'][i]=PRA01_MAR_1H['T'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v34.csv')


# #PRA-MAR
# i=0
# j=0
# while j<len(PRA01_MAR_1H):
#     while i< len(DADOS_1H):
#         if PRA01_MAR_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['Td_PRA01'][i]=PRA01_MAR_1H['Td'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v35.csv')


# #PRA-MAR
# i=0
# j=0
# while j<len(PRA01_MAR_1H):
#     while i< len(DADOS_1H):
#         if PRA01_MAR_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['UR_PRA01'][i]=PRA01_MAR_1H['UR'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v36.csv')


# #PRA-MAI
# i=0
# j=0
# while j<len(PRA01_MAI_1H):
#     while i< len(DADOS_1H):
#         if PRA01_MAI_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['mag10_PRA01'][i]=PRA01_MAI_1H['mag'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v37.csv')

# #PRA-MAI
# i=0
# j=0
# while j<len(PRA01_MAI_1H):
#     while i< len(DADOS_1H):
#         if PRA01_MAI_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['dir10_PRA01'][i]=PRA01_MAI_1H['dir'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v38.csv')


# #PRA-MAI
# i=0
# j=0
# while j<len(PRA01_MAI_1H):
#     while i< len(DADOS_1H):
#         if PRA01_MAI_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['T_PRA01'][i]=PRA01_MAI_1H['T'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v39.csv')


# #PRA-MAI
# i=0
# j=0
# while j<len(PRA01_MAI_1H):
#     while i< len(DADOS_1H):
#         if PRA01_MAI_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['Td_PRA01'][i]=PRA01_MAI_1H['Td'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v40.csv')


# #PRA-MAI
# i=0
# j=0
# while j<len(PRA01_MAI_1H):
#     while i< len(DADOS_1H):
#         if PRA01_MAI_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['UR_PRA01'][i]=PRA01_MAI_1H['UR'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v41.csv')


# #PRA-JUL
# i=0
# j=0
# while j<len(PRA01_JUL_1H):
#     while i< len(DADOS_1H):
#         if PRA01_JUL_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['mag10_PRA01'][i]=PRA01_JUL_1H['mag'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v42.csv')

# #PRA-JUL
# i=0
# j=0
# while j<len(PRA01_JUL_1H):
#     while i< len(DADOS_1H):
#         if PRA01_JUL_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['dir10_PRA01'][i]=PRA01_JUL_1H['dir'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v43.csv')


# #PRA-JUL
# i=0
# j=0
# while j<len(PRA01_JUL_1H):
#     while i< len(DADOS_1H):
#         if PRA01_JUL_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['T_PRA01'][i]=PRA01_JUL_1H['T'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v44.csv')


# #PRA-JUL
# i=0
# j=0
# while j<len(PRA01_JUL_1H):
#     while i< len(DADOS_1H):
#         if PRA01_JUL_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['Td_PRA01'][i]=PRA01_JUL_1H['Td'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v45.csv')


# #PRA-JUL
# i=0
# j=0
# while j<len(PRA01_JUL_1H):
#     while i< len(DADOS_1H):
#         if PRA01_JUL_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['UR_PRA01'][i]=PRA01_JUL_1H['UR'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v46.csv')


# #PRA-SET
# i=0
# j=0
# while j<len(PRA01_SET_1H):
#     while i< len(DADOS_1H):
#         if PRA01_SET_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['mag10_PRA01'][i]=PRA01_SET_1H['mag'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v47.csv')

# #PRA-SET
# i=0
# j=0
# while j<len(PRA01_SET_1H):
#     while i< len(DADOS_1H):
#         if PRA01_SET_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['dir10_PRA01'][i]=PRA01_SET_1H['dir'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v48.csv')


# #PRA-SET
# i=0
# j=0
# while j<len(PRA01_SET_1H):
#     while i< len(DADOS_1H):
#         if PRA01_SET_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['T_PRA01'][i]=PRA01_SET_1H['T'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v49.csv')


# #PRA-SET
# i=0
# j=0
# while j<len(PRA01_SET_1H):
#     while i< len(DADOS_1H):
#         if PRA01_SET_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['Td_PRA01'][i]=PRA01_SET_1H['Td'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v50.csv')


# #PRA-SET
# i=0
# j=0
# while j<len(PRA01_SET_1H):
#     while i< len(DADOS_1H):
#         if PRA01_SET_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['UR_PRA01'][i]=PRA01_SET_1H['UR'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v51.csv')


# #PRA-OUT
# i=0
# j=0
# while j<len(PRA01_OUT_1H):
#     while i< len(DADOS_1H):
#         if PRA01_OUT_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['mag10_PRA01'][i]=PRA01_OUT_1H['mag'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v52.csv')

# #PRA-OUT
# i=0
# j=0
# while j<len(PRA01_OUT_1H):
#     while i< len(DADOS_1H):
#         if PRA01_OUT_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['dir10_PRA01'][i]=PRA01_OUT_1H['dir'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v53.csv')


# #PRA-OUT
# i=0
# j=0
# while j<len(PRA01_OUT_1H):
#     while i< len(DADOS_1H):
#         if PRA01_OUT_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['T_PRA01'][i]=PRA01_OUT_1H['T'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v54.csv')


# #PRA-OUT
# i=0
# j=0
# while j<len(PRA01_OUT_1H):
#     while i< len(DADOS_1H):
#         if PRA01_OUT_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['Td_PRA01'][i]=PRA01_OUT_1H['Td'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

# #SALVANDO EM CSV - TROCAR VERSAO
# DADOS_1H.to_csv(pathname+'DADOS_1H_v55.csv')


# #PRA-OUT
# i=0
# j=0
# while j<len(PRA01_OUT_1H):
#     while i< len(DADOS_1H):
#         if PRA01_OUT_1H.index[j]==DADOS_1H.index[i]:
#             DADOS_1H['UR_PRA01'][i]=PRA01_OUT_1H['UR'][j]
#         else:
#             print("SCRIPT RODANDO - matriz ij - j: {} e i: {}".format(j,i))
        
#         i=i+1
#     i=0
#     j=j+1

#SALVANDO EM CSV - TROCAR VERSAO
DADOS_1H.to_csv(pathname+'DADOS_1H_v56.csv')

DADOS_1H = pd.read_csv(pathname+'DADOS_1H_v55.csv',sep=',')
DADOS_1H.index=DADOS_1H['time']
DADOS_1H.index=pd.to_datetime(DADOS_1H.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME
DADOS_1H_TESTE = DADOS_1H.replace(-9999.0, np.nan)
DADOS_1H_TESTE.to_csv(pathname+'DADOS_1H_vFINAL.csv')
DADOS_1H=DADOS_1H_TESTE

#CONTINUAR DAQUI
#DADOS_1H = pd.read_csv(pathname+'DADOS_1H_v56.csv',sep=',')
#DADOS_1H.index=DADOS_1H['time']
#DADOS_1H.index=pd.to_datetime(DADOS_1H.index, format="%Y-%m-%d %H:%M:%S") #CONVERTENDO PRA DATETIME


#RUGOSIDADES TABELAS DA DNV
zo_tab_menor=0.0001  #tabelaDNV - open sea without waves
zo_tab_maior=0.001  #tabelaDNV - Open sea with waves

##Funcoes para calcular os parametros de friccao e beta
def beta(h,zo):
    B=np.log(h/zo)
    return B

def ufric(h,uh,zo):
    ufric=(k*uh)/beta(h,zo)
    return ufric

#TODOS OS CALCULOS USAM WAVERYS, ERA5 E PRA01
#METODO 1 - DNV TABELADO
#METODO 2 - LEI DA POTENCIA
#METODO 3 - DONELAN 90 - NEUTRO
#METODO 4 - DONELAN 90 - ESTAVEL
#METODO 5 - DONELAN 93 - NEUTRO
#METODO 6 - DONELAN 93 - ESTAVEL
#METODO 7 - TAYLOR & YELLAND 01 - CARMO et al. 22
#METODO 8 - TAYLOR & YELLAND 01 - CARMO et al. 22 - ESTAVEL
#METODO 9 - CARMO ET AL - ADAPTADO PRA REGIAO
#METODO 10 - GARRETT 1980 ???


#metodo1- ERA5 - Z0 min
for i in niveis:
    DADOS_1H['mag10_ERA5_MET1']=DADOS_1H['mag10_ERA5']
    DADOS_1H['mag{}_ERA5_MET1'.format(i)]=(ufric(10,DADOS_1H['mag10_ERA5'],zo_tab_maior)/k)*(beta(i,zo_tab_maior))

#metodo1- WAVERYS - Z0 min
for i in niveis:
    DADOS_1H['mag10_c']=DADOS_1H['mag10m_WAVERYS']
    DADOS_1H['mag{}_WAVERYS_MET1'.format(i)]=(ufric(10,DADOS_1H['mag10m_WAVERYS'],zo_tab_maior)/k)*(beta(i,zo_tab_maior))

#metodo1- PRA01 - Z0 min
for i in niveis:
    DADOS_1H['mag10_PRA01_MET1']=DADOS_1H['mag10_PRA01']
    DADOS_1H['mag{}_PRA01_MET1'.format(i)]=(ufric(10,DADOS_1H['mag10_PRA01'],zo_tab_maior)/k)*(beta(i,zo_tab_maior))


#metodo2 - LEI DA POTENCIA - ERA5
for i in niveis:
    DADOS_1H['mag10_ERA5_MET2']=DADOS_1H['mag10_ERA5']
    DADOS_1H['mag{}_ERA5_MET2'.format(i)]=DADOS_1H['mag10_ERA5']*((i/10)**(alfa_oc))


#metodo2 - LEI DA POTENCIA - WAVERYS
for i in niveis:
    DADOS_1H['mag10_WAVERYS_MET2']=DADOS_1H['mag10m_WAVERYS']
    DADOS_1H['mag{}_WAVERYS_MET2'.format(i)]=DADOS_1H['mag10m_WAVERYS']*((i/10)**(alfa_oc))


#metodo2 - LEI DA POTENCIA - ERA5
for i in niveis:
    DADOS_1H['mag10_PRA01_MET2']=DADOS_1H['mag10_PRA01']
    DADOS_1H['mag{}_PRA01_MET2'.format(i)]=DADOS_1H['mag10_PRA01']*((i/10)**(alfa_oc))


#metodo3 - Zo Donelan 90 - Neutro 
DADOS_1H['zo_Do_90_ERA5']=0.033*(DADOS_1H['Hs_ERA5']/4)
DADOS_1H['zo_Do_90_WAVERYS']=0.033*(DADOS_1H['Hs_WAVERYS']/4)
DADOS_1H['zo_Do_90_P98']=0.033*(DADOS_1H['Hs_P98']/4)

for i in niveis:
    DADOS_1H['mag10_ERA5_MET3']=DADOS_1H['mag10_ERA5']
    DADOS_1H['mag{}_ERA5_MET3'.format(i)]=(ufric(10,DADOS_1H['mag10_ERA5'],DADOS_1H['zo_Do_90_ERA5'])/k)*(beta(i,DADOS_1H['zo_Do_90_ERA5']))

for i in niveis:
    DADOS_1H['mag10_WAVERYS_MET3']=DADOS_1H['mag10m_WAVERYS']
    DADOS_1H['mag{}_WAVERYS_MET3'.format(i)]=(ufric(10,DADOS_1H['mag10m_WAVERYS'],DADOS_1H['zo_Do_90_WAVERYS'])/k)*(beta(i,DADOS_1H['zo_Do_90_WAVERYS']))

for i in niveis:
    DADOS_1H['mag10_PRA01_MET3']=DADOS_1H['mag10_P98']
    DADOS_1H['mag{}_PRA01_MET3'.format(i)]=(ufric(10,DADOS_1H['mag10_P98'],DADOS_1H['zo_Do_90_P98'])/k)*(beta(i,DADOS_1H['zo_Do_90_P98']))


#metodo4 - Zo Donelan 90 - Estável 
for i in niveis:
    DADOS_1H['mag10_ERA5_MET4']=DADOS_1H['mag10_ERA5']
    DADOS_1H['mag{}_ERA5_MET4'.format(i)]=(ufric(10,DADOS_1H['mag10_ERA5'],DADOS_1H['zo_Do_90_ERA5'])/k)*((beta(i,DADOS_1H['zo_Do_90_ERA5']))+((4.7*(i))/80))

for i in niveis:
    DADOS_1H['mag10_WAVERYS_MET4']=DADOS_1H['mag10m_WAVERYS']
    DADOS_1H['mag{}_WAVERYS_MET4'.format(i)]=(ufric(10,DADOS_1H['mag10m_WAVERYS'],DADOS_1H['zo_Do_90_WAVERYS'])/k)*((beta(i,DADOS_1H['zo_Do_90_WAVERYS']))+((4.7*(i))/80))

for i in niveis:
    DADOS_1H['mag10_PRA01_MET4']=DADOS_1H['mag10_P98']
    DADOS_1H['mag{}_PRA01_MET4'.format(i)]=(ufric(10,DADOS_1H['mag10_P98'],DADOS_1H['zo_Do_90_P98'])/k)*((beta(i,DADOS_1H['zo_Do_90_P98']))+((4.7*(i))/80))


#metodo3 - Zo Donelan 90 - Neutro 

#metodo4 - Zo Donelan 93 - Neutro
DADOS_1H['zo_Do_93_ERA5']=(DADOS_1H['Hs_ERA5']/4)*(6.7*(10**-4))*((DADOS_1H['mag10_ERA5']/Cp)**2.6)
DADOS_1H['zo_Do_93_P98']=(DADOS_1H['Hs_P98']/4)*(6.7*(10**-4))*((DADOS_1H['mag10_P98']/Cp)**2.6)
DADOS_1H['zo_Do_93_WAVERYS']=(DADOS_1H['Hs_WAVERYS']/4)*(6.7*(10**-4))*((DADOS_1H['mag10m_WAVERYS']/Cp)**2.6)


for i in niveis:
    DADOS_1H['mag10_ERA5_MET5']=DADOS_1H['mag10_ERA5']
    DADOS_1H['mag{}_ERA5_MET5'.format(i)]=(ufric(10,DADOS_1H['mag10_ERA5'],DADOS_1H['zo_Do_93_ERA5'])/k)*(beta(i,DADOS_1H['zo_Do_93_ERA5']))

for i in niveis:
    DADOS_1H['mag10_WAVERYS_MET5']=DADOS_1H['mag10m_WAVERYS']
    DADOS_1H['mag{}_WAVERYS_MET5'.format(i)]=(ufric(10,DADOS_1H['mag10m_WAVERYS'],DADOS_1H['zo_Do_93_WAVERYS'])/k)*(beta(i,DADOS_1H['zo_Do_93_WAVERYS']))

for i in niveis:
    DADOS_1H['mag10_PRA01_MET5']=DADOS_1H['mag10_P98']
    DADOS_1H['mag{}_PRA01_MET5'.format(i)]=(ufric(10,DADOS_1H['mag10_P98'],DADOS_1H['zo_Do_93_P98'])/k)*(beta(i,DADOS_1H['zo_Do_93_P98']))


#metodo4 - Zo Donelan 90 - Estável 
for i in niveis:
    DADOS_1H['mag10_ERA5_MET6']=DADOS_1H['mag10_ERA5']
    DADOS_1H['mag{}_ERA5_MET6'.format(i)]=(ufric(10,DADOS_1H['mag10_ERA5'],DADOS_1H['zo_Do_93_ERA5'])/k)*((beta(i,DADOS_1H['zo_Do_93_ERA5']))+((4.7*(i))/80))

for i in niveis:
    DADOS_1H['mag10_WAVERYS_MET6']=DADOS_1H['mag10m_WAVERYS']
    DADOS_1H['mag{}_WAVERYS_MET6'.format(i)]=(ufric(10,DADOS_1H['mag10m_WAVERYS'],DADOS_1H['zo_Do_93_WAVERYS'])/k)*((beta(i,DADOS_1H['zo_Do_93_WAVERYS']))+((4.7*(i))/80))

for i in niveis:
    DADOS_1H['mag10_PRA01_MET6']=DADOS_1H['mag10_P98']
    DADOS_1H['mag{}_PRA01_MET6'.format(i)]=(ufric(10,DADOS_1H['mag10_P98'],DADOS_1H['zo_Do_93_P98'])/k)*((beta(i,DADOS_1H['zo_Do_93_P98']))+((4.7*(i))/80))



#metodo7 - Zo Taylor and Yelland 01 - Neutro 
DADOS_1H['Lp_ERA5']=(DADOS_1H['Tp_ERA5']**2*g)/(2*np.pi)
DADOS_1H['Lp_WAVERYS']=(DADOS_1H['Tp_WAVERYS']**2*g)/(2*np.pi)
DADOS_1H['Lp_P98']=(DADOS_1H['Tp_P98']**2*g)/(2*np.pi)
DADOS_1H['ustar_6_ERA5']=(k*DADOS_1H['mag10_ERA5'])/np.log(10/(mi*DADOS_1H['Hs_ERA5']*((DADOS_1H['Hs_ERA5']/DADOS_1H['Lp_ERA5'])**4.5)))
DADOS_1H['ustar_6_WAVERYS']=(k*DADOS_1H['mag10m_WAVERYS'])/np.log(10/(mi*DADOS_1H['Hs_WAVERYS']*((DADOS_1H['Hs_WAVERYS']/DADOS_1H['Lp_WAVERYS'])**4.5)))
DADOS_1H['ustar_6_P98']=(k*DADOS_1H['mag10_P98'])/np.log(10/(mi*DADOS_1H['Hs_P98']*((DADOS_1H['Hs_P98']/DADOS_1H['Lp_P98'])**4.5)))
DADOS_1H['zo_TaYe_ERA5']=(mi*DADOS_1H['Hs_ERA5']*((DADOS_1H['Hs_ERA5']/DADOS_1H['Lp_ERA5'])**4.5))
DADOS_1H['zo_TaYe_WAVERYS']=(mi*DADOS_1H['Hs_WAVERYS']*((DADOS_1H['Hs_WAVERYS']/DADOS_1H['Lp_WAVERYS'])**4.5))
DADOS_1H['zo_TaYe_P98']=(mi*DADOS_1H['Hs_P98']*((DADOS_1H['Hs_P98']/DADOS_1H['Lp_P98'])**4.5))

for i in niveis:
    DADOS_1H['mag10_ERA5_MET7']=DADOS_1H['mag10_ERA5']
    DADOS_1H['mag{}_ERA5_MET7'.format(i)]=(DADOS_1H['ustar_6_ERA5']/k)*((np.log((i)/DADOS_1H['zo_TaYe_ERA5'])))

for i in niveis:
    DADOS_1H['mag10_WAVERYS_MET7']=DADOS_1H['mag10m_WAVERYS']
    DADOS_1H['mag{}_WAVERYS_MET7'.format(i)]=(DADOS_1H['ustar_6_WAVERYS']/k)*((np.log((i)/DADOS_1H['zo_TaYe_WAVERYS'])))

for i in niveis:
    DADOS_1H['mag10_PRA01_MET7']=DADOS_1H['mag10_P98']
    DADOS_1H['mag{}_PRA01_MET7'.format(i)]=(DADOS_1H['ustar_6_P98']/k)*((np.log((i)/DADOS_1H['zo_TaYe_P98'])))


#metodo8 - Zo Taylor and Yelland 01 - Estavel 

for i in niveis:
    DADOS_1H['mag10_ERA5_MET8']=DADOS_1H['mag10_ERA5']
    DADOS_1H['mag{}_ERA5_MET8'.format(i)]=(DADOS_1H['ustar_6_ERA5']/k)*((np.log((i)/DADOS_1H['zo_TaYe_ERA5']))+((4.7*(i))/90))

for i in niveis:
    DADOS_1H['mag10_WAVERYS_MET8']=DADOS_1H['mag10m_WAVERYS']
    DADOS_1H['mag{}_WAVERYS_MET8'.format(i)]=(DADOS_1H['ustar_6_WAVERYS']/k)*((np.log((i)/DADOS_1H['zo_TaYe_WAVERYS']))+((4.7*(i))/90))

for i in niveis:
    DADOS_1H['mag10_PRA01_MET8']=DADOS_1H['mag10_P98']
    DADOS_1H['mag{}_PRA01_MET8'.format(i)]=(DADOS_1H['ustar_6_P98']/k)*((np.log((i)/DADOS_1H['zo_TaYe_P98']))+((4.7*(i))/90))


#metodo 9 - Zo CAR 23 - Estavel 

#metodo6A - Taylor Yelland neutro  - CORRIGIDO - TSM e T - ERA5
DADOS_1H['Lp_ERA5']=(DADOS_1H['Tp_ERA5']**2*g)/(2*np.pi)
DADOS_1H['Lp_WAVERYS']=(DADOS_1H['Tp_WAVERYS']**2*g)/(2*np.pi)
DADOS_1H['Lp_P98']=(DADOS_1H['Tp_P98']**2*g)/(2*np.pi)
DADOS_1H['ustar_6_ERA5']=(k*DADOS_1H['mag10_ERA5'])/np.log(10/(mi*DADOS_1H['Hs_ERA5']*((DADOS_1H['Hs_ERA5']/DADOS_1H['Lp_ERA5'])**4.5)))
DADOS_1H['ustar_6_WAVERYS']=(k*DADOS_1H['mag10m_WAVERYS'])/np.log(10/(mi*DADOS_1H['Hs_WAVERYS']*((DADOS_1H['Hs_WAVERYS']/DADOS_1H['Lp_WAVERYS'])**4.5)))
DADOS_1H['ustar_6_P98']=(k*DADOS_1H['mag10_P98'])/np.log(10/(mi*DADOS_1H['Hs_P98']*((DADOS_1H['Hs_P98']/DADOS_1H['Lp_P98'])**4.5)))
DADOS_1H['gama2_ERA5']=((np.abs(DADOS_1H['TSM_ERA5']*DADOS_1H['T_P98']*0.7))/epsilon) ## TERMO DE AJUSTE DA ESTABILIDADE PELA TSM E T
DADOS_1H['gama2_P98']=((np.abs(DADOS_1H['TSM_P98']*DADOS_1H['T_P98']*0.7))/epsilon) ## TERMO DE AJUSTE DA ESTABILIDADE PELA TSM E T
DADOS_1H['zo_TaYe_ERA5_NOVO']=(mi*DADOS_1H['Hs_ERA5']*((DADOS_1H['Hs_ERA5']/DADOS_1H['Lp_ERA5'])**4.5)*(DADOS_1H['gama2_ERA5']))
DADOS_1H['zo_TaYe_P98_NOVO']=(mi*DADOS_1H['Hs_P98']*((DADOS_1H['Hs_P98']/DADOS_1H['Lp_P98'])**4.5)*(DADOS_1H['gama2_P98']))
DADOS_1H['zo_TaYe_WAVERYS_NOVO']=(mi*DADOS_1H['Hs_WAVERYS']*((DADOS_1H['Hs_WAVERYS']/DADOS_1H['Lp_WAVERYS'])**4.5)*(DADOS_1H['gama2_ERA5']))

for i in niveis:
    DADOS_1H['mag10_ERA5_MET9']=DADOS_1H['mag10_ERA5']
    DADOS_1H['mag{}_ERA5_MET9'.format(i)]=(DADOS_1H['ustar_6_ERA5']/k)*((np.log((i)/DADOS_1H['zo_TaYe_ERA5_NOVO']))+((4.7*(i))/90))

for i in niveis:
    DADOS_1H['mag10_WAVERYS_MET9']=DADOS_1H['mag10m_WAVERYS']
    DADOS_1H['mag{}_WAVERYS_MET9'.format(i)]=(DADOS_1H['ustar_6_WAVERYS']/k)*((np.log((i)/DADOS_1H['zo_TaYe_WAVERYS_NOVO']))+((4.7*(i))/90))

for i in niveis:
    DADOS_1H['mag10_PRA01_MET9']=DADOS_1H['mag10_P98']
    DADOS_1H['mag{}_PRA01_MET9'.format(i)]=(DADOS_1H['ustar_6_P98']/k)*((np.log((i)/DADOS_1H['zo_TaYe_P98_NOVO']))+((4.7*(i))/90))


# criando df vazio para calcular os perfis - RS
df_perfis = pd.DataFrame(randn(17,1),index='10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170'.split(),
                  columns='z'.split())

df_perfis['z']=df_perfis.index

TIPOS=['ERA5','PRA01','WAVERYS']
for dados in TIPOS:
    df_perfis['{}_MET1'.format(dados)]=np.nan  #zo_tab2 = 0.01
    df_perfis['{}_MET2'.format(dados)]=np.nan  #zo_tab1 = 0.001
    df_perfis['{}_MET3'.format(dados)]=np.nan   #Lei da potencia alfa    
    df_perfis['{}_MET4'.format(dados)]=np.nan   #zo_Don1990
    df_perfis['{}_MET5'.format(dados)]=np.nan   #zo_Don1990 maior
    df_perfis['{}_MET6'.format(dados)]=np.nan   #zo_Don1990 maior
    df_perfis['{}_MET7'.format(dados)]=np.nan   #zo_Don1990 maior
    df_perfis['{}_MET8'.format(dados)]=np.nan   #zo_Don1990 maior
    df_perfis['{}_MET9'.format(dados)]=np.nan   #zo_Don1990 maior
    df_perfis['LIDAR'.format(dados)]=np.nan   #zo_Don1990 maior

metodologias=['ERA5_MET1','WAVERYS_MET1','PRA01_MET1',
              'ERA5_MET2','WAVERYS_MET2','PRA01_MET2',
              'ERA5_MET3','WAVERYS_MET3','PRA01_MET3',
              'ERA5_MET4','WAVERYS_MET4','PRA01_MET4',
              'ERA5_MET5','WAVERYS_MET5','PRA01_MET5',
              'ERA5_MET6','WAVERYS_MET6','PRA01_MET6',
              'ERA5_MET7','WAVERYS_MET7','PRA01_MET7',
              'ERA5_MET8','WAVERYS_MET8','PRA01_MET8',
              'ERA5_MET9','WAVERYS_MET9','PRA01_MET9']


i=0
while i < len(niveis):
    for metodos in metodologias:
        df_perfis['{}'.format(metodos)][i]=DADOS_1H['mag{}_{}'.format(df_perfis.index[i],metodos)].mean()
    i=i+1

##adicionando dados do lidar
df_perfis['LIDAR'][8]=DADOS_1H['mag90_LIDAR'].mean()
df_perfis['LIDAR'][9]=DADOS_1H['mag100_LIDAR'].mean()
df_perfis['LIDAR'][10]=DADOS_1H['mag110_LIDAR'].mean()
df_perfis['LIDAR'][11]=DADOS_1H['mag120_LIDAR'].mean()
df_perfis['LIDAR'][14]=DADOS_1H['mag150_LIDAR'].mean()
df_perfis['LIDAR'][16]=DADOS_1H['mag170_LIDAR'].mean()



# histograma direcional - FUNCAO PARA CALCULAR
def HistDir(P, D, arqname='', Pmax=None, MaxProb=None, par='hs', MyRange=[],
            dir16=True, interpolado=True, r_angle=None, porcentagem=False,data1=None,data2=None,path=None,localidade=None):
    '''
    Makes a directional histogram for a single point

    P = array/DataFrame(obrigatorio)
        intensidade
    D = array/DataFrame (obrigatorio)
        direção em graus
    arqname =  string (opcional default '')
        string com o nome do output desejado
    Pmax =  None/float(opcional default None)
        Caso queira fixar o eixo de velocidades entre com o valor de um
        inteiro, caso não use None e o valor sera calculado automaticamente
    MaxProb = None/int (opcional default None)
        Caso queira fixar o eixo de probabilidades entre com o valor de
        um float, caso não use None e o valor sera calculado automaticamente
    par = string (opcional default 'hs')
        defina o parametro que esta sendo analisado
    MyRange = list (opcional)
        entre com uma lista dos ranges de insidade, caso não entrar ele farão
        o calculo automático segundo o parametro
    dir16 = True/False (opcional)
        Se as direçõeses serão divididas em 16 ou 8 (default 16 (True))
    interpolado = True/False (opcional)
        Se desejar gerar um histograma interpolado ou rosa
        (default histograma interpolado (True))
    r_angle = None/float (opicional)
        Se desejar alterar a posiçãoo dos labels dos angulos entre com o
        float do angulo, caso nao None (default None)
    porcentagem True/False (opcional)
        Se desejar que o histograma calcule a porcentagem ou as ocorrencias
        (apenas na tabela)

    Detalhes :
        parametros - escolha dependendo da analise
                hs = Altura significativa
                tp = Período de onda
                corrente = corrente
                vento = vento
                energia = energia de onda

    output:
        A png figure and excel table with directional histogram analysis
    '''

    # Transforma as entradas para numpy arrays
    P = asarray(P)
    D = asarray(D)
    # verificar possivel erro do usuario
    if len(D.shape) > 1:
        if D.shape[1] == 1:
            D = D[:, 0]
        elif D.shape[0] == 1:
            D = D[0, :]
        else:
            print(u'Verifique as dimensões da direção')
            sys.pause('')
    if len(P.shape) > 1:
        if P.shape[1] == 1:
            P = P[:, 0]
        elif P.shape[0] == 1:
            P = P[0, :]
        else:
            print(u'Verifique as dimensões do parametro de entrada')
    # verifica se o valor de m�ximo esta definido
    if Pmax is None:
        Pmax = max(P)
    # verifica o parametro para estabelecer os dados de entrada
    if par == 'hs':
        # titulo do plot
        titulo2 = u'Altura de onda (m) - convenção meteorológica'
        # parte do cabeçalho da tabela de ocorrencia conjunta
        cabecalho = u'Altura (m)'
        # parte do nome do arquivo
        fname = arqname + '_altura'
        myfmt = '0.0'
        # define o numero de divisões do parametro na tabela de ocorrencia
        # verifica se há de 5 a 10 classes
        P_bins = arange(0, round(Pmax) + 0.5, 0.5)
        if len(P_bins) > 11:
            P_bins = arange(0, round(Pmax) + 1, 1)
    elif par == 'tp':
        titulo2 = u'Período de onda (s' + u') - convenção meteorológica'
        cabecalho = u'Período (s)'
        fname = arqname + '_periodo'
        myfmt = '0.0'
        P_bins = arange(0, round(Pmax) + 2, 2)
        if len(P_bins) > 11:
            P_bins = arange(0, round(Pmax) + 3, 3)
        elif len(P_bins) < 6:
            P_bins = arange(0, round(Pmax) + 1, 1)
    elif par == 'vento':
        titulo2 = u'ventos (m/s) ' + localidade +' ({} a {})'.format(data1,data2)
        cabecalho = u'Vento(m/s)'
        fname = arqname + '_vento'
        myfmt = '0.0'
        print(Pmax)
        P_bins = arange(0, Pmax + 2, 2)
        if len(P_bins) > 11:
            if len(P_bins) > 11:
                P_bins = arange(0, Pmax + 2.5, 2.5)
                if len(P_bins) > 11:
                    P_bins = arange(0, Pmax + 5, 5)
        elif len(P_bins) < 6:
            P_bins = arange(0, 10. + 1, 1)
            if len(P_bins) < 6:
                P_bins = arange(0, Pmax + 0.5, 0.5)
    elif par == 'energia':
        titulo2 = u'Energia de onda (kJ/m²) - convenção meteorológica'
        cabecalho = u'Energia (kJ/m²)'
        fname = arqname + '_energia'
        myfmt = '0.0'
        P_bins = arange(0, round(Pmax / 10.) * 10. + 2, 2)
        if len(P_bins) > 11:
            P_bins = arange(0, round(Pmax / 10.) * 10. + 3, 3)
        if len(P_bins) > 11:
            P_bins = arange(0, round(Pmax / 10.) * 10. + 5, 5)
        elif len(P_bins) < 6:
            P_bins = arange(0, round(Pmax / 10.) * 10. + 1, 1)
    elif par == 'corrente':
        titulo2 = u'Intensidade da Corrente (m/s) - convenção  oceanográfica'
        cabecalho = u'Corrente(m/s)'
        fname = arqname + '_corrente'
        myfmt = '0.00'
        P_bins = arange(0, round(Pmax * 10.) / 10. + 0.1, 0.1)
        if len(P_bins) > 11:
            if len(P_bins) > 11:
                P_bins = arange(0, round(Pmax * 10.) / 10. + 0.2, 0.2)
                if len(P_bins) > 11:
                    P_bins = arange(0, round(Pmax * 10.) / 10. + 0.25, 0.25)
        elif len(P_bins) < 6:
            P_bins = arange(0, round(Pmax * 10) / 10. + 0.05, 0.05)
            if len(P_bins) < 6:
                P_bins = arange(0, round(Pmax * 10) / 10. + 0.03, 0.03)
    else:
        print(u'defina um dos parâmetro \n hs = Altura significativa \n \
        tp = Período de onda \n corrente = corrente \n vento = vento \n \
        energia = energia de onda')

    # caso o usuario tenha definido o range
    if MyRange != []:
        del P_bins
        P_bins = arange(0, Pmax + MyRange, MyRange)
    # verifica se serão 16 ou 8 direções (se aplica pra tabela)
    if dir16 is True:
        # subtrair da direção
        dirdiff = 11.25
        # numero de bins da direção
        dirrg = 17
        # cabeçalho de direções da tabela
        head = [cabecalho, 'N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE',
                'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW', u'(%)']
    else:
        dirdiff = 22.5
        dirrg = 9
        head = [cabecalho, 'N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW', u'(%)']
    # plota o polar
    fig = plt.figure(figsize=(4.5,4))
    ax = plt.subplot(1, 1, 1, polar=True)
    # define que o zero sera no norte
    ax.set_theta_zero_location("N")
    # roda o grafico azimuital
    ax.set_theta_direction(-1)
    ax.yaxis.grid(color='#898989')
    ax.xaxis.grid(color='#898989')
    if interpolado:
        dir_bins = linspace(-dirdiff, 360 - dirdiff, dirrg)
        if np.any(D < 0):
            D[D < 0] += 360
        D[D > (360 - dirdiff)] -= 360
        P_bins_hist = linspace(0, Pmax, 10)
        # calcular o histograma 2d baseado nos bins
        table = (histogram2d(x=P,
                             y=D + dirdiff,
                             bins=[P_bins_hist, dir_bins],
                             normed=True))
        binarea = P_bins_hist[1] * dir_bins[1] * 200
        # grade para plot
        x, y = meshgrid(table[2][:], table[1][:])
        # calcular a porcentagem
        z = table[0] * binarea
        lowlev = z.max() * 0.05
        z[z < lowlev] = -lowlev * 0.9
        z = insert(z, z.shape[1], (z[:, 0]), axis=1)
        # Duplica ultima coluna para completar o diag polar
        z = insert(z, z.shape[0], (z[-1]), axis=0)
        # adiciona os zeros para que a interpolação não estoure o gráfico
        # converte para radianos
        x = radians(x)

        # interpola todos os parametros
        z = ndimage.zoom(z, 5)
        x = ndimage.zoom(x, 5)
        y = ndimage.zoom(y, 5)
        # calcular a porcentagem
        # definir a prob m5ax
        if MaxProb is None:
            MaxProb = z.max() * 1.1
        # plota dados,define os limites do grafico, palheta de cores, a origem
        cax = ax.contourf(x,
                          y,
                          z,
                          linspace(0, ceil(MaxProb), 20),
                          cmap=plt.cm.hot,
                          origin='lower',
                          antialiased=False)
        # numeros e a distância do grafico polar
        ff = ceil(MaxProb / 10) + 1
        cb = plt.colorbar(cax,
                          pad=.075,
                          shrink=0.8,
                          format='%i',
                          ticks=arange(0, ceil(MaxProb) + ff, ff))
        cb.ax.set_title(u'(%)\n')
        cb.set_clim(0, ceil(MaxProb))
        # define a cor em zero
        cb.vmin = 0
        # tipo
        tipo = u'_histograma_direcional.png'
        # titulo
        titulo1 = u'Histograma direcional -'
    else:

        try:
            P = ma.masked_equal(P, 0)
            D = ma.masked_array(D, P.mask)
            P = ma.MaskedArray.compressed(P)
            D = ma.MaskedArray.compressed(D)
        except BaseException:
            pass

        dir_bins = linspace(-dirdiff, 360 - dirdiff, dirrg)
        if np.any(D < 0):
            D[D < 0] += 360
        if np.any(D > 360):
            D[D > 360] -= 360
        if np.any(D > (360 - dirdiff)):
            D[D > (360 - dirdiff)] -= 360
        # calcular o histograma 2d baseado nos bins
        table = (histogram2d(x=P,
                             y=D,
                             bins=[P_bins, dir_bins],
                             normed=False)
                 )
        theta = radians(table[2][:])[:-1]

        stat = cumsum(table[0], axis=0) / table[0].sum() * 100.
        legenda = []

        windcolors = flipud(
            plt.cm.hot(
                list(map(
                    int, list(
                        round(
                            linspace(0, 250, len(P_bins)
                                     )
                        )
                    )
                )
                )
            )
        )
        for k in flipud(range(size(stat, 0))):
            ax.bar(theta + abs(np.min(theta)),
                   stat[k, :],
                   width=radians(dirdiff * 2),
                   bottom=0.0,
                   color=windcolors[k],
                   edgecolor="k")
            legenda.append('-'.join([str(table[1][k]), str(table[1][k + 1])]))
        if MaxProb is not None:
            ax.set_rmax(MaxProb)
        ax.tick_params(direction='out', length=6, color='r', zorder=10)
        legenda[-1] = u'>' + str(table[1][k])
        plt.legend(legenda, bbox_to_anchor=(1.55, 1.0), title=cabecalho)
        tipo = u'_rosa_direcional.png'
        titulo1 = 'Rosa de '
    # ajustar eixo automaticamente
    if r_angle is None:
        r_angle = dir_bins[argmin(sum(table[0], axis=0))]
    ax.set_rlabel_position(r_angle)
    # carregar o logo
    a = plt.axes([.001, .001, .2, .2], facecolor='None')
    im = plt.imshow(array(Image.open(r'/home/ladsin/Documents/ladsin/Petrobras/Perfis_Vento/ladsin_1.png')))
    plt.axis('off')
    plt.setp(a, xticks=[], yticks=[])
    #logo = plt.imread('logo_atmosmarine_10kb.jpg')
    #ax.figure.figimage(logo, 80, 10, alpha=1, zorder=1)
    # titulo do gráfico
    ax.set_title(titulo1 + titulo2, fontsize=9, y=1.1)
    if arqname == '':
        plt.show()
    else:
        # nome da figura
        plt.savefig(path + fname + tipo, format='png', dpi=300, bbox_inches='tight')
    # limpa a figura
    plt.clf()
    plt.cla()
    plt.close()

    if arqname != '':
        D[D > 360] -= 360
        D[D > (360 - dirdiff)] -= 360
        dir_bins = linspace(0, 360, dirrg) - dirdiff
        # calcular o histograma 2d baseado nos bins
        table = (histogram2d(x=P,
                             y=D,
                             bins=[P_bins, dir_bins],
                             normed=porcentagem)
                 )
        if porcentagem:
            binarea = P_bins[1] * dir_bins[1] * 200
        else:
            binarea = 1
        # escrever xlsx de sa�da
        workbook = xlsxwriter.Workbook(path + fname + u'_ocorrencia_conjunta.xlsx')
        # da o nome do ponto para a tabela
        worksheet = workbook.add_worksheet()
        # informa��es para formata��o
        # tamanho das colunas em cm
        worksheet.set_column('A:B', 10)
        worksheet.set_column('B:S', 5)
        # criar formatos (ja inserindo negrito)
        format1 = workbook.add_format({'bold': True,
                                       'font_name': 'Arial', 'font_size': 10,
                                       'align': 'center',
                                       'bg_color': '#C9C9C9'})
        format2 = workbook.add_format({'bold': False,
                                       'font_name': 'Arial', 'font_size': 10,
                                       'align': 'center',
                                       'bg_color': '#FFFFFF'})
        format3 = workbook.add_format({'bold': False,
                                       'font_name': 'Arial', 'font_size': 10,
                                       'align': 'center',
                                       'bg_color': '#FFFFFF'})
        format4 = workbook.add_format({'bold': False,
                                       'font_name': 'Arial', 'font_size': 10,
                                       'align': 'center',
                                       'bg_color': '#C9C9C9'})
        format5 = workbook.add_format({'bold': False,
                                       'font_name': 'Arial', 'font_size': 10,
                                       'align': 'center',
                                       'bg_color': '#C9C9C9'})
        format6 = workbook.add_format({'bold': True,
                                       'font_name': 'Arial', 'font_size': 10,
                                       'align': 'center',
                                       'bg_color': '#C9C9C9'})
        format7 = workbook.add_format({'bold': False,
                                       'font_name': 'Arial', 'font_size': 10,
                                       'align': 'center',
                                       'bg_color': '#C9C9C9'})
        # formata��o das casas decimais
        if porcentagem:
            format2.set_num_format('0.00')
        else:
            format2.set_num_format('0')
        format3.set_num_format('0.0')
        format5.set_num_format(myfmt)
        format7.set_num_format(myfmt)
        # inserir linhas de divis�o da c�lula
        format4.set_top(1)
        format6.set_bottom(1)
        format7.set_bottom(1)

        # insere o cabeçalho no arquivo
        for k, hd in enumerate(head):
            worksheet.write(0, k, hd, format6)
        # escreve as linhas de ocorrencia
        for j in range((len(table[1]) - 1)):
            worksheet.write(
                j + 1,
                0,
                '-'.join([
                    str(table[1][j]),
                    str(table[1][j + 1])]
                ).replace('.', ','),
                format1)

            for i in range(len(head) - 2):
                worksheet.write(
                    j + 1, i + 1, table[0][j, i] * binarea, format2)
            if porcentagem:
                worksheet.write(j +
                                1, i +
                                2, np.sum(table[0][j, :]) *
                                binarea, format3)
            else:
                worksheet.write(
                    j + 1,
                    i + 2,
                    np.sum(table[0][j, :]) / np.sum(table[0]) * 100,
                    format3)
        # escreve o total e a porcentagem de valores por direção
        worksheet.write(j + 2, 0, u'(%)', format1)
        totais = np.sum(table[0], axis=0) * binarea
        if porcentagem is False:
            totais = totais / np.sum(totais) * 100
        for i, total in enumerate(totais):
            worksheet.write(j + 2, i + 1, total, format5)

        worksheet.write(j + 2, i + 3, '', format5)
        worksheet.write(j + 2, i + 2, '', format5)
        worksheet.write(j + 3, i + 3, '', format5)
        worksheet.write(j + 3, i + 2, '', format5)
        worksheet.write(j + 4, i + 3, '', format7)
        worksheet.write(j + 4, i + 2, '', format7)

        # escreve as medias e os maximos de cada direção
        worksheet.write(j + 3, 0, u'Media', format5)
        worksheet.write(j + 4, 0, u'Max.', format7)
        if np.ma.isMaskedArray(P) is False:
            P = np.ma.masked_object(P, -99999)
        for l in range(len(table[2]) - 1):
            if len(np.ma.MaskedArray.compressed(
                    P[(D > table[2][l]) & (D < table[2][l + 1])])) == 0:
                worksheet.write(j + 4, l + 1, 0, format7)
                worksheet.write(j + 3, l + 1, 0, format5)
            else:
                worksheet.write(
                    j + 4, l + 1, np.nanmax(
                        P[(D > table[2][l]) & (D < table[2][l + 1])]
                        ),
                    format7)
                worksheet.write(
                    j + 3, l + 1, np.nanmean(
                        P[(D > table[2][l]) & (D < table[2][l + 1])]
                        ),
                    format5)

        # encerra o arquivo
        workbook.close()



#plotando a rosa dos ventos - mag10_P18
ax, fig = plt.subplots(figsize=(8,4))
data1=(str(DADOS_1H.index))[16:26]
data2=(str(DADOS_1H.index))[-97:-87]
HistDir(DADOS_1H['mag10_ERA5'].dropna(),
        DADOS_1H['dir10_ERA5'].dropna(),
        path=pathname+'/imagens/',
        par='vento',
        arqname = 'windrose_ERA5_10m',
        interpolado = False,
        MaxProb = 16,
        localidade='ERA5',
        data1=data1,
        data2=data2)


#plotando a rosa dos ventos - mag10_P18
ax, fig = plt.subplots(figsize=(8,4))
data1=(str(DADOS_1H.index))[16:26]
data2=(str(DADOS_1H.index))[-97:-87]
HistDir(DADOS_1H['mag10m_WAVERYS'].dropna(),
        DADOS_1H['dir10m_WAVERYS'].dropna(),
        path=pathname+'/imagens/',
        par='vento',
        arqname = 'windrose_WAVERYS_10m',
        interpolado = False,
        MaxProb = 16,
        localidade='WAVERYS',
        data1=data1,
        data2=data2)

#plotando a rosa dos ventos - mag10_P18
ax, fig = plt.subplots(figsize=(8,4))
data1=(str(DADOS_1H.index))[16:26]
data2=(str(DADOS_1H.index))[-97:-87]
HistDir(DADOS_1H['mag90_LIDAR'].dropna(),
        DADOS_1H['dir90_LIDAR'].dropna(),
        path=pathname+'/imagens/',
        par='vento',
        arqname = 'windrose_LIDAR90m',
        interpolado = False,
        MaxProb = 16,
        localidade='LIDAR - 90m',
        data1=data1,
        data2=data2)


#plotando a rosa dos ventos - mag10_P18
ax, fig = plt.subplots(figsize=(8,4))
data1=(str(DADOS_1H.index))[16:26]
data2=(str(DADOS_1H.index))[-97:-87]
HistDir(DADOS_1H['mag100_LIDAR'].dropna(),
        DADOS_1H['dir100_LIDAR'].dropna(),
        path=pathname+'/imagens/',
        par='vento',
        arqname = 'windrose_LIDAR100m',
        interpolado = False,
        MaxProb = 16,
        localidade='LIDAR - 100m',
        data1=data1,
        data2=data2)


#plotando a rosa dos ventos - mag10_P18
ax, fig = plt.subplots(figsize=(8,4))
data1=(str(DADOS_1H.index))[16:26]
data2=(str(DADOS_1H.index))[-97:-87]
HistDir(DADOS_1H['mag110_LIDAR'].dropna(),
        DADOS_1H['dir110_LIDAR'].dropna(),
        path=pathname+'/imagens/',
        par='vento',
        arqname = 'windrose_LIDAR110m',
        interpolado = False,
        MaxProb = 16,
        localidade='LIDAR - 110m',
        data1=data1,
        data2=data2)


#plotando a rosa dos ventos - mag10_P18
ax, fig = plt.subplots(figsize=(8,4))
data1=(str(DADOS_1H.index))[16:26]
data2=(str(DADOS_1H.index))[-97:-87]
HistDir(DADOS_1H['mag120_LIDAR'].dropna(),
        DADOS_1H['dir120_LIDAR'].dropna(),
        path=pathname+'/imagens/',
        par='vento',
        arqname = 'windrose_LIDAR120m',
        interpolado = False,
        MaxProb = 16,
        localidade='LIDAR - 120m',
        data1=data1,
        data2=data2)


#plotando a rosa dos ventos - mag10_P18
ax, fig = plt.subplots(figsize=(8,4))
data1=(str(DADOS_1H.index))[16:26]
data2=(str(DADOS_1H.index))[-97:-87]
HistDir(DADOS_1H['mag150_LIDAR'].dropna(),
        DADOS_1H['dir150_LIDAR'].dropna(),
        path=pathname+'/imagens/',
        par='vento',
        arqname = 'windrose_LIDAR150m',
        interpolado = False,
        MaxProb = 16,
        localidade='LIDAR - 150m',
        data1=data1,
        data2=data2)


#plotando a rosa dos ventos - mag10_P18
ax, fig = plt.subplots(figsize=(8,4))
data1=(str(DADOS_1H.index))[16:26]
data2=(str(DADOS_1H.index))[-97:-87]
HistDir(DADOS_1H['mag170_LIDAR'].dropna(),
        DADOS_1H['dir170_LIDAR'].dropna(),
        path=pathname+'/imagens/',
        par='vento',
        arqname = 'windrose_LIDAR170m',
        interpolado = False,
        MaxProb = 16,
        localidade='LIDAR - 170m',
        data1=data1,
        data2=data2)


metodologias=['ERA5_MET1','WAVERYS_MET1','PRA01_MET1',
              'ERA5_MET2','WAVERYS_MET2','PRA01_MET2',
              'ERA5_MET3','WAVERYS_MET3','PRA01_MET3',
              'ERA5_MET4','WAVERYS_MET4','PRA01_MET4',
              'ERA5_MET5','WAVERYS_MET5','PRA01_MET5',
              'ERA5_MET6','WAVERYS_MET6','PRA01_MET6',
              'ERA5_MET7','WAVERYS_MET7','PRA01_MET7',
              'ERA5_MET8','WAVERYS_MET8','PRA01_MET8',
              'ERA5_MET9','WAVERYS_MET9','PRA01_MET9']

#TODOS OS CALCULOS USAM WAVERYS, ERA5 E PRA01
#METODO 1 - DNV TABELADO
#METODO 2 - LEI DA POTENCIA
#METODO 3 - DONELAN 90 - NEUTRO
#METODO 4 - DONELAN 90 - ESTAVEL
#METODO 5 - DONELAN 93 - NEUTRO
#METODO 6 - DONELAN 93 - ESTAVEL
#METODO 7 - TAYLOR & YELLAND 01 - NEUTRO
#METODO 8 - TAYLOR & YELLAND 01 - ESTAVEL
#METODO 9 - CARMO ET AL - ADAPTADO PRA REGIAO

#PERFIS DE VENTO - TODOS
fig,ax=plt.subplots(figsize=(12,15))
line1= ax.plot(df_perfis['ERA5_MET1'],df_perfis.index, color='red',label='DNV', linestyle='--')
line2= ax.plot(df_perfis['ERA5_MET2'],df_perfis.index, color='lightgray', label='LEI DA POTÊNCIA', linestyle='--')
line3= ax.plot(df_perfis['ERA5_MET3'],df_perfis.index, color='turquoise',label='DONELAN 90 - NEUTRO', linestyle='--')
#line4= ax.plot(df_perfis['ERA5_MET4'],df_perfis.index, color='turquoise',label='DONELAN 90 - ESTÁVEL', linestyle=':')
line5= ax.plot(df_perfis['ERA5_MET5'],df_perfis.index, color='lawngreen',label='DONELAN ET AL 93 - NEUTRO', linestyle='--')
#line6= ax.plot(df_perfis['ERA5_MET6'],df_perfis.index, color='lawngreen',label='DONELAN ET AL 93 - ESTÁVEL', linestyle=':')
line7= ax.plot(df_perfis['ERA5_MET7'],df_perfis.index, color='indigo',label='CARMO ET AL 22 - NEUTRO', linestyle='--')
#line8= ax.plot(df_perfis['ERA5_MET8'],df_perfis.index, color='indigo',label='CARMO ET AL 22 - ESTÁVEL', linestyle=':')
line9= ax.plot(df_perfis['ERA5_MET9'],df_perfis.index, color='green',label='NOVO MÉTODO', linestyle='--')
line10=plt.scatter(df_perfis['LIDAR'],df_perfis.index, marker="*",s=200, color='black', label='LIDAR')

ax.set_xlim(2,17)
ax.legend()
#criando o for para fazer o vetor - Maranhão
for k in range (0,np.size(df_perfis.index)):
    qax = plt.quiver(0,df_perfis.index[k],df_perfis['ERA5_MET7'][k],0, units="dots",
             width=2.5,
             headwidth = 2,
             headlength = 4,
             headaxislength = 5,
             alpha=0.5,
             color='gray',
             scale=0.0225)
ax.set_title('Perfis de vento médios - ERA5 - 2020', fontsize=20)
ax.set_ylabel('Altura (m)', fontsize=20)
ax.set_xlabel('Velocidade (m/s)', fontsize=20)
plt.savefig(pathname+'/imagens/perfis_medios_ERA5.png') #coloca 2 digitos
plt.clf() #fecha o arquivo para nao plotar por cima


#PERFIS DE VENTO - TODOS
fig,ax=plt.subplots(figsize=(12,15))
line1= ax.plot(df_perfis['WAVERYS_MET1'],df_perfis.index, color='red',label='DNV', linestyle='--')
line2= ax.plot(df_perfis['WAVERYS_MET2'],df_perfis.index, color='lightgray', label='LEI DA POTÊNCIA', linestyle='--')
line3= ax.plot(df_perfis['WAVERYS_MET3'],df_perfis.index, color='turquoise',label='DONELAN 90 - NEUTRO', linestyle='--')
#line4= ax.plot(df_perfis['ERA5_MET4'],df_perfis.index, color='turquoise',label='DONELAN 90 - ESTÁVEL', linestyle=':')
line5= ax.plot(df_perfis['WAVERYS_MET5'],df_perfis.index, color='lawngreen',label='DONELAN ET AL 93 - NEUTRO', linestyle='--')
#line6= ax.plot(df_perfis['ERA5_MET6'],df_perfis.index, color='lawngreen',label='DONELAN ET AL 93 - ESTÁVEL', linestyle=':')
line7= ax.plot(df_perfis['WAVERYS_MET7'],df_perfis.index, color='indigo',label='CARMO ET AL 22 - NEUTRO', linestyle='--')
#line8= ax.plot(df_perfis['ERA5_MET8'],df_perfis.index, color='indigo',label='CARMO ET AL 22 - ESTÁVEL', linestyle=':')
line9= ax.plot(df_perfis['WAVERYS_MET9'],df_perfis.index, color='green',label='NOVO MÉTODO', linestyle='--')
line10=plt.scatter(df_perfis['LIDAR'],df_perfis.index, marker="*",s=200, color='black', label='LIDAR')

ax.set_xlim(2,17)
ax.legend()
#criando o for para fazer o vetor - Maranhão
for k in range (0,np.size(df_perfis.index)):
    qax = plt.quiver(0,df_perfis.index[k],df_perfis['WAVERYS_MET7'][k],0, units="dots",
             width=2.5,
             headwidth = 2,
             headlength = 4,
             headaxislength = 5,
             alpha=0.5,
             color='gray',
             scale=0.0225)
ax.set_title('Perfis de vento médios - WAVERYS - 2020', fontsize=20)
ax.set_ylabel('Altura (m)', fontsize=20)
ax.set_xlabel('Velocidade (m/s)', fontsize=20)
plt.savefig(pathname+'/imagens/perfis_medios_WAVERYS.png') #coloca 2 digitos
plt.clf() #fecha o arquivo para nao plotar por cima




#PERFIS DE VENTO - TODOS
fig,ax=plt.subplots(figsize=(12,15))
line1= ax.plot(df_perfis['PRA01_MET1'],df_perfis.index, color='red',label='DNV', linestyle='--')
line2= ax.plot(df_perfis['PRA01_MET2'],df_perfis.index, color='lightgray', label='LEI DA POTÊNCIA', linestyle='--')
line3= ax.plot(df_perfis['PRA01_MET3'],df_perfis.index, color='turquoise',label='DONELAN 90 - NEUTRO', linestyle='--')
#line4= ax.plot(df_perfis['ERA5_MET4'],df_perfis.index, color='turquoise',label='DONELAN 90 - ESTÁVEL', linestyle=':')
line5= ax.plot(df_perfis['PRA01_MET5'],df_perfis.index, color='lawngreen',label='DONELAN ET AL 93 - NEUTRO', linestyle='--')
#line6= ax.plot(df_perfis['ERA5_MET6'],df_perfis.index, color='lawngreen',label='DONELAN ET AL 93 - ESTÁVEL', linestyle=':')
line7= ax.plot(df_perfis['PRA01_MET7'],df_perfis.index, color='indigo',label='CARMO ET AL 22 - NEUTRO', linestyle='--')
#line8= ax.plot(df_perfis['ERA5_MET8'],df_perfis.index, color='indigo',label='CARMO ET AL 22 - ESTÁVEL', linestyle=':')
line9= ax.plot(df_perfis['PRA01_MET9'],df_perfis.index, color='green',label='NOVO MÉTODO', linestyle='--')
line10=plt.scatter(df_perfis['LIDAR'],df_perfis.index, marker="*",s=200, color='black', label='LIDAR')

ax.set_xlim(2,17)
ax.legend()
#criando o for para fazer o vetor - Maranhão
for k in range (0,np.size(df_perfis.index)):
    qax = plt.quiver(0,df_perfis.index[k],df_perfis['PRA01_MET7'][k],0, units="dots",
             width=2.5,
             headwidth = 2,
             headlength = 4,
             headaxislength = 5,
             alpha=0.5,
             color='gray',
             scale=0.0225)
ax.set_title('Perfis de vento médios - PRA01 - 2020', fontsize=20)
ax.set_ylabel('Altura (m)', fontsize=20)
ax.set_xlabel('Velocidade (m/s)', fontsize=20)
plt.savefig(pathname+'/imagens/perfis_medios_PRA01.png') #coloca 2 digitos
plt.clf() #fecha o arquivo para nao plotar por cima



#GRAFICOS DE DISPERSAO
#ESTIMATIVA X BOIA P18
DADOS_1H['linhazero']=0
fig, ax=plt.subplots(2,1,sharex=True,figsize=(20,8))
fig.suptitle('MÉTODO DA DNV - VENTO A 90 METROS', fontsize=20)
lineMET= ax[0].plot(DADOS_1H.index,DADOS_1H['mag90_ERA5_MET1'],color='blue',label='ERA5', linestyle=':')
lineMET2= ax[0].plot(DADOS_1H.index,DADOS_1H['mag90_WAVERYS_MET1'],color='green',label='WAVERYS', linestyle=':')
lineP18= ax[0].plot(DADOS_1H.index,DADOS_1H['mag90_LIDAR'],color='red',label='LIDAR', linestyle=':')

#diferenças
difMET= ax[1].plot(DADOS_1H.index,DADOS_1H['mag90_ERA5_MET1']-DADOS_1H['mag90_LIDAR'], color='darkgoldenrod',label='DIF = ERA5 - LIDAR', linestyle=':')
difMET2= ax[1].plot(DADOS_1H.index,DADOS_1H['mag90_WAVERYS_MET1']-DADOS_1H['mag90_LIDAR'], color='deepskyblue',label='DIF = WAVERYS - LIDAR', linestyle=':')
eixo0= ax[1].plot(DADOS_1H.index, DADOS_1H['linhazero'], color='gray', linestyle='--')

ax[0].set_title('Magnitude do vento (m/s)', fontsize=14)
ax[1].set_title('DIFERENÇA ENTRE ESTIMATIVAS E LIDAR', fontsize=14)
ax[0].set_ylabel('m/s', fontsize=15)
ax[1].set_ylabel('m/s', fontsize=15)
ax[1].set_xlabel('data', fontsize=15)
ax[0].legend()
ax[1].legend()
ax[0].set_ylim(0,30)
ax[1].set_ylim(-10,10)
plt.savefig(pathname+'/imagens/DISP_ERA5_WAV_DNV_90m.png') #coloca 2 digitos
plt.clf() #fecha o arquivo para nao plotar por cima


#GRAFICOS DE DISPERSAO
#ESTIMATIVA X BOIA P18
DADOS_1H['linhazero']=0
fig, ax=plt.subplots(2,1,sharex=True,figsize=(20,8))
fig.suptitle('MÉTODO DA DNV - VENTO A 100 METROS', fontsize=20)
lineMET= ax[0].plot(DADOS_1H.index,DADOS_1H['mag100_ERA5_MET1'],color='blue',label='ERA5', linestyle=':')
lineMET2= ax[0].plot(DADOS_1H.index,DADOS_1H['mag100_WAVERYS_MET1'],color='green',label='WAVERYS', linestyle=':')
lineP18= ax[0].plot(DADOS_1H.index,DADOS_1H['mag100_LIDAR'],color='red',label='LIDAR', linestyle=':')

#diferenças
difMET= ax[1].plot(DADOS_1H.index,DADOS_1H['mag100_ERA5_MET1']-DADOS_1H['mag100_LIDAR'], color='darkgoldenrod',label='DIF = ERA5 - LIDAR', linestyle=':')
difMET2= ax[1].plot(DADOS_1H.index,DADOS_1H['mag100_WAVERYS_MET1']-DADOS_1H['mag100_LIDAR'], color='deepskyblue',label='DIF = WAVERYS - LIDAR', linestyle=':')
eixo0= ax[1].plot(DADOS_1H.index, DADOS_1H['linhazero'], color='gray', linestyle='--')

ax[0].set_title('Magnitude do vento (m/s)', fontsize=14)
ax[1].set_title('DIFERENÇA ENTRE ESTIMATIVAS E LIDAR', fontsize=14)
ax[0].set_ylabel('m/s', fontsize=15)
ax[1].set_ylabel('m/s', fontsize=15)
ax[1].set_xlabel('data', fontsize=15)
ax[0].legend()
ax[1].legend()
ax[0].set_ylim(0,30)
ax[1].set_ylim(-10,10)
plt.savefig(pathname+'/imagens/DISP_ERA5_WAV_DNV_100m.png') #coloca 2 digitos
plt.clf() #fecha o arquivo para nao plotar por cima


#GRAFICOS DE DISPERSAO
#ESTIMATIVA X BOIA P18
DADOS_1H['linhazero']=0
fig, ax=plt.subplots(2,1,sharex=True,figsize=(20,8))
fig.suptitle('MÉTODO DA DNV - VENTO A 110 METROS', fontsize=20)
lineMET= ax[0].plot(DADOS_1H.index,DADOS_1H['mag110_ERA5_MET1'],color='blue',label='ERA5', linestyle=':')
lineMET2= ax[0].plot(DADOS_1H.index,DADOS_1H['mag110_WAVERYS_MET1'],color='green',label='WAVERYS', linestyle=':')
lineP18= ax[0].plot(DADOS_1H.index,DADOS_1H['mag110_LIDAR'],color='red',label='LIDAR', linestyle=':')

#diferenças
difMET= ax[1].plot(DADOS_1H.index,DADOS_1H['mag110_ERA5_MET1']-DADOS_1H['mag110_LIDAR'], color='darkgoldenrod',label='DIF = ERA5 - LIDAR', linestyle=':')
difMET2= ax[1].plot(DADOS_1H.index,DADOS_1H['mag110_WAVERYS_MET1']-DADOS_1H['mag110_LIDAR'], color='deepskyblue',label='DIF = WAVERYS - LIDAR', linestyle=':')
eixo0= ax[1].plot(DADOS_1H.index, DADOS_1H['linhazero'], color='gray', linestyle='--')

ax[0].set_title('Magnitude do vento (m/s)', fontsize=14)
ax[1].set_title('DIFERENÇA ENTRE ESTIMATIVAS E LIDAR', fontsize=14)
ax[0].set_ylabel('m/s', fontsize=15)
ax[1].set_ylabel('m/s', fontsize=15)
ax[1].set_xlabel('data', fontsize=15)
ax[0].legend()
ax[1].legend()
ax[0].set_ylim(0,30)
ax[1].set_ylim(-10,10)
plt.savefig(pathname+'/imagens/DISP_ERA5_WAV_DNV_110m.png') #coloca 2 digitos
plt.clf() #fecha o arquivo para nao plotar por cima



#GRAFICOS DE DISPERSAO
#ESTIMATIVA X BOIA P18
DADOS_1H['linhazero']=0
fig, ax=plt.subplots(2,1,sharex=True,figsize=(20,8))
fig.suptitle('MÉTODO DA DNV - VENTO A 120 METROS', fontsize=20)
lineMET= ax[0].plot(DADOS_1H.index,DADOS_1H['mag120_ERA5_MET1'],color='blue',label='ERA5', linestyle=':')
lineMET2= ax[0].plot(DADOS_1H.index,DADOS_1H['mag120_WAVERYS_MET1'],color='green',label='WAVERYS', linestyle=':')
lineP18= ax[0].plot(DADOS_1H.index,DADOS_1H['mag120_LIDAR'],color='red',label='LIDAR', linestyle=':')

#diferenças
difMET= ax[1].plot(DADOS_1H.index,DADOS_1H['mag120_ERA5_MET1']-DADOS_1H['mag120_LIDAR'], color='darkgoldenrod',label='DIF = ERA5 - LIDAR', linestyle=':')
difMET2= ax[1].plot(DADOS_1H.index,DADOS_1H['mag120_WAVERYS_MET1']-DADOS_1H['mag120_LIDAR'], color='deepskyblue',label='DIF = WAVERYS - LIDAR', linestyle=':')
eixo0= ax[1].plot(DADOS_1H.index, DADOS_1H['linhazero'], color='gray', linestyle='--')

ax[0].set_title('Magnitude do vento (m/s)', fontsize=14)
ax[1].set_title('DIFERENÇA ENTRE ESTIMATIVAS E LIDAR', fontsize=14)
ax[0].set_ylabel('m/s', fontsize=15)
ax[1].set_ylabel('m/s', fontsize=15)
ax[1].set_xlabel('data', fontsize=15)
ax[0].legend()
ax[1].legend()
ax[0].set_ylim(0,30)
ax[1].set_ylim(-10,10)
plt.savefig(pathname+'/imagens/DISP_ERA5_WAV_DNV_120m.png') #coloca 2 digitos
plt.clf() #fecha o arquivo para nao plotar por cima



#import oceanspy as ospy
#z0 = ospy.rugosidade.garratt1980(Hs)


#######PROXIMOS PASSOS
#### REFAZER O PROCESSO PARA 10minutos
