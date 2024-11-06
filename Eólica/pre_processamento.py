#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 10:41:23 2024

@author: ladsin
"""

# Faz um NCzão juntando

import xarray as xr
import numpy as np
pathname2='/media/ladsin/Seagate Expansion Drive/'
#                                 ABRINDO INTERPPOLAÇÕES 
nc_wind=xr.open_mfdataset(pathname2+'/ERA5VENTO(NOVO)/ERA5_WIND_INTERPOLADO.nc') #abrindo os arquivos nc de vento
nc_wave=xr.open_mfdataset(pathname2+'/ERA5VENTO(NOVO)/ERA5_WAVE_INTERPOLADO.nc') #abrindo os arquivos nc de vento
nc_T=xr.open_mfdataset(pathname2+'/ERA5VENTO(NOVO)/ERA5_T_INTERPOLADO.nc') #abrindo os arquivos nc de vento
nc_waverys=xr.open_mfdataset(pathname2+'/ERA5VENTO(NOVO)/WAVERYS_INTERPOLADO_6H.nc') #abrindo os arquivos nc de vento
nc_oi=xr.open_mfdataset(pathname2+'/ERA5VENTO(NOVO)/OISST_INTERPOLADO_6H.nc') #abrindo os arquivos nc de vento
nc_mur=xr.open_mfdataset(pathname2+'/ERA5VENTO(NOVO)/TSM_MUR_INTERPOLADO_6H.nc') #abrindo os arquivos nc de vento



#                                 ABRINDO VARIÁVEIS 
nc_out = xr.Dataset()

u10=nc_wind.u10
nc_out['u10'] = xr.DataArray(u10, dims=('time','latitude','longitude'),
                         coords={'time': nc_wind.time, 'latitude': nc_wind.latitude,
                                 'longitude':nc_wind.longitude})
nc_out.u10.attrs['long_name'] = "componente u a 10m"
nc_out.u10.attrs['units'] = "m s**-1"


v10=nc_wind.v10
nc_out['v10'] = xr.DataArray(v10, dims=('time','latitude','longitude'),
                         coords={'time': nc_wind.time, 'latitude': nc_wind.latitude,
                                 'longitude':nc_wind.longitude})
nc_out.v10.attrs['long_name'] = "componente v a 10m"
nc_out.v10.attrs['units'] = "m s**-1"


mwp=nc_wave.mwp
nc_out['mwp'] = xr.DataArray(mwp, dims=('time','latitude','longitude'),
                         coords={'time': nc_wind.time, 'latitude': nc_wind.latitude,
                                 'longitude':nc_wind.longitude})
nc_out.mwp.attrs['long_name'] = "Mean wave Period"
nc_out.mwp.attrs['units'] = "s"  

mwd=nc_wave.mwd
nc_out['mwd'] = xr.DataArray(mwd, dims=('time','latitude','longitude'),
                         coords={'time': nc_wind.time, 'latitude': nc_wind.latitude,
                                 'longitude':nc_wind.longitude})
nc_out.mwd.attrs['long_name'] = "Mean Wave Direction"
nc_out.mwd.attrs['units'] = "Degree true"

pp1d=nc_wave.pp1d
nc_out['pp1d'] = xr.DataArray(pp1d, dims=('time','latitude','longitude'),
                         coords={'time': nc_wind.time, 'latitude': nc_wind.latitude,
                                 'longitude':nc_wind.longitude})
nc_out.pp1d.attrs['long_name'] = "Peak Wave Period"
nc_out.pp1d.attrs['units'] = "s"

swh=nc_wave.swh
nc_out['swh'] = xr.DataArray(swh, dims=('time','latitude','longitude'),
                         coords={'time': nc_wind.time, 'latitude': nc_wind.latitude,
                                 'longitude':nc_wind.longitude})
nc_out.swh.attrs['long_name'] = "Significant height of combined wind waves and swell"
nc_out.swh.attrs['units'] = "m"

shts=nc_wave.shts
nc_out['shts'] = xr.DataArray(shts, dims=('time','latitude','longitude'),
                         coords={'time': nc_wind.time, 'latitude': nc_wind.latitude,
                                 'longitude':nc_wind.longitude})
nc_out.shts.attrs['long_name'] = "Significant height of total swell"
nc_out.shts.attrs['units'] = "m"

shww=nc_wave.shww
nc_out['shww'] = xr.DataArray(shww, dims=('time','latitude','longitude'),
                         coords={'time': nc_wind.time, 'latitude': nc_wind.latitude,
                                 'longitude':nc_wind.longitude})
nc_out.shww.attrs['long_name'] = "Significant height of of wind waves"
nc_out.shww.attrs['units'] = "m"

hmax=nc_wave.hmax
nc_out['hmax'] = xr.DataArray(hmax, dims=('time','latitude','longitude'),
                         coords={'time': nc_wind.time, 'latitude': nc_wind.latitude,
                                 'longitude':nc_wind.longitude})
nc_out.hmax.attrs['long_name'] = "Maximum individual wave height"
nc_out.hmax.attrs['units'] = "m"

t2m=nc_T.t2m
nc_out['t2m'] = xr.DataArray(t2m, dims=('time','latitude','longitude'),
                         coords={'time': nc_wind.time, 'latitude': nc_wind.latitude,
                                 'longitude':nc_wind.longitude})
nc_out.t2m.attrs['long_name'] = "Temperatura a 2 metros"
nc_out.t2m.attrs['units'] = "K"

d2m=nc_T.d2m
nc_out['d2m'] = xr.DataArray(d2m, dims=('time','latitude','longitude'),
                         coords={'time': nc_wind.time, 'latitude': nc_wind.latitude,
                                 'longitude':nc_wind.longitude})
nc_out.d2m.attrs['long_name'] = "Temperatura do Ponto de Orvalho a 2 metros"
nc_out.d2m.attrs['units'] = "K"

sp=nc_T.sp
nc_out['sp'] = xr.DataArray(sp, dims=('time','latitude','longitude'),
                         coords={'time': nc_wind.time, 'latitude': nc_wind.latitude,
                                 'longitude':nc_wind.longitude})
nc_out.sp.attrs['long_name'] = "Pressão na Superfície"
nc_out.sp.attrs['units'] = "Pa"

sst=nc_T.sst
nc_out['sst'] = xr.DataArray(sst, dims=('time','latitude','longitude'),
                         coords={'time': nc_wind.time, 'latitude': nc_wind.latitude,
                                 'longitude':nc_wind.longitude})
nc_out.sst.attrs['long_name'] = "Temperatura da Superfície do Mar"
nc_out.sst.attrs['units'] = "K"

VHM0=nc_waverys.VHM0
nc_out['VHM0'] = xr.DataArray(VHM0, dims=('time','latitude','longitude'),
                         coords={'time': nc_wind.time, 'latitude': nc_wind.latitude,
                                 'longitude':nc_wind.longitude})
nc_out.VHM0.attrs['long_name'] = "Spectral significant wave height"
nc_out.VHM0.attrs['units'] = "m"

VTPK=nc_waverys.VTPK
nc_out['VTPK'] = xr.DataArray(VTPK, dims=('time','latitude','longitude'),
                         coords={'time': nc_wind.time, 'latitude': nc_wind.latitude,
                                 'longitude':nc_wind.longitude})
nc_out.VTPK.attrs['long_name'] = "Wave period at spectral peak / peak period (Tp)"
nc_out.VTPK.attrs['units'] = "s"

VTM01_WW=nc_waverys.VTM01_WW
nc_out['VTM01_WW'] = xr.DataArray(VTM01_WW, dims=('time','latitude','longitude'),
                         coords={'time': nc_wind.time, 'latitude': nc_wind.latitude,
                                 'longitude':nc_wind.longitude})
nc_out.VTM01_WW.attrs['long_name'] = "Sea surface wind wave mean period"
nc_out.VTM01_WW.attrs['units'] = "s"

VMDR=nc_waverys.VMDR
nc_out['VMDR'] = xr.DataArray(VMDR, dims=('time','latitude','longitude'),
                         coords={'time': nc_wind.time, 'latitude': nc_wind.latitude,
                                 'longitude':nc_wind.longitude})
nc_out.VMDR.attrs['long_name'] = "Mean wave direction from (Mdir)"
nc_out.VMDR.attrs['units'] = "degree"

sst_OI=nc_oi.sst_OI
nc_out['sst'] = xr.DataArray(sst_OI, dims=('time','latitude','longitude'),
                         coords={'time': nc_wind.time, 'latitude': nc_wind.latitude,
                                 'longitude':nc_wind.longitude})
nc_out.sst_OI.attrs['long_name'] = "Daily Sea Surface Temperature"
nc_out.sst_OI.attrs['units'] = "degC"

analysed_sst=nc_mur.analysed_sst
nc_out['analysed_sst'] = xr.DataArray(analysed_sst, dims=('time','latitude','longitude'),
                         coords={'time': nc_wind.time, 'latitude': nc_wind.latitude,
                                 'longitude':nc_wind.longitude})
nc_out.analysed_sst.attrs['long_name'] = "Analysed sea surface temperature"
nc_out.analysed_sst.attrs['units'] = "kelvin"

sst_anomaly=nc_mur.sst_anomaly
nc_out['analysed_sst'] = xr.DataArray(sst_anomaly, dims=('time','latitude','longitude'),
                         coords={'time': nc_wind.time, 'latitude': nc_wind.latitude,
                                 'longitude':nc_wind.longitude})
nc_out.sst_anomaly.attrs['long_name'] = "SST anomaly from a seasonal SST climatology based"
nc_out.sst_anomaly.attrs['units'] = "kelvin"

nc_out.to_netcdf(pathname2+'/ERA5VENTO(NOVO)/DATA_FINAL.nc', format='NETCDF4_CLASSIC' )