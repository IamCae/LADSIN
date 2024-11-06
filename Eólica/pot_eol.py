#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 12:32:27 2024

@author: ladsin
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

pathname2='/media/ladsin/Seagate Expansion Drive/'

# Carrega o arquivo NetCDF com os dados de vento
ds = xr.open_dataset(pathname2+'/ERA5VENTO(NOVO)/ERA5_WIND_INTERPOLADO.nc')  # substitua pelo caminho do seu arquivo
ds1 = ds.isel(time=0)
#ds1 = ds.sel(time='2018-12-31T18:00:00')

# Supondo que a variável de velocidade do vento seja 'wind_speed' (ajuste conforme necessário)
wind_speed = np.sqrt(ds1['u10']**2+ds1['v10']**2)

# Calcula o potecial eólico (potência por unidade de área)
rho = 1.225  # densidade do ar em kg/m³
wind_power_density = 0.5 * rho * wind_speed**3

# Define a projeção para o mapa e o layout da figura
fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={'projection': ccrs.PlateCarree()})
#ax.set_extent([-60, -30, -35, 4.5], crs=ccrs.PlateCarree())  # Limites para a costa do Brasil

# Plota o potencial eólico
ax = plt.pcolormesh(ds['lon'], ds['lat'], wind_power_density,
                           cmap='viridis', shading='auto')
cbar = plt.colorbar(orientation='vertical', pad=0.02)
cbar.set_label('Potencial Eólico (W/m²)')

# Adiciona recursos ao mapa
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.LAND, edgecolor='black')
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.LAKES, alpha=0.5)
ax.add_feature(cfeature.RIVERS)

# Título e labels
plt.title('Potencial Eólico ao longo da Costa do Brasil')
plt.show()
