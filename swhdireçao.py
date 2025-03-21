#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 15:05:37 2024

@author: ladsin
"""

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.units import units
import numpy as np
import xarray as xr
import cartopy.io.shapereader as shpreader  # Import shapefiles
from datetime import datetime, timedelta  # Basic dates and time types
import cmocean
import matplotlib.colors as mcolors
import os  # Importando módulo os para manipulação de diretórios

# Abrindo os dados do ERA5 
altura = xr.open_dataset('/home/ladsin/Chalo/.nc/isaoceano.nc').metpy.parse_cf()
direçao = xr.open_dataset('/home/ladsin/Chalo/.nc/isaoceano.nc').metpy.parse_cf()

# Recorte da latitude e longitude
lon_slice = slice(-100, -10)
lat_slice = slice(10, -60)

lon_slice2 = slice(-100, -10, 5)
lat_slice2 = slice(10, -60, 5)

# Puxar os dados do arquivo já com o recorte
lons = altura.longitude.sel(longitude=lon_slice).values
lats = altura.latitude.sel(latitude=lat_slice).values
lons2 = direçao.longitude.sel(longitude=lon_slice2).values
lats2 = direçao.latitude.sel(latitude=lat_slice2).values

# Define o diretório onde as figuras serão salvas
output_dir = '/home/ladsin/Chalo/IMAGENS/isa'
# Cria o diretório se ele não existir
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Definir cada variável 
for i in range(len(altura.variables['valid_time'])):  # Alterado para valid_time
    
    direçao['mwd'].attrs['units'] = 'degrees'
    d = direçao['mwd'].metpy.sel(
    valid_time=direçao.valid_time[i],
    latitude=lat_slice2,
    longitude=lon_slice2
    ).metpy.unit_array.squeeze().to('degrees').magnitude  # .magnitude para converter para um array NumPy

        
    swh = altura['swh'].metpy.sel(
        valid_time=altura.valid_time[i],  # Usar valid_time aqui
        latitude=lat_slice,
        longitude=lon_slice
    ).metpy.unit_array.squeeze()

    vtime = altura.valid_time.data[i].astype('datetime64[ms]')

    plt.figure(figsize=(10, 10))
    
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black', linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=0.5)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), color='gray', alpha=1.0, linestyle='--', linewidth=0.25, 
                      xlocs=np.arange(-100, -10, 5), ylocs=np.arange(-60, 10, 5), draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    ax.add_feature(cfeature.LAND)
    
    colors = ["#F8F8FF", "#00BFFF", "#1E90FF", "#228B22", "#7CFC00", 
              "#ADFF2F", "#F0E68C", "#FFFF00", "#FFD700", "#FFA500", 
              "#FF0000", "#800000", '#A020F0']

    cmap2 = mcolors.LinearSegmentedColormap.from_list("", colors)
    cmap2.set_over('#A020F0')
    cmap2.set_under('#F8F8FF')
    levels = np.arange(2.5, 9, .5)
    contour = ax.contourf(lons, lats, swh,
                          cmap=cmap2, levels=levels, extend='both')
    u = -np.sin(np.deg2rad(d))  # Conversão para componente U do vento
    v = -np.cos(np.deg2rad(d))  # Conversão para componente V do vento
    
    ax.quiver(lons2,
              lats2,
              u,
              v,
              scale=None,
              headwidth=5,
              headlength=5,
              headaxislength=2.5,
              minshaft=1,
              minlength=1,
              units='width',
              scale_units=None,
              angles='uv',
              width=0.002,
              pivot='tail',
              color='black')
    
    shapefile = list(shpreader.Reader('/home/ladsin/BR_UF_2019.shp').geometries())
    ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='black', facecolor='none', linewidth=0.8)
    
    cbar = plt.colorbar(contour, ax=ax, pad=0.06, fraction=0.04)
    cbar.set_label(label='Altura Significativa (m)', size=12)
    cbar.ax.tick_params(labelsize=12)
    
    # Defina o título
    title_str = f'Altura Significativa/Direçao Media: {vtime}'
    plt.title(title_str, fontsize=16, loc='right')

    # Salva a figura em um arquivo com o mesmo nome do título
    # Remove caracteres indesejados do título para criar um nome de arquivo válido
    file_name = f"Altura_Significativa_{vtime.astype('datetime64[s]')}.png"
    
    plt.savefig(os.path.join(output_dir, file_name), dpi=300, bbox_inches='tight')  # Salvar a figura
    plt.close()  # Fecha a figura para liberar memória
