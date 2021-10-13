# -*- coding: utf-8 -*-
# Example plotting for the MBARI-WEC 2021-10 dataset. For more information, 
# please see accompanying README.md file.

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import os
import getpass
from mhkit.wave.resource import energy_flux, energy_period
from mhkit.wave.performance import capture_length_matrix
import seaborn as sns
from palettable.colorbrewer.sequential import YlOrRd_3 as my_colors01
from palettable.colorbrewer.qualitative import Set2_3 as qual_colors
import scipy.io as sio
from datetime import datetime
import cartopy.crs as ccrs
import geopy.distance

file_name = "SPOT-0589_2021-09-03_2021-09-29_download.csv"
base_name = os.path.splitext(file_name)[0]
date_parser = lambda epoch: pd.to_datetime(epoch, unit='s')

rho = 1025      # density [kg/m^3]
h = 100         # depth [m]
g = 9.81        # gravity [m/s^2]

#%% Read bulk parameters

dat = pd.read_csv(file_name,
	index_col=3,
    usecols=np.insert(np.arange(13),-1,[364,365,366]),
    parse_dates=['Epoch Time'],
    date_parser=date_parser,
    )

b = dat.to_xarray()

print(b)

#%% Frequency array

dat1 = pd.read_csv(file_name,
                   index_col=[],
                   usecols=np.arange(13,13+38+1),
                   )
freq_array = dat1.iloc[0].to_xarray()
freq_array.name = 'Frequency'

dat2 = pd.read_csv(file_name,
                   index_col=[],
                   usecols=np.arange(13+38+1,13+2*(38+1)),
                   )
df_array = dat2.iloc[0].to_xarray().assign_coords(dict(index=freq_array.values)).rename(dict(index='Frequency'))
df_array.name = 'df'
df_array.attrs['long_name'] = 'Frequency spacing'
df_array.attrs['units'] = 'Hz'

#%% a and b parameters

names = ['a1','b1','a2','b2']
tmp_list = []
for idx, name in enumerate(names):
    dat_tmp = pd.read_csv(file_name,
                       index_col=[0],
                       usecols=np.insert(np.arange(13+(2+idx)*(38+1),13+(3+idx)*(38+1)),0,3),
                       date_parser=date_parser,
                       )
    tmp_da = dat_tmp.to_xarray().to_array(dim='Frequency', name=name)
    tmp_da = tmp_da.assign_coords({'Frequency':freq_array.values})
    tmp_list.append(tmp_da)

ab_ds = xr.merge(tmp_list)

#%% Spectral density, spreading, etc.

dat_S = pd.read_csv(file_name,
                   index_col=[0],
                   usecols=np.insert(np.arange(13+6*(38+1),13+7*(38+1)),0,3),
                   date_parser=date_parser,
                   )
S = dat_S.to_xarray().to_array(dim='Frequency', name='Variance density')
S = S.assign_coords({'Frequency':freq_array.values})

dat_dir = pd.read_csv(file_name,
                   index_col=[0],
                   usecols=np.insert(np.arange(13+7*(38+1),13+8*(38+1)),0,3),
                   date_parser=date_parser,
                   )
Dir = dat_dir.to_xarray().to_array(dim='Frequency', name='Direction')
Dir = Dir.assign_coords({'Frequency':freq_array.values})

dat_spread = pd.read_csv(file_name,
                   index_col=[0],
                   usecols=np.insert(np.arange(13+8*(38+1),13+9*(38+1)),0,3),
                   date_parser=date_parser,
                   )
spread = dat_spread.to_xarray().to_array(dim='Frequency', name='Directional spread')
spread = spread.assign_coords({'Frequency':freq_array.values})

#%% Combine, clean up, and save to netcdf

da = xr.merge([b,ab_ds,S,Dir,spread,df_array])
da['Battery Voltage (V)'].attrs['units'] = 'V'
da['Battery Voltage (V)'].attrs['long_name'] = 'Battery voltage'

da['Power (W)'].attrs['units'] = 'W'
da['Power (W)'].attrs['long_name'] = 'Battery power'

da['Humidity (%rel)'].attrs['units'] = '1'
da['Humidity (%rel)'].attrs['standard_name'] = 'relative_humidity'
da['Humidity (%rel)'].attrs['long_name'] = 'Relative humidity'

da['Significant Wave Height (m)'].attrs['units'] = 'm'
da['Significant Wave Height (m)'].attrs['standard_name'] = 'sea_surface_wave_significant_height'
da['Significant Wave Height (m)'].attrs['long_name'] = 'Significant wave height'

da['Direction'].attrs['units'] = 'degree'
da['Direction'].attrs['long_name'] = ''

da['Peak Period (s)'].attrs['units'] = 's'
da['Peak Period (s)'].attrs['standard_name'] = 'sea_surface_wave_period_at_variance_spectral_density_maximum'
da['Peak Period (s)'].attrs['long_name'] = 'Peak period'

da['Mean Period (s)'].attrs['units'] = 's'
da['Mean Period (s)'].attrs['standard_name'] = 'sea_surface_wave_zero_upcrossing_period'
da['Mean Period (s)'].attrs['long_name'] = 'Mean period'

da['Peak Direction (deg)'].attrs['units'] = 'degree'
da['Peak Direction (deg)'].attrs['standard_name'] = 'sea_surface_wave_from_direction_at_variance_spectral_density_maximum'
da['Peak Direction (deg)'].attrs['long_name'] = 'Peak direction'

da['Peak Directional Spread (deg)'].attrs['units'] = 'degree'
da['Peak Directional Spread (deg)'].attrs['standard_name'] = 'sea_surface_wave_directional_spread_at_variance_spectral_density_maximum'
da['Peak Directional Spread (deg)'].attrs['long_name'] = 'Peak directional spread'

da['Mean Direction (deg)'].attrs['units'] = 'degree'
da['Mean Direction (deg)'].attrs['standard_name'] = 'sea_surface_wave_from_direction'
da['Mean Direction (deg)'].attrs['long_name'] = 'Mean direction'

da['Mean Directional Spread (deg)'].attrs['units'] = 'degree'
da['Mean Directional Spread (deg)'].attrs['long_name'] = 'Mean directional spread'

da['Latitude (deg)'].attrs['units'] = 'degree_north'
da['Latitude (deg)'].attrs['standard_name'] = 'latitude'
da['Latitude (deg)'].attrs['long_name'] = 'Latitude'

da['Longitude (deg)'].attrs['units'] = 'degree_east'
da['Longitude (deg)'].attrs['standard_name'] = 'longitude'
da['Longitude (deg)'].attrs['long_name'] = 'Longitude'

da['Wind Speed (m/s)'].attrs['units'] = 'm/s'
da['Wind Speed (m/s)'].attrs['standard_name'] = 'wind_speed'
da['Wind Speed (m/s)'].attrs['long_name'] = 'Wind speed'

da['Wind Direction (deg)'].attrs['units'] = 'degree'
da['Wind Direction (deg)'].attrs['standard_name'] = 'wind_from_direction'
da['Wind Direction (deg)'].attrs['long_name'] = 'Wind direction'

da['Surface Temperature (°C)'] = 274.15*da['Surface Temperature (°C)']
da['Surface Temperature (°C)'].attrs['units'] = 'K'
da['Surface Temperature (°C)'].attrs['standard_name'] = 'sea_surface_temperature'
da['Surface Temperature (°C)'].attrs['long_name'] = 'Surface temperature'

da['Frequency'].attrs['units'] = 'Hz'
da['Frequency'].attrs['standard_name'] = 'wave_frequency'
da['Frequency'].attrs['long_name'] = 'Frequency'

da['Variance density'].attrs['units'] = 'm^2/Hz'
da['Variance density'].attrs['standard_name'] = 'sea_surface_wave_variance_spectral_density'
da['Variance density'].attrs['long_name'] = 'Spectral density'

da['Directional spread'].attrs['units'] = 'degree'
da['Directional spread'].attrs['standard_name'] = 'sea_surface_wave_directional_spread'
da['Directional spread'].attrs['long_name'] = 'Directional spreading'


da = da.rename({'Epoch Time':'time',
                'Frequency':'freq',
                'Battery Voltage (V)':'batter_voltage',
                'Variance density':'S',
                'Direction':'wave_dir',
               'Power (W)':'battery_power',
               'Humidity (%rel)':'humidity',
               'Significant Wave Height (m)':'Hm0',
               'Peak Period (s)':'Tp',
               'Mean Period (s)':'Tm',
               'Peak Direction (deg)':'peak_dir',
               'Peak Directional Spread (deg)':'peak_spread',
               'Mean Direction (deg)':'mean_dir',
               'Mean Directional Spread (deg)':'mean_spread',
               'Directional spread':'spread',
               'Latitude (deg)':'spot_lat',
               'Longitude (deg)':'spot_lon',
               'Wind Speed (m/s)':'wind_speed',
               'Wind Direction (deg)':'wind_dir',
               'Surface Temperature (°C)':'temperature'})

J = xr.DataArray(np.array([energy_flux(da.isel(time=idx)['S'].to_pandas(),h=h,rho=rho).values[0][0] 
                           for idx in range(len(da.time))]),
                 dims='time',
                 name='J').assign_coords(dict(time=da.time.values))
J.attrs['units'] = 'W/m'
J.attrs['long_name'] = 'Energy flux'

Te = xr.DataArray(np.array([energy_period(da.isel(time=idx)['S'].to_pandas(),
                                          ) for idx in range(len(da.time))]).squeeze(),
                  dims='time',
                  name='Te').assign_coords(dict(time=da.time.values))
Te.attrs['units'] = 's'
Te.attrs['long_name'] = 'Energy period'

da = xr.merge([da, J, Te])

da.time.attrs['long_name'] = 'Epoch time'

da.attrs['institution'] = 'Sandia National Laboratories and Monterey Bay Aquarium Research Institute'
da.attrs['Conventions'] = 'CF-1.8'
da.attrs['title'] = base_name
da.attrs['source'] = 'Sofar spotter buoy'
da.attrs['history'] = 'generated {:} by {:}'.format(datetime.now().strftime('%Y-%m-%d @ %H:%M:%S'),
                                  getpass.getuser())
da.attrs['references'] = 'https://content.sofarocean.com/hubfs/Technical_Reference_Manual.pdf'
da = da.sortby('time')
da = da.drop_isel(time=0) # first sample appears anomalous
print(da)
da.to_netcdf(base_name + '.nc',
             encoding={'time':{'units': "seconds since 2000-01-01 00:00:00"}})
ds = da

#%% WEC data
# see https://doi.org/10.1007/s40722-021-00197-9 for variable descriptions

mat = sio.loadmat('HourlyAverages_2021_09_28.mat')
wec_power = mat['HourlyData']['PC_Power'][0][0].squeeze()
wec_damping = mat['HourlyData']['PC_Scale'][0][0].squeeze()
wec_date = pd.to_datetime(mat['HourlyData']['Time_DateNum'][0][0].squeeze() \
                          -719529, unit='D')
    
wec_lat = mat['HourlyData']['XB_Lat'][0][0].squeeze()
wec_lon = mat['HourlyData']['XB_Lon'][0][0].squeeze()

wec = xr.Dataset(data_vars=dict(power=('time',wec_power), 
                                damping=('time',wec_damping),
                                wec_lat=('time',wec_lat),
                                wec_lon=('time',wec_lon),
                                ),
                 coords=dict(time=wec_date))

wec['power'].attrs['long_name'] = 'WEC power'
wec['power'].attrs['units'] = 'W'

wec['damping'].attrs['long_name'] = 'Damping factor'
wec['damping'].attrs['units'] = ' '

wec['wec_lat'].attrs['long_name'] = 'Latitude'
wec['wec_lat'].attrs['standard_name'] = 'latitude'
wec['wec_lat'].attrs['units'] = 'degree_north'

wec['wec_lon'].attrs['long_name'] = 'Longitude'
wec['wec_lon'].attrs['standard_name'] = 'longitude'
wec['wec_lon'].attrs['units'] = 'degree_east'


wec['cw'] = wec.power / ds.J.interp_like(wec.power, method='nearest')
wec['cw'].attrs['long_name'] = 'Capture width'
wec['cw'].attrs['units'] = 'm'

print(wec)

#%% Combine datasets

ds = xr.merge([da, wec.interp_like(da)])

dist = [geopy.distance.distance((ds.spot_lat.values[idx], ds.spot_lon.values[idx]), 
               (ds.wec_lat.values[idx], ds.wec_lon.values[idx])).km \
        for idx in range(len(ds.spot_lat.values))]
dist_da = xr.DataArray(data=dist, 
                       dims='time', 
                       coords=dict(time=ds.time.values),
                       name='dist',
                       attrs=dict(
                           units='km',
                           long_name='Distance',
                           ))

ds = xr.merge([ds, dist_da])

#%%

Te_lims = [5,13]
Hm0_lims = [0.5, 2.5]

fig, ax = plt.subplots()
ds.S.sel(freq=slice(1/4)).dropna('time').plot.contourf(ax=ax,
                    levels=12,
                    cmap=my_colors01.mpl_colormap)
fig.tight_layout()
ax.autoscale(enable=True, axis='x', tight=True)
fig.savefig(base_name + '_spectrogram.pdf')

#%% 

vars_to_plot = ['J',
                'Hm0',
                'Te',
                'mean_dir',
                'power',
                'cw',
                ]
fig, ax = plt.subplots(nrows=len(vars_to_plot),
                       figsize=(8,12),
                       sharex=True)

for axi, var in zip(ax, vars_to_plot):
    axi.set_prop_cycle('color', qual_colors.mpl_colors)
    qual_colors
    ds[var].plot(ax=axi,
                 marker='.',
                 )
    axi.label_outer()
    axi.spines['right'].set_visible(False)
    axi.spines['top'].set_visible(False)
    axi.autoscale(enable=True, axis='x', tight=True)

ax[3].fill_between(ds.time.values, ds.mean_dir, ds.mean_dir + ds.mean_spread,
                   color=qual_colors.mpl_colors[0],
                   alpha=0.25,
                   )
p1 = ax[3].fill_between(ds.time.values, ds.mean_dir, ds.mean_dir - ds.mean_spread,
                   color=qual_colors.mpl_colors[0],
                   alpha=0.25,
                   )

p1.set_label('Mean directional spread')
ax[3].legend(fontsize='small')
for axi in ax:
    for item in ([axi.title, axi.xaxis.label, axi.yaxis.label] + axi.get_xticklabels() + axi.get_yticklabels()):
        item.set_fontsize('x-small')

fig.savefig(base_name + '_time_history.pdf')

#%% 

fig, ax = plt.subplots()
ds.plot.scatter(x="Te", y="Hm0", 
    hue="J", 
    ax=ax,
    # s=4,
    alpha=0.5)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
fig.savefig(base_name + '_energy_flux_scatter.pdf')

#%% 

fig, ax = plt.subplots()
ds.plot.scatter(x='Te', y='Hm0', hue='cw',
                ax=ax,
                alpha=0.5)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

ax.set_xlim(Te_lims)
ax.set_ylim(Hm0_lims)

fig.savefig(base_name + '_capture_width_scatter.pdf')

#%%

def my_bins(var, spacing=1.0):
    bins = np.arange(np.floor(var.min()), np.ceil(var.max() + spacing), spacing)   
    return bins

Hm0_bins = my_bins(ds.dropna(dim='time', subset=['cw']).Hm0, spacing=0.5)
Te_bins = my_bins(ds.dropna(dim='time', subset=['cw']).Te, spacing=0.5)

fig = plt.figure()

ax = sns.histplot(x=ds.Te, 
                 y=ds.Hm0, 
                 # color=my_colors.mpl_colors[0],
                  bins=(Te_bins, Hm0_bins),
                 stat='count',
                 cbar=True,
                 cmap=my_colors01.mpl_colormap,
                  cbar_kws={'label': 'Number of observations'},
                 alpha=0.5,
                 )
ds.plot.scatter(x="Te", y="Hm0", 
    c='k',
    ax=ax,
    s=2,
    alpha=0.5)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlim(Te_lims)
ax.set_ylim(Hm0_lims)
fig.savefig(base_name + '_hist2d.pdf')

#%%

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
ds.plot.scatter(x='spot_lon', y='spot_lat',
                ax=ax,
                subplot_kws=dict(projection=ccrs.Orthographic(-80, 35), facecolor="gray"),
                transform=ccrs.PlateCarree(),
                label='Spotter',
                )
ds.plot.scatter(x='wec_lon', y='wec_lat',
                ax=ax,
                subplot_kws=dict(projection=ccrs.Orthographic(-80, 35), facecolor="gray"),
                transform=ccrs.PlateCarree(),
                label='WEC',
                )

ax.legend()
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, 
                  color='gray', 
                  alpha=0.5, 
                  linestyle='--',
                  zorder=0,
                  )
gl.xlabels_top = False
gl.ylabels_left = False

fig.savefig(base_name + '_geolocation.pdf')

#%% 

m = capture_length_matrix(ds.dropna(dim='time', subset=['cw']).Hm0.values, 
                          ds.dropna(dim='time', subset=['cw']).Te.values, 
                          ds.dropna(dim='time', subset=['cw']).cw.values, 
                          'mean', 
                          Hm0_bins[:-1] + np.diff(Hm0_bins)/2, 
                          Te_bins[:-1] + np.diff(Te_bins)/2)
cwm = xr.DataArray(m, dims=('Hm0','Te'),
                   name='Capture width',
                   attrs=dict(units='m'))

fig, ax = plt.subplots()
cwm.plot.pcolormesh(ax=ax,
                    alpha=0.5,
                    cmap=my_colors01.mpl_colormap)
ds.dropna(dim='time', subset=['power']).plot.scatter(x="Te", y="Hm0",
                ax=ax,
                c='k',
                s=2,
                )
ax.set_xlim(Te_lims)
ax.set_ylim(Hm0_lims)

fig.savefig(base_name + '_captureWidthMatrix.pdf')
