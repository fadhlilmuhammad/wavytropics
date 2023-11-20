
import numpy as np
import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker
from cartopy.util import add_cyclic_point
import pandas as pd
from matplotlib import pyplot as plt
import xarray as xr
# import iris 

#adapted from 
#https://scitools-iris.readthedocs.io/en/v3.0.1/generated/gallery/general/plot_SOI_filtering.html

#Also, thanks to Lingwei Meng for the ideas in applying the Lanczos bandpass filter; see her GitHub (https://github.com/liv0505/Lanczos-Filter)

def low_pass_weights(window, cutoff):
    # """Calculate weights for a low pass Lanczos filter.

    # Args:

    # window: int
    #     The length of the filter window.

    # cutoff: float
    #     The cutoff frequency in inverse time steps.

    # """
    order = ((window - 1) // 2) + 1
    nwts = 2 * order + 1
    w = np.zeros([nwts])
    n = nwts // 2
    w[n] = 2 * cutoff
    k = np.arange(1.0, n)
    sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
    firstfactor = np.sin(2.0 * np.pi * cutoff * k) / (np.pi * k)
    w[n - 1 : 0 : -1] = firstfactor * sigma
    w[n + 1 : -1] = firstfactor * sigma
    return w[1:-1]


def band_pass(data, lofr, hifr, nwts):
    #  """Applying band pass to a time series

    # Args:

    # nwts: int
    #     The length of the filter window.

    # lofr: float
    #     low frequency time-cutoff (1/days).
        
    # hifr: float
    #     High frequency time-cutoff (1/days).

    # """
    wgts_lo = low_pass_weights(nwts, lofr)
    wgts_hi = low_pass_weights(nwts, hifr)

    da_wgtlo = xr.DataArray(wgts_lo, dims = ['window'])
    da_wgthi = xr.DataArray(wgts_hi, dims = ['window'])

    data_lo = data.rolling(time = len(wgts_lo), center = True).construct('window').dot(da_wgtlo)   
    data_hi = data.rolling(time = len(wgts_hi), center = True).construct('window').dot(da_wgthi)

    #the bandpass is the difference between high and low
    data_bp = data_hi - data_lo

    return(data_bp)

    


#Plotting function
def plot_composites(composites,
                    namefig='composite_plot',
                    cmap='RdYlBu_r', 
                    levels=np.arange(-16,17,2), 
                    latN=25,latS=-25,
                    central_longitude=120,
                    lonL=-180,LonR=180,
                    lonticks = np.arange(-180,181,60),
                    latticks = np.arange(-25,26,15),
                    add_cyclic = True,
                    coord_name = ['lat','lon'],
                    title = 'MJO',
                    figsize = (15,15),
                   ):

    
    import matplotlib.pyplot as plt
    levels = levels
    fig = plt.figure(figsize=figsize, dpi=300)
    
    # Set the axes using the specified map projection
    ax1 = fig.add_subplot(4,2,1,projection=ccrs.PlateCarree(central_longitude=central_longitude))
    ax2 = fig.add_subplot(4,2,3,projection=ccrs.PlateCarree(central_longitude=central_longitude),sharex=ax1)
    ax3 = fig.add_subplot(4,2,5,projection=ccrs.PlateCarree(central_longitude=central_longitude),sharex=ax2)
    ax4 = fig.add_subplot(4,2,7,projection=ccrs.PlateCarree(central_longitude=central_longitude),sharex=ax2)
    ax5 = fig.add_subplot(4,2,2,projection=ccrs.PlateCarree(central_longitude=central_longitude),sharex=ax2)
    ax6 = fig.add_subplot(4,2,4,projection=ccrs.PlateCarree(central_longitude=central_longitude),sharex=ax2)
    ax7 = fig.add_subplot(4,2,6,projection=ccrs.PlateCarree(central_longitude=central_longitude),sharex=ax2)
    ax8 = fig.add_subplot(4,2,8,projection=ccrs.PlateCarree(central_longitude=central_longitude),sharex=ax2)
    

    if add_cyclic:
    # Add cyclic point to data
        data1=composites["phase-1"]
        data1, lons = add_cyclic_point(data1, coord=data1[coord_name[1]])
        data2=composites["phase-2"]
        data2, lons = add_cyclic_point(data2, coord=data2[coord_name[1]])
        data3=composites["phase-3"]
        data3, lons = add_cyclic_point(data3, coord=data3[coord_name[1]])
        data4=composites["phase-4"]
        data4, lons = add_cyclic_point(data4, coord=data4[coord_name[1]])
        data5=composites["phase-5"]
        data5, lons = add_cyclic_point(data5, coord=data5[coord_name[1]])
        data6=composites["phase-6"]
        data6, lons = add_cyclic_point(data6, coord=data6[coord_name[1]])
        data7=composites["phase-7"]
        data7, lons = add_cyclic_point(data7, coord=data7[coord_name[1]])
        data8=composites["phase-8"]
        data8, lons = add_cyclic_point(data8, coord=data8[coord_name[1]])
    else:
        data1=composites["phase-1"]
        data2=composites["phase-2"]
        data3=composites["phase-3"]
        data4=composites["phase-4"]
        data5=composites["phase-5"]
        data6=composites["phase-6"]
        data7=composites["phase-7"]
        data8=composites["phase-8"]
    

    coord = composites['phase-1']
    # Make a filled contour plot
    cs = ax1.contourf(lons, coord[coord_name[0]], data1,
                projection = ccrs.PlateCarree(central_longitude=central_longitude),
                transform = ccrs.PlateCarree(),cmap=cmap,extend='both',levels=levels)
    ax2.contourf(lons, coord[coord_name[0]], data2,
                projection = ccrs.PlateCarree(central_longitude=central_longitude),
                transform = ccrs.PlateCarree(),cmap=cmap,extend='both',levels=levels)
    ax3.contourf(lons, coord[coord_name[0]], data3,
                projection = ccrs.PlateCarree(central_longitude=central_longitude),
                transform = ccrs.PlateCarree(),cmap=cmap,extend='both',levels=levels)
    ax4.contourf(lons, coord[coord_name[0]], data4,
                projection = ccrs.PlateCarree(central_longitude=central_longitude),
                transform = ccrs.PlateCarree(),cmap=cmap,extend='both',levels=levels)
    ax5.contourf(lons, coord[coord_name[0]], data5,
                projection = ccrs.PlateCarree(central_longitude=central_longitude),
                transform = ccrs.PlateCarree(),cmap=cmap,extend='both',levels=levels)
    ax6.contourf(lons, coord[coord_name[0]], data6,
                projection = ccrs.PlateCarree(central_longitude=central_longitude),
                transform = ccrs.PlateCarree(),cmap=cmap,extend='both',levels=levels)
    ax7.contourf(lons, coord[coord_name[0]], data7,
                projection = ccrs.PlateCarree(central_longitude=central_longitude),
                transform = ccrs.PlateCarree(),cmap=cmap,extend='both',levels=levels)
    ax8.contourf(lons, coord[coord_name[0]], data8,
                projection = ccrs.PlateCarree(central_longitude=central_longitude),
                transform = ccrs.PlateCarree(),cmap=cmap,extend='both',levels=levels)
    
    # Add coastlines
    ax1.coastlines()
    ax2.coastlines()
    ax3.coastlines()
    ax4.coastlines()
    ax5.coastlines()
    ax6.coastlines()
    ax7.coastlines()
    ax8.coastlines()
    
    # Define the xticks for longitude
    ax4.set_xticks(lonticks, crs=ccrs.PlateCarree(central_longitude=central_longitude))
    ax4.set_xticklabels(lonticks,fontsize=12)  
    ax8.set_xticks(lonticks, crs=ccrs.PlateCarree(central_longitude=central_longitude))
    ax8.set_xticklabels(lonticks, fontsize=12)  
    
    lon_formatter = cticker.LongitudeFormatter()
    ax4.xaxis.set_major_formatter(lon_formatter)
    ax8.xaxis.set_major_formatter(lon_formatter)
    
    # Define the yticks for latitude
    ax1.set_yticks(latticks, crs=ccrs.PlateCarree())
    ax2.set_yticks(latticks, crs=ccrs.PlateCarree())
    ax3.set_yticks(latticks, crs=ccrs.PlateCarree())
    ax4.set_yticks(latticks, crs=ccrs.PlateCarree())
    lat_formatter = cticker.LatitudeFormatter()
    ax1.yaxis.set_major_formatter(lat_formatter)
    ax2.yaxis.set_major_formatter(lat_formatter)
    ax3.yaxis.set_major_formatter(lat_formatter)
    ax4.yaxis.set_major_formatter(lat_formatter)
    
    
    ax1.text(0.013, 0.75, 'Phase-1', transform=ax1.transAxes, 
                size=10, weight='bold').set_bbox(dict(facecolor='w', alpha=1.0, edgecolor='k'))
    ax2.text(0.013, -0.62, 'Phase-2', transform=ax1.transAxes, 
                size=10, weight='bold').set_bbox(dict(facecolor='w', alpha=1.0, edgecolor='k'))
    ax3.text(0.013, -2.01, 'Phase-3', transform=ax1.transAxes, 
                size=10, weight='bold').set_bbox(dict(facecolor='w', alpha=1.0, edgecolor='k'))
    ax4.text(0.013, -3.38, 'Phase-4', transform=ax1.transAxes, 
                size=10, weight='bold').set_bbox(dict(facecolor='w', alpha=1.0, edgecolor='k'))
    ax5.text(1.094, 0.75, 'Phase-5', transform=ax1.transAxes, 
                size=10, weight='bold').set_bbox(dict(facecolor='w', alpha=1.0, edgecolor='k'))
    ax6.text(1.094, -0.62, 'Phase-6', transform=ax1.transAxes, 
                size=10, weight='bold').set_bbox(dict(facecolor='w', alpha=1.0, edgecolor='k'))
    ax7.text(1.094, -2.01, 'Phase-7', transform=ax1.transAxes, 
                size=10, weight='bold').set_bbox(dict(facecolor='w', alpha=1.0, edgecolor='k'))
    ax8.text(1.094, -3.38, 'Phase-8', transform=ax1.transAxes, 
                size=10, weight='bold').set_bbox(dict(facecolor='w', alpha=1.0, edgecolor='k'))
    
    plt.subplots_adjust(left=0.,
                        bottom=0., 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.08, 
                        hspace=-0.9)
    
    axs = (ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8)
    # Add colorbar
    cbar = plt.colorbar(cs, ax = axs,shrink=.3,pad=0.02)
    
    plt.suptitle(title, y=0.6,x=0.375)
    # plt.tight_layout()
    plt.savefig(namefig, format = 'png',dpi=fig.dpi)
    plt.show()