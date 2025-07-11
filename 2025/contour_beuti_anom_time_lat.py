import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import calendar
from matplotlib import gridspec
from matplotlib import rcParams
from matplotlib import interactive
interactive(True)
# pylint: disable=C0103



# plot paramaters
params = {
    'text.latex.preamble': '\\usepackage{gensymb}',
    'image.origin': 'lower',
    'image.interpolation': 'nearest',
    'image.cmap': 'RdYlBu_r',
    'axes.grid': False,
    'savefig.dpi': 300,  # to adjust notebook inline plot size
    'xtick.top':        False,  # shold the top and bottom have tick marks
    'xtick.bottom':     True,
    'xtick.major.size': 2.5,
    'ytick.major.size': 2.5,
    'ytick.direction': 'out',
    'xtick.direction': 'out',
    'axes.labelsize': 12,  # fontsize for x and y labels 
    'axes.titlesize': 12,
    'font.size': 12,  # was 10
    'legend.fontsize': 12,  # was 10
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'figure.figsize': [8.5, 11],
    'font.family': 'STIXGeneral',
    'toolbar': 'None',
    'savefig.bbox': 'tight',
    'axes.spines.top': True,
    'axes.spines.bottom': True,
    'axes.spines.left': True,
    'axes.spines.right': True,
    'font.family': 'Arial',
}
rcParams.update(params)


class nf(float):
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()


# -------------------------------------------------------
# -- Input variables, change these
# -------------------------------------------------------

# clim year perio
yr_clim_bgn = 1982
yr_clim_end = 2010

# end year
yr_end = 2025

# dir_out for plots
dir_out = './'

# lats of the sst data
lat_wnt = [31, 47]

# distance want
dis_wnt = [0, 75]
# dis_wnt = [0, 150]


# file name of the BEUTI data
dir_in = './'
fn_beuti = '{}BEUTI_monthly.nc'.format(dir_in)

# roms variable name in the xr.ds
var_roms = ['sst_oi']

# beuti variable name in the xr.ds
var_beuti = ['ui_mon']

# roms dimensions labels
dim_roms = ['time', 'latitude', 'longitude']

# beuti dimensions labels
dim_beuti = ['lat', 'time']

# dim of the final anomaly xr.da
dim_da = ['lat', 'time']
dim_da_sst = ['latitude', 'time']

# figure usually has 1 columns, 2 rows
num_clmn = 1
num_row = 2

# ylim and ticks
ylm = [31, 47]
y_tck = np.arange(ylm[0], ylm[1]+3, 3)

# CalCOFI figures usually focus discussion over the last 3 years, change
# if more years are wanted
num_xyrs = 5
# num_xyrs = 28

# color for the xyrs
color_xyrs = ['green', 'orange', 'dodgerblue', 'red']

# figure size
fig_wdth = 8.5
fig_hght = 6

# nlevels
d_beuti = 3
dmin_beuti = -18
dmax_beuti = 18
nlvl_beuti = np.arange(dmin_beuti-6*d_beuti, dmax_beuti+7*d_beuti, d_beuti)
nlvl_beuti1 = np.arange(dmin_beuti-1*d_beuti, dmax_beuti+2*d_beuti, 2*d_beuti)

# colorbar labelsb
clrbr_lbl_beuti = 'BEUTI Anoms (mmol s^-1 m^-1)'

# month begin and end
month_bgn = 1
month_end = 12

# -------------------------------------------------------
# -- END: Input variables, change these
# -------------------------------------------------------

# beuti xr.ds
ds_beuti = xr.open_dataset(fn_beuti)

# get xr.da with sst and beuti data
da_beuti = ds_beuti[var_beuti[0]]

# beuti anom
in_clim_beuti = np.logical_and(
    da_beuti.time.dt.year >= yr_clim_bgn, da_beuti.time.dt.year <= yr_clim_end)
beuti_clim = da_beuti[:, in_clim_beuti].groupby('time.month').mean('time')
beuti_anom = da_beuti.groupby('time.month') - beuti_clim

# x limit
yr_bgn = yr_end - num_xyrs + 1
x_bgn = '{}-{:02d}'.format(yr_bgn, month_bgn)
x_end = '{}-{:02d}'.format(yr_end, month_end)
xlm = [np.datetime64(x_bgn), np.datetime64(x_end)]

# vlines
vln = np.zeros(num_xyrs, dtype='datetime64[M]')
for i in range(0, num_xyrs-1):
    vln[i] = '{}-01'.format(yr_bgn+1+i)

# setup subplots, spacing and figure size
plt.close()
gs1 = gridspec.GridSpec(num_row, num_clmn)
gs1.update(left=0.05, right=0.85, bottom=0.05, top=0.9, wspace=0.1, hspace=0.1)
fig = plt.figure(figsize=(fig_wdth, fig_hght))

# ------------------------------------------------------------
# contour beuti anom
# ------------------------------------------------------------
# setup subplots, spacing and figure size
plt.close()
gs1 = gridspec.GridSpec(num_row, num_clmn)
gs1.update(left=0.05, right=0.85, bottom=0.05, top=0.9, wspace=0.1, hspace=0.1)
fig = plt.figure(figsize=(fig_wdth, fig_hght))



ax = fig.add_subplot(gs1[0])
x_beuti = beuti_anom[dim_da[1]].data.astype('datetime64[M]')
y_beuti = beuti_anom[dim_da[0]].data

# contour
CS = plt.contour(x_beuti, y_beuti, beuti_anom.data, nlvl_beuti1,
                 vmin=dmin_beuti, vmax=dmax_beuti, colors='gray', linewidths=0.5)
# Define a class that forces representation of float to look a certain way
# This remove trailing zero so '1.0' becomes '1'

# Recast levels to new class
CS.levels = [nf(val) for val in CS.levels]

labels1 = plt.clabel(CS, CS.levels, inline=False, fmt='%r',
                     fontsize=7, colors='k', rightside_up=True)

# rotate and change integers (ie 1.0) to whole number (ie 1)
# for l in labels1:
#     txt1 = l.get_text()
#     if float(txt1).is_integer():
#         l.set_text(txt1.split('.')[0])
#     l.set_rotation(0)

# contourf
dmin = np.ceil(np.nanmin(beuti_anom.data))
dmax = np.ceil(np.nanmax(beuti_anom.data))
lvl1 = np.arange(dmin, dmax+0.5, 0.5)

chck_min = 0
chck_max = 0
if dmin < dmin_beuti:
    chck_min = 1
    extnd1 = 'min'
if dmax > dmax_beuti:
    chck_max = 1
    extnd1 = 'max'
if chck_max+chck_min==2:
    extnd1 = 'both'

# dminmax = np.min([np.abs(dmin), dmax])
# dmin = -1*dminmax
# dmax = dminmax


plt.contourf(x_beuti, y_beuti, beuti_anom.data, lvl1,
             vmin=dmin_beuti, vmax=dmax_beuti, cmap='bwr')

# xtick labels
mon_lbl = list()
for i in x_beuti:
    dti = pd.to_datetime(i)
    moni = calendar.month_name[dti.month]
    if moni == 'January':
        mon_lbl.append('{}\n      {}'.format(moni[0], dti.year))
    else:
        mon_lbl.append(moni[0])
plt.xticks(x_beuti, mon_lbl, fontsize=7)

# ytick labesl
lat_lbl = list()
for i in y_tck:
    lat_lbl.append('{}$\degree$N'.format(i))
plt.yticks(y_tck, lat_lbl, fontsize=7)

# x, y limits
plt.xlim(xlm)
plt.ylim(ylm)

# vlines
plt.vlines(vln, ylm[0], ylm[1], colors='red', linestyles='dashed', linewidth=2)

# colorbar
ax_pos = ax.get_position()
x_cb = ax_pos.x0 + ax_pos.width + ax_pos.width/40.0
y_cb = ax_pos.y0
y_hght = ax_pos.height
cbaxes = plt.gcf().add_axes([x_cb, y_cb+(y_hght*0.1)/2.0, 0.005, y_hght*0.9])

m = plt.cm.ScalarMappable(cmap='bwr')
m.set_array(beuti_anom.data)
m.set_clim(dmin_beuti, dmax_beuti)
plt.colorbar(m, cax=cbaxes, label=clrbr_lbl_beuti,
             format='%4.1f', extend=extnd1,
             boundaries=np.arange(dmin_beuti, dmax_beuti+d_beuti, d_beuti))


# save figure
fn_fig = '{}contour_beuti_anom_lat_{}_{}.png'.format(dir_out, lat_wnt[0], lat_wnt[1])

plt.savefig(fn_fig)
