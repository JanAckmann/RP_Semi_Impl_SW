# -*- coding: utf-8 -*-
#module load python/2.7-ve3


import os as os
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
#from pathlib import Path
from scipy.linalg import norm


from matplotlib.colors import LinearSegmentedColormap
import matplotlib

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)


path=  'data_ADI_Precon/'

# res4



#path_true='../data/data_ADI_LinfPhiEXIT1M10_Dp_refInst_dt200_res8/'
#path_ref='../data/data_ADI_LinfPhiEXIT1M10_Dp_refInst_dt200_res4/'
pathlist= ['../data/sheusp_DP_L2Exit_1M3_dt200_res4/'] # '../data/data_ADI_NOTgcr_D6_2M4_RP_res4/'
#pathlist= ['../data/sheusp_IMPR_SP_L2Exit_1M3_dt200_res4_speed/'] # '../data/data_ADI_NOTgcr_D6_2M4_RP_res4/'



#path_refRP='../data/data_ADI_NOTgcr_D2_3M4_RP_Ext_res2/'
#path_refRPD0='../data/data_ADI_NOTgcr_D0_3M4_RP_Ext_res2/'
#path_true='../data/data_ADI_NOTgcr_EXITcond1M10_exp2_res8/'
#path_ref='../data/data_ADI_NOTgcr_EXITcond3M4_Dp_res2/'#'../data/data_ADI_NOTgcr_D6_2M4_RP_res4/'##'data_ADI_Precon/' #'data_ADI_Precon/Precon' ##'../#explicit/data/' '../data/data_ADI_NOTgcr_D6_2M4_RP_res4/'#
#pathlist= ['../data/data_ADI_NOTgcr_EXITcond10M10_Dp_res2/'] # '../data/data_ADI_NOTgcr_D6_2M4_RP_res4/'
#res_diff=4
#latitude=0 #119

#path_refRP='../data/data_ADI_NOTgcr_D0_1M3_RP_Ext_res1/'
#path_refRPD0='../data/data_ADI_NOTgcr_D0_1M3_RP_Ext_res1/'
#path_true='../data/data_ADI_NOTgcr_EXITcond1M10_exp2_res8/'
#path_ref='../data/data_ADI_NOTgcr_EXITcond1M3_Dp_res1/'#'../data/data_ADI_NOTgcr_D6_2M4_RP_res4/'##'data_ADI_Precon/' #'data_ADI_Precon/Precon' ##'../#explicit/data/' '../data/data_ADI_NOTgcr_D6_2M4_RP_res4/'#
#pathlist= ['../data/data_ADI_NOTgcr_EXITcond10M10_Dp_res1/'] # '../data/data_ADI_NOTgcr_D6_2M4_RP_res4/'
#res_diff=8
#latitude=0 #59

varlist=['R']
var='R'
exp='3'
codesQ='F'
codesD='F'

res=4
# timestep 200 seconds
times=[0]
tprint= 797.0*400.0 # in seconds
# timestep 270 seconds
#times=[0, 44, 88, 132, 177, 221, 265, 309, 354]
#tprint= 797.0*(200.0/270.0)*270.0 # in seconds

cmap = plt.get_cmap('nipy_spectral')
new_cmap = truncate_colormap(cmap, 0.15, 1.0)
  
############## MAIN PROGRAMM ###################### 
for path in pathlist: 

    for Precon in range(33,34,5): 

        pdf_pages = PdfPages(path+'P0_'+'Precon'+str(Precon)+'_exp'+str(exp)+'_codes_'+codesQ+codesD+'DP.pdf')
        for P_steps in range(len(times)):
          time=times[P_steps]
          plotnum=0
          fig = plt.figure(figsize=(8.27, 11.69), dpi=100)
          plt.suptitle('TIME in days: '+str(P_steps*tprint/(3600.0*24.0)))
          for var in varlist: 
            plotnum=plotnum+1

            filename =path+'Precon'+str(Precon)+'_' + var+'_exp'+str(exp)+'_time'+str(time)+'_iter_0_codes_'+codesQ+codesD+'_bits'+str(23)+'.txt'
            #filename_ref =path+'Precon'+str(Precon)+'_' + var+'_exp'+str(exp)+'_time'+str(0)+'_codes_'+codesQ+codesD+'_bits'+str(23)+'.txt'
            #filenameU =path+'Precon'+str(Precon)+'_U_exp'+str(exp)+'_time'+str(time)+'_codes_'+codesQ+codesD+'_bits'+str(23)+'.txt'
            #filenameV =path+'Precon'+str(Precon)+'_V_exp'+str(exp)+'_time'+str(time)+'_codes_'+codesQ+codesD+'_bits'+str(23)+'.txt'

            xcoord, ycoord, mean_var=np.loadtxt( filename, usecols=(0,1,2), unpack=True)
            #xcoord, ycoord, mean_varU=np.loadtxt( filenameU, usecols=(0,1,2), unpack=True)
            #xcoord, ycoord, mean_varV=np.loadtxt( filenameV, usecols=(0,1,2), unpack=True)
            #xcoord, ycoord, mean_var_ref=np.loadtxt( filename_ref, usecols=(0,1,2), unpack=True)

            ncols, nrows = len(set(xcoord)), len(set(ycoord)) 
       
            grid_var = (mean_var.reshape((nrows, ncols), order='F'))
            #grid_varU = (mean_varU.reshape((nrows, ncols), order='F'))
            #grid_varV = (mean_varV.reshape((nrows, ncols), order='F'))
            #grid_var_ref = (mean_var_ref.reshape((nrows, ncols), order='F'))

            grid_xcoord= (xcoord.reshape((nrows, ncols), order='F'))
            grid_ycoord= (ycoord.reshape((nrows, ncols), order='F'))


            #if (exp == '1'):
            #  grid_var=(grid_var-8000.)/8000.
            #elif (exp == '3'):
            #  grid_var=(grid_var-grid_var_ref)/8000.
            #  print(np.amin(grid_var),np.amax(grid_var))

            plotid='31'+str(plotnum)
            plt.subplot(plotid)
            plt.title(var, loc='left')
            plt.xlabel('lon', fontsize=18)
            plt.ylabel('lat', fontsize=18)
            #if (var == 'V'):
            plt.imshow(grid_var, extent=(xcoord.min(), xcoord.max(), 
               ycoord.min(), ycoord.max()),
               interpolation='None', cmap=new_cmap, vmin=np.amin((grid_var)), vmax=np.amax(abs(grid_var)))
          plt.tight_layout()
          pdf_pages.savefig(fig)
          plt.close()
        pdf_pages.close()
        print(path+'Precon'+str(Precon)+'_' + var+'_exp'+str(exp)+'_time'+str(time)+'_codes_'+codesQ+codesD+'_bits'+str(23)+'.pdf')



