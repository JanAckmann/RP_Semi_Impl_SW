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
        'size'   : 13}

matplotlib.rc('font', **font)


path=  'data_ADI_Precon/'

# res4



#path_true='../data/data_ADI_LinfPhiEXIT1M10_Dp_refInst_dt200_res8/'
#path_ref='../data/data_ADI_LinfPhiEXIT1M10_Dp_refInst_dt200_res4/'
#pathlist= ['../data/sheusp_DP_L2Exit_1M3_dt200_res4/'] # '../data/data_ADI_NOTgcr_D6_2M4_RP_res4/'
pathlist= ['../data_RP_Semi_Impl_SW/sheusp_SP_L2Exit_1M3_dt200_res4/'] # '../data/data_ADI_NOTgcr_D6_2M4_RP_res4/'
codesQ='F'
codesD='F'
#pathlist= ['../data/sheusp_SP_1iteration_L2Exit_1M3_dt200_res4/'] # '../data/data_ADI_NOTgcr_D6_2M4_RP_res4/'
#pathlist= ['../data/sheusp_IMPR_SP_L2Exit_1M3_dt200_res4/'] # '../data/data_ADI_NOTgcr_D6_2M4_RP_res4/'
path_ref='../data_RP_Semi_Impl_SW/sheusp_DP_L2Exit_1M3_dt200_res4/'
#codesQ='F'
#codesD='F'


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

varlist=['H']
var='H'
exp_list=['1', '3']
Precon='7'
res=4
# timestep 200 seconds
times=[354]# # 354]
tprint= 797.0*200.0 # in seconds
# timestep 270 seconds
#times=[0, 44, 88, 132, 177, 221, 265, 309, 354]
#tprint= 797.0*(200.0/270.0)*270.0 # in seconds

cmap = plt.get_cmap('nipy_spectral')
new_cmap = truncate_colormap(cmap, 0.15, 1.0)
  
def L2_norm2D(field,n1, n2):
    L22D=0.0
    for i in range(n1):
        for j in range(n2):
          L22D+=((field[j,i]))
    return np.sqrt(L22D)

def min2D(field, n1, n2): 
    min2D=field[0,0]
    for i in range(n1):
        if(min2D>min((field[:,i]))):
            min2D=min((field[:,i]))
    return min2D

def max2D(field, n1, n2): 
    max2D=0.0
    for i in range(n1):
        if(max2D<max((field[:,i]))):
            max2D=max((field[:,i]))
    return max2D

def Linf_norm2D(field, n1, n2): 
    Linf2D=0.0
    for i in range(n1):
        if(Linf2D<max(abs(field[:,i]))):
            Linf2D=max(abs(field[:,i]))
    return Linf2D

def L2_norm(field):
   
   L2=np.sum((field))
   return np.sqrt(L2)

def Linf_norm(field): 
   Linf=max(abs(field))
   return Linf

def Pot_Energy(fieldH,fieldU, fieldV,lats, n1, n2):
    pot_energy=np.zeros((n2, n1))
    for i in range(n1):
        for j in range(n2):
          R=6371.22E+03
          x_len=R*np.cos(lats[j])
          #print(x_len, R, n1, lats[j], np.cos(lats[j]))
          y_len=R
          Volume=x_len*y_len
          #Mass=1.0 #Volume*1.225
          #print(Mass)
          pot_energy[j,i]=Volume*(fieldH[j,i])**2 
    return pot_energy

def Kinetic_Energy_noVol(fieldH,fieldU, fieldV,lats, n1, n2):
    kin_energy=np.zeros((n2, n1))
    for i in range(n1):
        for j in range(n2):
          kin_energy[j,i]=((fieldU[j,i])**2 + (fieldV[j,i])**2)
    return kin_energy

def Kinetic_Energy(fieldH,fieldU, fieldV,lats, n1, n2):
    kin_energy=np.zeros((n2, n1))
    for i in range(n1):
        for j in range(n2):
          R=6371.22E+03
          x_len=R*np.cos(lats[j])
          #print(x_len, R, n1, lats[j], np.cos(lats[j]))
          y_len=R
          Volume=x_len*y_len
          #Mass=1.0 #Volume*1.225
          #print(Mass)
          kin_energy[j,i]=Volume*((fieldU[j,i])**2 + (fieldV[j,i])**2)**2
    return kin_energy
############## MAIN PROGRAMM ###################### 
codesQ_save=codesQ
codesD_save=codesD
for path in pathlist: 
 for time in times:
  for exp in exp_list:
    codesQ=codesQ_save
    codesD=codesD_save
    print(exp)
    filename =path+'Precon'+str(Precon)+'_' + 'H_exp'+str(exp)+'_time'+str(time)+'_codes_'+codesQ+codesD+'_bits'+str(23)+'.txt'
    filename_gen =path+'Precon'+str(Precon)+'_' + 'H_exp'+str(exp)+'_time'+str(0)+'_codes_'+codesQ+codesD+'_bits'+str(23)+'.txt'
    filenameU =path+'Precon'+str(Precon)+'_U_exp'+str(exp)+'_time'+str(time)+'_codes_'+codesQ+codesD+'_bits'+str(23)+'.txt'
    filenameU_gen =path+'Precon'+str(Precon)+'_' + 'U_exp'+str(exp)+'_time'+str(0)+'_codes_'+codesQ+codesD+'_bits'+str(23)+'.txt'
    filenameV =path+'Precon'+str(Precon)+'_V_exp'+str(exp)+'_time'+str(time)+'_codes_'+codesQ+codesD+'_bits'+str(23)+'.txt'
    filenameV_gen =path+'Precon'+str(Precon)+'_' + 'V_exp'+str(exp)+'_time'+str(0)+'_codes_'+codesQ+codesD+'_bits'+str(23)+'.txt'
  
    xcoord, ycoord, H=np.loadtxt( filename, usecols=(0,1,2), unpack=True)
    xcoord, ycoord, H_gen=np.loadtxt( filename_gen, usecols=(0,1,2), unpack=True)
  
    xcoord, ycoord, U=np.loadtxt( filenameU, usecols=(0,1,2), unpack=True)
    xcoord, ycoord, U_gen=np.loadtxt( filenameU_gen, usecols=(0,1,2), unpack=True)
    
    xcoord, ycoord, V=np.loadtxt( filenameV, usecols=(0,1,2), unpack=True)
    xcoord, ycoord, V_gen=np.loadtxt( filenameV_gen, usecols=(0,1,2), unpack=True)
    ncols, nrows = len(set(xcoord)), len(set(ycoord))
    #print(sorted(set(ycoord)))

    grid_varH = (H.reshape((nrows, ncols), order='F'))   
    grid_varH_gen = (H_gen.reshape((nrows, ncols), order='F'))   

    grid_varU = (U.reshape((nrows, ncols), order='F'))   
    grid_varU_gen = (U_gen.reshape((nrows, ncols), order='F'))   

    grid_varV = (V.reshape((nrows, ncols), order='F'))   
    grid_varV_gen = (V_gen.reshape((nrows, ncols), order='F'))   

    #pot=Pot_Energy(grid_varH-grid_varH_gen, grid_varU-grid_varU_gen, grid_varV-grid_varV_gen, sorted(set(ycoord)), ncols, nrows)
    #pot_gen=Pot_Energy(grid_varH_gen, grid_varU-grid_varU_gen, grid_varV-grid_varV_gen, sorted(set(ycoord)), ncols, nrows)
    #print('pot/pot_gen')
    #print((L2_norm2D(pot, ncols, nrows)/L2_norm2D(pot_gen, ncols, nrows)))# ,Linf_norm(kinetic-kinetic_gen)/ Linf_norm(kinetic_gen))
    #print('inf')
    #print(Linf_norm(H-H_gen)/Linf_norm(H_gen))# ,Linf_norm(kinetic-kinetic_gen)/ Linf_norm(kinetic_gen))

    #print('kinetic energy')
    kinetic=Kinetic_Energy_noVol(abs(grid_varH-grid_varH_gen), grid_varU-grid_varU_gen, grid_varV-grid_varV_gen, sorted(set(ycoord)), ncols, nrows)
    kinetic_gen=Kinetic_Energy_noVol(grid_varH_gen, grid_varU_gen, grid_varV_gen, sorted(set(ycoord)), ncols, nrows)
    kinetic=np.sqrt(kinetic)
    kinetic_gen=np.sqrt(kinetic_gen)
    maximum=(max2D(kinetic, ncols, nrows)) #, max2D(kinetic_gen, ncols, nrows))    
    minimum=(min2D(kinetic, ncols, nrows)) #, min2D(kinetic_gen, ncols, nrows))    

    xcoord[:]=xcoord[:]/np.pi
    ycoord[:]=ycoord[:]/np.pi
    plt.xlabel('x/Pi', fontsize=15)
    plt.ylabel('y/Pi', fontsize=15)
    plt.imshow(np.flipud(kinetic),extent=(xcoord.min(), xcoord.max(),
               ycoord.min(), ycoord.max()), vmin=0, vmax=maximum,
              interpolation='nearest', cmap=cm.jet)
    clb=plt.colorbar(shrink=0.7)
    clb.set_label(r'$\frac{m}{s}$',fontsize=18, rotation=0)
    plt.savefig('kinetic'+str(exp)+'_SPpertDP.png')
    plt.close()

    


    codesQ='F'
    codesD='F'

    filename =path_ref+'Precon'+str(Precon)+'_' + 'H_exp'+str(exp)+'_time'+str(time)+'_codes_'+codesQ+codesD+'_bits'+str(23)+'.txt'
    filename_gen =path_ref+'Precon'+str(Precon)+'_' + 'H_exp'+str(exp)+'_time'+str(0)+'_codes_'+codesQ+codesD+'_bits'+str(23)+'.txt'
    filenameU =path_ref+'Precon'+str(Precon)+'_U_exp'+str(exp)+'_time'+str(time)+'_codes_'+codesQ+codesD+'_bits'+str(23)+'.txt'
    filenameU_gen =path_ref+'Precon'+str(Precon)+'_' + 'U_exp'+str(exp)+'_time'+str(0)+'_codes_'+codesQ+codesD+'_bits'+str(23)+'.txt'
    filenameV =path_ref+'Precon'+str(Precon)+'_V_exp'+str(exp)+'_time'+str(time)+'_codes_'+codesQ+codesD+'_bits'+str(23)+'.txt'
    filenameV_gen =path_ref+'Precon'+str(Precon)+'_' + 'V_exp'+str(exp)+'_time'+str(0)+'_codes_'+codesQ+codesD+'_bits'+str(23)+'.txt'
  
    xcoord, ycoord, H=np.loadtxt( filename, usecols=(0,1,2), unpack=True)
    xcoord, ycoord, H_gen=np.loadtxt( filename_gen, usecols=(0,1,2), unpack=True)
  
    xcoord, ycoord, U=np.loadtxt( filenameU, usecols=(0,1,2), unpack=True)
    xcoord, ycoord, U_gen=np.loadtxt( filenameU_gen, usecols=(0,1,2), unpack=True)
    
    xcoord, ycoord, V=np.loadtxt( filenameV, usecols=(0,1,2), unpack=True)
    xcoord, ycoord, V_gen=np.loadtxt( filenameV_gen, usecols=(0,1,2), unpack=True)
    ncols, nrows = len(set(xcoord)), len(set(ycoord))
    #print(sorted(set(ycoord)))

    grid_varH = (H.reshape((nrows, ncols), order='F'))   
    grid_varH_gen = (H_gen.reshape((nrows, ncols), order='F'))   

    grid_varU = (U.reshape((nrows, ncols), order='F'))   
    grid_varU_gen = (U_gen.reshape((nrows, ncols), order='F'))   

    grid_varV = (V.reshape((nrows, ncols), order='F'))   
    grid_varV_gen = (V_gen.reshape((nrows, ncols), order='F'))   

    #pot=Pot_Energy(grid_varH-grid_varH_gen, grid_varU-grid_varU_gen, grid_varV-grid_varV_gen, sorted(set(ycoord)), ncols, nrows)
    #pot_gen=Pot_Energy(grid_varH_gen, grid_varU-grid_varU_gen, grid_varV-grid_varV_gen, sorted(set(ycoord)), ncols, nrows)
    #print('pot/pot_gen')
    #print((L2_norm2D(pot, ncols, nrows)/L2_norm2D(pot_gen, ncols, nrows)))# ,Linf_norm(kinetic-kinetic_gen)/ Linf_norm(kinetic_gen))
    #print('inf')
    #print(Linf_norm(H-H_gen)/Linf_norm(H_gen))# ,Linf_norm(kinetic-kinetic_gen)/ Linf_norm(kinetic_gen))

    #print('kinetic energy')
    kinetic_ref=Kinetic_Energy_noVol(abs(grid_varH-grid_varH_gen), grid_varU-grid_varU_gen, grid_varV-grid_varV_gen, sorted(set(ycoord)), ncols, nrows)
    kinetic_gen_ref=Kinetic_Energy_noVol(grid_varH_gen, grid_varU_gen, grid_varV_gen, sorted(set(ycoord)), ncols, nrows)

    xcoord[:]=xcoord[:]/np.pi
    ycoord[:]=ycoord[:]/np.pi
    kinetic_ref=abs(kinetic-np.sqrt(kinetic_ref))
    kinetic_gen_ref=abs(kinetic_ref-np.sqrt(kinetic_gen_ref))
    if (str(exp)=='3'):
        maximum=100#(max2D(kinetic_ref, ncols, nrows)) #, max2D(kinetic_gen, ncols, nrows))
    else:
        maximum=1
    plt.xlabel('x/Pi', fontsize=15)
    plt.ylabel('y/Pi', fontsize=15)

    plt.imshow(np.flipud(kinetic_ref),extent=((xcoord.min(), xcoord.max(),
               ycoord.min(), ycoord.max())), vmin=maximum*10**(-4), vmax=maximum,
              interpolation='nearest', norm=matplotlib.colors.LogNorm(),  cmap=cm.jet)
    clb=plt.colorbar(shrink=0.7)
    clb.set_label(r'$\frac{m}{s}$',fontsize=18, rotation=0)
    plt.savefig('kinetic_diff'+str(exp)+'_SPpertDP.png')
    plt.close()

