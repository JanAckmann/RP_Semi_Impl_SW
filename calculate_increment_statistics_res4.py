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
#pathlist= ['../data/sheusp_SP_L2Exit_1M3_dt200_res4/'] # '../data/data_ADI_NOTgcr_D6_2M4_RP_res4/'
#pathlist= ['../data/sheusp_IMPR_SP_L2Exit_1M3_dt200_res4/'] # '../data/data_ADI_NOTgcr_D6_2M4_RP_res4/'

codesQ='F'
codesD='F'


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
exp_list=['1', '3']
Precon='2'
iteration=4
res=4
# timestep 200 seconds
#times=[0,48,96,144,192,240,288,336]
times=[i for i in range(0,360,12)]
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

def Linf_norm2D(field, n1, n2): 
    Linf2D=0.0
    for i in range(n1):
        if(Linf2D<max(abs(field[:,i]))):
            Linf2D=max(abs(field[:,i]))
    return Linf2D

def L2_norm(field):
   
   L2=np.mean(abs(field))
   return L2

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
for path in pathlist: 
 for exp in exp_list:
   for time in times:
    print(exp)
    filename =path+'Precon'+str(Precon)+'_' + 'R_exp'+str(exp)+'_time'+str(time)+'_iter_'+str(iteration)+'_codes_'+codesQ+codesD+'_bits'+str(23)+'.txt'
    filename2 =path+'Precon'+str(Precon)+'_' + 'R_exp'+str(exp)+'_time'+str(time)+'_iter_'+str(iteration-1)+'_codes_'+codesQ+codesD+'_bits'+str(23)+'.txt'
    filename3 =path+'Precon'+str(Precon)+'_' + 'R_exp'+str(exp)+'_time'+str(time)+'_iter_'+str(iteration-2)+'_codes_'+codesQ+codesD+'_bits'+str(23)+'.txt'
    filename4 =path+'Precon'+str(Precon)+'_' + 'R_exp'+str(exp)+'_time'+str(time)+'_iter_'+str(iteration-3)+'_codes_'+codesQ+codesD+'_bits'+str(23)+'.txt'
    filename5 =path+'Precon'+str(Precon)+'_' + 'R_exp'+str(exp)+'_time'+str(time)+'_iter_'+str(iteration-4)+'_codes_'+codesQ+codesD+'_bits'+str(23)+'.txt'
    print(filename)  
    print(filename5)  
    xcoord, ycoord, R=np.loadtxt( filename, usecols=(0,1,2), unpack=True)
    xcoord, ycoord, R2=np.loadtxt( filename2, usecols=(0,1,2), unpack=True)
    xcoord, ycoord, R3=np.loadtxt( filename3, usecols=(0,1,2), unpack=True)
    xcoord, ycoord, R4=np.loadtxt( filename4, usecols=(0,1,2), unpack=True)
    xcoord, ycoord, R5=np.loadtxt( filename5, usecols=(0,1,2), unpack=True)
    print(R)
    ncols, nrows = len(set(xcoord)), len(set(ycoord))
    #print(sorted(set(ycoord)))

    grid_varR = (R.reshape((nrows, ncols), order='F'))   
    
    print(L2_norm(R/9.81))# ,Linf_norm(kinetic-kinetic_gen)/ Linf_norm(kinetic_gen))
    print(Linf_norm(R/9.81))# ,Linf_norm(kinetic-kinetic_gen)/ Linf_norm(kinetic_gen))
    print('solver increment-1')
    print(L2_norm((R-R2)/9.81))# ,Linf_norm(kinetic-kinetic_gen)/ Linf_norm(kinetic_gen))
    print(Linf_norm((R-R2)/9.81))# ,Linf_norm(kinetic-kinetic_gen)/ Linf_norm(kinetic_gen))
    print('solver increment-2')
    print(L2_norm((R2-R3)/9.81))# ,Linf_norm(kinetic-kinetic_gen)/ Linf_norm(kinetic_gen))
    print(Linf_norm((R2-R3)/9.81))# ,Linf_norm(kinetic-kinetic_gen)/ Linf_norm(kinetic_gen))
    print('solver increment-3')
    print(L2_norm((R3-R4)/9.81))# ,Linf_norm(kinetic-kinetic_gen)/ Linf_norm(kinetic_gen))
    print(Linf_norm((R3-R4)/9.81))# ,Linf_norm(kinetic-kinetic_gen)/ Linf_norm(kinetic_gen))
    print('solver increment-4')
    print(L2_norm((R4-R5)/9.81))# ,Linf_norm(kinetic-kinetic_gen)/ Linf_norm(kinetic_gen))
    print(Linf_norm((R4-R5)/9.81))# ,Linf_norm(kinetic-kinetic_gen)/ Linf_norm(kinetic_gen))
