# -*- coding: utf-8 -*-
#module load python/2.7-ve3


import os as os
import os.path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
from matplotlib import rcParams

from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator

#path='../data/sheusp_IMPR_SP_L2Exit_1M3_dt200_res4/'
#path='../data/sheusp_IMPR_SP_plusLat_L2Exit_1M3_dt200_res4/'
#path='../data/sheusp_IMPR_SP_RPSolver_DP3_L2Exit_1M3_dt200_res4/'
#path='../data/sheusp_IMPR_SP_RPGCR_Prec_DP128_L2Exit_1M3_dt200_res4/'
path='../data/sheusp_IMPR_SP_FRPPreAll_FRPLap_FullSPr0_DP2_L2Exit_1M3_dt200_res4/'
#path='../data/sheusp_DP_L2Exit_1M3_dt200_res4/'
#path='../data/sheusp_SP_L2Exit_1M3_dt200_res4/'
#path='../data/sheusp_SP_L2Exit_1M3_dt200_noabsorb_res4/'
#path='../data/sheusp_SP_L2Exit_1M3_dt200_LAbsorb_res4/'
#path='../data/sheusp_SP_L2Exit_1M3_dt200_HAbsorb_res4/'
#path='../data/sheusp_DP_L2Exit_1M3_dt200_NOAbsorb_res4/'
#path='../data/sheusp_DP_L2Exit_1M3_dt200_LAbsorb_res4/'
#path='../data/sheusp_DP_L2Exit_1M3_dt200_HAbsorb_res4/'
#path_ref= '../data/sheusp_DP_L2Exit_1M3_dt200_HAbsorb_res4/'
plotpath= path

varlist=['R']
explist=['1' , '3']
if (path=='../data/sheusp_IMPR_SP_L2Exit_1M3_dt200_res4/' or path=='../data/sheusp_IMPR_SP_plusLat_L2Exit_1M3_dt200_res4/' ):
  Precon='9'
  codesD='T'
elif (path=='../data/sheusp_IMPR_SP_FRPPreAll_FRPLap_FullSPr0_DP2_L2Exit_1M3_dt200_res4/'):
  Precon='23'
  codesD='T'
else:
  Precon='7'
  codesD='F'

codesQ='F'

iteration=0
solver='implicit_HP_implPD'

start_time= 0
end_time=  360#648#336#648

jump=48


font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)
rcParams.update({'figure.autolayout': True})

import numpy as np


def calc_average(in_file, start, end, jump):
  file_timeslice= in_file+str(start)
  counter=1
  x,y,average = np.loadtxt(file_timeslice, usecols=(0,1,2), unpack=True)

  for timestep in range(start,end+1, jump):
    print( timestep)
    file_timeslice= filename+str(timestep)
    counter+=1
    x,y,temp = np.loadtxt(file_timeslice, usecols=(0,1,2), unpack=True)
    average+= temp
  
  return x, y, average/float(counter)
  
def calc_std(in_file, average, start, end, jump):
  counter=0
  std=np.zeros(len(average), dtype=float)

  for timestep in range(start,end+1, jump):
    file_timeslice= filename+str(timestep)
    counter+=1
    x,y,temp = np.loadtxt(file_timeslice, usecols=(0,1,2), unpack=True)
    std+= (temp-average)**2.0
  
  return x, y, np.sqrt(std/float(counter-1))
TIMESTEPS=[0,48,96,144,192,240,288,336]

#TIMESTEPS=[0,24 ,48,72,96,120,144 ,168 ,192,216,240,264,288,312,336]
############## MAIN PROGRAMM ######################  

number=len(TIMESTEPS) #(end_time-start_time)/11
cmap = plt.get_cmap('viridis')
colors = [cmap(i) for i in np.linspace(1, 0,int( number)+1)]
for exp in explist:
 mean_full_residuum_norm=0.0
 if (iteration==4):
  if(exp=='1'):  
    y_min = 10**(-8)
    y_max = 10**(0)
  elif(exp=='3'):
    y_min = 10**(-6)
    y_max = 10**(-1)
 if (iteration==0):
  if(exp=='1'):  
    y_min = 10**(-12)
    y_max = 10**(1)
  elif(exp=='3'):
    y_min = 10**(-12)
    y_max = 10**(1)
 for var in varlist: 
  for bits in range(23,22,-2):
    for time in TIMESTEPS:# range(start_time, end_time+1,jump):
      print(time)
      residuum_norm= np.zeros((256))
      residuum_norm_sp= np.zeros((256))
      number_of_iters=np.zeros((256))#np.zeros((iteration))
      for iteration in range(iteration,iteration+1,1):

        filenameT = path+'Precon'+Precon+'_' + var+'_exp'+exp+'_time'+str(time)+'_iter_'+str(iteration)+'_codes_'+str(codesD)+str(codesQ)+'_bits'+str(bits)+'.txt'
        filenameT_sp = path+'Precon9_' + var+'_exp'+exp+'_time'+str(time)+'_iter_'+str(iteration)+'_codes_'+str(codesD)+str(codesQ)+'_bits'+str(bits)+'.txt'
        #filenameT_re = path_ref+'Precon'+Precon+'_' + var+'_exp'+exp+'_time'+str(time)+'_iter_'+str(iteration)+'_codes_'+str(codesD)+str(codesQ)+'_bits'+str(bits)+'.txt'
        #print 'calc_average'
        print(filenameT)

        xcoord, ycoord, Residuum, exitcon   =np.loadtxt( filenameT, usecols=(0,1,2,3), unpack=True)
        xcoord, ycoord, Residuum_sp, exitcon   =np.loadtxt( filenameT_sp, usecols=(0,1,2,3), unpack=True)
        ncols, nrows = len(set(xcoord)), len(set(ycoord))

        grid_residuum = (Residuum.reshape((nrows, ncols), order='F'))
        grid_residuum_sp = (Residuum_sp.reshape((nrows, ncols), order='F'))
        full_residuum_norm=np.sqrt(np.sum (grid_residuum[1,:]**2)) 
        full_residuum_norm_sp=np.sqrt(np.sum (grid_residuum_sp[1,:]**2))
        print(abs(full_residuum_norm-full_residuum_norm_sp)/abs(full_residuum_norm))
        mean_full_residuum_norm=mean_full_residuum_norm+(1.0/float(len(TIMESTEPS)))*abs(full_residuum_norm-full_residuum_norm_sp)/abs(full_residuum_norm)
        print('MEAN error in residuum',mean_full_residuum_norm)
        print (nrows)
        for lat in range(nrows):
          residuum_norm[lat]= np.sqrt(np.sum (grid_residuum[lat,:]**2)) #np.mean (np.abs(Residuum))#np.linalg.norm (Residuum)
          residuum_norm_sp[lat]= np.sqrt(np.sum (grid_residuum_sp[lat,:]**2)) #np.mean (np.abs(Residuum))#np.linalg.norm (Residuum)
        print(abs(residuum_norm- residuum_norm_sp)/abs(residuum_norm))
        number_of_iters[:]=range(0,256) #np.zeros((iteration))
        number_of_iters=((number_of_iters)-128.0)/256.0 #np.zeros((iteration))
        plt.plot(number_of_iters, abs(residuum_norm[:]-residuum_norm_sp[:])/abs(residuum_norm[:]), color=colors[int((time-start_time)/jump)], label ='day'+str(1+int(time/24)))
    plt.xlabel('Y/Pi')
    plt.ylabel(r'Rel. Error of $r_0$ in L2 norm ')
    #plt.plot((0.0,float(max_iteration) ), (1.0,1.0), 'k-')
    print(path)
    print(path=='../data/sheusp_SP_L2Exit_1M3_dt200_res4/')
    if (path=='../data/sheusp_SP_L2Exit_1M3_dt200_res4/'):
     print(path)
     plt.legend(loc='center', prop={'size':14}, ncol=2, bbox_to_anchor=(0.6,0.91), framealpha=1)
    else:
     plt.legend(loc='center', prop={'size':14}, ncol=2, bbox_to_anchor=(0.6,0.91), framealpha=1)
    #plt.legend(loc='upper left', prop={'size':6}, bbox_to_anchor=(1,1)) #(loc='upper right')  #['day0'], loc='upper right', 'day1','day2', 'day3','day4', 'day5','day6', 'day7', 'day8', 'day9','day10', 'day11','day12', 'day13','day14', 'day15', 'day16'
    #plt.title(r'Evolution of Residuals')
    plt.yscale('log')
    plt.ylim(ymin = y_min,ymax = y_max)
    plt.axis
    #plt.show()
    #plt.savefig(plotpath +'L2_Precon'+Precon+'_' + var+'_exp'+exp+'_codes_FF'+'_bits'+str(bits)+'.pdf', bbox_inches=0)
    #plt.close()
  
    #print( plotpath +'L2_Precon'+Precon+'_' + var+'_exp'+exp+'_codes_FF'+'_bits'+str(bits)+'.pdf')

    plt.savefig(plotpath +'NumErr1_Lat_Precon'+Precon+'_' + var+'_exp'+exp+'_codes_'+str(codesD)+str(codesQ)+'_bits'+str(bits)+'_iter'+str(iteration)+'L2norm.pdf', bbox_inches=0)
    plt.close()
  
    print( plotpath +'NumErr1_Lat_Precon'+Precon+'_' + var+'_exp'+exp+'_codes_'+str(codesD)+str(codesQ)+'_bits'+str(bits)+'_iter'+str(iteration)+'L2norm.pdf')


#compare_histograms_normal('h', 31, sigma_h) #, -0.0526478751975,0.0526478751975 )
#compare_histograms_talone('t', 31, sigma_t, -0.0280507825342,0.0280507825342 )
#compare_histograms('vn', 31, sigma_vn,-0.012620211929 , 0.012620211929)

#call_truncation_errors('h',n_cells, 1, num_timesteps)
#call_truncation_errors('t',n_cells, 32, num_timesteps)
#call_truncation_errors('vn',n_edges, 32, num_timesteps)


