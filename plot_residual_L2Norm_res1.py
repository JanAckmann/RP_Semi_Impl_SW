# -*- coding: utf-8 -*-
#module load python/2.7-ve3


import os as os
import os.path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib

from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator

#path='../data/sheusp_IMPR_SP_L2Exit_1M3_dt200_res4/'
#path='../data/sheusp_DP_L2Exit_1M3_dt200_1Preconiter_res4/'
#path='../data/sheusp_IMPR_SP_L2Exit_1M3_dt200_res4/'
#path='../data/sheusp_IMPR_SP_RPGCR_Prec_DP16_L2Exit_1M3_dt200_res4/'
#path='../data/sheusp_IMPR_SP_RPPreAll_RPLap_FullSPr0_DP0_L2Exit_1M3_dt200_res4/'
#path='../data/sheusp_IMPR_SP_FRPPreAll_FRPLap_SPr0_DP2_L2Exit_1M3_dt200_res4/' # no single precision left in precon
path='../data/sheusp_IMPR_SP_FRPPreAll_FRPLap_SPr0V2_RPxAx_DP3_L2Exit_1M3_dt800_res1/' # + reduced precisision in conjugacy 
#path='../data/sheusp_IMPR_SP_FRPPreAll_FRPLap_SPr0V2_DP4_L2Exit_1M3_dt200_res4/' # + reduced precisision in conjugacy 
#path='../data/sheusp_IMPR_SP_plusLat_L2Exit_1M3_dt200_res4/'
Precon='23'
codesD='T'
codesQ='F'

#path='../data/sheusp_DP_L2Exit_1M3_dt800_res1/'
#path='../data/degraded_sheusp_DP_L2Exit_1M3_dt200_res4/'
#path='../data/sheusp_SP_L2Exit_1M3_dt200_res4/'
#path='../data/CrazyAbsorb_sheusp_SP_L2Exit_1M3_dt200_res4/'
#path='../data/sheusp_SP_L2Exit_1M3_dt200_noabsorb_res4/'
#path='../data/sheusp_SP_L2Exit_1M3_dt200_LAbsorb_res4/'
#path='../data/sheusp_SP_L2Exit_1M3_dt200_HAbsorb_res4/'
#path='../data/sheusp_DP_L2Exit_1M3_dt200_NOAbsorb_res4/'
#path='../data/sheusp_DP_L2Exit_1M3_dt200_LAbsorb_res4/'
#path='../data/sheusp_DP_L2Exit_1M3_dt200_HAbsorb_res4/'
#path_ref= '../data/sheusp_DP_L2Exit_1M3_dt200_HAbsorb_res4/'
#Precon='7'
#codesD='F'
#codesQ='F'

plotpath= path

varlist=['R']
#explist=['1' , '3']
explist=['1', '3']
#Precon='7'
#codesD='F'
#codesQ='F'

solver='implicit_HP_implPD'

start_time= 0
end_time=  360#648#336#648

jump=48


font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)


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

#TIMESTEPS=[0,24,48,72,96,120,144,168,192,216,240,264,288,312,336]
TIMESTEPS=[0,48,96,144,192,240,288,336]
  
#TIMESTEPS=[i for i in range(360,12)]
############## MAIN PROGRAMM ######################  

number=len(TIMESTEPS) #(end_time-start_time)/11
cmap = plt.get_cmap('viridis')
colors = [cmap(i) for i in np.linspace(1, 0,int( number)+1)]

for exp in explist:
 for var in varlist: 
  for bits in range(23,22,-2):
    #plt.subplot(2, 1, 1)
    errors= np.zeros((3, 17))
    count=0
    max_iteration=0
    for time in TIMESTEPS:# range(start_time, end_time+1,jump):
      plotted= False
      print(time)
      residuum_norm= np.zeros((100))
      residuum_infnorm= np.zeros((100))
      exit= np.zeros((100))
      for iteration in range(0,100,1):

        filenameT = path+'Precon'+Precon+'_' + var+'_exp'+exp+'_time'+str(time)+'_iter_'+str(iteration)+'_codes_'+str(codesD)+str(codesQ)+'_bits'+str(bits)+'.txt'
        #filenameT_re = path_ref+'Precon'+Precon+'_' + var+'_exp'+exp+'_time'+str(time)+'_iter_'+str(iteration)+'_codes_'+str(codesD)+str(codesQ)+'_bits'+str(bits)+'.txt'
        #print 'calc_average'
        #print(filenameT)

        #xcoord, ycoord, mean_var= calc_average(filename, start_time, end_time, jump)
        #print 'done'
        #print 'savetxt mean'
        #np.savetxt( filename+'avg'+str(start_time)+'to'+str(end_time), np.transpose([xcoord,ycoord,mean_var])) 

  

        #print 'calc_standard dev'


        #xcoord, ycoord, std_var= calc_std(filename, mean_var, start_time, end_time, 5000)
        #print 'done'
        #print 'savetxt std'
        #np.savetxt( filename+'std'+str(start_time)+'to'+str(end_time), np.transpose([xcoord,ycoord,std_var])) 


        ### plot std and mean


        #print( 'load residual')
        #xcoord, ycoord, mean_var_refF =np.loadtxt( file_refF, usecols=(0,1,2), unpack=True)
        #xcoord, ycoord, mean_var_refT =np.loadtxt( file_refT, usecols=(0,1,2), unpack=True)
        if (os.path.exists(filenameT)):
          xcoord, ycoord, Residuum, exitcon   =np.loadtxt( filenameT, usecols=(0,1,2,3), unpack=True)
          #xcoord, ycoord, Residuum_re, exitcon   =np.loadtxt( filenameT_re, usecols=(0,1,2,3), unpack=True)
          #xcoord, ycoord, mean_varF     =np.loadtxt( filenameF, usecols=(0,1,2), unpack=True)

          residuum_norm[iteration]= np.sqrt(np.sum (Residuum**2)) #np.mean (np.abs(Residuum))#np.linalg.norm (Residuum)
          print(residuum_norm[iteration]/residuum_norm[0])
          residuum_infnorm[iteration]=np.max (np.abs(Residuum)) #np.linalg.norm (Residuum , np.inf)
          #print(residuum_norm[iteration], np.sqrt(np.mean (Residuum_re**2)), (np.sqrt(np.mean (Residuum_re**2))-residuum_norm[iteration] )/residuum_norm[iteration])
          #exit[iteration]=np.linalg.norm (exitcon , np.inf)
          max_iteration=max(iteration,max_iteration)
          #print(max_iteration)
        else:
          if(plotted== False ):
            number_of_iters=range(0,iteration,1)#np.zeros((iteration))
            #number_of_iters=(range[0,iteration,1])
            print(residuum_infnorm[0:iteration]/residuum_infnorm[0])
            
            plt.plot(number_of_iters, residuum_norm[0:iteration], color=colors[int((time-start_time)/jump)], label ='day'+str(1+int(time/24)))
            #plt.plot(number_of_iters, residuum_norm[0:iteration], label ='day'+str(int(time/24)) )
            plotted= True

        #errors[0,count]= bits
        #errors[1,count]= errorT
        #errors[2,count]= errorF
        #print errors
        #count=count+1
    number_of_iters=range(0,6,1)
    plt.xticks(number_of_iters)
    plt.xlabel('Solver Iterations')
    plt.ylabel('L2 Residual')
    #plt.plot((0.0,float(max_iteration) ), (1.0,1.0), 'k-')
    plt.legend(loc='upper left', prop={'size':14}, ncol=2, bbox_to_anchor=(0.50,1.04), framealpha=1)
    #plt.legend(loc='upper left', prop={'size':6}, bbox_to_anchor=(1,1)) #(loc='upper right')  #['day0'], loc='upper right', 'day1','day2', 'day3','day4', 'day5','day6', 'day7', 'day8', 'day9','day10', 'day11','day12', 'day13','day14', 'day15', 'day16'
    #plt.title(r'Evolution of Residuals')
    plt.yscale('log')
    plt.ylim(ymin = 5*10**(-5),ymax = 10**(6))
    plt.axis
    #plt.show()
    #plt.savefig(plotpath +'L2_Precon'+Precon+'_' + var+'_exp'+exp+'_codes_FF'+'_bits'+str(bits)+'.pdf', bbox_inches=0)
    #plt.close()
  
    #print( plotpath +'L2_Precon'+Precon+'_' + var+'_exp'+exp+'_codes_FF'+'_bits'+str(bits)+'.pdf')

    plt.savefig(plotpath +'Norm_true_Errors_Precon'+Precon+'_' + var+'_exp'+exp+'_codes_'+str(codesD)+str(codesQ)+'_bits'+str(bits)+'L2norm.pdf', bbox_inches=0)
    plt.close()
  
    print( plotpath +'Norm_Errors_Precon'+Precon+'_' + var+'_exp'+exp+'_codes_'+str(codesD)+str(codesQ)+'_bits'+str(bits)+'L2norm.pdf')


#compare_histograms_normal('h', 31, sigma_h) #, -0.0526478751975,0.0526478751975 )
#compare_histograms_talone('t', 31, sigma_t, -0.0280507825342,0.0280507825342 )
#compare_histograms('vn', 31, sigma_vn,-0.012620211929 , 0.012620211929)

#call_truncation_errors('h',n_cells, 1, num_timesteps)
#call_truncation_errors('t',n_cells, 32, num_timesteps)
#call_truncation_errors('vn',n_edges, 32, num_timesteps)


