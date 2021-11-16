import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 13}

matplotlib.rc('font', **font)

print('scaling version 3')
resolutions, rp_times=np.loadtxt('Speed_test_allres_preconLapl.o6674', usecols=(1, 2),unpack=True) # scaling of ax

dp_times=np.loadtxt('Speed_test_allres_preconLapl_dp.o6659', usecols=(2))
print(resolutions[:60]**2*32*64, dp_times[:60]/rp_times[:60])
print(resolutions[:223]**2*32*64, dp_times[:223]/rp_times[:223])

plt.plot(resolutions[:len(dp_times)]**2*32*64, dp_times[:], label='DP')
plt.plot(resolutions[:]**2*32*64, rp_times[:], label='MP')
plt.legend(loc='upper left', prop={'size':14}, ncol=2, bbox_to_anchor=(0.2,1.0), framealpha=1)
plt.xlabel('#Grid Points')
plt.ylabel('Measured Computing Time [s]')
plt.xscale('log')
plt.yscale('log')

plt.savefig('Time_measurement.png', bbox_inches=0)
plt.close()

plt.plot(resolutions[:len(dp_times)]**2*32*64, dp_times[:]/rp_times[:len(dp_times)])
plt.axhline(y=4.0, color='r', linestyle='-')
plt.axhline(y=2.0, color='g', linestyle='-')

#plt.legend(loc='upper left', prop={'size':14}, ncol=2, bbox_to_anchor=(0.2,1.0), framealpha=1)
plt.xlabel('#Grid Points')
plt.ylabel('Measured Speed-up')
plt.xscale('log')
plt.ylim(ymin = 0,ymax = 8)
plt.savefig('Speed_up.png', bbox_inches=0)
plt.close()
