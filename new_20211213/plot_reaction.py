import math
import numpy as np
import matplotlib.pyplot as plt
import linecache
from matplotlib.ticker import AutoMinorLocator

##########################
readin_file='/home/jzl/plot_parkinson_paper/new_20211213/vtgcm.VEN1.JCO2.dat'

new_height=np.zeros(55)
read_height=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+7)
  read_height[i,:]=a.split()
a=linecache.getline(readin_file,9+7)
read_height[9,0:1]=a.split()
for i in range(0,55):
  new_height[i]=read_height[i//6,i%6]

height_readin_file='/home/jzl/plot_parkinson_paper/field7/vtgcm.VEN1.ZP.FLDS7.dat'
matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(height_readin_file,i+32)
  matrix_2[i,:]=a.split()
a=linecache.getline(height_readin_file,839+32)
matrix_2[839,0:3]=a.split()
height=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    height[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
height[73,:]=height[0,:]
height_12=np.zeros(69)
height_24=np.zeros(69)
for i in range(0,69):
  height_12[i]=height[36,i]
for i in range(0,69):
  height_24[i]=height[1,i]
  
def forward_12(xx):
  return np.interp(xx,new_height,height_12[0:55])

def inverse_12(xx):
  return np.interp(xx,height_12[0:55],new_height)
  
def forward_24(xx):
  return np.interp(xx,new_height,height_24[0:55])

def inverse_24(xx):
  return np.interp(xx,height_24[0:55],new_height)
  
def forward2(xxx):
  return (5.0*(10.0**-6)/(np.exp(xxx)))

def inverse2(xxx):
  with np.errstate(divide='ignore'):
    return (np.log(5.0*(10.0**-6)/xxx))

RJIS=np.zeros(55)
read_RJIS=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+17)
  read_RJIS[i,:]=a.split()
a=linecache.getline(readin_file,9+17)
read_RJIS[9,0:1]=a.split()
for i in range(0,55):
  RJIS[i]=read_RJIS[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(RJIS,new_height,linewidth=2,linestyle='-',color='black',label='JIS')

plt.legend(fontsize=7)
plt.xlabel('RJIS')
#plt.xlim([-100,3600])
plt.ylim([-16,11])
plt.ylabel('Ln(P$_{{0}}$/P)')

secax_y = ax.secondary_yaxis('right', functions=(forward_12,inverse_12))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,110,150,190,230])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_RJIS.png', dpi=300)

plt.savefig('Figure_RJIS.eps', format='eps')