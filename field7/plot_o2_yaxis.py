import math
import numpy as np
import matplotlib.pyplot as plt
import linecache
from matplotlib.ticker import AutoMinorLocator

readin_file='vtgcm.VEN1.ZP.FLDS7.dat'

degree=np.zeros(74)
read_degree=np.zeros(shape=(13,6))
for i in range(0,12):
  a=linecache.getline(readin_file,i+6)
  read_degree[i,:]=a.split()
a=linecache.getline(readin_file,18)
read_degree[12,0]=a
for i in range(0,73):
  degree[i]=read_degree[i//6,i%6]
for i in range(37,73):
  degree[i]=360+degree[i]
degree[73]=360

new_height=np.zeros(69)
read_height=np.zeros(shape=(12,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+20)
  read_height[i,:]=a.split()
a=linecache.getline(readin_file,11+20)
read_height[11,0:3]=a.split()
for i in range(0,69):
  new_height[i]=read_height[i//6,i%6]
  

matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+32)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+32)
matrix_2[839,0:3]=a.split()
height=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    height[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
height[73,:]=height[0,:]

average_height=np.zeros(69)
for i in range(0,69):
  average_height[i]=np.mean(height[:,i])
  


def forward(xx):
  return np.interp(xx,new_height,average_height)

def inverse(xx):
  return np.interp(xx,average_height,new_height)

def forward2(xxx):
  return (5.0*(10.0**-6)/(np.exp(xxx)))

def inverse2(xxx):
  with np.errstate(divide='ignore'):
    return (np.log(5.0*(10.0**-6)/xxx))

def forwardx(xx):
  return (xx/180.0*12.0)

def inversex(xx):
  return (xx*180.0/12.0)







matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+4387)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+4387)
matrix_2[839,0:3]=a.split()
o2=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    o2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
o2[73,:]=o2[0,:]
o2[:,68]=o2[:,67]


o2_12=np.zeros(66)
o2_12=o2[36,:]
for i in range(0,66):
  o2_12[i]=10**(o2_12[i])

fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

plt.semilogx(o2_12,new_height,linewidth=2)
plt.xlabel('O$_{{2}}$ Mixing Ratio')
plt.xlim([10**(-7.5),10**(-3.0)])

plt.ylabel('Ln(P$_{{0}}$/P)')



secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,110,150,190,230,270])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_a_o2_plot_yaxis.png', dpi=300)

plt.savefig('Figure_a_o2_plot_yaxis.eps', format='eps')
plt.close()

#############################################

#####################################################

readin_file='vtgcm.VEN2.ZP.FLDS7.dat'

degree=np.zeros(74)
read_degree=np.zeros(shape=(13,6))
for i in range(0,12):
  a=linecache.getline(readin_file,i+6)
  read_degree[i,:]=a.split()
a=linecache.getline(readin_file,18)
read_degree[12,0]=a
for i in range(0,73):
  degree[i]=read_degree[i//6,i%6]
for i in range(37,73):
  degree[i]=360+degree[i]
degree[73]=360

new_height=np.zeros(69)
read_height=np.zeros(shape=(12,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+20)
  read_height[i,:]=a.split()
a=linecache.getline(readin_file,11+20)
read_height[11,0:3]=a.split()
for i in range(0,69):
  new_height[i]=read_height[i//6,i%6]
  

matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+32)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+32)
matrix_2[839,0:3]=a.split()
height=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    height[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
height[73,:]=height[0,:]

average_height=np.zeros(69)
for i in range(0,69):
  average_height[i]=np.mean(height[:,i])
  


def forward(xx):
  return np.interp(xx,new_height,average_height)

def inverse(xx):
  return np.interp(xx,average_height,new_height)

matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+4387)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+4387)
matrix_2[839,0:3]=a.split()
o2=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    o2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
o2[73,:]=o2[0,:]
o2[:,68]=o2[:,67]


o2_12=np.zeros(66)
o2_12=o2[36,:]
for i in range(0,66):
  o2_12[i]=10**(o2_12[i])

fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

plt.semilogx(o2_12,new_height,linewidth=2)
plt.xlabel('O$_{{2}}$ Mixing Ratio')
plt.xlim([10**(-7.5),10**(-3.0)])

plt.ylabel('Ln(P$_{{0}}$/P)')



secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,110,150,190,230,270])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_b_o2_plot_yaxis.png', dpi=300)

plt.savefig('Figure_b_o2_plot_yaxis.eps', format='eps')
plt.close()
#############################################



#########################################
readin_file='vtgcm.VEN3.ZP.FLDS7.dat'

degree=np.zeros(74)
read_degree=np.zeros(shape=(13,6))
for i in range(0,12):
  a=linecache.getline(readin_file,i+6)
  read_degree[i,:]=a.split()
a=linecache.getline(readin_file,18)
read_degree[12,0]=a
for i in range(0,73):
  degree[i]=read_degree[i//6,i%6]
for i in range(37,73):
  degree[i]=360+degree[i]
degree[73]=360

new_height=np.zeros(69)
read_height=np.zeros(shape=(12,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+20)
  read_height[i,:]=a.split()
a=linecache.getline(readin_file,11+20)
read_height[11,0:3]=a.split()
for i in range(0,69):
  new_height[i]=read_height[i//6,i%6]
  

matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+32)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+32)
matrix_2[839,0:3]=a.split()
height=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    height[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
height[73,:]=height[0,:]

average_height=np.zeros(69)
for i in range(0,69):
  average_height[i]=np.mean(height[:,i])
  


def forward(xx):
  return np.interp(xx,new_height,average_height)

def inverse(xx):
  return np.interp(xx,average_height,new_height)

matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+4387)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+4387)
matrix_2[839,0:3]=a.split()
o2=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    o2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
o2[73,:]=o2[0,:]
o2[:,68]=o2[:,67]


o2_12=np.zeros(66)
o2_12=o2[36,:]
for i in range(0,66):
  o2_12[i]=10**(o2_12[i])

fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

plt.semilogx(o2_12,new_height,linewidth=2)
plt.xlabel('O$_{{2}}$ Mixing Ratio')
plt.xlim([10**(-7.5),10**(-3.0)])

plt.ylabel('Ln(P$_{{0}}$/P)')



secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,110,150,190,230,270])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_c_o2_plot_yaxis.png', dpi=300)

plt.savefig('Figure_c_o2_plot_yaxis.eps', format='eps')
plt.close()

#########################################


readin_file='vtgcm.VEN4.ZP.FLDS7.dat'

degree=np.zeros(74)
read_degree=np.zeros(shape=(13,6))
for i in range(0,12):
  a=linecache.getline(readin_file,i+6)
  read_degree[i,:]=a.split()
a=linecache.getline(readin_file,18)
read_degree[12,0]=a
for i in range(0,73):
  degree[i]=read_degree[i//6,i%6]
for i in range(37,73):
  degree[i]=360+degree[i]
degree[73]=360

new_height=np.zeros(69)
read_height=np.zeros(shape=(12,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+20)
  read_height[i,:]=a.split()
a=linecache.getline(readin_file,11+20)
read_height[11,0:3]=a.split()
for i in range(0,69):
  new_height[i]=read_height[i//6,i%6]
  

matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+32)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+32)
matrix_2[839,0:3]=a.split()
height=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    height[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
height[73,:]=height[0,:]

average_height=np.zeros(69)
for i in range(0,69):
  average_height[i]=np.mean(height[:,i])
  


def forward(xx):
  return np.interp(xx,new_height,average_height)

def inverse(xx):
  return np.interp(xx,average_height,new_height)

matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+4387)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+4387)
matrix_2[839,0:3]=a.split()
o2=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    o2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
o2[73,:]=o2[0,:]
o2[:,68]=o2[:,67]


o2_12=np.zeros(66)
o2_12=o2[36,:]
for i in range(0,66):
  o2_12[i]=10**(o2_12[i])

fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

plt.semilogx(o2_12,new_height,linewidth=2)
plt.xlabel('O$_{{2}}$ Mixing Ratio')
plt.xlim([10**(-7.5),10**(-3.0)])

plt.ylabel('Ln(P$_{{0}}$/P)')



secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,140,210,280,350,420])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_d_o2_plot_yaxis.png', dpi=300)

plt.savefig('Figure_d_o2_plot_yaxis.eps', format='eps')
plt.close()


#############################################

