import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import linecache
import math

readin_file='ven1.flds2.baseline.dat'

degree=np.zeros(74)
read_degree=np.zeros(shape=(13,6))
for i in range(0,12):
  a=linecache.getline(readin_file,i+6+833)
  read_degree[i,:]=a.split()
a=linecache.getline(readin_file,18+833)
read_degree[12,0]=a
for i in range(0,73):
  degree[i]=read_degree[i//6,i%6]
for i in range(37,73):
  degree[i]=360+degree[i]
degree[73]=360

height=np.zeros(66)
matrix_2=np.zeros(shape=(814,6))
for i in range(0,814):
  a=linecache.getline(readin_file,i+20+833)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  height[i]=matrix_2[i//6,i%6]
  
o2=np.zeros(shape=(74,66))
for i in range(0,66):
  for j in range(0,73):
    o2[j,i]=matrix_2[(i*73+j)//6+11,(i*73+j)%6]
o2[73,:]=o2[0,:]    
X,Y=np.meshgrid(degree,height)
plt.figure()
#plt.contourf(X,Y,temperature.T,15,cmap=plt.cm.nipy_spectral)
plt.contourf(X,Y,o2.T,18,vmin=-6.8,vmax=-3.2,cmap=plt.cm.jet)
C=plt.contour(X,Y,o2.T,18,linestyles='solid',vmin=-6.8,vmax=-3.2,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,fmt='%.2f',fontsize=6)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)

plt.xlabel('Longitude (degree)')
plt.ylabel('Height (km)')
my_x_ticks = np.arange(0, 360, 60)
plt.xticks(my_x_ticks)
plt.savefig('Figure_o2_a.png', dpi=300)

plt.savefig('Figure_o2_a.eps', format='eps')

o2_12=np.zeros(66)
plt.figure()
o2_12=o2[36,:]
for i in range(0,66):
  o2_12[i]=10**(o2_12[i])
plt.semilogx(o2_12,height,linewidth=2)
plt.xlabel('O$_{{2}}$ Mixing Ratio')
plt.xlim([10**(-5.25),10**(-3.0)])
#plt.xticks([-5,-4,-3],[r'10$^{{-5}}$',r'10$^{{-4}}$',r'10$^{{-3}}$'])
plt.ylabel('Height (km)')
plt.savefig('Figure_o2_2a.png', dpi=300)

plt.savefig('Figure_o2_2a.eps', format='eps')
#########################################################
readin_file='ven2.flds2.imaginus1.dat'

degree=np.zeros(74)
read_degree=np.zeros(shape=(13,6))
for i in range(0,12):
  a=linecache.getline(readin_file,i+6+833)
  read_degree[i,:]=a.split()
a=linecache.getline(readin_file,18+833)
read_degree[12,0]=a
for i in range(0,73):
  degree[i]=read_degree[i//6,i%6]
for i in range(37,73):
  degree[i]=360+degree[i]
degree[73]=360

height=np.zeros(66)
matrix_2=np.zeros(shape=(814,6))
for i in range(0,814):
  a=linecache.getline(readin_file,i+20+833)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  height[i]=matrix_2[i//6,i%6]
  
o2=np.zeros(shape=(74,66))
for i in range(0,66):
  for j in range(0,73):
    o2[j,i]=matrix_2[(i*73+j)//6+11,(i*73+j)%6]
o2[73,:]=o2[0,:]    
X,Y=np.meshgrid(degree,height)
plt.figure()
#plt.contourf(X,Y,temperature.T,15,cmap=plt.cm.nipy_spectral)
plt.contourf(X,Y,o2.T,13,vmin=-6.8,vmax=-3.2,cmap=plt.cm.jet)
C=plt.contour(X,Y,o2.T,13,linestyles='solid',vmin=-6.8,vmax=-3.2,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,fmt='%.2f',fontsize=6)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)

plt.xlabel('Longitude (degree)')
plt.ylabel('Height (km)')
my_x_ticks = np.arange(0, 360, 60)
plt.xticks(my_x_ticks)
plt.savefig('Figure_o2_b.png', dpi=300)

plt.savefig('Figure_o2_b.eps', format='eps')

plt.figure()
o2_12=o2[36,:]
for i in range(0,66):
  o2_12[i]=10**(o2_12[i])
plt.semilogx(o2_12,height,linewidth=2)
plt.xlabel('O$_{{2}}$ Mixing Ratio')
plt.xlim([10**(-5.25),10**(-3.0)])
#plt.xticks([-5,-4,-3],[r'10$^{{-5}}$',r'10$^{{-4}}$',r'10$^{{-3}}$'])
plt.ylabel('Height (km)')
plt.savefig('Figure_o2_2b.png', dpi=300)

plt.savefig('Figure_o2_2b.eps', format='eps')
##############################################################
readin_file='ven3.flds2.imaginus2.dat'

degree=np.zeros(74)
read_degree=np.zeros(shape=(13,6))
for i in range(0,12):
  a=linecache.getline(readin_file,i+6+833)
  read_degree[i,:]=a.split()
a=linecache.getline(readin_file,18+833)
read_degree[12,0]=a
for i in range(0,73):
  degree[i]=read_degree[i//6,i%6]
for i in range(37,73):
  degree[i]=360+degree[i]
degree[73]=360

height=np.zeros(66)
matrix_2=np.zeros(shape=(814,6))
for i in range(0,814):
  a=linecache.getline(readin_file,i+20+833)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  height[i]=matrix_2[i//6,i%6]
  
o2=np.zeros(shape=(74,66))
for i in range(0,66):
  for j in range(0,73):
    o2[j,i]=matrix_2[(i*73+j)//6+11,(i*73+j)%6]
o2[73,:]=o2[0,:]    
X,Y=np.meshgrid(degree,height)
plt.figure()
#plt.contourf(X,Y,temperature.T,15,cmap=plt.cm.nipy_spectral)
plt.contourf(X,Y,o2.T,15,vmin=-6.8,vmax=-3.2,cmap=plt.cm.jet)
C=plt.contour(X,Y,o2.T,15,linestyles='solid',vmin=-6.8,vmax=-3.2,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,fmt='%.2f',fontsize=6)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)

plt.xlabel('Longitude (degree)')
plt.ylabel('Height (km)')
my_x_ticks = np.arange(0, 360, 60)
plt.xticks(my_x_ticks)
plt.savefig('Figure_o2_c.png', dpi=300)

plt.savefig('Figure_o2_c.eps', format='eps')

plt.figure()
o2_12=o2[36,:]
for i in range(0,66):
  o2_12[i]=10**(o2_12[i])
plt.semilogx(o2_12,height,linewidth=2)
plt.xlabel('O$_{{2}}$ Mixing Ratio')
plt.xlim([10**(-5.25),10**(-3.0)])
#plt.xticks([-5,-4,-3],[r'10$^{{-5}}$',r'10$^{{-4}}$',r'10$^{{-3}}$'])
plt.ylabel('Height (km)')
plt.savefig('Figure_o2_2c.png', dpi=300)

plt.savefig('Figure_o2_2c.eps', format='eps')
##########################################################
readin_file='ven4.flds2.imaginus5.dat'

degree=np.zeros(74)
read_degree=np.zeros(shape=(13,6))
for i in range(0,12):
  a=linecache.getline(readin_file,i+6+833)
  read_degree[i,:]=a.split()
a=linecache.getline(readin_file,18+833)
read_degree[12,0]=a
for i in range(0,73):
  degree[i]=read_degree[i//6,i%6]
for i in range(37,73):
  degree[i]=360+degree[i]
degree[73]=360

height=np.zeros(66)
matrix_2=np.zeros(shape=(814,6))
for i in range(0,814):
  a=linecache.getline(readin_file,i+20+833)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  height[i]=matrix_2[i//6,i%6]
  
o2=np.zeros(shape=(74,66))
for i in range(0,66):
  for j in range(0,73):
    o2[j,i]=matrix_2[(i*73+j)//6+11,(i*73+j)%6]
o2[73,:]=o2[0,:]    
X,Y=np.meshgrid(degree,height)
plt.figure()
#plt.contourf(X,Y,temperature.T,15,cmap=plt.cm.nipy_spectral)
plt.contourf(X,Y,o2.T,18,vmin=-6.8,vmax=-3.2,cmap=plt.cm.jet)
C=plt.contour(X,Y,o2.T,18,linestyles='solid',vmin=-6.8,vmax=-3.2,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,fmt='%.2f',fontsize=6)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)

plt.xlabel('Longitude (degree)')
plt.ylabel('Height (km)')
my_x_ticks = np.arange(0, 360, 60)
plt.xticks(my_x_ticks)
plt.savefig('Figure_o2_d.png', dpi=300)

plt.savefig('Figure_o2_d.eps', format='eps')

plt.figure()
o2_12=o2[36,:]
for i in range(0,66):
  o2_12[i]=10**(o2_12[i])
plt.semilogx(o2_12,height,linewidth=2)
plt.xlabel('O$_{{2}}$ Mixing Ratio')
plt.xlim([10**(-5.25),10**(-3.0)])
#plt.xticks([-5,-4,-3],[r'10$^{{-5}}$',r'10$^{{-4}}$',r'10$^{{-3}}$'])
plt.ylabel('Height (km)')
plt.savefig('Figure_o2_2d.png', dpi=300)

plt.savefig('Figure_o2_2d.eps', format='eps')