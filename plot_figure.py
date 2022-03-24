import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import linecache

readin_file='ven1.flds2.baseline.dat'

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

height=np.zeros(66)
matrix_2=np.zeros(shape=(814,6))
for i in range(0,814):
  a=linecache.getline(readin_file,i+20)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  height[i]=matrix_2[i//6,i%6]
  
temperature=np.zeros(shape=(74,66))
for i in range(0,66):
  for j in range(0,73):
    temperature[j,i]=matrix_2[(i*73+j)//6+11,(i*73+j)%6]
temperature[73,:]=temperature[0,:]    
X,Y=np.meshgrid(degree,height)
plt.figure()
#plt.contourf(X,Y,temperature.T,15,cmap=plt.cm.nipy_spectral)
plt.contourf(X,Y,temperature.T,17,vmin=110,vmax=260,cmap=plt.cm.jet)
C=plt.contour(X,Y,temperature.T,17,vmin=110,vmax=260,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,colors='k',fmt='%i',fontsize=7)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)

plt.xlabel('Longitude (degree)')
plt.ylabel('Height (km)')
my_x_ticks = np.arange(0, 360, 60)
plt.xticks(my_x_ticks)
plt.savefig('Figure_a.png', dpi=300)

plt.savefig('Figure_a.eps', format='eps')

#plt.figure()
#plt.contourf(X,Y,temperature.T,30,cmap=plt.cm.nipy_spectral)
#cbar=plt.colorbar()
#cbar.set_label('Temperature (K)')
#plt.xlabel('Longitude (degree)')
#plt.ylabel('Height (km)')
#plt.savefig('test_baseline_3.png', dpi=300)
#
#plt.figure()
#plt.contourf(X,Y,temperature.T,30,cmap=plt.cm.jet)
#cbar=plt.colorbar()
#cbar.set_label('Temperature (K)')
#plt.xlabel('Longitude (degree)')
#plt.ylabel('Height (km)')
#plt.savefig('test_baseline_4.png', dpi=300)
#
#plt.figure()
#plt.contourf(X,Y,temperature.T,15,cmap=plt.cm.nipy_spectral)
#C=plt.contour(X,Y,temperature.T,8,colors='black',linewidth=1)
#plt.clabel(C,inline=True,fmt='%i',fontsize=7)
#
#plt.xlabel('Longitude (degree)')
#plt.ylabel('Height (km)')
#plt.savefig('test_baseline_5.png', dpi=300)
#
#plt.figure()
#plt.contourf(X,Y,temperature.T,15,cmap=plt.cm.jet)
#C=plt.contour(X,Y,temperature.T,8,colors='black',linewidth=1)
#plt.clabel(C,inline=True,fmt='%i',fontsize=7)
#
#plt.xlabel('Longitude (degree)')
#plt.ylabel('Height (km)')
#plt.savefig('test_baseline_6.png', dpi=300)
readin_file='ven2.flds2.imaginus1.dat'

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

height=np.zeros(66)
matrix_2=np.zeros(shape=(814,6))
for i in range(0,814):
  a=linecache.getline(readin_file,i+20)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  height[i]=matrix_2[i//6,i%6]
  
temperature=np.zeros(shape=(74,66))
for i in range(0,66):
  for j in range(0,73):
    temperature[j,i]=matrix_2[(i*73+j)//6+11,(i*73+j)%6]
temperature[73,:]=temperature[0,:]    
X,Y=np.meshgrid(degree,height)
plt.figure()
#plt.contourf(X,Y,temperature.T,15,cmap=plt.cm.nipy_spectral)
plt.contourf(X,Y,temperature.T,15,vmin=110,vmax=260,cmap=plt.cm.jet)
C=plt.contour(X,Y,temperature.T,15,vmin=110,vmax=260,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,fmt='%i',fontsize=7)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)

plt.xlabel('Longitude (degree)')
plt.ylabel('Height (km)')
my_x_ticks = np.arange(0, 360, 60)
plt.xticks(my_x_ticks)
plt.savefig('Figure_b.png', dpi=300)

plt.savefig('Figure_b.eps', format='eps')
#########################################################################################
readin_file='ven3.flds2.imaginus2.dat'

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

height=np.zeros(66)
matrix_2=np.zeros(shape=(814,6))
for i in range(0,814):
  a=linecache.getline(readin_file,i+20)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  height[i]=matrix_2[i//6,i%6]
  
temperature=np.zeros(shape=(74,66))
for i in range(0,66):
  for j in range(0,73):
    temperature[j,i]=matrix_2[(i*73+j)//6+11,(i*73+j)%6]
temperature[73,:]=temperature[0,:]    
X,Y=np.meshgrid(degree,height)
plt.figure()
#plt.contourf(X,Y,temperature.T,15,cmap=plt.cm.nipy_spectral)
plt.contourf(X,Y,temperature.T,17,vmin=110,vmax=260,cmap=plt.cm.jet)
C=plt.contour(X,Y,temperature.T,17,vmin=110,vmax=260,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,fmt='%i',fontsize=7)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)

plt.xlabel('Longitude (degree)')
plt.ylabel('Height (km)')
my_x_ticks = np.arange(0, 360, 60)
plt.xticks(my_x_ticks)
plt.savefig('Figure_c.png', dpi=300)

plt.savefig('Figure_c.eps', format='eps')
################################################################################################
readin_file='ven4.flds2.imaginus5.dat'

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

height=np.zeros(66)
matrix_2=np.zeros(shape=(814,6))
for i in range(0,814):
  a=linecache.getline(readin_file,i+20)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  height[i]=matrix_2[i//6,i%6]
  
temperature=np.zeros(shape=(74,66))
for i in range(0,66):
  for j in range(0,73):
    temperature[j,i]=matrix_2[(i*73+j)//6+11,(i*73+j)%6]
temperature[73,:]=temperature[0,:]    
X,Y=np.meshgrid(degree,height)
plt.figure()
#plt.contourf(X,Y,temperature.T,15,cmap=plt.cm.nipy_spectral)
lines_lv=[110,120,130,140,150,160,170,180,190,200,210,220,230,240,260,280,300,320,340,360,380,400,420,440,460,480,500,520,540,560,580]
cf_lv=[0,110,120,130,140,150,160,170,180,190,200,210,220,230,240,260,280,300,320,340,360,380,400,420,440,460,480,500,520,540,560,580,600]
plt.contourf(X,Y,temperature.T,levels=cf_lv,vmin=100,vmax=600,cmap=plt.cm.jet)
C=plt.contour(X,Y,temperature.T,levels=lines_lv,vmin=100,vmax=600,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,fmt='%i',fontsize=7)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)

plt.xlabel('Longitude (degree)')
plt.ylabel('Height (km)')

my_x_ticks = np.arange(0, 360, 60)
plt.xticks(my_x_ticks)
plt.savefig('Figure_d.png', dpi=300)

plt.savefig('Figure_d.eps', format='eps')

plt.figure()
#plt.contourf(X,Y,temperature.T,15,cmap=plt.cm.nipy_spectral)
lines_lv=[110,120,140,160,180,200,220,240,280,310,340,370,400,430,460,490,520,550,580]
cf_lv=[0,110,120,140,160,180,200,220,240,280,310,340,370,400,430,460,490,520,550,580,600]
plt.contourf(X,Y,temperature.T,levels=cf_lv,vmin=100,vmax=600,cmap=plt.cm.jet)
C=plt.contour(X,Y,temperature.T,levels=lines_lv,vmin=100,vmax=600,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,fmt='%i',fontsize=7)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)

plt.xlabel('Longitude (degree)')
plt.ylabel('Height (km)')

my_x_ticks = np.arange(0, 360, 60)
plt.xticks(my_x_ticks)
plt.savefig('Figure_d2.png', dpi=300)

plt.savefig('Figure_d2.eps', format='eps')
