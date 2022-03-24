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


X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
#plt.contourf(X,Y,temperature.T,15,cmap=plt.cm.nipy_spectral)
lines_lv=[80,100,120,140,160,180,200,220,240,260,280,300,320]
cf_lv=[0,80,100,120,140,160,180,200,220,240,260,280,300,320,360]
plt.contourf(X,Y,height.T,levels=cf_lv,vmin=70,vmax=500,cmap=plt.cm.jet)
C=plt.contour(X,Y,height.T,levels=lines_lv,linestyles='solid',vmin=70,vmax=500,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('white','white','white','black','black','black','black','black','black','black','black','black','black'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)
plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,110,150,190,230,270])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_a_height_yaxis.png', dpi=300)

plt.savefig('Figure_a_height_yaxis.eps', format='eps')
plt.close()

####################################


matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+903)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+903)
matrix_2[839,0:3]=a.split()

temperature=np.zeros(shape=(74,69))
temperature_log=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    temperature[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
temperature[73,:]=temperature[0,:]
temperature[:,68]=temperature[:,67]
for i in range(0,69):
  for j in range(0,74):
    temperature_log[j,i]=math.log10(temperature[j,i])
    

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
#plt.contourf(X,Y,temperature.T,15,cmap=plt.cm.nipy_spectral)
lines_lv=[110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260]
plt.contourf(X,Y,temperature_log.T,100,vmin=math.log10(100),vmax=math.log10(600),cmap=plt.cm.jet)
C=plt.contour(X,Y,temperature.T,levels=lines_lv,vmin=110,vmax=260,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('white','white','white','white','black','black','black','black','black','black','black','black','black','black','black','black'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)

plt.xlabel('Longitude (degrees)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,110,150,190,230,270])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')

plt.savefig('Figure_a_temperature_yaxis.png', dpi=300)

plt.savefig('Figure_a_temperature_yaxis.eps', format='eps')
plt.close()

###########################################

matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+1774)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+1774)
matrix_2[839,0:3]=a.split()
o1=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    o1[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
o1[73,:]=o1[0,:]

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-9.5,-9.0,-8.5,-8.0,-7.5,-7.0,-6.5,-6.0,-5.5,-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5]
cf_lv=[-20.0,-9.5,-9.0,-8.5,-8.0,-7.5,-7.0,-6.5,-6.0,-5.5,-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0]
plt.contourf(X,Y,o1.T,levels=cf_lv,vmin=-10.0,vmax=-0.0,cmap=plt.cm.jet)
C=plt.contour(X,Y,o1.T,levels=lines_lv,linestyles='solid',vmin=-10.0,vmax=-0.0,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('white','white','white','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','white'),fmt='%.2f',fontsize=5.5)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,110,150,190,230,270])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_a_o1_yaxis.png', dpi=300)

plt.savefig('Figure_a_o1_yaxis.eps', format='eps')
plt.close()

#############################################

matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+2645)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+2645)
matrix_2[839,0:3]=a.split()
co=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    co[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
co[73,:]=co[0,:]

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-5.5,-5.25,-5.0,-4.75,-4.5,-4.25,-4.0,-3.75,-3.5,-3.25,-3.0,-2.75,-2.5,-2.25,-2.0,-1.75,-1.5,-1.25,-1.0,-0.75]
cf_lv=[-10.0,-5.5,-5.25,-5.0,-4.75,-4.5,-4.25,-4.0,-3.75,-3.5,-3.25,-3.0,-2.75,-2.5,-2.25,-2.0,-1.75,-1.5,-1.25,-1.0,-0.75,0.0]
plt.contourf(X,Y,co.T,levels=cf_lv,vmin=-5.25,vmax=-0.5,cmap=plt.cm.jet)
C=plt.contour(X,Y,co.T,levels=lines_lv,linestyles='solid',vmin=-5.25,vmax=-0.5,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('white','white','white','white','black','black','black','black','black','black','black','black','black','black','black','black','black','black','white','white'),fmt='%.2f',fontsize=5.5)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,110,150,190,230,270])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_a_co_yaxis.png', dpi=300)

plt.savefig('Figure_a_co_yaxis.eps', format='eps')
plt.close()
#############################################

matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+3516)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+3516)
matrix_2[839,0:3]=a.split()
co2=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    co2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
co2[73,:]=co2[0,:]

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5]
cf_lv=[-10.0,-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0]
plt.contourf(X,Y,co2.T,levels=cf_lv,vmin=-5.0,vmax=-0.0,cmap=plt.cm.jet)
C=plt.contour(X,Y,co2.T,levels=lines_lv,linestyles='solid',vmin=-5.0,vmax=-0.0,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('white','white','white','black','black','black','black','black','black','white'),fmt='%.2f',fontsize=5.5)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,110,150,190,230,270])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_a_co2_yaxis.png', dpi=300)

plt.savefig('Figure_a_co2_yaxis.eps', format='eps')
plt.close()

#############################################

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

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-7.5,-7.25,-7.0,-6.75,-6.5,-6.25,-6.0,-5.75,-5.5,-5.25,-5.0,-4.75,-4.5,-4.25,-4.0,-3.75,-3.5,-3.25]
cf_lv=[-9.0,-7.5,-7.25,-7.0,-6.75,-6.5,-6.25,-6.0,-5.75,-5.5,-5.25,-5.0,-4.75,-4.5,-4.25,-4.0,-3.75,-3.5,-3.25,-3.0]
plt.contourf(X,Y,o2.T,levels=cf_lv,vmin=-7.5,vmax=-3.25,cmap=plt.cm.jet)
C=plt.contour(X,Y,o2.T,levels=lines_lv,linestyles='solid',vmin=-7.5,vmax=-3.25,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('white','white','white','white','black','black','black','black','black','black','black','black','black','black','black','black','white','white'),fmt='%.2f',fontsize=5.5)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,110,150,190,230,270])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_a_o2_yaxis.png', dpi=300)

plt.savefig('Figure_a_o2_yaxis.eps', format='eps')
plt.close()

#############################################

matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+5258)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+5258)
matrix_2[839,0:3]=a.split()
wind=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    wind[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
wind[73,:]=wind[0,:]
wind[:,68]=wind[:,67]

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-150,-125,-100,-75,-50,-25,0,25,50,75,100,125,150,175,200,225,250,275]
cf_lv=[-175,-150,-125,-100,-75,-50,-25,0,25,50,75,100,125,150,175,200,225,250,275,300]
plt.contourf(X,Y,wind.T,levels=cf_lv,vmin=-350,vmax=450,cmap=plt.cm.jet)
C=plt.contour(X,Y,wind.T,levels=lines_lv,linestyles='solid',vmin=-350,vmax=450,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black'),fmt='%i',fontsize=5.5)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,110,150,190,230,270])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_a_wind_yaxis.png', dpi=300)

plt.savefig('Figure_a_wind_yaxis.eps', format='eps')
plt.close()

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

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
#plt.contourf(X,Y,temperature.T,15,cmap=plt.cm.nipy_spectral)
lines_lv=[80,100,120,140,160,180,200,220,240,260,280,300,320]
cf_lv=[0,80,100,120,140,160,180,200,220,240,260,280,300,320,340]
plt.contourf(X,Y,height.T,levels=cf_lv,vmin=70,vmax=500,cmap=plt.cm.jet)
C=plt.contour(X,Y,height.T,levels=lines_lv,linestyles='solid',vmin=70,vmax=500,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('white','white','white','black','black','black','black','black','black','black','black','black','black'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)

plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,110,150,190,230,270])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_b_height_yaxis.png', dpi=300)

plt.savefig('Figure_b_height_yaxis.eps', format='eps')
plt.close()

####################################


matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+903)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+903)
matrix_2[839,0:3]=a.split()

temperature=np.zeros(shape=(74,69))
temperature_log=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    temperature[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
temperature[73,:]=temperature[0,:]
temperature[:,68]=temperature[:,67]
for i in range(0,69):
  for j in range(0,74):
    temperature_log[j,i]=math.log10(temperature[j,i])
    

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
#plt.contourf(X,Y,temperature.T,15,cmap=plt.cm.nipy_spectral)
lines_lv=[110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260]
plt.contourf(X,Y,temperature_log.T,100,vmin=math.log10(100),vmax=math.log10(600),cmap=plt.cm.jet)
C=plt.contour(X,Y,temperature.T,levels=lines_lv,vmin=110,vmax=260,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('white','white','white','white','black','black','black','black','black','black','black','black','black','black','black','black'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)

plt.xlabel('Longitude (degrees)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,110,150,190,230,270])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_b_temperature_yaxis.png', dpi=300)

plt.savefig('Figure_b_temperature_yaxis.eps', format='eps')
plt.close()
###########################################
matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+1774)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+1774)
matrix_2[839,0:3]=a.split()
o1=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    o1[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
o1[73,:]=o1[0,:]

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-9.5,-9.0,-8.5,-8.0,-7.5,-7.0,-6.5,-6.0,-5.5,-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5]
cf_lv=[-20.0,-9.5,-9.0,-8.5,-8.0,-7.5,-7.0,-6.5,-6.0,-5.5,-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0]
plt.contourf(X,Y,o1.T,levels=cf_lv,vmin=-10.0,vmax=-0.0,cmap=plt.cm.jet)
C=plt.contour(X,Y,o1.T,levels=lines_lv,linestyles='solid',vmin=-10.0,vmax=-0.0,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('white','white','white','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','white'),fmt='%.2f',fontsize=5.5)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,110,150,190,230,270])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_b_o1_axis.png', dpi=300)

plt.savefig('Figure_b_o1_axis.eps', format='eps')
plt.close()
#############################################

matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+2645)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+2645)
matrix_2[839,0:3]=a.split()
co=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    co[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
co[73,:]=co[0,:]

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-5.5,-5.25,-5.0,-4.75,-4.5,-4.25,-4.0,-3.75,-3.5,-3.25,-3.0,-2.75,-2.5,-2.25,-2.0,-1.75,-1.5,-1.25,-1.0,-0.75]
cf_lv=[-10.0,-5.5,-5.25,-5.0,-4.75,-4.5,-4.25,-4.0,-3.75,-3.5,-3.25,-3.0,-2.75,-2.5,-2.25,-2.0,-1.75,-1.5,-1.25,-1.0,-0.75,0.0]
plt.contourf(X,Y,co.T,levels=cf_lv,vmin=-5.25,vmax=-0.5,cmap=plt.cm.jet)
C=plt.contour(X,Y,co.T,levels=lines_lv,linestyles='solid',vmin=-5.25,vmax=-0.5,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('white','white','white','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','white','white'),fmt='%.2f',fontsize=5.5)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,110,150,190,230,270])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_b_co_yaxis.png', dpi=300)

plt.savefig('Figure_b_co_yaxis.eps', format='eps')
plt.close()

#############################################

matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+3516)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+3516)
matrix_2[839,0:3]=a.split()
co2=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    co2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
co2[73,:]=co2[0,:]

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5]
cf_lv=[-10.0,-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0]
plt.contourf(X,Y,co2.T,levels=cf_lv,vmin=-5.0,vmax=-0.0,cmap=plt.cm.jet)
C=plt.contour(X,Y,co2.T,levels=lines_lv,linestyles='solid',vmin=-5.0,vmax=-0.0,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('white','white','white','black','black','black','black','black','black','white'),fmt='%.2f',fontsize=6)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,110,150,190,230,270])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_b_co2_yaxis.png', dpi=300)

plt.savefig('Figure_b_co2_yaxis.eps', format='eps')
plt.close()

#############################################

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

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-7.5,-7.25,-7.0,-6.75,-6.5,-6.25,-6.0,-5.75,-5.5,-5.25,-5.0,-4.75,-4.5,-4.25,-4.0,-3.75,-3.5,-3.25]
cf_lv=[-9.0,-7.5,-7.25,-7.0,-6.75,-6.5,-6.25,-6.0,-5.75,-5.5,-5.25,-5.0,-4.75,-4.5,-4.25,-4.0,-3.75,-3.5,-3.25,-3.0]
plt.contourf(X,Y,o2.T,levels=cf_lv,vmin=-7.5,vmax=-3.25,cmap=plt.cm.jet)
C=plt.contour(X,Y,o2.T,levels=lines_lv,linestyles='solid',vmin=-7.5,vmax=-3.25,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('white','white','white','white','black','black','black','black','black','black','black','black','black','black','black','black','white','white'),fmt='%.2f',fontsize=5.5)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,110,150,190,230,270])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_b_o2_yaxis.png', dpi=300)

plt.savefig('Figure_b_o2_yaxis.eps', format='eps')
plt.close()
#############################################

matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+5258)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+5258)
matrix_2[839,0:3]=a.split()
wind=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    wind[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
wind[73,:]=wind[0,:]
wind[:,68]=wind[:,67]

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-150,-125,-100,-75,-50,-25,0,25,50,75,100,125,150,175,200,225,250,275]
cf_lv=[-175,-150,-125,-100,-75,-50,-25,0,25,50,75,100,125,150,175,200,225,250,275,300]
plt.contourf(X,Y,wind.T,levels=cf_lv,vmin=-350,vmax=450,cmap=plt.cm.jet)
C=plt.contour(X,Y,wind.T,levels=lines_lv,linestyles='solid',vmin=-350,vmax=450,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black'),fmt='%i',fontsize=5.5)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,110,150,190,230,270])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_b_wind_yaxis.png', dpi=300)

plt.savefig('Figure_b_wind_yaxis.eps', format='eps')
plt.close()

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

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)


lines_lv=[80,100,120,140,160,180,200,220,240,260,280,300,320]
cf_lv=[0,80,100,120,140,160,180,200,220,240,260,280,300,320,360]
plt.contourf(X,Y,height.T,levels=cf_lv,vmin=70,vmax=500,cmap=plt.cm.jet)
C=plt.contour(X,Y,height.T,levels=lines_lv,linestyles='solid',vmin=70,vmax=500,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('white','white','white','black','black','black','black','black','black','black','black','black','black'),fmt='%i',fontsize=5.5)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,110,150,190,230,270])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_c_height_yaxis.png', dpi=300)

plt.savefig('Figure_c_height_yaxis.eps', format='eps')
plt.close()
####################################


matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+903)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+903)
matrix_2[839,0:3]=a.split()

temperature=np.zeros(shape=(74,69))
temperature_log=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    temperature[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
temperature[73,:]=temperature[0,:]
temperature[:,68]=temperature[:,67]
for i in range(0,69):
  for j in range(0,74):
    temperature_log[j,i]=math.log10(temperature[j,i])
    

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
#plt.contourf(X,Y,temperature.T,15,cmap=plt.cm.nipy_spectral)
lines_lv=[110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270]
plt.contourf(X,Y,temperature_log.T,100,vmin=math.log10(100),vmax=math.log10(600),cmap=plt.cm.jet)
C=plt.contour(X,Y,temperature.T,levels=lines_lv,vmin=110,vmax=270,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('white','white','white','white','black','black','black','black','black','black','black','black','black','black','black','black','black'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)

plt.xlabel('Longitude (degrees)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,110,150,190,230,270])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_c_temperature_yaxis.png', dpi=300)

plt.savefig('Figure_c_temperature_yaxis.eps', format='eps')
plt.close()
###########################################

matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+1774)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+1774)
matrix_2[839,0:3]=a.split()
o1=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    o1[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
o1[73,:]=o1[0,:]

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-9.5,-9.0,-8.5,-8.0,-7.5,-7.0,-6.5,-6.0,-5.5,-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5]
cf_lv=[-20.0,-9.5,-9.0,-8.5,-8.0,-7.5,-7.0,-6.5,-6.0,-5.5,-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0]
plt.contourf(X,Y,o1.T,levels=cf_lv,vmin=-10.0,vmax=-0.0,cmap=plt.cm.jet)
C=plt.contour(X,Y,o1.T,levels=lines_lv,linestyles='solid',vmin=-10.0,vmax=-0.0,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('white','white','white','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','white'),fmt='%.2f',fontsize=5.5)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,110,150,190,230,270])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_c_o1_yaxis.png', dpi=300)

plt.savefig('Figure_c_o1_yaxis.eps', format='eps')
plt.close()
#############################################

matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+2645)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+2645)
matrix_2[839,0:3]=a.split()
co=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    co[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
co[73,:]=co[0,:]

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-5.5,-5.25,-5.0,-4.75,-4.5,-4.25,-4.0,-3.75,-3.5,-3.25,-3.0,-2.75,-2.5,-2.25,-2.0,-1.75,-1.5,-1.25,-1.0,-0.75]
cf_lv=[-10.0,-5.5,-5.25,-5.0,-4.75,-4.5,-4.25,-4.0,-3.75,-3.5,-3.25,-3.0,-2.75,-2.5,-2.25,-2.0,-1.75,-1.5,-1.25,-1.0,-0.75,0.0]
plt.contourf(X,Y,co.T,levels=cf_lv,vmin=-5.25,vmax=-0.5,cmap=plt.cm.jet)
C=plt.contour(X,Y,co.T,levels=lines_lv,linestyles='solid',vmin=-5.25,vmax=-0.5,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('white','white','white','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','white','white'),fmt='%.2f',fontsize=5.5)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,110,150,190,230,270])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_c_co_yaxis.png', dpi=300)

plt.savefig('Figure_c_co_yaxis.eps', format='eps')
plt.close()
#############################################

matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+3516)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+3516)
matrix_2[839,0:3]=a.split()
co2=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    co2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
co2[73,:]=co2[0,:]

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5]
cf_lv=[-10.0,-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0]
plt.contourf(X,Y,co2.T,levels=cf_lv,vmin=-5.0,vmax=-0.0,cmap=plt.cm.jet)
C=plt.contour(X,Y,co2.T,levels=lines_lv,linestyles='solid',vmin=-5.0,vmax=-0.0,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('white','white','white','black','black','black','black','black','black','white'),fmt='%.2f',fontsize=5.5)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,110,150,190,230,270])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_c_co2_yaxis.png', dpi=300)

plt.savefig('Figure_c_co2_yaxis.eps', format='eps')
plt.close()
#############################################

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

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-7.5,-7.25,-7.0,-6.75,-6.5,-6.25,-6.0,-5.75,-5.5,-5.25,-5.0,-4.75,-4.5,-4.25,-4.0,-3.75,-3.5,-3.25]
cf_lv=[-9.0,-7.5,-7.25,-7.0,-6.75,-6.5,-6.25,-6.0,-5.75,-5.5,-5.25,-5.0,-4.75,-4.5,-4.25,-4.0,-3.75,-3.5,-3.25,-3.0]
plt.contourf(X,Y,o2.T,levels=cf_lv,vmin=-7.5,vmax=-3.25,cmap=plt.cm.jet)
C=plt.contour(X,Y,o2.T,levels=lines_lv,linestyles='solid',vmin=-7.5,vmax=-3.25,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('white','white','white','white','black','black','black','black','black','black','black','black','black','black','black','black','white','white'),fmt='%.2f',fontsize=5.5)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,110,150,190,230,270])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_c_o2_yaxis.png', dpi=300)

plt.savefig('Figure_c_o2_yaxis.eps', format='eps')
plt.close()
#############################################

matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+5258)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+5258)
matrix_2[839,0:3]=a.split()
wind=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    wind[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
wind[73,:]=wind[0,:]
wind[:,68]=wind[:,67]

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-150,-125,-100,-75,-50,-25,0,25,50,75,100,125,150,175,200,225,250,275]
cf_lv=[-175,-150,-125,-100,-75,-50,-25,0,25,50,75,100,125,150,175,200,225,250,275,300]
plt.contourf(X,Y,wind.T,levels=cf_lv,vmin=-350,vmax=450,cmap=plt.cm.jet)
C=plt.contour(X,Y,wind.T,levels=lines_lv,linestyles='solid',vmin=-350,vmax=450,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black'),fmt='%i',fontsize=5.5)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,110,150,190,230,270])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_c_wind_yaxis.png', dpi=300)

plt.savefig('Figure_c_wind_yaxis.eps', format='eps')
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

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
#plt.contourf(X,Y,temperature.T,15,cmap=plt.cm.nipy_spectral)
lines_lv=[80,100,120,140,160,180,200,220,240,260,280,300,350,400,450,500,550]
cf_lv=[0,80,100,120,140,160,180,200,220,240,260,280,300,350,400,450,500,550,700]
plt.contourf(X,Y,height.T,levels=cf_lv,vmin=70,vmax=500,cmap=plt.cm.jet)
C=plt.contour(X,Y,height.T,levels=lines_lv,linestyles='solid',vmin=70,vmax=500,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('white','white','black','black','black','black','black','black','black','black','black','black','black','black','black','white','white'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)

plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,140,210,280,350,420])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_d_height_yaxis.png', dpi=300)

plt.savefig('Figure_d_height_yaxis.eps', format='eps')
plt.close()
####################################


matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+903)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+903)
matrix_2[839,0:3]=a.split()

temperature=np.zeros(shape=(74,69))
temperature_log=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    temperature[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
temperature[73,:]=temperature[0,:]
temperature[:,68]=temperature[:,67]
for i in range(0,69):
  for j in range(0,74):
    temperature_log[j,i]=math.log10(temperature[j,i])
    

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
#plt.contourf(X,Y,temperature.T,15,cmap=plt.cm.nipy_spectral)
lines_lv=[110,120,140,160,180,200,220,240,280,310,340,370,400,430,460,490,520,550,580,610,640]
plt.contourf(X,Y,temperature_log.T,100,vmin=math.log10(100),vmax=math.log10(600),cmap=plt.cm.jet)
C=plt.contour(X,Y,temperature.T,levels=lines_lv,vmin=100,vmax=600,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('white','white','white','black','black','black','black','black','black','black','black','black','black','black','black','black','white','white','white','white','white'),fmt='%i',fontsize=5.5)


#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)

plt.xlabel('Longitude (degrees)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,140,210,280,350,420])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_d_temperature_yaxis.png', dpi=300)

plt.savefig('Figure_d_temperature_yaxis.eps', format='eps')
plt.close()
###########################################

matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+1774)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+1774)
matrix_2[839,0:3]=a.split()
o1=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    o1[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
o1[73,:]=o1[0,:]

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-9.5,-9.0,-8.5,-8.0,-7.5,-7.0,-6.5,-6.0,-5.5,-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5]
cf_lv=[-20.0,-9.5,-9.0,-8.5,-8.0,-7.5,-7.0,-6.5,-6.0,-5.5,-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0]
plt.contourf(X,Y,o1.T,levels=cf_lv,vmin=-10.0,vmax=-0.0,cmap=plt.cm.jet)
C=plt.contour(X,Y,o1.T,levels=lines_lv,linestyles='solid',vmin=-10.0,vmax=-0.0,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('white','white','white','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','white'),fmt='%.2f',fontsize=5.5)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,140,210,280,350,420])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_d_o1_yaxis.png', dpi=300)

plt.savefig('Figure_d_o1_yaxis.eps', format='eps')
plt.close()

#############################################

matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+2645)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+2645)
matrix_2[839,0:3]=a.split()
co=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    co[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
co[73,:]=co[0,:]

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-5.5,-5.25,-5.0,-4.75,-4.5,-4.25,-4.0,-3.75,-3.5,-3.25,-3.0,-2.75,-2.5,-2.25,-2.0,-1.75,-1.5,-1.25,-1.0,-0.75]
cf_lv=[-10.0,-5.5,-5.25,-5.0,-4.75,-4.5,-4.25,-4.0,-3.75,-3.5,-3.25,-3.0,-2.75,-2.5,-2.25,-2.0,-1.75,-1.5,-1.25,-1.0,-0.75,0.0]
plt.contourf(X,Y,co.T,levels=cf_lv,vmin=-5.25,vmax=-0.5,cmap=plt.cm.jet)
C=plt.contour(X,Y,co.T,levels=lines_lv,linestyles='solid',vmin=-5.25,vmax=-0.5,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('white','white','white','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','white','white'),fmt='%.2f',fontsize=5.5)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,140,210,280,350,420])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_d_co_yaxis.png', dpi=300)

plt.savefig('Figure_d_co_yaxis.eps', format='eps')
plt.close()
#############################################

matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+3516)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+3516)
matrix_2[839,0:3]=a.split()
co2=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    co2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
co2[73,:]=co2[0,:]

for i in range(57,69):
  for j in range(15,62):
    co2[j,i]=-5.5

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5]
cf_lv=[-10.0,-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0]
plt.contourf(X,Y,co2.T,levels=cf_lv,vmin=-5.0,vmax=-0.0,cmap=plt.cm.jet)
C=plt.contour(X,Y,co2.T,levels=lines_lv,linestyles='solid',vmin=-5.0,vmax=-0.0,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('white','white','white','black','black','black','black','black','black','white'),fmt='%.2f',fontsize=5.5)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,140,210,280,350,420])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_d_co2_yaxis.png', dpi=300)

plt.savefig('Figure_d_co2_yaxis.eps', format='eps')
plt.close()
#############################################

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

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-7.5,-7.25,-7.0,-6.75,-6.5,-6.25,-6.0,-5.75,-5.5,-5.25,-5.0,-4.75,-4.5,-4.25,-4.0,-3.75,-3.5,-3.25]
cf_lv=[-9.0,-7.5,-7.25,-7.0,-6.75,-6.5,-6.25,-6.0,-5.75,-5.5,-5.25,-5.0,-4.75,-4.5,-4.25,-4.0,-3.75,-3.5,-3.25,-3.0]
plt.contourf(X,Y,o2.T,levels=cf_lv,vmin=-7.5,vmax=-3.25,cmap=plt.cm.jet)
C=plt.contour(X,Y,o2.T,levels=lines_lv,linestyles='solid',vmin=-7.5,vmax=-3.25,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('white','white','white','white','black','black','black','black','black','black','black','black','black','black','black','black','white','white'),fmt='%.2f',fontsize=5.5)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,140,210,280,350,420])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_d_o2_yaxis.png', dpi=300)

plt.savefig('Figure_d_o2_yaxis.eps', format='eps')
plt.close()
#############################################

matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+5258)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+5258)
matrix_2[839,0:3]=a.split()
wind=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    wind[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
wind[73,:]=wind[0,:]
wind[:,68]=wind[:,67]

X,Y=np.meshgrid(degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-350,-300,-250,-200,-150,-100,-50,0,50,100,150,200,250,300,350,400,450]
cf_lv=[-400,-350,-300,-250,-200,-150,-100,-50,0,50,100,150,200,250,300,350,400,450,500]
plt.contourf(X,Y,wind.T,levels=cf_lv,vmin=-350,vmax=450,cmap=plt.cm.jet)
C=plt.contour(X,Y,wind.T,levels=lines_lv,linestyles='solid',vmin=-350,vmax=450,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=8,colors=('white','white','black','black','black','black','black','black','black','black','black','black','black','black','black','white','white'),fmt='%i',fontsize=5.5)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','-120','-60','0'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([0.0,6.0,12.0,18.0,24.0])

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,140,210,280,350,420])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')
plt.savefig('Figure_d_wind_yaxis.png', dpi=300)

plt.savefig('Figure_d_wind_yaxis.eps', format='eps')
plt.close()