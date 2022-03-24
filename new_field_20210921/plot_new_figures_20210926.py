import math
import numpy as np
import matplotlib.pyplot as plt
import linecache
from matplotlib.ticker import AutoMinorLocator

##########################
readin_file='/home/jzl/plot_parkinson_paper/field6/vtgcm.VEN2.Flds6.dat'

degree=np.zeros(73)
read_degree=np.zeros(shape=(13,6))
for i in range(0,12):
  a=linecache.getline(readin_file,i+6)
  read_degree[i,:]=a.split()
a=linecache.getline(readin_file,18)
read_degree[12,0]=a
for i in range(0,73):
  degree[i]=read_degree[i//6,i%6]
new_degree=np.zeros(73)
new_degree[0:37]=degree[36:73]
new_degree[37:73]=degree[1:37]+360

new_height=np.zeros(68)
read_height=np.zeros(shape=(12,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+20)
  read_height[i,:]=a.split()
a=linecache.getline(readin_file,11+20)
read_height[11,0:2]=a.split()
for i in range(0,68):
  new_height[i]=read_height[i//6,i%6]
##########################
readin_file='/home/jzl/plot_parkinson_paper/field7/vtgcm.VEN4.ZP.FLDS7.dat'
matrix_2=np.zeros(shape=(840,6))
for i in range(0,839):
  a=linecache.getline(readin_file,i+903)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,839+903)
matrix_2[839,0:3]=a.split()
temperature=np.zeros(shape=(74,69))
for i in range(0,69):
  for j in range(0,73):
    temperature[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
temperature[73,:]=temperature[0,:]
temperature_old=np.zeros(shape=(73,68))
temperature_old[0:36,:]=temperature[37:73,0:68]
temperature_old[36:73,:]=temperature[0:37,0:68]

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
wind_old=np.zeros(shape=(73,68))
wind_old[0:36,:]=wind[37:73,0:68]
wind_old[36:73,:]=wind[0:37,0:68]

readin_file='vtgcm.VEN4.v2e.QNIR0p7.pFlds7.dat'

height_readin_file='vtgcm.VEN4.v2e.QNIR0p7.pFlds7.dat'
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(height_readin_file,i+32)
  matrix_2[i,:]=a.split()
a=linecache.getline(height_readin_file,827+32)
matrix_2[827,0:2]=a.split()
height=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    height[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
average_height=np.zeros(68)
for i in range(0,68):
  average_height[i]=np.mean(height[:,i])
  
def forward(xx):
  return np.interp(xx,new_height,average_height[0:68])

def inverse(xx):
  return np.interp(xx,average_height[0:68],new_height)

def forward2(xxx):
  return (5.0*(10.0**-6)/(np.exp(xxx)))

def inverse2(xxx):
  with np.errstate(divide='ignore'):
    return (np.log(5.0*(10.0**-6)/xxx))

def forwardx(xx):
  return (xx/180.0*12.0)

def inversex(xx):
  return (xx*180.0/12.0)

matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+891)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+891)
matrix_2[827,0:2]=a.split()
temperature=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    temperature[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
temperature_log=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    temperature_log[j,i]=math.log10(temperature[j,i])
new_temperature=np.zeros(shape=(73,68))
new_temperature[0:37,:]=temperature[36:73,:]
new_temperature[37:73,:]=temperature[1:37,:]
new_temperature_log=np.zeros(shape=(73,68))
new_temperature_log[0:37,:]=temperature_log[36:73,:]
new_temperature_log[37:73,:]=temperature_log[1:37,:]

X,Y=np.meshgrid(new_degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
lines_lv=[110,120,140,160,180,200,220,240,280,310,340,370,400,430,460,490,520,550,580,610,640]
plt.contourf(X,Y,new_temperature_log.T,100,vmin=math.log10(100),vmax=math.log10(650),cmap=plt.cm.jet)
C=plt.contour(X,Y,new_temperature.T,levels=lines_lv,vmin=100,vmax=600,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=5,colors=('white','white','white','black','black','black','black','black','black','black','black','black','black','black','black','black','white','white','white','white','white'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)

plt.xlabel('Longitude (degrees)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','240','300','360'])

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

plt.savefig('Figure_ven4_0p7_temperature_yaxis.png', dpi=300)

plt.savefig('Figure_ven4_0p7_temperature_yaxis.eps', format='eps')
plt.close() 

############################
readin_file='vtgcm.VEN4.v3e.QNIR1p3.pFlds7.dat'

height_readin_file='vtgcm.VEN4.v3e.QNIR1p3.pFlds7.dat'
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(height_readin_file,i+32)
  matrix_2[i,:]=a.split()
a=linecache.getline(height_readin_file,827+32)
matrix_2[827,0:2]=a.split()
height=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    height[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
average_height=np.zeros(68)
for i in range(0,68):
  average_height[i]=np.mean(height[:,i])
  
def forward(xx):
  return np.interp(xx,new_height,average_height[0:68])

def inverse(xx):
  return np.interp(xx,average_height[0:68],new_height)

def forward2(xxx):
  return (5.0*(10.0**-6)/(np.exp(xxx)))

def inverse2(xxx):
  with np.errstate(divide='ignore'):
    return (np.log(5.0*(10.0**-6)/xxx))

def forwardx(xx):
  return (xx/180.0*12.0)

def inversex(xx):
  return (xx*180.0/12.0)

matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+891)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+891)
matrix_2[827,0:2]=a.split()
temperature=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    temperature[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
temperature_log=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    temperature_log[j,i]=math.log10(temperature[j,i])
new_temperature=np.zeros(shape=(73,68))
new_temperature[0:37,:]=temperature[36:73,:]
new_temperature[37:73,:]=temperature[1:37,:]
new_temperature_log=np.zeros(shape=(73,68))
new_temperature_log[0:37,:]=temperature_log[36:73,:]
new_temperature_log[37:73,:]=temperature_log[1:37,:]
    
X,Y=np.meshgrid(new_degree,new_height)
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
lines_lv=[110,120,140,160,180,200,220,240,280,310,340,370,400,430,460,490,520,550,580,610,640]
plt.contourf(X,Y,new_temperature_log.T,100,vmin=math.log10(100),vmax=math.log10(650),cmap=plt.cm.jet)
C=plt.contour(X,Y,new_temperature.T,levels=lines_lv,vmin=100,vmax=600,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=5,colors=('white','white','white','black','black','black','black','black','black','black','black','black','black','black','black','black','white','white','white','white','white'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)

plt.xlabel('Longitude (degrees)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','240','300','360'])

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

plt.savefig('Figure_ven4_1p3_temperature_yaxis.png', dpi=300)

plt.savefig('Figure_ven4_1p3_temperature_yaxis.eps', format='eps')
plt.close() 