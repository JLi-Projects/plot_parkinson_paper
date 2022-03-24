import math
import numpy as np
import matplotlib.pyplot as plt
import linecache
from matplotlib.ticker import AutoMinorLocator

readin_file='vtgcm.VEN1.ZP.5SECFLDS.dat'

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
average_height=np.zeros(69)
for i in range(0,69):
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
  a=linecache.getline(readin_file,i+32)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+32)
matrix_2[827,0:2]=a.split()
o2ir=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    o2ir[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
o2ir_log=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    if(o2ir[j,i]<0.00001):
      o2ir_log[j,i]=math.log10(0.010000001)
    else:
      o2ir_log[j,i]=math.log10(o2ir[j,i])

X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[2.5,5,7.5,10,15,20,25,30,35,40]
cf_lv=[-0.1,2.5,5,7.5,10,15,20,25,30,35,40,50]
plt.contourf(X,Y,o2ir_log.T,levels=100,vmin=math.log10(0.01),vmax=math.log10(50),cmap=plt.cm.jet)
C=plt.contour(X,Y,o2ir.T,levels=lines_lv,linestyles='solid',vmin=math.log10(0.01),vmax=math.log10(50),colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black','black','white'),fmt='%.1f',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.ylim(-11,-5.5)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])
ax.yaxis.set_minor_locator(AutoMinorLocator())
secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([90,100,110,120])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([10**-0,8*10**-1,6*10**-1,4*10**-1,2*10**-1,10**-1,8*10**-2,6*10**-2,4*10**-2,2*10**-2,10**-2,8*10**-3,6*10**-3,4*10**-3,2*10**-3,10**-3])
secax_y2.set_yticklabels(('10$^{{-0}}$','','','','','10$^{{-1}}$','','','','','10$^{{-2}}$','','','','','10$^{{-3}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator(5))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_a_o2ir_yaxis.png', dpi=300)

plt.savefig('Figure_a_o2ir.eps', format='eps')
plt.close()

####################################
  
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+891)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+891)
matrix_2[827,0:2]=a.split()
ohir=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    ohir[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
ohir_log=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    if(ohir[j,i]<0.00001):
      ohir_log[j,i]=math.log10(0.0010000001)
    else:
      ohir_log[j,i]=math.log10(ohir[j,i])

X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
cf_lv=[-0.1,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.5]
plt.contourf(X,Y,ohir_log.T,levels=300,vmin=math.log10(0.001),vmax=math.log10(1.1),cmap=plt.cm.jet)
C=plt.contour(X,Y,ohir.T,levels=lines_lv,linestyles='solid',vmin=math.log10(0.001),vmax=math.log10(1.1),colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black','black','black','black','black','white'),fmt='%.2f',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.ylim(-11,-5.5)
plt.xlim(-90,90)
my_x_ticks = np.arange(-90, 91, 30)
plt.xticks(my_x_ticks,['-90','-60','-30','0','30','60','90'])
ax.yaxis.set_minor_locator(AutoMinorLocator())
secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-4.0,0.0,4.0])
secax_x.set_xticklabels(('20','0','4'))

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([90,100,110,120])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([10**-0,8*10**-1,6*10**-1,4*10**-1,2*10**-1,10**-1,8*10**-2,6*10**-2,4*10**-2,2*10**-2,10**-2,8*10**-3,6*10**-3,4*10**-3,2*10**-3,10**-3])
secax_y2.set_yticklabels(('10$^{{-0}}$','','','','','10$^{{-1}}$','','','','','10$^{{-2}}$','','','','','10$^{{-3}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_a_ohir_yaxis.png', dpi=300)

plt.savefig('Figure_a_ohir.eps', format='eps')
plt.close()  

####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+1750)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+1750)
matrix_2[827,0:2]=a.split()
o3=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    o3[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]

X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4]
cf_lv=[-22,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,0.1]
plt.contourf(X,Y,o3.T,levels=cf_lv,vmin=-20,vmax=-4,cmap=plt.cm.jet)
C=plt.contour(X,Y,o3.T,levels=lines_lv,linestyles='solid',vmin=-20,vmax=-4,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('white','white','black','black','black','black','black','black','black','black','black','black','black','black','white'),fmt='%.1f',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.xlim(-69,69)
my_x_ticks = np.arange(-60, 61, 30)
plt.xticks(my_x_ticks,['-60','-30','0','30','60'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-4.0,0.0,4.0])
secax_x.set_xticklabels(('20','0','4'))

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


plt.savefig('Figure_a_o3_yaxis.png', dpi=300)

plt.savefig('Figure_a_o3.eps', format='eps')
plt.close()  

####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+2609)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+2609)
matrix_2[827,0:2]=a.split()
so=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    so[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
    
new_so=np.zeros(shape=(73,68))
new_so[0:37,:]=so[36:73,:]
new_so[37:73,:]=so[1:37,:]

X,Y=np.meshgrid(new_degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6]
cf_lv=[-20,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5]
plt.contourf(X,Y,new_so.T,levels=cf_lv,vmin=-18,vmax=-4,cmap=plt.cm.jet)
C=plt.contour(X,Y,new_so.T,levels=lines_lv,linestyles='solid',vmin=-18,vmax=-4,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('white','black','black','black','black','black','black','black','black','black','black','black','white'),fmt='%.1f',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(-180, 181, 60)
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','240','300','360'])

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


plt.savefig('Figure_a_so_yaxis.png', dpi=300)

plt.savefig('Figure_a_so.eps', format='eps')
plt.close()  
  
 ####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+3468)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+3468)
matrix_2[827,0:2]=a.split()
so2=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    so2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]

new_so2=np.zeros(shape=(73,68))
new_so2[0:37,:]=so2[36:73,:]
new_so2[37:73,:]=so2[1:37,:]

X,Y=np.meshgrid(new_degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-24,-22,-20,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5]
cf_lv=[-25,-24,-22,-20,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4]
plt.contourf(X,Y,new_so2.T,levels=cf_lv,vmin=-24,vmax=-5,cmap=plt.cm.jet)
C=plt.contour(X,Y,new_so2.T,levels=lines_lv,linestyles='solid',vmin=-24,vmax=-5,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('white','white','black','black','black','black','black','black','black','black','black','black','black','black','black','white','white'),fmt='%.1f',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','240','300','360'])

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


plt.savefig('Figure_a_so2_yaxis.png', dpi=300)

plt.savefig('Figure_a_so2.eps', format='eps')
plt.close()   
  
####################################################################################  
readin_file='vtgcm.VEN2.ZP.5SECFLDS.dat'

degree=np.zeros(73)
read_degree=np.zeros(shape=(13,6))
for i in range(0,12):
  a=linecache.getline(readin_file,i+6)
  read_degree[i,:]=a.split()
a=linecache.getline(readin_file,18)
read_degree[12,0]=a
for i in range(0,73):
  degree[i]=read_degree[i//6,i%6]


new_height=np.zeros(68)
read_height=np.zeros(shape=(12,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+20)
  read_height[i,:]=a.split()
a=linecache.getline(readin_file,11+20)
read_height[11,0:2]=a.split()
for i in range(0,68):
  new_height[i]=read_height[i//6,i%6]
  
    
height_readin_file='/home/jzl/plot_parkinson_paper/field7/vtgcm.VEN2.ZP.FLDS7.dat'
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
average_height=np.zeros(69)
for i in range(0,69):
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
  a=linecache.getline(readin_file,i+32)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+32)
matrix_2[827,0:2]=a.split()
o2ir=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    o2ir[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
o2ir_log=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    if(o2ir[j,i]<0.00001):
      o2ir_log[j,i]=math.log10(0.010000001)
    else:
      o2ir_log[j,i]=math.log10(o2ir[j,i])

X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
cf_lv=[-0.1,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.5]
plt.contourf(X,Y,o2ir_log.T,levels=100,vmin=math.log10(0.01),vmax=math.log10(50),cmap=plt.cm.jet)
C=plt.contour(X,Y,o2ir.T,levels=lines_lv,linestyles='solid',vmin=math.log10(0.01),vmax=math.log10(50),colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('white','black','black','black','black','black','black','black','black','black','black','black'),fmt='%.2f',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.ylim(-11,-5.5)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])
ax.yaxis.set_minor_locator(AutoMinorLocator())
secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([90,100,110,120])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([10**-0,8*10**-1,6*10**-1,4*10**-1,2*10**-1,10**-1,8*10**-2,6*10**-2,4*10**-2,2*10**-2,10**-2,8*10**-3,6*10**-3,4*10**-3,2*10**-3,10**-3])
secax_y2.set_yticklabels(('10$^{{-0}}$','','','','','10$^{{-1}}$','','','','','10$^{{-2}}$','','','','','10$^{{-3}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_b_o2ir_yaxis.png', dpi=300)

plt.savefig('Figure_b_o2ir.eps', format='eps')
plt.close()

####################################
  
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+891)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+891)
matrix_2[827,0:2]=a.split()
ohir=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    ohir[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
ohir_log=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    if(ohir[j,i]<0.00001):
      ohir_log[j,i]=math.log10(0.0010000001)
    else:
      ohir_log[j,i]=math.log10(ohir[j,i])

X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[0.005,0.01,0.015,0.02,0.025,0.03,0.04,0.05,0.06,0.07,0.08]

plt.contourf(X,Y,ohir_log.T,levels=300,vmin=math.log10(0.001),vmax=math.log10(1.1),cmap=plt.cm.jet)
C=plt.contour(X,Y,ohir.T,levels=lines_lv,linestyles='solid',vmin=math.log10(0.001),vmax=math.log10(1.1),colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('white','black','black','black','black','black','black','black','black','black','black'),fmt='%.2f',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.ylim(-11,-5.5)
plt.xlim(-90,90)
my_x_ticks = np.arange(-90, 91, 30)
plt.xticks(my_x_ticks,['-90','-60','-30','0','30','60','90'])
ax.yaxis.set_minor_locator(AutoMinorLocator())
secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-4.0,0.0,4.0])
secax_x.set_xticklabels(('20','0','4'))

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([90,100,110,120])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([10**-0,8*10**-1,6*10**-1,4*10**-1,2*10**-1,10**-1,8*10**-2,6*10**-2,4*10**-2,2*10**-2,10**-2,8*10**-3,6*10**-3,4*10**-3,2*10**-3,10**-3])
secax_y2.set_yticklabels(('10$^{{-0}}$','','','','','10$^{{-1}}$','','','','','10$^{{-2}}$','','','','','10$^{{-3}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_b_ohir_yaxis.png', dpi=300)

plt.savefig('Figure_b_ohir.eps', format='eps')
plt.close()  

####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+1750)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+1750)
matrix_2[827,0:2]=a.split()
o3=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    o3[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]

X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4]
cf_lv=[-22,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,0.1]
plt.contourf(X,Y,o3.T,levels=cf_lv,vmin=-20,vmax=-4,cmap=plt.cm.jet)
C=plt.contour(X,Y,o3.T,levels=lines_lv,linestyles='solid',vmin=-20,vmax=-4,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('white','white','black','black','black','black','black','black','black','black','black','black','black','black','white'),fmt='%.1f',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.xlim(-69,69)
my_x_ticks = np.arange(-60, 61, 30)
plt.xticks(my_x_ticks,['-60','-30','0','30','60'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-4.0,0.0,4.0])
secax_x.set_xticklabels(('20','0','4'))

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


plt.savefig('Figure_b_o3_yaxis.png', dpi=300)

plt.savefig('Figure_b_o3.eps', format='eps')
plt.close()  

####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+2609)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+2609)
matrix_2[827,0:2]=a.split()
so=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    so[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]

new_so=np.zeros(shape=(73,68))
new_so[0:37,:]=so[36:73,:]
new_so[37:73,:]=so[1:37,:]

X,Y=np.meshgrid(new_degree,new_height)

fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6]
cf_lv=[-20,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5]
plt.contourf(X,Y,new_so.T,levels=cf_lv,vmin=-18,vmax=-4,cmap=plt.cm.jet)
C=plt.contour(X,Y,new_so.T,levels=lines_lv,linestyles='solid',vmin=-18,vmax=-4,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('white','black','black','black','black','black','black','black','black','black','black','black','white'),fmt='%.1f',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','240','300','360'])

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


plt.savefig('Figure_b_so_yaxis.png', dpi=300)

plt.savefig('Figure_b_so.eps', format='eps')
plt.close()  
  
 ####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+3468)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+3468)
matrix_2[827,0:2]=a.split()
so2=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    so2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]

new_so2=np.zeros(shape=(73,68))
new_so2[0:37,:]=so2[36:73,:]
new_so2[37:73,:]=so2[1:37,:]

X,Y=np.meshgrid(new_degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-24,-22,-20,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5]
cf_lv=[-25,-24,-22,-20,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4]
plt.contourf(X,Y,new_so2.T,levels=cf_lv,vmin=-24,vmax=-5,cmap=plt.cm.jet)
C=plt.contour(X,Y,new_so2.T,levels=lines_lv,linestyles='solid',vmin=-24,vmax=-5,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('white','white','black','black','black','black','black','black','black','black','black','black','black','black','black','white','white'),fmt='%.1f',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','240','300','360'])

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


plt.savefig('Figure_b_so2_yaxis.png', dpi=300)

plt.savefig('Figure_b_so2.eps', format='eps')
plt.close()     
  
####################################################################################  
readin_file='vtgcm.VEN3.ZP.5SECFLDS.dat'

degree=np.zeros(73)
read_degree=np.zeros(shape=(13,6))
for i in range(0,12):
  a=linecache.getline(readin_file,i+6)
  read_degree[i,:]=a.split()
a=linecache.getline(readin_file,18)
read_degree[12,0]=a
for i in range(0,73):
  degree[i]=read_degree[i//6,i%6]


new_height=np.zeros(68)
read_height=np.zeros(shape=(12,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+20)
  read_height[i,:]=a.split()
a=linecache.getline(readin_file,11+20)
read_height[11,0:2]=a.split()
for i in range(0,68):
  new_height[i]=read_height[i//6,i%6]
  
    
height_readin_file='/home/jzl/plot_parkinson_paper/field7/vtgcm.VEN3.ZP.FLDS7.dat'
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
average_height=np.zeros(69)
for i in range(0,69):
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
  a=linecache.getline(readin_file,i+32)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+32)
matrix_2[827,0:2]=a.split()
o2ir=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    o2ir[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
o2ir_log=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    if(o2ir[j,i]<0.00001):
      o2ir_log[j,i]=math.log10(0.010000001)
    else:
      o2ir_log[j,i]=math.log10(o2ir[j,i])

X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0]
plt.contourf(X,Y,o2ir_log.T,levels=100,vmin=math.log10(0.01),vmax=math.log10(50),cmap=plt.cm.jet)
C=plt.contour(X,Y,o2ir.T,levels=lines_lv,linestyles='solid',vmin=math.log10(0.01),vmax=math.log10(50),colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black','black'),fmt='%.1f',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.ylim(-11,-5.5)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])
ax.yaxis.set_minor_locator(AutoMinorLocator())
secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([90,100,110,120])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([10**-0,8*10**-1,6*10**-1,4*10**-1,2*10**-1,10**-1,8*10**-2,6*10**-2,4*10**-2,2*10**-2,10**-2,8*10**-3,6*10**-3,4*10**-3,2*10**-3,10**-3])
secax_y2.set_yticklabels(('10$^{{-0}}$','','','','','10$^{{-1}}$','','','','','10$^{{-2}}$','','','','','10$^{{-3}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_c_o2ir_yaxis.png', dpi=300)

plt.savefig('Figure_c_o2ir.eps', format='eps')
plt.close()

####################################
  
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+891)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+891)
matrix_2[827,0:2]=a.split()
ohir=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    ohir[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
ohir_log=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    if(ohir[j,i]<0.00001):
      ohir_log[j,i]=math.log10(0.0010000001)
    else:
      ohir_log[j,i]=math.log10(ohir[j,i])

X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.20,0.22,0.24]

plt.contourf(X,Y,ohir_log.T,levels=300,vmin=math.log10(0.001),vmax=math.log10(1.1),cmap=plt.cm.jet)
C=plt.contour(X,Y,ohir.T,levels=lines_lv,linestyles='solid',vmin=math.log10(0.001),vmax=math.log10(1.1),colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black','black','black','black','black'),fmt='%.2f',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.ylim(-11,-5.5)
plt.xlim(-90,90)
my_x_ticks = np.arange(-90, 91, 30)
plt.xticks(my_x_ticks,['-90','-60','-30','0','30','60','90'])
ax.yaxis.set_minor_locator(AutoMinorLocator())
secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-4.0,0.0,4.0])
secax_x.set_xticklabels(('20','0','4'))

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([90,100,110,120])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([10**-0,8*10**-1,6*10**-1,4*10**-1,2*10**-1,10**-1,8*10**-2,6*10**-2,4*10**-2,2*10**-2,10**-2,8*10**-3,6*10**-3,4*10**-3,2*10**-3,10**-3])
secax_y2.set_yticklabels(('10$^{{-0}}$','','','','','10$^{{-1}}$','','','','','10$^{{-2}}$','','','','','10$^{{-3}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_c_ohir_yaxis.png', dpi=300)

plt.savefig('Figure_c_ohir.eps', format='eps')
plt.close()  

####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+1750)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+1750)
matrix_2[827,0:2]=a.split()
o3=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    o3[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]

X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4]
cf_lv=[-22,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,0.1]
plt.contourf(X,Y,o3.T,levels=cf_lv,vmin=-20,vmax=-4,cmap=plt.cm.jet)
C=plt.contour(X,Y,o3.T,levels=lines_lv,linestyles='solid',vmin=-20,vmax=-4,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('white','white','black','black','black','black','black','black','black','black','black','black','black','black','white'),fmt='%.1f',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.xlim(-69,69)
my_x_ticks = np.arange(-60, 61, 30)
plt.xticks(my_x_ticks,['-60','-30','0','30','60'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-4.0,0.0,4.0])
secax_x.set_xticklabels(('20','0','4'))

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


plt.savefig('Figure_c_o3_yaxis.png', dpi=300)

plt.savefig('Figure_c_o3.eps', format='eps')
plt.close()  

####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+2609)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+2609)
matrix_2[827,0:2]=a.split()
so=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    so[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]

new_so=np.zeros(shape=(73,68))
new_so[0:37,:]=so[36:73,:]
new_so[37:73,:]=so[1:37,:]

X,Y=np.meshgrid(new_degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6]
cf_lv=[-20,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5]
plt.contourf(X,Y,new_so.T,levels=cf_lv,vmin=-18,vmax=-4,cmap=plt.cm.jet)
C=plt.contour(X,Y,new_so.T,levels=lines_lv,linestyles='solid',vmin=-18,vmax=-4,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('white','black','black','black','black','black','black','black','black','black','black','black','white'),fmt='%.1f',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','240','300','360'])

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


plt.savefig('Figure_c_so_yaxis.png', dpi=300)

plt.savefig('Figure_c_so.eps', format='eps')
plt.close()  
  
 ####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+3468)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+3468)
matrix_2[827,0:2]=a.split()
so2=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    so2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]

new_so2=np.zeros(shape=(73,68))
new_so2[0:37,:]=so2[36:73,:]
new_so2[37:73,:]=so2[1:37,:]

X,Y=np.meshgrid(new_degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-24,-22,-20,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5]
cf_lv=[-25,-24,-22,-20,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4]
plt.contourf(X,Y,new_so2.T,levels=cf_lv,vmin=-24,vmax=-5,cmap=plt.cm.jet)
C=plt.contour(X,Y,new_so2.T,levels=lines_lv,linestyles='solid',vmin=-24,vmax=-5,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('white','white','black','black','black','black','black','black','black','black','black','black','black','black','black','white','white'),fmt='%.1f',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
my_x_ticks = np.arange(0, 361, 60)

plt.xticks(my_x_ticks,['0','60','120','180','240','300','360'])

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


plt.savefig('Figure_c_so2_yaxis.png', dpi=300)

plt.savefig('Figure_c_so2.eps', format='eps')
plt.close()      
  
##########################################
readin_file='vtgcm.VEN4.ZP.5SECFLDS.dat'

degree=np.zeros(73)
read_degree=np.zeros(shape=(13,6))
for i in range(0,12):
  a=linecache.getline(readin_file,i+6)
  read_degree[i,:]=a.split()
a=linecache.getline(readin_file,18)
read_degree[12,0]=a
for i in range(0,73):
  degree[i]=read_degree[i//6,i%6]


new_height=np.zeros(68)
read_height=np.zeros(shape=(12,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+20)
  read_height[i,:]=a.split()
a=linecache.getline(readin_file,11+20)
read_height[11,0:2]=a.split()
for i in range(0,68):
  new_height[i]=read_height[i//6,i%6]
  
    
height_readin_file='/home/jzl/plot_parkinson_paper/field7/vtgcm.VEN4.ZP.FLDS7.dat'
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
average_height=np.zeros(69)
for i in range(0,69):
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
  a=linecache.getline(readin_file,i+32)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+32)
matrix_2[827,0:2]=a.split()
o2ir=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    o2ir[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
o2ir_log=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    if(o2ir[j,i]<0.00001):
      o2ir_log[j,i]=math.log10(0.010000001)
    else:
      o2ir_log[j,i]=math.log10(o2ir[j,i])

X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[4,8,12,16,20,24,28,32,36,40,44]

plt.contourf(X,Y,o2ir_log.T,levels=100,vmin=math.log10(0.01),vmax=math.log10(50),cmap=plt.cm.jet)
C=plt.contour(X,Y,o2ir.T,levels=lines_lv,linestyles='solid',vmin=math.log10(0.01),vmax=math.log10(50),colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black','black','white','white'),fmt='%.1f',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.ylim(-11,-5.5)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])
ax.yaxis.set_minor_locator(AutoMinorLocator())
secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([90,100,110,120])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([10**-0,8*10**-1,6*10**-1,4*10**-1,2*10**-1,10**-1,8*10**-2,6*10**-2,4*10**-2,2*10**-2,10**-2,8*10**-3,6*10**-3,4*10**-3,2*10**-3,10**-3])
secax_y2.set_yticklabels(('10$^{{-0}}$','','','','','10$^{{-1}}$','','','','','10$^{{-2}}$','','','','','10$^{{-3}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_d_o2ir_yaxis.png', dpi=300)

plt.savefig('Figure_d_o2ir.eps', format='eps')
plt.close()

####################################
  
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+891)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+891)
matrix_2[827,0:2]=a.split()
ohir=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    ohir[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
ohir_log=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    if(ohir[j,i]<0.00001):
      ohir_log[j,i]=math.log10(0.0010000001)
    else:
      ohir_log[j,i]=math.log10(ohir[j,i])

X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
cf_lv=[-0.1,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.5]
plt.contourf(X,Y,ohir_log.T,levels=300,vmin=math.log10(0.001),vmax=math.log10(1.1),cmap=plt.cm.jet)
C=plt.contour(X,Y,ohir.T,levels=lines_lv,linestyles='solid',vmin=math.log10(0.001),vmax=math.log10(1.1),colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black','black','black','black','black','white'),fmt='%.2f',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.ylim(-11,-5.5)
plt.xlim(-90,90)
my_x_ticks = np.arange(-90, 91, 30)
plt.xticks(my_x_ticks,['-90','-60','-30','0','30','60','90'])
ax.yaxis.set_minor_locator(AutoMinorLocator())
secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-4.0,0.0,4.0])
secax_x.set_xticklabels(('20','0','4'))

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([90,100,110,120])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([10**-0,8*10**-1,6*10**-1,4*10**-1,2*10**-1,10**-1,8*10**-2,6*10**-2,4*10**-2,2*10**-2,10**-2,8*10**-3,6*10**-3,4*10**-3,2*10**-3,10**-3])
secax_y2.set_yticklabels(('10$^{{-0}}$','','','','','10$^{{-1}}$','','','','','10$^{{-2}}$','','','','','10$^{{-3}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_d_ohir_yaxis.png', dpi=300)

plt.savefig('Figure_d_ohir.eps', format='eps')
plt.close()  

####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+1750)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+1750)
matrix_2[827,0:2]=a.split()
o3=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    o3[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]

X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4]
cf_lv=[-22,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,0.1]
plt.contourf(X,Y,o3.T,levels=cf_lv,vmin=-20,vmax=-4,cmap=plt.cm.jet)
C=plt.contour(X,Y,o3.T,levels=lines_lv,linestyles='solid',vmin=-20,vmax=-4,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('white','white','black','black','black','black','black','black','black','black','black','black','black','black','white'),fmt='%.1f',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.xlim(-69,69)
my_x_ticks = np.arange(-60, 61, 30)
plt.xticks(my_x_ticks,['-60','-30','0','30','60'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-4.0,0.0,4.0])
secax_x.set_xticklabels(('20','0','4'))

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,140,210,280,350,420])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_d_o3_yaxis.png', dpi=300)

plt.savefig('Figure_d_o3.eps', format='eps')
plt.close()  

####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+2609)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+2609)
matrix_2[827,0:2]=a.split()
so=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    so[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]

new_so=np.zeros(shape=(73,68))
new_so[0:37,:]=so[36:73,:]
new_so[37:73,:]=so[1:37,:]

X,Y=np.meshgrid(new_degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6]
cf_lv=[-20,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5]
plt.contourf(X,Y,new_so.T,levels=cf_lv,vmin=-18,vmax=-4,cmap=plt.cm.jet)
C=plt.contour(X,Y,new_so.T,levels=lines_lv,linestyles='solid',vmin=-18,vmax=-4,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('white','black','black','black','black','black','black','black','black','black','black','black','white'),fmt='%.1f',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
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
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_d_so_yaxis.png', dpi=300)

plt.savefig('Figure_d_so.eps', format='eps')
plt.close()  
  
 ####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+3468)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+3468)
matrix_2[827,0:2]=a.split()
so2=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    so2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]

new_so2=np.zeros(shape=(73,68))
new_so2[0:37,:]=so2[36:73,:]
new_so2[37:73,:]=so2[1:37,:]

X,Y=np.meshgrid(new_degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-24,-22,-20,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5]
cf_lv=[-25,-24,-22,-20,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4]
plt.contourf(X,Y,new_so2.T,levels=cf_lv,vmin=-24,vmax=-5,cmap=plt.cm.jet)
C=plt.contour(X,Y,new_so2.T,levels=lines_lv,linestyles='solid',vmin=-24,vmax=-5,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('white','white','black','black','black','black','black','black','black','black','black','black','black','black','black','white','white'),fmt='%.1f',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
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
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_d_so2_yaxis.png', dpi=300)

plt.savefig('Figure_d_so2.eps', format='eps')
plt.close()     
  
  
  
  
  
  