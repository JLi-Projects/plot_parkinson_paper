import math
import numpy as np
import matplotlib.pyplot as plt
import linecache
from matplotlib.ticker import AutoMinorLocator

readin_file='vtgcm.VEN1.Flds6.dat'

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
  a=linecache.getline(readin_file,i+4327)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+4327)
matrix_2[827,0:2]=a.split()
noir=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    noir[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
noir_log=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    if(noir[j,i]<0.000001):
      noir_log[j,i]=math.log10(0.000001)
    else:
      noir_log[j,i]=math.log10(noir[j,i])

X,Y=np.meshgrid(degree,new_height)

fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[0.0025,0.0050,0.0075,0.01,0.0125,0.015,0.0175,0.02,0.025,0.03,0.035]
cf_lv=[-0.1,2.5,5,7.5,10,15,20,25,30,35,40,50]
plt.contourf(X,Y,noir_log.T,levels=100,vmin=math.log10(0.0001),vmax=math.log10(0.4),cmap=plt.cm.jet)
C=plt.contour(X,Y,noir.T,levels=lines_lv,linestyles='solid',vmin=math.log10(0.0001),vmax=math.log10(0.4),colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black','black','black','black'),fmt='%.4f',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.ylim(-11,0)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])
ax.yaxis.set_minor_locator(AutoMinorLocator())
secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([90,100,110,120,130])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([8*10**-1,6*10**-1,4*10**-1,2*10**-1,10**-1,8*10**-2,6*10**-2,4*10**-2,2*10**-2,10**-2,8*10**-3,6*10**-3,4*10**-3,2*10**-3,10**-3,8*10**-4,6*10**-4,4*10**-4,2*10**-4,10**-4,8*10**-5,6*10**-5,4*10**-5,2*10**-5,10**-5])
secax_y2.set_yticklabels(('','','','','10$^{{-1}}$','','','','','10$^{{-2}}$','','','','','10$^{{-3}}$','','','','','10$^{{-4}}$','','','','','10$^{{-5}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator(5))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_a_noir_yaxis.png', dpi=300)

plt.savefig('Figure_a_noir.eps', format='eps')
plt.close()

####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+32)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+32)
matrix_2[827,0:2]=a.split()
qeuv2=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    qeuv2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
qeuv2_log=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    if(qeuv2[j,i]<1.0):
      qeuv2_log[j,i]=math.log10(1.0)
    else:
      qeuv2_log[j,i]=math.log10(qeuv2[j,i])

X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[100,300,500,700,900,1100,1300,1500,1700,1900,2100]
cf_lv=[-22,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,0.1]
plt.contourf(X,Y,qeuv2_log.T,levels=100,vmin=math.log10(10),vmax=math.log10(10000),cmap=plt.cm.jet)
C=plt.contour(X,Y,qeuv2.T,levels=lines_lv,linestyles='solid',vmin=math.log10(10),vmax=math.log10(10000),colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black','black','black','black'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.xlim(-69,69)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

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


plt.savefig('Figure_a_qeuv2_yaxis.png', dpi=300)

plt.savefig('Figure_a_qeuv2.eps', format='eps')
plt.close()  

####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+891)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+891)
matrix_2[827,0:2]=a.split()
qnir2=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    qnir2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]


X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400]
cf_lv=[-0.1,200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2800]
plt.contourf(X,Y,qnir2.T,levels=cf_lv,vmin=200,vmax=2600,cmap=plt.cm.jet)
C=plt.contour(X,Y,qnir2.T,levels=lines_lv,linestyles='solid',vmin=200,vmax=2600,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('white','black','black','black','black','black','black','black','black','black','black','white'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.xlim(-69,69)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

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


plt.savefig('Figure_a_qnir2_yaxis.png', dpi=300)

plt.savefig('Figure_a_qnir2.eps', format='eps')
plt.close()  
####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+1750)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+1750)
matrix_2[827,0:2]=a.split()
qir2=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    qir2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]


X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-3000,-2750,-2500,-2250,-2000,-1750,-1500,-1250,-1000,-750,-500,-250]
cf_lv=[-4000,-3000,-2750,-2500,-2250,-2000,-1750,-1500,-1250,-1000,-750,-500,-250,0]
plt.contourf(X,Y,qir2.T,levels=cf_lv,vmin=-11000,vmax=0,cmap=plt.cm.jet)
C=plt.contour(X,Y,qir2.T,levels=lines_lv,linestyles='solid',vmin=-11000,vmax=0,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black','black','black','white','white'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.xlim(-69,69)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

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


plt.savefig('Figure_a_qir2_yaxis.png', dpi=300)

plt.savefig('Figure_a_qir2.eps', format='eps')
plt.close()

####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+2609)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+2609)
matrix_2[827,0:2]=a.split()
cond2=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    cond2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]


X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-4000,-3500,-3000,-2500,-2000,-1500,-1000,-500,-250,0]
cf_lv=[-4500,-4000,-3500,-3000,-2500,-2000,-1500,-1000,-500,-250,0,200]
plt.contourf(X,Y,cond2.T,levels=cf_lv,vmin=-15000,vmax=1000,cmap=plt.cm.jet)
C=plt.contour(X,Y,cond2.T,levels=lines_lv,linestyles='solid',vmin=-15000,vmax=1000,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black','black','black'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.xlim(-69,69)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

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


plt.savefig('Figure_a_cond2_yaxis.png', dpi=300)

plt.savefig('Figure_a_cond2.eps', format='eps')
plt.close()

####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+3468)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+3468)
matrix_2[827,0:2]=a.split()
totdynm=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    totdynm[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]


X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-300,-200,-100,0,100,200,300,400,500,700,900]
cf_lv=[-400,-300,-200,-100,0,100,200,300,400,500,600,700,800,900,1100]
plt.contourf(X,Y,totdynm.T,levels=cf_lv,vmin=-1000,vmax=6000,cmap=plt.cm.jet)
C=plt.contour(X,Y,totdynm.T,levels=lines_lv,linestyles='solid',vmin=-1000,vmax=6000,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black','black','black','white'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.xlim(-69,69)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

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


plt.savefig('Figure_a_totdynm_yaxis.png', dpi=300)

plt.savefig('Figure_a_totdynm.eps', format='eps')
plt.close()

################################################################################################################
readin_file='vtgcm.VEN2.Flds6.dat'

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
  a=linecache.getline(readin_file,i+4327)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+4327)
matrix_2[827,0:2]=a.split()
noir=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    noir[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
noir_log=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    if(noir[j,i]<0.000001):
      noir_log[j,i]=math.log10(0.000001)
    else:
      noir_log[j,i]=math.log10(noir[j,i])

X,Y=np.meshgrid(degree,new_height)

fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[0.0005,0.001,0.0015,0.002,0.003,0.004,0.005,0.006,0.007,0.008]
#cf_lv=[-0.1,2.5,5,7.5,10,15,20,25,30,35,40,50]
plt.contourf(X,Y,noir_log.T,levels=100,vmin=math.log10(0.0001),vmax=math.log10(0.4),cmap=plt.cm.jet)
C=plt.contour(X,Y,noir.T,levels=lines_lv,linestyles='solid',vmin=math.log10(0.0001),vmax=math.log10(0.4),colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black','black','black'),fmt='%.4f',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.ylim(-11,0)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])
ax.yaxis.set_minor_locator(AutoMinorLocator())
secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([90,100,110,120,130])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([8*10**-1,6*10**-1,4*10**-1,2*10**-1,10**-1,8*10**-2,6*10**-2,4*10**-2,2*10**-2,10**-2,8*10**-3,6*10**-3,4*10**-3,2*10**-3,10**-3,8*10**-4,6*10**-4,4*10**-4,2*10**-4,10**-4,8*10**-5,6*10**-5,4*10**-5,2*10**-5,10**-5])
secax_y2.set_yticklabels(('','','','','10$^{{-1}}$','','','','','10$^{{-2}}$','','','','','10$^{{-3}}$','','','','','10$^{{-4}}$','','','','','10$^{{-5}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator(5))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_b_noir_yaxis.png', dpi=300)

plt.savefig('Figure_b_noir.eps', format='eps')
plt.close()

####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+32)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+32)
matrix_2[827,0:2]=a.split()
qeuv2=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    qeuv2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
qeuv2_log=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    if(qeuv2[j,i]<1.0):
      qeuv2_log[j,i]=math.log10(1.0)
    else:
      qeuv2_log[j,i]=math.log10(qeuv2[j,i])

X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[100,200,300,400,500,600,700,800]
cf_lv=[-22,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,0.1]
plt.contourf(X,Y,qeuv2_log.T,levels=100,vmin=math.log10(10),vmax=math.log10(10000),cmap=plt.cm.jet)
C=plt.contour(X,Y,qeuv2.T,levels=lines_lv,linestyles='solid',vmin=math.log10(10),vmax=math.log10(10000),colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.xlim(-69,69)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

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


plt.savefig('Figure_b_qeuv2_yaxis.png', dpi=300)

plt.savefig('Figure_b_qeuv2.eps', format='eps')
plt.close()  

####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+891)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+891)
matrix_2[827,0:2]=a.split()
qnir2=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    qnir2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]


X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400]
cf_lv=[-0.1,200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2800]
plt.contourf(X,Y,qnir2.T,levels=cf_lv,vmin=200,vmax=2600,cmap=plt.cm.jet)
C=plt.contour(X,Y,qnir2.T,levels=lines_lv,linestyles='solid',vmin=200,vmax=2600,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('white','black','black','black','black','black','black','black','black','black','black','white'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.xlim(-69,69)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

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


plt.savefig('Figure_b_qnir2_yaxis.png', dpi=300)

plt.savefig('Figure_b_qnir2.eps', format='eps')
plt.close()  
####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+1750)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+1750)
matrix_2[827,0:2]=a.split()
qir2=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    qir2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]


X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-2400,-2200,-2000,-1800,-1600,-1400,-1200,-1000,-800,-600,-400,-200]
cf_lv=[-2600,-2400,-2200,-2000,-1800,-1600,-1400,-1200,-1000,-800,-600,-400,-200,0]
plt.contourf(X,Y,qir2.T,levels=cf_lv,vmin=-11000,vmax=0,cmap=plt.cm.jet)
C=plt.contour(X,Y,qir2.T,levels=lines_lv,linestyles='solid',vmin=-11000,vmax=0,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black','black','black','white','white'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.xlim(-69,69)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

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


plt.savefig('Figure_b_qir2_yaxis.png', dpi=300)

plt.savefig('Figure_b_qir2.eps', format='eps')
plt.close()

####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+2609)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+2609)
matrix_2[827,0:2]=a.split()
cond2=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    cond2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]


X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-1400,-1200,-1000,-900,-800,-700,-600,-500,-400,-300,-200,-100,0]
cf_lv=[-2400,-1200,-1000,-900,-800,-700,-600,-500,-400,-300,-200,-100,0,100]
plt.contourf(X,Y,cond2.T,levels=cf_lv,vmin=-15000,vmax=1000,cmap=plt.cm.jet)
C=plt.contour(X,Y,cond2.T,levels=lines_lv,linestyles='solid',vmin=-15000,vmax=1000,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black','black','black','black','black','white'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.xlim(-69,69)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

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


plt.savefig('Figure_b_cond2_yaxis.png', dpi=300)

plt.savefig('Figure_b_cond2.eps', format='eps')
plt.close()

####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+3468)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+3468)
matrix_2[827,0:2]=a.split()
totdynm=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    totdynm[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]


X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-300,-200,-100,0,100,200,300,400,500,600]
cf_lv=[-400,-300,-200,-100,0,100,200,300,400,500,600,700]
plt.contourf(X,Y,totdynm.T,levels=cf_lv,vmin=-1000,vmax=6000,cmap=plt.cm.jet)
C=plt.contour(X,Y,totdynm.T,levels=lines_lv,linestyles='solid',vmin=-1000,vmax=6000,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black','black','black'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.xlim(-69,69)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

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


plt.savefig('Figure_b_totdynm_yaxis.png', dpi=300)

plt.savefig('Figure_b_totdynm.eps', format='eps')
plt.close()


################################################################################################################
readin_file='vtgcm.VEN3.Flds6.dat'

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
  a=linecache.getline(readin_file,i+4327)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+4327)
matrix_2[827,0:2]=a.split()
noir=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    noir[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
noir_log=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    if(noir[j,i]<0.000001):
      noir_log[j,i]=math.log10(0.000001)
    else:
      noir_log[j,i]=math.log10(noir[j,i])

X,Y=np.meshgrid(degree,new_height)

fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[0.002,0.004,0.006,0.008,0.010,0.012,0.014,0.018,0.022,0.026]
#cf_lv=[-0.1,2.5,5,7.5,10,15,20,25,30,35,40,50]
plt.contourf(X,Y,noir_log.T,levels=100,vmin=math.log10(0.0001),vmax=math.log10(0.4),cmap=plt.cm.jet)
C=plt.contour(X,Y,noir.T,levels=lines_lv,linestyles='solid',vmin=math.log10(0.0001),vmax=math.log10(0.4),colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black','black','black'),fmt='%.3f',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.ylim(-11,0)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])
ax.yaxis.set_minor_locator(AutoMinorLocator())
secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([90,100,110,120,130])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([8*10**-1,6*10**-1,4*10**-1,2*10**-1,10**-1,8*10**-2,6*10**-2,4*10**-2,2*10**-2,10**-2,8*10**-3,6*10**-3,4*10**-3,2*10**-3,10**-3,8*10**-4,6*10**-4,4*10**-4,2*10**-4,10**-4,8*10**-5,6*10**-5,4*10**-5,2*10**-5,10**-5])
secax_y2.set_yticklabels(('','','','','10$^{{-1}}$','','','','','10$^{{-2}}$','','','','','10$^{{-3}}$','','','','','10$^{{-4}}$','','','','','10$^{{-5}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator(5))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_c_noir_yaxis.png', dpi=300)

plt.savefig('Figure_c_noir.eps', format='eps')
plt.close()

####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+32)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+32)
matrix_2[827,0:2]=a.split()
qeuv2=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    qeuv2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
qeuv2_log=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    if(qeuv2[j,i]<1.0):
      qeuv2_log[j,i]=math.log10(1.0)
    else:
      qeuv2_log[j,i]=math.log10(qeuv2[j,i])

X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2600]
cf_lv=[-22,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,0.1]
plt.contourf(X,Y,qeuv2_log.T,levels=100,vmin=math.log10(10),vmax=math.log10(10000),cmap=plt.cm.jet)
C=plt.contour(X,Y,qeuv2.T,levels=lines_lv,linestyles='solid',vmin=math.log10(10),vmax=math.log10(10000),colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.xlim(-69,69)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

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


plt.savefig('Figure_c_qeuv2_yaxis.png', dpi=300)

plt.savefig('Figure_c_qeuv2.eps', format='eps')
plt.close()  

####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+891)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+891)
matrix_2[827,0:2]=a.split()
qnir2=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    qnir2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]


X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400]
cf_lv=[-0.1,200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2800]
plt.contourf(X,Y,qnir2.T,levels=cf_lv,vmin=200,vmax=2600,cmap=plt.cm.jet)
C=plt.contour(X,Y,qnir2.T,levels=lines_lv,linestyles='solid',vmin=200,vmax=2600,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('white','black','black','black','black','black','black','black','black','black','black','white'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.xlim(-69,69)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

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


plt.savefig('Figure_c_qnir2_yaxis.png', dpi=300)

plt.savefig('Figure_c_qnir2.eps', format='eps')
plt.close()  
####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+1750)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+1750)
matrix_2[827,0:2]=a.split()
qir2=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    qir2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]


X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-3250,-3000,-2750,-2500,-2250,-2000,-1750,-1500,-1250,-1000,-750,-500,-250]
cf_lv=[-3500,-3250,-3000,-2750,-2500,-2250,-2000,-1750,-1500,-1250,-1000,-750,-500,-250,0]
plt.contourf(X,Y,qir2.T,levels=cf_lv,vmin=-11000,vmax=0,cmap=plt.cm.jet)
C=plt.contour(X,Y,qir2.T,levels=lines_lv,linestyles='solid',vmin=-11000,vmax=0,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black','black','black','white','white'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.xlim(-69,69)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

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


plt.savefig('Figure_c_qir2_yaxis.png', dpi=300)

plt.savefig('Figure_c_qir2.eps', format='eps')
plt.close()

####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+2609)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+2609)
matrix_2[827,0:2]=a.split()
cond2=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    cond2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]


X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-5000,-4500,-4000,-3500,-3000,-2500,-2000,-1500,-1000,-500,0]
cf_lv=[-6000,-5000,-4500,-4000,-3500,-3000,-2500,-2000,-1500,-1000,-500,0,500]
plt.contourf(X,Y,cond2.T,levels=cf_lv,vmin=-15000,vmax=1000,cmap=plt.cm.jet)
C=plt.contour(X,Y,cond2.T,levels=lines_lv,linestyles='solid',vmin=-15000,vmax=1000,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black','black','black','black'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.xlim(-69,69)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

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


plt.savefig('Figure_c_cond2_yaxis.png', dpi=300)

plt.savefig('Figure_c_cond2.eps', format='eps')
plt.close()

####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+3468)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+3468)
matrix_2[827,0:2]=a.split()
totdynm=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    totdynm[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]


X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-400,-300,-200,-100,0,100,200,300,400,500,600,700,900,1100]
cf_lv=[-500,-400,-300,-200,-100,0,100,200,300,400,500,600,700,900,1100,1300]
plt.contourf(X,Y,totdynm.T,levels=cf_lv,vmin=-1000,vmax=6000,cmap=plt.cm.jet)
C=plt.contour(X,Y,totdynm.T,levels=lines_lv,linestyles='solid',vmin=-1000,vmax=6000,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black','black','black','white','white'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.xlim(-69,69)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

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


plt.savefig('Figure_c_totdynm_yaxis.png', dpi=300)

plt.savefig('Figure_c_totdynm.eps', format='eps')
plt.close()

################################################################################################################
readin_file='vtgcm.VEN4.Flds6.dat'

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
  a=linecache.getline(readin_file,i+4327)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+4327)
matrix_2[827,0:2]=a.split()
noir=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    noir[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
noir_log=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    if(noir[j,i]<0.000001):
      noir_log[j,i]=math.log10(0.000001)
    else:
      noir_log[j,i]=math.log10(noir[j,i])

X,Y=np.meshgrid(degree,new_height)

fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.25,0.3,0.35,0.4]
#cf_lv=[-0.1,2.5,5,7.5,10,15,20,25,30,35,40,50]
plt.contourf(X,Y,noir_log.T,levels=100,vmin=math.log10(0.0001),vmax=math.log10(0.4),cmap=plt.cm.jet)
C=plt.contour(X,Y,noir.T,levels=lines_lv,linestyles='solid',vmin=math.log10(0.0001),vmax=math.log10(0.4),colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black','black','black','white','white'),fmt='%.3f',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.ylim(-11,0)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])
ax.yaxis.set_minor_locator(AutoMinorLocator())
secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

secax_y = ax.secondary_yaxis('right', functions=(forward,inverse))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([90,100,110,120,130])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([8*10**-1,6*10**-1,4*10**-1,2*10**-1,10**-1,8*10**-2,6*10**-2,4*10**-2,2*10**-2,10**-2,8*10**-3,6*10**-3,4*10**-3,2*10**-3,10**-3,8*10**-4,6*10**-4,4*10**-4,2*10**-4,10**-4,8*10**-5,6*10**-5,4*10**-5,2*10**-5,10**-5])
secax_y2.set_yticklabels(('','','','','10$^{{-1}}$','','','','','10$^{{-2}}$','','','','','10$^{{-3}}$','','','','','10$^{{-4}}$','','','','','10$^{{-5}}$'))
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator(5))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_d_noir_yaxis.png', dpi=300)

plt.savefig('Figure_d_noir.eps', format='eps')
plt.close()

####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+32)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+32)
matrix_2[827,0:2]=a.split()
qeuv2=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    qeuv2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
qeuv2_log=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    if(qeuv2[j,i]<1.0):
      qeuv2_log[j,i]=math.log10(1.0)
    else:
      qeuv2_log[j,i]=math.log10(qeuv2[j,i])

X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[500,1500,2500,3500,4500,5500,6500,7500,8500,9000,9500,10000,10500]
cf_lv=[-22,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,0.1]
plt.contourf(X,Y,qeuv2_log.T,levels=100,vmin=math.log10(10),vmax=math.log10(10000),cmap=plt.cm.jet)
C=plt.contour(X,Y,qeuv2.T,levels=lines_lv,linestyles='solid',vmin=math.log10(10),vmax=math.log10(10000),colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black','black','black','black','black','black'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.xlim(-69,69)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

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


plt.savefig('Figure_d_qeuv2_yaxis.png', dpi=300)

plt.savefig('Figure_d_qeuv2.eps', format='eps')
plt.close()  

####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+891)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+891)
matrix_2[827,0:2]=a.split()
qnir2=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    qnir2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]


X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400]
cf_lv=[-0.1,200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2800]
plt.contourf(X,Y,qnir2.T,levels=cf_lv,vmin=200,vmax=2600,cmap=plt.cm.jet)
C=plt.contour(X,Y,qnir2.T,levels=lines_lv,linestyles='solid',vmin=200,vmax=2600,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('white','black','black','black','black','black','black','black','black','black','black','white'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.xlim(-69,69)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

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


plt.savefig('Figure_d_qnir2_yaxis.png', dpi=300)

plt.savefig('Figure_d_qnir2.eps', format='eps')
plt.close()  
####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+1750)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+1750)
matrix_2[827,0:2]=a.split()
qir2=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    qir2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]


X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-10000,-9000,-8000,-7000,-6000,-5000,-4000,-3000,-2500,-2000,-1500,-1000,-500]
cf_lv=[-12000,-10000,-9000,-8000,-7000,-6000,-5000,-4000,-3000,-2500,-2000,-1500,-1000,-500,0]
plt.contourf(X,Y,qir2.T,levels=cf_lv,vmin=-11000,vmax=0,cmap=plt.cm.jet)
C=plt.contour(X,Y,qir2.T,levels=lines_lv,linestyles='solid',vmin=-11000,vmax=0,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black','black','black','white','white'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.xlim(-69,69)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

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


plt.savefig('Figure_d_qir2_yaxis.png', dpi=300)

plt.savefig('Figure_d_qir2.eps', format='eps')
plt.close()

####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+2609)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+2609)
matrix_2[827,0:2]=a.split()
cond2=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    cond2[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]


X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-16000,-14000,-12000,-10000,-8000,-6000,-4000,-2000,0]
cf_lv=[-30000,-16000,-14000,-12000,-10000,-8000,-6000,-4000,-2000,0,1200]
plt.contourf(X,Y,cond2.T,levels=cf_lv,vmin=-15000,vmax=1000,cmap=plt.cm.jet)
C=plt.contour(X,Y,cond2.T,levels=lines_lv,linestyles='solid',vmin=-15000,vmax=1000,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black','black'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.xlim(-69,69)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

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


plt.savefig('Figure_d_cond2_yaxis.png', dpi=300)

plt.savefig('Figure_d_cond2.eps', format='eps')
plt.close()

####################################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+3468)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+3468)
matrix_2[827,0:2]=a.split()
totdynm=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    totdynm[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]


X,Y=np.meshgrid(degree,new_height)


fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)

lines_lv=[-1000,-750,-500,-250,0,500,1000,2000,3000,4000,5000,6000]
cf_lv=[-1500,-1000,-750,-500,-250,0,500,1000,2000,3000,4000,5000,6000,7000]
plt.contourf(X,Y,totdynm.T,levels=cf_lv,vmin=-1000,vmax=6000,cmap=plt.cm.jet)
C=plt.contour(X,Y,totdynm.T,levels=lines_lv,linestyles='solid',vmin=-1000,vmax=6000,colors='black',linewidths=0.7)
plt.clabel(C,inline=True,inline_spacing=3,colors=('black','black','black','black','black','black','black','black','black','black','white','white'),fmt='%i',fontsize=5.5)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.jet)
#plt.pcolor(X,Y,temperature.T,cmap=plt.cm.nipy_spectral)


plt.xlabel('Longitude (degree)')
plt.ylabel('Ln(P$_{{0}}$/P)')
plt.xlim(-69,69)
my_x_ticks = np.arange(-180, 181, 60)
plt.xticks(my_x_ticks,['-180','-120','-60','0','60','120','180'])

secax_x = ax.secondary_xaxis('top',functions=(forwardx,inversex))
secax_x.set_xlabel('Local Time (hrs)')
secax_x.set_xticks([-12.0,-8.0,-4.0,0.0,4.0,8.0,12.0])
secax_x.set_xticklabels(('12','16','20','0','4','8','12'))

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


plt.savefig('Figure_d_totdynm_yaxis.png', dpi=300)

plt.savefig('Figure_d_totdynm.eps', format='eps')
plt.close()
