import math
import numpy as np
import matplotlib.pyplot as plt
import linecache
from matplotlib.ticker import AutoMinorLocator
####################################################
readin_file='/home/jzl/plot_parkinson_paper/field6/vtgcm.VEN3.Flds6.dat'

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

######################
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+891)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+891)
matrix_2[827,0:2]=a.split()
qnir3=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    qnir3[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+32)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+32)
matrix_2[827,0:2]=a.split()
qeuv3=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    qeuv3[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+1750)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+1750)
matrix_2[827,0:2]=a.split()
qir3=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    qir3[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+2609)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+2609)
matrix_2[827,0:2]=a.split()
cond3=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    cond3[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
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
    
#######################
readin_file='vtgcm.VEN3.v2c.QNIR0p7.sFlds11.dat'
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+891)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+891)
matrix_2[827,0:2]=a.split()
qnir3_0p7=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    qnir3_0p7[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]

matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+32)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+32)
matrix_2[827,0:2]=a.split()
qeuv3_0p7=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    qeuv3_0p7[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+1750)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+1750)
matrix_2[827,0:2]=a.split()
qir3_0p7=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    qir3_0p7[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+2609)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+2609)
matrix_2[827,0:2]=a.split()
cond3_0p7=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    cond3_0p7[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+3468)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+3468)
matrix_2[827,0:2]=a.split()
totdynm_0p7=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    totdynm_0p7[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
###############
readin_file='vtgcm.VEN3.v3c.QNIR1p3.sFlds11.dat'
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+891)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+891)
matrix_2[827,0:2]=a.split()
qnir3_1p3=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    qnir3_1p3[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]

matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+32)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+32)
matrix_2[827,0:2]=a.split()
qeuv3_1p3=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    qeuv3_1p3[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+1750)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+1750)
matrix_2[827,0:2]=a.split()
qir3_1p3=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    qir3_1p3[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+2609)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+2609)
matrix_2[827,0:2]=a.split()
cond3_1p3=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    cond3_1p3[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]
matrix_2=np.zeros(shape=(828,6))
for i in range(0,827):
  a=linecache.getline(readin_file,i+3468)
  matrix_2[i,:]=a.split()
a=linecache.getline(readin_file,827+3468)
matrix_2[827,0:2]=a.split()
totdynm_1p3=np.zeros(shape=(73,68))
for i in range(0,68):
  for j in range(0,73):
    totdynm_1p3[j,i]=matrix_2[(i*73+j)//6,(i*73+j)%6]

################

fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(qnir3[0,:],new_height,linewidth=2,linestyle='-',color='green',label='Near IR heating standard')
plt.plot(qnir3_0p7[0,:],new_height,linewidth=2,linestyle='--',color='red',label='Near IR heating X0.7')
plt.plot(qnir3_1p3[0,:],new_height,linewidth=2,linestyle=':',color='blue',label='Near IR heating X1.3')
plt.legend(fontsize=7)
plt.xlabel('[K/Day]')
plt.xlim([-100,3600])
plt.ylim([-16,16])
plt.ylabel('Ln(P$_{{0}}$/P)')

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

plt.savefig('Figure_nir_VEN3.png', dpi=300)

plt.savefig('Figure_nir_VEN3.eps', format='eps')

fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(qeuv3[0,:],new_height,linewidth=2,linestyle='-',color='green',label='EUV-UV heating standard')
plt.plot(qeuv3_0p7[0,:],new_height,linewidth=2,linestyle='--',color='red',label='EUV-UV heating X0.7')
plt.plot(qeuv3_1p3[0,:],new_height,linewidth=2,linestyle=':',color='blue',label='EUV-UV heating X1.3')
plt.legend(fontsize=7)
plt.xlabel('[K/Day]')
#plt.xlim([-100,3600])
plt.ylim([-16,16])
plt.ylabel('Ln(P$_{{0}}$/P)')

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


plt.savefig('Figure_euv_VEN3.png', dpi=300)

plt.savefig('Figure_euv_VEN3.eps', format='eps')

fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(qir3[0,:],new_height,linewidth=2,linestyle='-',color='green',label='CO$_{{2}}$ 15 ${\mu}$m cooling standard')
plt.plot(qir3_0p7[0,:],new_height,linewidth=2,linestyle='--',color='red',label='CO$_{{2}}$ 15 ${\mu}$m cooling X0.7')
plt.plot(qir3_1p3[0,:],new_height,linewidth=2,linestyle=':',color='blue',label='CO$_{{2}}$ 15 ${\mu}$m cooling X1.3')
plt.legend(fontsize=7)
plt.xlabel('[K/Day]')
#plt.xlim([-100,3600])
plt.ylim([-16,16])
plt.ylabel('Ln(P$_{{0}}$/P)')

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


plt.savefig('Figure_qir_VEN3.png', dpi=300)

plt.savefig('Figure_qir_VEN3.eps', format='eps')

fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(cond3[0,:],new_height,linewidth=2,linestyle='-',color='green',label='Conduction standard')
plt.plot(cond3_0p7[0,:],new_height,linewidth=2,linestyle='--',color='red',label='Conduction X0.7')
plt.plot(cond3_1p3[0,:],new_height,linewidth=2,linestyle=':',color='blue',label='Conduction X1.3')
plt.legend(fontsize=7)
plt.xlabel('[K/Day]')
#plt.xlim([-100,3600])
plt.ylim([-16,16])
plt.ylabel('Ln(P$_{{0}}$/P)')

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


plt.savefig('Figure_cond3_VEN3.png', dpi=300)

plt.savefig('Figure_cond3_VEN3.eps', format='eps')

fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(totdynm[0,:],new_height,linewidth=2,linestyle='-',color='green',label='Total advection standard')
plt.plot(totdynm_0p7[0,:],new_height,linewidth=2,linestyle='--',color='red',label='Total advection X0.7')
plt.plot(totdynm_1p3[0,:],new_height,linewidth=2,linestyle=':',color='blue',label='Total advection X1.3')
plt.legend(fontsize=7)
plt.xlabel('[K/Day]')
#plt.xlim([-100,3600])
plt.ylim([-16,16])
plt.ylabel('Ln(P$_{{0}}$/P)')

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


plt.savefig('Figure_totdynm_VEN3.png', dpi=300)

plt.savefig('Figure_totdynm_VEN3.eps', format='eps')

