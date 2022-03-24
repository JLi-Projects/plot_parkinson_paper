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
###############################
readin_file='/home/jzl/plot_parkinson_paper/new_20211213/vtgcm.VEN1.JO2JSO2aJSO2b.dat'

RJO2=np.zeros(55)
read_RJO2=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+18)
  read_RJO2[i,:]=a.split()
a=linecache.getline(readin_file,9+18)
read_RJO2[9,:]=a.split()
for i in range(0,55):
  RJO2[i]=read_RJO2[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(RJO2,new_height,linewidth=2,linestyle='-',color='black',label='JO$_{{2}}$')

plt.legend(fontsize=7)
plt.xlabel('JO$_{{2}}$')
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


plt.savefig('Figure_JO2.png', dpi=300)

plt.savefig('Figure_JO2.eps', format='eps')

RJSO2a=np.zeros(55)
read_RJSO2a=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+46)
  read_RJSO2a[i,:]=a.split()
a=linecache.getline(readin_file,9+46)
read_RJSO2a[9,:]=a.split()
for i in range(0,55):
  RJSO2a[i]=read_RJSO2a[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(RJSO2a,new_height,linewidth=2,linestyle='-',color='black',label='JSO$_{{2}}$a')

plt.legend(fontsize=7)
plt.xlabel('JSO$_{{2}}$a')
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


plt.savefig('Figure_JSO2a.png', dpi=300)

plt.savefig('Figure_JSO2a.eps', format='eps')

RJSO2b=np.zeros(55)
read_RJSO2b=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+74)
  read_RJSO2b[i,:]=a.split()
a=linecache.getline(readin_file,9+74)
read_RJSO2b[9,:]=a.split()
for i in range(0,55):
  RJSO2b[i]=read_RJSO2b[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(RJSO2b,new_height,linewidth=2,linestyle='-',color='black',label='JSO$_{{2}}$b')

plt.legend(fontsize=7)
plt.xlabel('JSO$_{{2}}$b')
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


plt.savefig('Figure_JSO2b.png', dpi=300)

plt.savefig('Figure_JSO2b.eps', format='eps')



###############################
readin_file='/home/jzl/plot_parkinson_paper/new_20211213/vtgcm.VEN1.O3.LT24.dat'

O3SRCA=np.zeros(55)
read_O3SRCA=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+18)
  read_O3SRCA[i,:]=a.split()
a=linecache.getline(readin_file,9+18)
read_O3SRCA[9,:]=a.split()
for i in range(0,55):
  O3SRCA[i]=read_O3SRCA[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(O3SRCA,new_height,linewidth=2,linestyle='-',color='black',label='O$_{{3}}$SRCA')

plt.legend(fontsize=7)
plt.xlabel('O$_{{3}}$SRCA')
#plt.xlim([-100,3600])
plt.ylim([-16,11])
plt.ylabel('Ln(P$_{{0}}$/P)')

secax_y = ax.secondary_yaxis('right', functions=(forward_24,inverse_24))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,100,130,160])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_O3SRCA.png', dpi=300)

plt.savefig('Figure_O3SRCA.eps', format='eps')

O3SRCB=np.zeros(55)
read_O3SRCB=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+46)
  read_O3SRCB[i,:]=a.split()
a=linecache.getline(readin_file,9+46)
read_O3SRCB[9,:]=a.split()
for i in range(0,55):
  O3SRCB[i]=read_O3SRCB[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(O3SRCB,new_height,linewidth=2,linestyle='-',color='black',label='O$_{{3}}$SRCB')

plt.legend(fontsize=7)
plt.xlabel('O$_{{3}}$SRCB')
#plt.xlim([-100,3600])
plt.ylim([-16,11])
plt.ylabel('Ln(P$_{{0}}$/P)')

secax_y = ax.secondary_yaxis('right', functions=(forward_24,inverse_24))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,100,130,160])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_O3SRCB.png', dpi=300)

plt.savefig('Figure_O3SRCB.eps', format='eps')

O3SRCC=np.zeros(55)
read_O3SRCC=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+74)
  read_O3SRCC[i,:]=a.split()
a=linecache.getline(readin_file,9+74)
read_O3SRCC[9,:]=a.split()
for i in range(0,55):
  O3SRCC[i]=read_O3SRCC[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(O3SRCC,new_height,linewidth=2,linestyle='-',color='black',label='O$_{{3}}$SRCC')

plt.legend(fontsize=7)
plt.xlabel('O$_{{3}}$SRCC')
#plt.xlim([-100,3600])
plt.ylim([-16,11])
plt.ylabel('Ln(P$_{{0}}$/P)')

secax_y = ax.secondary_yaxis('right', functions=(forward_24,inverse_24))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,100,130,160])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_O3SRCC.png', dpi=300)

plt.savefig('Figure_O3SRCC.eps', format='eps')

O3SRCD=np.zeros(55)
read_O3SRCD=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+102)
  read_O3SRCD[i,:]=a.split()
a=linecache.getline(readin_file,9+102)
read_O3SRCD[9,:]=a.split()
for i in range(0,55):
  O3SRCD[i]=read_O3SRCD[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(O3SRCD,new_height,linewidth=2,linestyle='-',color='black',label='O$_{{3}}$SRCD')

plt.legend(fontsize=7)
plt.xlabel('O$_{{3}}$SRCD')
#plt.xlim([-100,3600])
plt.ylim([-16,11])
plt.ylabel('Ln(P$_{{0}}$/P)')

secax_y = ax.secondary_yaxis('right', functions=(forward_24,inverse_24))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,100,130,160])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_O3SRCD.png', dpi=300)

plt.savefig('Figure_O3SRCD.eps', format='eps')


O3SRCE=np.zeros(55)
read_O3SRCE=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+130)
  read_O3SRCE[i,:]=a.split()
a=linecache.getline(readin_file,9+130)
read_O3SRCE[9,:]=a.split()
for i in range(0,55):
  O3SRCE[i]=read_O3SRCE[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(O3SRCE,new_height,linewidth=2,linestyle='-',color='black',label='O$_{{3}}$SRCE')

plt.legend(fontsize=7)
plt.xlabel('O$_{{3}}$SRCE')
#plt.xlim([-100,3600])
plt.ylim([-16,11])
plt.ylabel('Ln(P$_{{0}}$/P)')

secax_y = ax.secondary_yaxis('right', functions=(forward_24,inverse_24))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,100,130,160])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_O3SRCE.png', dpi=300)

plt.savefig('Figure_O3SRCE.eps', format='eps')

O3SRCT=np.zeros(55)
read_O3SRCT=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+158)
  read_O3SRCT[i,:]=a.split()
a=linecache.getline(readin_file,9+158)
read_O3SRCT[9,:]=a.split()
for i in range(0,55):
  O3SRCT[i]=read_O3SRCT[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(O3SRCT,new_height,linewidth=2,linestyle='-',color='black',label='O$_{{3}}$SRCT')

plt.legend(fontsize=7)
plt.xlabel('O$_{{3}}$SRCT')
#plt.xlim([-100,3600])
plt.ylim([-16,11])
plt.ylabel('Ln(P$_{{0}}$/P)')

secax_y = ax.secondary_yaxis('right', functions=(forward_24,inverse_24))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,100,130,160])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_O3SRCT.png', dpi=300)

plt.savefig('Figure_O3SRCT.eps', format='eps')

O3LOSSA=np.zeros(55)
read_O3LOSSA=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+186)
  read_O3LOSSA[i,:]=a.split()
a=linecache.getline(readin_file,9+186)
read_O3LOSSA[9,:]=a.split()
for i in range(0,55):
  O3LOSSA[i]=read_O3LOSSA[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(O3LOSSA,new_height,linewidth=2,linestyle='-',color='black',label='O$_{{3}}$LOSSA')

plt.legend(fontsize=7)
plt.xlabel('O$_{{3}}$LOSSA')
#plt.xlim([-100,3600])
plt.ylim([-16,11])
plt.ylabel('Ln(P$_{{0}}$/P)')

secax_y = ax.secondary_yaxis('right', functions=(forward_24,inverse_24))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,100,130,160])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_O3LOSSA.png', dpi=300)

plt.savefig('Figure_O3LOSSA.eps', format='eps')

O3LOSSB=np.zeros(55)
read_O3LOSSB=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+214)
  read_O3LOSSB[i,:]=a.split()
a=linecache.getline(readin_file,9+214)
read_O3LOSSB[9,:]=a.split()
for i in range(0,55):
  O3LOSSB[i]=read_O3LOSSB[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(O3LOSSB,new_height,linewidth=2,linestyle='-',color='black',label='O$_{{3}}$LOSSB')

plt.legend(fontsize=7)
plt.xlabel('O$_{{3}}$LOSSB')
#plt.xlim([-100,3600])
plt.ylim([-16,11])
plt.ylabel('Ln(P$_{{0}}$/P)')

secax_y = ax.secondary_yaxis('right', functions=(forward_24,inverse_24))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,100,130,160])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_O3LOSSB.png', dpi=300)

plt.savefig('Figure_O3LOSSB.eps', format='eps')

O3LOSSC=np.zeros(55)
read_O3LOSSC=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+242)
  read_O3LOSSC[i,:]=a.split()
a=linecache.getline(readin_file,9+242)
read_O3LOSSC[9,:]=a.split()
for i in range(0,55):
  O3LOSSC[i]=read_O3LOSSC[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(O3LOSSC,new_height,linewidth=2,linestyle='-',color='black',label='O$_{{3}}$LOSSC')

plt.legend(fontsize=7)
plt.xlabel('O$_{{3}}$LOSSC')
#plt.xlim([-100,3600])
plt.ylim([-16,11])
plt.ylabel('Ln(P$_{{0}}$/P)')

secax_y = ax.secondary_yaxis('right', functions=(forward_24,inverse_24))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,100,130,160])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_O3LOSSC.png', dpi=300)

plt.savefig('Figure_O3LOSSC.eps', format='eps')

O3LOSSD=np.zeros(55)
read_O3LOSSD=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+270)
  read_O3LOSSD[i,:]=a.split()
a=linecache.getline(readin_file,9+270)
read_O3LOSSD[9,:]=a.split()
for i in range(0,55):
  O3LOSSD[i]=read_O3LOSSD[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(O3LOSSD,new_height,linewidth=2,linestyle='-',color='black',label='O$_{{3}}$LOSSD')

plt.legend(fontsize=7)
plt.xlabel('O$_{{3}}$LOSSD')
#plt.xlim([-100,3600])
plt.ylim([-16,11])
plt.ylabel('Ln(P$_{{0}}$/P)')

secax_y = ax.secondary_yaxis('right', functions=(forward_24,inverse_24))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,100,130,160])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_O3LOSSD.png', dpi=300)

plt.savefig('Figure_O3LOSSD.eps', format='eps')


O3LOSSE=np.zeros(55)
read_O3LOSSE=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+298)
  read_O3LOSSE[i,:]=a.split()
a=linecache.getline(readin_file,9+298)
read_O3LOSSE[9,:]=a.split()
for i in range(0,55):
  O3LOSSE[i]=read_O3LOSSE[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(O3LOSSE,new_height,linewidth=2,linestyle='-',color='black',label='O$_{{3}}$LOSSE')

plt.legend(fontsize=7)
plt.xlabel('O$_{{3}}$LOSSE')
#plt.xlim([-100,3600])
plt.ylim([-16,11])
plt.ylabel('Ln(P$_{{0}}$/P)')

secax_y = ax.secondary_yaxis('right', functions=(forward_24,inverse_24))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,100,130,160])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_O3LOSSE.png', dpi=300)

plt.savefig('Figure_O3LOSSE.eps', format='eps')

O3LOSSG=np.zeros(55)
read_O3LOSSG=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+326)
  read_O3LOSSG[i,:]=a.split()
a=linecache.getline(readin_file,9+326)
read_O3LOSSG[9,:]=a.split()
for i in range(0,55):
  O3LOSSG[i]=read_O3LOSSG[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(O3LOSSG,new_height,linewidth=2,linestyle='-',color='black',label='O$_{{3}}$LOSSG')

plt.legend(fontsize=7)
plt.xlabel('O$_{{3}}$LOSSG')
#plt.xlim([-100,3600])
plt.ylim([-16,11])
plt.ylabel('Ln(P$_{{0}}$/P)')

secax_y = ax.secondary_yaxis('right', functions=(forward_24,inverse_24))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,100,130,160])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_O3LOSSG.png', dpi=300)

plt.savefig('Figure_O3LOSSG.eps', format='eps')

O3LOSST=np.zeros(55)
read_O3LOSST=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+354)
  read_O3LOSST[i,:]=a.split()
a=linecache.getline(readin_file,9+354)
read_O3LOSST[9,:]=a.split()
for i in range(0,55):
  O3LOSST[i]=read_O3LOSST[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(O3LOSST,new_height,linewidth=2,linestyle='-',color='black',label='O$_{{3}}$LOSST')

plt.legend(fontsize=7)
plt.xlabel('O$_{{3}}$LOSST')
#plt.xlim([-100,3600])
plt.ylim([-16,11])
plt.ylabel('Ln(P$_{{0}}$/P)')

secax_y = ax.secondary_yaxis('right', functions=(forward_24,inverse_24))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,100,130,160])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_O3LOSST.png', dpi=300)

plt.savefig('Figure_O3LOSST.eps', format='eps')

###############################
readin_file='/home/jzl/plot_parkinson_paper/new_20211213/vtgcm.VEN1.O2.LT12&24.dat'

O2SRCA=np.zeros(55)
read_O2SRCA=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+18)
  read_O2SRCA[i,:]=a.split()
a=linecache.getline(readin_file,9+18)
read_O2SRCA[9,:]=a.split()
for i in range(0,55):
  O2SRCA[i]=read_O2SRCA[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(O2SRCA,new_height,linewidth=2,linestyle='-',color='black',label='O$_{{2}}$SRCA')

plt.legend(fontsize=7)
plt.xlabel('O$_{{2}}$SRCA')
#plt.xlim([-100,3600])
plt.ylim([-16,11])
plt.ylabel('Ln(P$_{{0}}$/P)')

secax_y = ax.secondary_yaxis('right', functions=(forward_24,inverse_24))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,100,130,160])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_O2SRCA_24.png', dpi=300)

plt.savefig('Figure_O2SRCA_24.eps', format='eps')

O2SRCA=np.zeros(55)
read_O2SRCA=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+46)
  read_O2SRCA[i,:]=a.split()
a=linecache.getline(readin_file,9+46)
read_O2SRCA[9,:]=a.split()
for i in range(0,55):
  O2SRCA[i]=read_O2SRCA[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(O2SRCA,new_height,linewidth=2,linestyle='-',color='black',label='O$_{{2}}$SRCA')

plt.legend(fontsize=7)
plt.xlabel('O$_{{2}}$SRCA')
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


plt.savefig('Figure_O2SRCA_12.png', dpi=300)

plt.savefig('Figure_O2SRCA_12.eps', format='eps')

O2SRCB=np.zeros(55)
read_O2SRCB=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+74)
  read_O2SRCB[i,:]=a.split()
a=linecache.getline(readin_file,9+74)
read_O2SRCB[9,:]=a.split()
for i in range(0,55):
  O2SRCB[i]=read_O2SRCB[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(O2SRCB,new_height,linewidth=2,linestyle='-',color='black',label='O$_{{2}}$SRCB')

plt.legend(fontsize=7)
plt.xlabel('O$_{{2}}$SRCB')
#plt.xlim([-100,3600])
plt.ylim([-16,11])
plt.ylabel('Ln(P$_{{0}}$/P)')

secax_y = ax.secondary_yaxis('right', functions=(forward_24,inverse_24))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,100,130,160])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_O2SRCB_24.png', dpi=300)

plt.savefig('Figure_O2SRCB_24.eps', format='eps')

O2SRCB=np.zeros(55)
read_O2SRCB=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+102)
  read_O2SRCB[i,:]=a.split()
a=linecache.getline(readin_file,9+102)
read_O2SRCB[9,:]=a.split()
for i in range(0,55):
  O2SRCB[i]=read_O2SRCB[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(O2SRCB,new_height,linewidth=2,linestyle='-',color='black',label='O$_{{2}}$SRCB')

plt.legend(fontsize=7)
plt.xlabel('O$_{{2}}$SRCB')
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


plt.savefig('Figure_O2SRCB_12.png', dpi=300)

plt.savefig('Figure_O2SRCB_12.eps', format='eps')

O2SRCC=np.zeros(55)
read_O2SRCC=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+130)
  read_O2SRCC[i,:]=a.split()
a=linecache.getline(readin_file,9+130)
read_O2SRCC[9,:]=a.split()
for i in range(0,55):
  O2SRCC[i]=read_O2SRCC[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(O2SRCC,new_height,linewidth=2,linestyle='-',color='black',label='O$_{{2}}$SRCC')

plt.legend(fontsize=7)
plt.xlabel('O$_{{2}}$SRCC')
#plt.xlim([-100,3600])
plt.ylim([-16,11])
plt.ylabel('Ln(P$_{{0}}$/P)')

secax_y = ax.secondary_yaxis('right', functions=(forward_24,inverse_24))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,100,130,160])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_O2SRCC_24.png', dpi=300)

plt.savefig('Figure_O2SRCC_24.eps', format='eps')

O2SRCC=np.zeros(55)
read_O2SRCC=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+158)
  read_O2SRCC[i,:]=a.split()
a=linecache.getline(readin_file,9+158)
read_O2SRCC[9,:]=a.split()
for i in range(0,55):
  O2SRCC[i]=read_O2SRCC[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(O2SRCC,new_height,linewidth=2,linestyle='-',color='black',label='O$_{{2}}$SRCC')

plt.legend(fontsize=7)
plt.xlabel('O$_{{2}}$SRCC')
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


plt.savefig('Figure_O2SRCC_12.png', dpi=300)

plt.savefig('Figure_O2SRCC_12.eps', format='eps')

O2SRCT=np.zeros(55)
read_O2SRCT=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+186)
  read_O2SRCT[i,:]=a.split()
a=linecache.getline(readin_file,9+186)
read_O2SRCT[9,:]=a.split()
for i in range(0,55):
  O2SRCT[i]=read_O2SRCT[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(O2SRCT,new_height,linewidth=2,linestyle='-',color='black',label='O$_{{2}}$SRCT')

plt.legend(fontsize=7)
plt.xlabel('O$_{{2}}$SRCT')
#plt.xlim([-100,3600])
plt.ylim([-16,11])
plt.ylabel('Ln(P$_{{0}}$/P)')

secax_y = ax.secondary_yaxis('right', functions=(forward_24,inverse_24))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,100,130,160])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_O2SRCT_24.png', dpi=300)

plt.savefig('Figure_O2SRCT_24.eps', format='eps')

O2SRCT=np.zeros(55)
read_O2SRCT=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+214)
  read_O2SRCT[i,:]=a.split()
a=linecache.getline(readin_file,9+214)
read_O2SRCT[9,:]=a.split()
for i in range(0,55):
  O2SRCT[i]=read_O2SRCT[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(O2SRCT,new_height,linewidth=2,linestyle='-',color='black',label='O$_{{2}}$SRCT')

plt.legend(fontsize=7)
plt.xlabel('O$_{{2}}$SRCT')
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


plt.savefig('Figure_O2SRCT_12.png', dpi=300)

plt.savefig('Figure_O2SRCT_12.eps', format='eps')

O2LOSSA=np.zeros(55)
read_O2LOSSA=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+242)
  read_O2LOSSA[i,:]=a.split()
a=linecache.getline(readin_file,9+242)
read_O2LOSSA[9,:]=a.split()
for i in range(0,55):
  O2LOSSA[i]=read_O2LOSSA[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(O2LOSSA,new_height,linewidth=2,linestyle='-',color='black',label='O$_{{2}}$LOSSA')

plt.legend(fontsize=7)
plt.xlabel('O$_{{2}}$LOSSA')
#plt.xlim([-100,3600])
plt.ylim([-16,11])
plt.ylabel('Ln(P$_{{0}}$/P)')

secax_y = ax.secondary_yaxis('right', functions=(forward_24,inverse_24))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,100,130,160])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_O2LOSSA_24.png', dpi=300)

plt.savefig('Figure_O2LOSSA_24.eps', format='eps')

O2LOSSA=np.zeros(55)
read_O2LOSSA=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+270)
  read_O2LOSSA[i,:]=a.split()
a=linecache.getline(readin_file,9+270)
read_O2LOSSA[9,:]=a.split()
for i in range(0,55):
  O2LOSSA[i]=read_O2LOSSA[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(O2LOSSA,new_height,linewidth=2,linestyle='-',color='black',label='O$_{{2}}$LOSSA')

plt.legend(fontsize=7)
plt.xlabel('O$_{{2}}$LOSSA')
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


plt.savefig('Figure_O2LOSSA_12.png', dpi=300)

plt.savefig('Figure_O2LOSSA_12.eps', format='eps')

O2LOSSB=np.zeros(55)
read_O2LOSSB=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+298)
  read_O2LOSSB[i,:]=a.split()
a=linecache.getline(readin_file,9+298)
read_O2LOSSB[9,:]=a.split()
for i in range(0,55):
  O2LOSSB[i]=read_O2LOSSB[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(O2LOSSB,new_height,linewidth=2,linestyle='-',color='black',label='O$_{{2}}$LOSSB')

plt.legend(fontsize=7)
plt.xlabel('O$_{{2}}$LOSSB')
#plt.xlim([-100,3600])
plt.ylim([-16,11])
plt.ylabel('Ln(P$_{{0}}$/P)')

secax_y = ax.secondary_yaxis('right', functions=(forward_24,inverse_24))


secax_y.set_ylabel(r'Average height (km)')
secax_y.set_yticks([70,100,130,160])
#secax_y.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2 = ax.secondary_yaxis(1.2, functions=(forward2,inverse2))
secax_y2.set_yticks([1,10**-2,10**-4,10**-6,10**-8,10**-10,10**-12])
secax_y2.set_yticklabels(('10$^{{0}}$','10$^{{-2}}$','10$^{{-4}}$','10$^{{-6}}$','10$^{{-8}}$','10$^{{-10}}$','10$^{{-12}}$'))
#secax_y2.set_yscale("log")
#secax_y2.yaxis.set_minor_locator(AutoMinorLocator())
secax_y2.set_ylabel(r'Pressure (Mb)')


plt.savefig('Figure_O2LOSSB_24.png', dpi=300)

plt.savefig('Figure_O2LOSSB_24.eps', format='eps')

O2LOSSB=np.zeros(55)
read_O2LOSSB=np.zeros(shape=(10,6))
for i in range(0,9):
  a=linecache.getline(readin_file,i+326)
  read_O2LOSSB[i,:]=a.split()
a=linecache.getline(readin_file,9+326)
read_O2LOSSB[9,:]=a.split()
for i in range(0,55):
  O2LOSSB[i]=read_O2LOSSB[i//6,i%6]
  
fig=plt.figure(figsize=(10,4))
ax=fig.add_subplot(1,2,1)
plt.plot(O2LOSSB,new_height,linewidth=2,linestyle='-',color='black',label='O$_{{2}}$LOSSB')

plt.legend(fontsize=7)
plt.xlabel('O$_{{2}}$LOSSB')
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


plt.savefig('Figure_O2LOSSB_12.png', dpi=300)

plt.savefig('Figure_O2LOSSB_12.eps', format='eps')
