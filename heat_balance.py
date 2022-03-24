import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import linecache

readin_file='vtgcm.HTBAL.VEN1.dat'

height=np.zeros(66)
matrix_2=np.zeros(shape=(11,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+7)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  height[i]=matrix_2[i//6,i%6]
  
QEUV=np.zeros(66)
matrix_2=np.zeros(shape=(11,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+18)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  QEUV[i]=matrix_2[i//6,i%6]

QNIR=np.zeros(66)
matrix_2=np.zeros(shape=(11,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+46)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  QNIR[i]=matrix_2[i//6,i%6]

QIR=np.zeros(66)
matrix_2=np.zeros(shape=(11,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+74)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  QIR[i]=matrix_2[i//6,i%6]

COND=np.zeros(66)
matrix_2=np.zeros(shape=(11,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+102)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  COND[i]=matrix_2[i//6,i%6]

TOT=np.zeros(66)
matrix_2=np.zeros(shape=(11,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+130)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  TOT[i]=matrix_2[i//6,i%6]
  
  
plt.figure()
plt.plot(QEUV,height,linewidth=2,color='blue',label='EUV-UV heating')
plt.plot(QNIR,height,linewidth=2,color='green',label='Near IR heating')
plt.plot(QIR,height,linewidth=2,color='red',label='CO$_{{2}}$ 15 ${\mu}$m cooling')
plt.plot(COND,height,linewidth=2,color='orange',label='Conduction')
plt.plot(TOT,height,linewidth=2,color='purple',label='Total advection')
plt.legend(fontsize=8)
plt.xlabel('[K/Day]')
plt.ylabel('Height [km]')
plt.xlim([-3500,3500])
plt.ylim([70,200])

plt.savefig('Figure_heat_a.png', dpi=300)

plt.savefig('Figure_heat_a.eps', format='eps')

#####################################
readin_file='vtgcm.HTBAL.VEN2.dat'

height=np.zeros(66)
matrix_2=np.zeros(shape=(11,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+7)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  height[i]=matrix_2[i//6,i%6]
  
QEUV=np.zeros(66)
matrix_2=np.zeros(shape=(11,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+18)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  QEUV[i]=matrix_2[i//6,i%6]

QNIR=np.zeros(66)
matrix_2=np.zeros(shape=(11,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+46)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  QNIR[i]=matrix_2[i//6,i%6]

QIR=np.zeros(66)
matrix_2=np.zeros(shape=(11,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+74)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  QIR[i]=matrix_2[i//6,i%6]

COND=np.zeros(66)
matrix_2=np.zeros(shape=(11,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+102)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  COND[i]=matrix_2[i//6,i%6]

TOT=np.zeros(66)
matrix_2=np.zeros(shape=(11,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+130)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  TOT[i]=matrix_2[i//6,i%6]
  
  
plt.figure()
plt.plot(QEUV,height,linewidth=2,color='blue',label='EUV-UV heating')
plt.plot(QNIR,height,linewidth=2,color='green',label='Near IR heating')
plt.plot(QIR,height,linewidth=2,color='red',label='CO$_{{2}}$ 15 ${\mu}$m cooling')
plt.plot(COND,height,linewidth=2,color='orange',label='Conduction')
plt.plot(TOT,height,linewidth=2,color='purple',label='Total advection')
plt.legend(fontsize=8)
plt.xlabel('[K/Day]')
plt.ylabel('Height [km]')
plt.xlim([-3500,3500])
plt.ylim([70,200])

plt.savefig('Figure_heat_b.png', dpi=300)

plt.savefig('Figure_heat_b.eps', format='eps')

##########################################
readin_file='vtgcm.HTBAL.VEN3.dat'

height=np.zeros(66)
matrix_2=np.zeros(shape=(11,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+7)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  height[i]=matrix_2[i//6,i%6]
  
QEUV=np.zeros(66)
matrix_2=np.zeros(shape=(11,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+18)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  QEUV[i]=matrix_2[i//6,i%6]

QNIR=np.zeros(66)
matrix_2=np.zeros(shape=(11,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+46)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  QNIR[i]=matrix_2[i//6,i%6]

QIR=np.zeros(66)
matrix_2=np.zeros(shape=(11,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+74)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  QIR[i]=matrix_2[i//6,i%6]

COND=np.zeros(66)
matrix_2=np.zeros(shape=(11,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+102)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  COND[i]=matrix_2[i//6,i%6]

TOT=np.zeros(66)
matrix_2=np.zeros(shape=(11,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+130)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  TOT[i]=matrix_2[i//6,i%6]
  
  
plt.figure()
plt.plot(QEUV,height,linewidth=2,color='blue',label='EUV-UV heating')
plt.plot(QNIR,height,linewidth=2,color='green',label='Near IR heating')
plt.plot(QIR,height,linewidth=2,color='red',label='CO$_{{2}}$ 15 ${\mu}$m cooling')
plt.plot(COND,height,linewidth=2,color='orange',label='Conduction')
plt.plot(TOT,height,linewidth=2,color='purple',label='Total advection')
plt.legend(fontsize=8)
plt.xlabel('[K/Day]')
plt.ylabel('Height [km]')
plt.xlim([-3500,3500])
plt.ylim([70,200])

plt.savefig('Figure_heat_c.png', dpi=300)

plt.savefig('Figure_heat_c.eps', format='eps')

###################################
readin_file='vtgcm.HTBAL.VEN4.dat'

height=np.zeros(66)
matrix_2=np.zeros(shape=(11,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+7)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  height[i]=matrix_2[i//6,i%6]
  
QEUV=np.zeros(66)
matrix_2=np.zeros(shape=(11,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+18)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  QEUV[i]=matrix_2[i//6,i%6]

QNIR=np.zeros(66)
matrix_2=np.zeros(shape=(11,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+46)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  QNIR[i]=matrix_2[i//6,i%6]

QIR=np.zeros(66)
matrix_2=np.zeros(shape=(11,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+74)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  QIR[i]=matrix_2[i//6,i%6]

COND=np.zeros(66)
matrix_2=np.zeros(shape=(11,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+102)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  COND[i]=matrix_2[i//6,i%6]

TOT=np.zeros(66)
matrix_2=np.zeros(shape=(11,6))
for i in range(0,11):
  a=linecache.getline(readin_file,i+130)
  matrix_2[i,:]=a.split()
for i in range(0,66):
  TOT[i]=matrix_2[i//6,i%6]
  
  
plt.figure()
plt.plot(QEUV,height,linewidth=2,color='blue',label='EUV-UV heating')
plt.plot(QNIR,height,linewidth=2,color='green',label='Near IR heating')
plt.plot(QIR,height,linewidth=2,color='red',label='CO$_{{2}}$ 15 ${\mu}$m cooling')
plt.plot(COND,height,linewidth=2,color='orange',label='Conduction')
plt.plot(TOT,height,linewidth=2,color='purple',label='Total advection')
plt.legend(fontsize=8)
plt.xlabel('[K/Day]')
plt.ylabel('Height [km]')
plt.xlim([-13000,13000])
plt.ylim([70,200])

plt.savefig('Figure_heat_d.png', dpi=300)

plt.savefig('Figure_heat_d.eps', format='eps')



  
  
