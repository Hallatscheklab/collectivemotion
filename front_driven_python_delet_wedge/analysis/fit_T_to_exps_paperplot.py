import math
import numpy as np
import matplotlib.pyplot as plt
 
d=5.
off=230.
binsize=8.

minw=20.
maxw=330.

# import data - FDM

lfA=[]
wfA=[]
f=open('../data/bound_L600.0_D2.5_s-0.06_diff0.0_seed1.dat','r')
for line in f.readlines():
    linesplit=line.split()
    lfA.append((d*float(linesplit[3])+d*float(linesplit[1]))/2.)
    wfA.append(d*float(linesplit[2])-d*float(linesplit[0]))

lfB=[]
wfB=[]
f=open('../data/bound_L600.0_D3.0_s-0.06_diff0.0_seed1.dat','r')
for line in f.readlines():
    linesplit=line.split()
    lfB.append((d*float(linesplit[3])+d*float(linesplit[1]))/2.)
    wfB.append(d*float(linesplit[2])-d*float(linesplit[0]))

lfC=[]
wfC=[]
f=open('../data/bound_L600.0_D4.0_s-0.06_diff0.0_seed1.dat','r')
for line in f.readlines():
    linesplit=line.split()
    lfC.append((d*float(linesplit[3])+d*float(linesplit[1]))/2.)
    wfC.append(d*float(linesplit[2])-d*float(linesplit[0]))

lfD=[]
wfD=[]
f=open('../data/bound_L600.0_D5.0_s-0.06_diff0.0_seed1.dat','r')
for line in f.readlines():
    linesplit=line.split()
    lfD.append((d*float(linesplit[3])+d*float(linesplit[1]))/2.)
    wfD.append(d*float(linesplit[2])-d*float(linesplit[0]))

lfE=[]
wfE=[]
f=open('../data/bound_L600.0_D6.0_s-0.06_diff0.0_seed1.dat','r')
for line in f.readlines():
    linesplit=line.split()
    lfE.append((d*float(linesplit[3])+d*float(linesplit[1]))/2.)
    wfE.append(d*float(linesplit[2])-d*float(linesplit[0]))

lfF=[]
wfF=[]
f=open('../data/bound_L600.0_D7.0_s-0.06_diff0.0_seed1.dat','r')
for line in f.readlines():
    linesplit=line.split()
    lfF.append((d*float(linesplit[3])+d*float(linesplit[1]))/2.)
    wfF.append(d*float(linesplit[2])-d*float(linesplit[0]))

lfG=[]
wfG=[]
f=open('../data/bound_L600.0_D8.0_s-0.06_diff0.0_seed1.dat','r')
for line in f.readlines():
    linesplit=line.split()
    lfG.append((d*float(linesplit[3])+d*float(linesplit[1]))/2.)
    wfG.append(d*float(linesplit[2])-d*float(linesplit[0]))

lfH=[]
wfH=[]
f=open('../data/bound_L600.0_D9.0_s-0.06_diff0.0_seed1.dat','r')
for line in f.readlines():
    linesplit=line.split()
    lfH.append((d*float(linesplit[3])+d*float(linesplit[1]))/2.)
    wfH.append(d*float(linesplit[2])-d*float(linesplit[0]))

# import data - exps

le=[]
we1=[]
we2=[]
we3=[]
we4=[]
we5=[]
we6=[]
we7=[]
we8=[]
we9=[]
f=open('../data/width_vs_length_exp.dat','r')
for line in f.readlines():
    linesplit=line.split()
    le.append(float(linesplit[0]))
    we1.append(float(linesplit[1]))
    we2.append(float(linesplit[2]))
    we3.append(float(linesplit[3]))
    we4.append(float(linesplit[4]))
    we5.append(float(linesplit[5]))
    we6.append(float(linesplit[6]))
    we7.append(float(linesplit[7]))
    we8.append(float(linesplit[8]))
    we9.append(float(linesplit[9]))

# align data - FDM

lfAoff=np.zeros(len(lfA))
for i in range(0,len(lfA)):
    if(wfA[i]>off):
        ifA=i
for i in range(0,len(lfA)):
    lfAoff[i]=lfA[i]-lfA[ifA]

lfBoff=np.zeros(len(lfB))
for i in range(0,len(lfB)):
    if(wfB[i]>off):
        ifB=i
for i in range(0,len(lfB)):
    lfBoff[i]=lfB[i]-lfB[ifB]

lfCoff=np.zeros(len(lfC))
for i in range(0,len(lfC)):
    if(wfC[i]>off):
        ifC=i
for i in range(0,len(lfC)):
    lfCoff[i]=lfC[i]-lfC[ifC]

lfDoff=np.zeros(len(lfD))
for i in range(0,len(lfD)):
    if(wfD[i]>off):
        ifD=i
for i in range(0,len(lfD)):
    lfDoff[i]=lfD[i]-lfD[ifD]

lfEoff=np.zeros(len(lfE))
for i in range(0,len(lfE)):
    if(wfE[i]>off):
        ifE=i
for i in range(0,len(lfE)):
    lfEoff[i]=lfE[i]-lfE[ifE]

lfFoff=np.zeros(len(lfF))
for i in range(0,len(lfF)):
    if(wfF[i]>off):
        ifF=i
for i in range(0,len(lfF)):
    lfFoff[i]=lfF[i]-lfF[ifF]

lfGoff=np.zeros(len(lfG))
for i in range(0,len(lfG)):
    if(wfG[i]>off):
        ifG=i
for i in range(0,len(lfG)):
    lfGoff[i]=lfG[i]-lfG[ifG]

lfHoff=np.zeros(len(lfH))
for i in range(0,len(lfH)):
    if(wfH[i]>off):
        ifH=i
for i in range(0,len(lfH)):
    lfHoff[i]=lfH[i]-lfH[ifH]

# align data - exps

le1off=np.zeros(len(le))
for i in range(0,len(le)):
    if(we1[i]>off):
        ie=i
for i in range(0,len(le)):
    le1off[i]=le[i]-le[ie]

le2off=np.zeros(len(le))
for i in range(0,len(le)):
    if(we2[i]>off):
        ie=i
for i in range(0,len(le)):
    le2off[i]=le[i]-le[ie]

le3off=np.zeros(len(le))
for i in range(0,len(le)):
    if(we3[i]>off):
        ie=i
for i in range(0,len(le)):
    le3off[i]=le[i]-le[ie]

le4off=np.zeros(len(le))
for i in range(0,len(le)):
    if(we4[i]>off):
        ie=i
for i in range(0,len(le)):
    le4off[i]=le[i]-le[ie]

le5off=np.zeros(len(le))
for i in range(0,len(le)):
    if(we5[i]>off):
        ie=i
for i in range(0,len(le)):
    le5off[i]=le[i]-le[ie]

le6off=np.zeros(len(le))
for i in range(0,len(le)):
    if(we6[i]>off):
        ie=i
for i in range(0,len(le)):
    le6off[i]=le[i]-le[ie]

le7off=np.zeros(len(le))
for i in range(0,len(le)):
    if(we7[i]>off):
        ie=i
for i in range(0,len(le)):
    le7off[i]=le[i]-le[ie]

le8off=np.zeros(len(le))
for i in range(0,len(le)):
    if(we8[i]>off):
        ie=i
for i in range(0,len(le)):
    le8off[i]=le[i]-le[ie]

le9off=np.zeros(len(le))
for i in range(0,len(le)):
    if(we9[i]>off):
        ie=i
for i in range(0,len(le)):
    le9off[i]=le[i]-le[ie]

mindis=-300.
maxdis=1500.
numbins=int((maxdis-mindis)/binsize)+1
minbin=int(mindis/binsize)

lemean=np.zeros(numbins)
wemean=np.zeros(numbins)
wesq=np.zeros(numbins)
count=np.zeros(numbins)
for i in range(0,len(le)):
    bin=int(le1off[i]/binsize)-minbin
    if(bin>=0 and bin<numbins):
        wemean[bin]+=we1[i]
        wesq[bin]+=we1[i]*we1[i]
        count[bin]+=1.
    bin=int(le2off[i]/binsize)-minbin
    if(bin>=0 and bin<numbins):
        wemean[bin]+=we2[i]
        wesq[bin]+=we2[i]*we2[i]
        count[bin]+=1.
    bin=int(le3off[i]/binsize)-minbin
    if(bin>=0 and bin<numbins):
        wemean[bin]+=we3[i]
        wesq[bin]+=we3[i]*we3[i]
        count[bin]+=1.
    bin=int(le4off[i]/binsize)-minbin
    if(bin>=0 and bin<numbins):
        wemean[bin]+=we4[i]
        wesq[bin]+=we4[i]*we4[i]
        count[bin]+=1.
    bin=int(le5off[i]/binsize)-minbin
    if(bin>=0 and bin<numbins):
        wemean[bin]+=we5[i]
        wesq[bin]+=we5[i]*we5[i]
        count[bin]+=1.
    bin=int(le6off[i]/binsize)-minbin
    if(bin>=0 and bin<numbins):
        wemean[bin]+=we6[i]
        wesq[bin]+=we6[i]*we6[i]
        count[bin]+=1.
    bin=int(le7off[i]/binsize)-minbin
    if(bin>=0 and bin<numbins):
        wemean[bin]+=we7[i]
        wesq[bin]+=we7[i]*we7[i]
        count[bin]+=1.
    bin=int(le8off[i]/binsize)-minbin
    if(bin>=0 and bin<numbins):
        wemean[bin]+=we8[i]
        wesq[bin]+=we8[i]*we8[i]
        count[bin]+=1.
    bin=int(le9off[i]/binsize)-minbin
    if(bin>=0 and bin<numbins):
        wemean[bin]+=we9[i]
        wesq[bin]+=we9[i]*we9[i]
        count[bin]+=1.
for bin in range(0,numbins):
    if(count[bin]>0):
        lemean[bin]=float(bin+minbin)*binsize
        wemean[bin]=wemean[bin]/count[bin]
        wesq[bin]=wesq[bin]/count[bin]
westd=np.zeros(numbins)
wem1std=np.zeros(numbins)
wep1std=np.zeros(numbins)
for bin in range(0,numbins):
    westd[bin]=math.sqrt(wesq[bin]-wemean[bin]**2)
    wem1std[bin]=wemean[bin]-westd[bin]
    wep1std[bin]=wemean[bin]+westd[bin]

lfAbin=np.zeros(numbins)
wfAbin=np.zeros(numbins)
count=np.zeros(numbins)
for i in range(0,len(lfA)):
    bin=int(lfAoff[i]/binsize)-minbin
    if(bin>=0 and bin<numbins):
        wfAbin[bin]+=wfA[i]
        count[bin]+=1.
for bin in range(0,numbins):
    if(count[bin]>0):
        lfAbin[bin]=float(bin+minbin)*binsize
        wfAbin[bin]=wfAbin[bin]/count[bin]
lfBbin=np.zeros(numbins)
wfBbin=np.zeros(numbins)
count=np.zeros(numbins)
for i in range(0,len(lfB)):
    bin=int(lfBoff[i]/binsize)-minbin
    if(bin>=0 and bin<numbins):
        wfBbin[bin]+=wfB[i]
        count[bin]+=1.
for bin in range(0,numbins):
    if(count[bin]>0):
        lfBbin[bin]=float(bin+minbin)*binsize
        wfBbin[bin]=wfBbin[bin]/count[bin]
lfCbin=np.zeros(numbins)
wfCbin=np.zeros(numbins)
count=np.zeros(numbins)
for i in range(0,len(lfC)):
    bin=int(lfCoff[i]/binsize)-minbin
    if(bin>=0 and bin<numbins):
        wfCbin[bin]+=wfC[i]
        count[bin]+=1.
for bin in range(0,numbins):
    if(count[bin]>0):
        lfCbin[bin]=float(bin+minbin)*binsize
        wfCbin[bin]=wfCbin[bin]/count[bin]
lfDbin=np.zeros(numbins)
wfDbin=np.zeros(numbins)
count=np.zeros(numbins)
for i in range(0,len(lfD)):
    bin=int(lfDoff[i]/binsize)-minbin
    if(bin>=0 and bin<numbins):
        wfDbin[bin]+=wfD[i]
        count[bin]+=1.
for bin in range(0,numbins):
    if(count[bin]>0):
        lfDbin[bin]=float(bin+minbin)*binsize
        wfDbin[bin]=wfDbin[bin]/count[bin]
lfEbin=np.zeros(numbins)
wfEbin=np.zeros(numbins)
count=np.zeros(numbins)
for i in range(0,len(lfE)):
    bin=int(lfEoff[i]/binsize)-minbin
    if(bin>=0 and bin<numbins):
        wfEbin[bin]+=wfE[i]
        count[bin]+=1.
for bin in range(0,numbins):
    if(count[bin]>0):
        lfEbin[bin]=float(bin+minbin)*binsize
        wfEbin[bin]=wfEbin[bin]/count[bin]
lfFbin=np.zeros(numbins)
wfFbin=np.zeros(numbins)
count=np.zeros(numbins)
for i in range(0,len(lfF)):
    bin=int(lfFoff[i]/binsize)-minbin
    if(bin>=0 and bin<numbins):
        wfFbin[bin]+=wfF[i]
        count[bin]+=1.
for bin in range(0,numbins):
    if(count[bin]>0):
        lfFbin[bin]=float(bin+minbin)*binsize
        wfFbin[bin]=wfFbin[bin]/count[bin]
lfGbin=np.zeros(numbins)
wfGbin=np.zeros(numbins)
count=np.zeros(numbins)
for i in range(0,len(lfG)):
    bin=int(lfGoff[i]/binsize)-minbin
    if(bin>=0 and bin<numbins):
        wfGbin[bin]+=wfG[i]
        count[bin]+=1.
for bin in range(0,numbins):
    if(count[bin]>0):
        lfGbin[bin]=float(bin+minbin)*binsize
        wfGbin[bin]=wfGbin[bin]/count[bin]
lfHbin=np.zeros(numbins)
wfHbin=np.zeros(numbins)
count=np.zeros(numbins)
for i in range(0,len(lfH)):
    bin=int(lfHoff[i]/binsize)-minbin
    if(bin>=0 and bin<numbins):
        wfHbin[bin]+=wfH[i]
        count[bin]+=1.
for bin in range(0,numbins):
    if(count[bin]>0):
        lfHbin[bin]=float(bin+minbin)*binsize
        wfHbin[bin]=wfHbin[bin]/count[bin]

offETA=off/(2.*math.sqrt(0.06*(2.-0.06))/(1.-0.06))

#plt.plot(le1off,we1,'gray')
#plt.plot(le2off,we2,'gray')
#plt.plot(le3off,we3,'gray')
#plt.plot(le4off,we4,'gray')
#plt.plot(le5off,we5,'gray')
#plt.plot(le6off,we6,'gray')
#plt.plot(le7off,we7,'gray')
#plt.plot(le8off,we8,'gray')
#plt.plot(le9off,we9,'gray')
#plt.plot(lemean,wemean,'k',linewidth=2.5)
#plt.plot(lfDoff,wfD,'r',linewidth=2.5)
#plt.plot([-offETA,offETA],[2.*off,0],'--k',linewidth=2.5)
#plt.xlim([-offETA,800.])
#plt.ylim([0,2.*off])
#plt.xlabel('length ($\mu$m)')
#plt.ylabel('width ($\mu$m)')
#plt.savefig('width_exp_fdm_lines.pdf')
#plt.show()

plt.fill_between(lemean,wem1std,wep1std,color='gray')
plt.plot(lemean,wemean,'k',linewidth=2.5,label='Experiment')
plt.plot(lfBoff,wfB,'--b',linewidth=2,label='Simulation, T=15$\mu$m')
plt.plot(lfDoff,wfD,'--r',linewidth=2,label='Simulation, T=25$\mu$m') # best fit
plt.plot(lfFoff,wfF,'--g',linewidth=2,label='Simulation, T=35$\mu$m')
plt.plot(lfHoff,wfH,'--m',linewidth=2,label='Simulation, T=45$\mu$m')
plt.plot([-offETA,offETA],[2.*off,0],':k',linewidth=2.5,label='Without surface tension')
plt.xlim([-150.,500.])
plt.ylim([0,350.])
plt.xlabel('Distance ($\mu$m)',fontsize=16)
plt.ylabel('Width ($\mu$m)',fontsize=16)
plt.tick_params(labelsize=16)
plt.legend(loc=1,fontsize=12)
plt.tight_layout()
plt.savefig('width_exp_fdm_shaded_5LT.pdf')
plt.show()

delA=0.
delB=0.
delC=0.
delD=0.
delE=0.
delF=0.
delG=0.
delH=0.
countA=0.
countB=0.
countC=0.
countD=0.
countE=0.
countF=0.
countG=0.
countH=0.
for bin in range(0,numbins):
    if(wemean[bin]>minw and wemean[bin]<maxw):
        if(wfAbin[bin]>minw and wfAbin[bin]<maxw):
            delA+=(wemean[bin]-wfAbin[bin])**2
            countA+=1.
        if(wfBbin[bin]>minw and wfBbin[bin]<maxw):
            delB+=(wemean[bin]-wfBbin[bin])**2
            countB+=1.
        if(wfCbin[bin]>minw and wfCbin[bin]<maxw):
            delC+=(wemean[bin]-wfCbin[bin])**2
            countC+=1.
        if(wfDbin[bin]>minw and wfDbin[bin]<maxw):
            delD+=(wemean[bin]-wfDbin[bin])**2
            countD+=1.
        if(wfEbin[bin]>minw and wfEbin[bin]<maxw):
            delE+=(wemean[bin]-wfEbin[bin])**2
            countE+=1.
        if(wfFbin[bin]>minw and wfFbin[bin]<maxw):
            delF+=(wemean[bin]-wfFbin[bin])**2
            countF+=1.
        if(wfGbin[bin]>minw and wfGbin[bin]<maxw):
            delG+=(wemean[bin]-wfGbin[bin])**2
            countG+=1.
        if(wfHbin[bin]>minw and wfHbin[bin]<maxw):
            delH+=(wemean[bin]-wfHbin[bin])**2
            countH+=1.
delA=delA/countA
delB=delB/countB
delC=delC/countC
delD=delD/countD
delE=delE/countE
delF=delF/countF
delG=delG/countG
delH=delH/countH

Tvals=[2.5*d,3.*d,4.*d,5.*d,6.*d,7.*d,8.*d,9.*d]
delvals=[delA,delB,delC,delD,delE,delF,delG,delH]
plt.plot(Tvals,delvals,'-ok')
plt.xlabel('T ($\mu$m)',fontsize=16)
plt.ylabel('Exp-Sim Deviation',fontsize=16)
plt.xlim([0,10.*d])
plt.ylim([0,600.])
plt.tick_params(labelsize=16)
plt.tight_layout()
plt.savefig('deviation_vs_T_paper.pdf')
plt.show()

f=open('../data/width_vs_length_mean_exp.dat','w')
for i in range(1,len(lemean)):
    print>>f,lemean[i],wemean[i],westd[i]

f=open('../data/width_vs_length_fdm_T2.5.dat','w')
for i in range(0,len(lfAoff)):
    print>>f,lfAoff[i],wfA[i]
f=open('../data/width_vs_length_fdm_T3.dat','w')
for i in range(0,len(lfBoff)):
    print>>f,lfBoff[i],wfB[i]
f=open('../data/width_vs_length_fdm_T4.dat','w')
for i in range(0,len(lfCoff)):
    print>>f,lfCoff[i],wfC[i]
f=open('../data/width_vs_length_fdm_T5.dat','w')
for i in range(0,len(lfDoff)):
    print>>f,lfDoff[i],wfD[i]
f=open('../data/width_vs_length_fdm_T6.dat','w')
for i in range(0,len(lfEoff)):
    print>>f,lfEoff[i],wfE[i]
f=open('../data/width_vs_length_fdm_T7.dat','w')
for i in range(0,len(lfFoff)):
    print>>f,lfFoff[i],wfF[i]
f=open('../data/width_vs_length_fdm_T8.dat','w')
for i in range(0,len(lfGoff)):
    print>>f,lfGoff[i],wfG[i]
f=open('../data/width_vs_length_fdm_T9.dat','w')
for i in range(0,len(lfHoff)):
    print>>f,lfHoff[i],wfH[i]

n=2
lfDsmooth1=[np.mean(lfDoff[i:i+n]) for i in range(0,len(lfDoff)-n)]
wfDsmooth1=[np.mean(wfD[i:i+n]) for i in range(0,len(wfD)-n)]
dwdlD1=[-np.polyfit(lfDoff[i:i+n],wfD[i:i+n],1)[0] for i in range(0,len(wfD)-n)]

n=40
lfDsmooth2=[np.mean(lfDoff[i:i+n]) for i in range(0,len(lfDoff)-n)]
wfDsmooth2=[np.mean(wfD[i:i+n]) for i in range(0,len(wfD)-n)]
dwdlD2=[-np.polyfit(lfDoff[i:i+n],wfD[i:i+n],1)[0] for i in range(0,len(wfD)-n)]

f=open('../data/width_deriv_fdm_T5.dat','w')
for i in range(0,len(wfDsmooth2)):
    print>>f,wfDsmooth2[i],dwdlD2[i]
