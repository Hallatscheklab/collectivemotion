import numpy as np
import matplotlib.pyplot as plt

d=5.
layer=8.*d

n=25
nskip=15

# import data - growing outward

t=[]
r=[]
f=open('/global/cscratch1/sd/cschreck/calc_vel_budding_damped_annulus/radius_radial_layer8.0_desync0.4_b4e3_att0.0_seed-1_axial.dat','r')
for line in f.readlines():
    linesplit=line.split()
    t.append(d*float(linesplit[0]))
    r.append(d*float(linesplit[1]))
    #r.append(d*(float(linesplit[1])-0.5))


tave=[np.mean(t[i:i+n]) for i in range(0,len(t)-n,nskip)]
rave=[np.mean(r[i:i+n]) for i in range(0,len(t)-n,nskip)]
curv=[-1./np.mean(r[i:i+n]) for i in range(0,len(t)-n,nskip)]
dldt=[np.polyfit(t[i:i+n],r[i:i+n],1)[0] for i in range(0,len(t)-n,nskip)]

# import data - growing outward

skip=5
t2=[]
r2=[]
f=open('/global/cscratch1/sd/cschreck/calc_vel_budding_damped_annulus/radius_reverse_radial_layer8.0_desync0.4_b4e3_att0.0_seed-1_axial.dat','r')
frame=0
for line in f.readlines():
    linesplit=line.split()
    frame+=1
    if(frame>skip):
        t2.append(d*float(linesplit[0]))
        r2.append(d*(float(linesplit[1])-1.))
        #r2.append(d*(float(linesplit[1])-0.5))

t2ave=[np.mean(t2[i:i+n]) for i in range(0,len(t2)-n,nskip)]
r2ave=[np.mean(r2[i:i+n]) for i in range(0,len(t2)-n,nskip)]
curv2=[1./np.mean(r2[i:i+n]) for i in range(0,len(t2)-n,nskip)]
dldt2=[-np.polyfit(t2[i:i+n],r2[i:i+n],1)[0] for i in range(0,len(t2)-n,nskip)]
 
# plot

plt.plot(curv,dldt,'-or',linewidth=2)
plt.plot(curv2,dldt2,'-om',linewidth=2)
plt.show()

mincurv=1./layer
mincurvfit=mincurv*0.8

print mincurvfit*1000.

curvfit=[]
dldtfit=[]
for i in range(0,len(curv)):
    if(abs(curv[i])<mincurvfit):
        curvfit.append(curv[i])
        dldtfit.append(dldt[i])
for i in range(0,len(curv2)):
    if(abs(curv2[i])<mincurvfit):
        curvfit.append(curv2[i])
        dldtfit.append(dldt2[i])
fitparams=np.polyfit(curvfit,dldtfit,1)
vflat=fitparams[1]
T=fitparams[0]/fitparams[1]

plt.plot(curv,dldt,'-or',linewidth=2)
plt.plot(curv2,dldt2,'-om',linewidth=2)
plt.plot([-mincurv,mincurv],[vflat*(1.-mincurv*T),vflat*(1.+mincurv*T)],'k',linewidth=2.5)
plt.plot([0.,0.],[-5.,20.],'k',linewidth=2.5)
plt.plot([0.0035,0.0035],[-5.,20.],'k',linewidth=2.5)
plt.plot([-0.3,0.3],[vflat,vflat],'k',linewidth=2.5)
plt.xlim([-0.02,0.02])
plt.ylim([-4,20])
plt.show()

print vflat

dldtscale=np.divide(dldt,vflat)
dldt2scale=np.divide(dldt2,vflat)
curvscale=np.multiply(curv,1000.)
curv2scale=np.multiply(curv2,1000.)
#layerscale=np.divide((layer+1.),1000.)
Tscale=T/1000.
mincurvscale=mincurv*1000.

print layer, T

for i in range(0,len(curvscale)):
    cR=min(1.,max(0.,curvscale[i]/(1000.*0.012)))
    cG=min(1.,max(0.,1.-2.*abs(curvscale[i]/(1000.*0.012))))
    cB=min(1.,max(0.,-curvscale[i]/(1000.*0.012)))
    plt.plot(curvscale[i],dldtscale[i],'o',markersize=10,color=(cR,cG,cB))
for i in range(0,len(curv2scale)):
    cR=min(1.,max(0.,curv2scale[i]/(1000.*0.012)))
    cG=min(1.,max(0.,1.-2.*abs(curv2scale[i]/(1000.*0.012))))
    cB=min(1.,max(0.,-curv2scale[i]/(1000.*0.012)))
    plt.plot(curv2scale[i],dldt2scale[i],'o',markersize=10,color=(cR,cG,cB))
plt.plot([-mincurvscale,mincurvscale],[1.-mincurvscale*Tscale,1.+mincurvscale*Tscale],'k',linewidth=2)
plt.xlim([-20.,20.])
plt.ylim([0.5,1.5])
plt.plot([-25.,25.],[1.,1.],'k',linewidth=2)
plt.plot([0.,0.],[0.,2.],'k',linewidth=2)
plt.xlabel('Curvature (mm$^{-1}$)',fontsize=16)
plt.ylabel('Normalized Velocity',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig('curv_vs_vel_color_lambda8.pdf')
plt.tight_layout()
plt.show()
