import time
import math
import numpy as np
import matplotlib.pyplot as plt

t1=time.time()

# import data
f=open('/global/cscratch1/sd/cschreck/calc_vel_budding_damped_annulus/prod_radial_layer8.0_desync0.4_b4e3_att0.0_seed-1_axial.dat','r')

# array sizes
nframes=100
maxcells=400000

# for calc
layer=8.

# for plottings
frameplot=96

# initialize arrays
n=np.zeros(nframes,dtype=np.int)
x=np.zeros((nframes,maxcells))
y=np.zeros((nframes,maxcells))
d=np.zeros((nframes,maxcells))
depth=np.zeros((nframes,maxcells))
type=np.zeros((nframes,maxcells))
radius=np.zeros(nframes)

# import data
frame=0
ncell=0
for line in f.readlines():
    linesplit=line.split()
    if len(linesplit)==1:
        frame=frame+1
        ncell=0
    elif(frame<nframes):
        ncell=ncell+1      
        n[frame-1]=ncell
        x[frame-1][ncell-1]=float(linesplit[0])
        y[frame-1][ncell-1]=float(linesplit[1])
        d[frame-1][ncell-1]=float(linesplit[2])
        depth[frame-1][ncell-1]=float(linesplit[3])
        type[frame-1][ncell-1]=float(linesplit[4])

for frame in range(0,nframes):
    if(n[frame]):
        count=0.
        for i in range(0,n[frame]):    
            if(depth[frame][i]<2.):
                radius[frame]+=x[frame][i]**2+y[frame][i]**2
                count+=1.
        radius[frame]=math.sqrt(radius[frame]/count)
        print frame, radius[frame]

# calculate boundaries
binsize=3.2
frontbins=int(2.*np.pi*radius[frameplot]/binsize)
binanglesize=2.*np.pi/np.float(frontbins)
xfront=np.zeros((nframes,frontbins+1))
yfront=np.zeros((nframes,frontbins+1))
xlayer=np.zeros((nframes,frontbins+1))+1e16
ylayer=np.zeros((nframes,frontbins+1))+1e16
for frame in range(0,nframes):
    # calc frontline
    for i in range(0,n[frame]): 
        if(depth[frame][i]<0.5):
            xi=x[frame][i]
            yi=y[frame][i]
            angle=np.arctan(yi/xi)
            if(xi<0. and yi>0.):
                angle+=np.pi
            elif(xi<0. and yi<0.):
                angle-=np.pi
            angle+=np.pi
            bin=int(angle/binanglesize)
            rsq=xi**2+yi**2
            xf=xfront[frame][bin]
            yf=yfront[frame][bin]
            rsqfront=xf**2+yf**2
            if(bin>=0 and bin<frontbins):
                if(rsq>rsqfront):
                    ri=math.sqrt(rsq)#+d[frame][i]/2.             
                    anglef=(np.float(bin)+0.5)*binanglesize-np.pi
                    xfront[frame][bin]=np.cos(anglef)*ri
                    yfront[frame][bin]=np.sin(anglef)*ri
        elif(depth[frame][i]<layer):
            xi=x[frame][i]
            yi=y[frame][i]
            angle=np.arctan(yi/xi)
            if(xi<0. and yi>0.):
                angle+=np.pi
            elif(xi<0. and yi<0.):
                angle-=np.pi
            angle+=np.pi
            bin=int(angle/binanglesize)
            rsq=xi**2+yi**2
            xf=xlayer[frame][bin]
            yf=ylayer[frame][bin]
            rsqlayer=xf**2+yf**2
            if(bin>=0 and bin<frontbins):
                if(rsq<rsqlayer):
                    ri=math.sqrt(rsq)#-d[frame][i]/2.         
                    anglef=(np.float(bin)+0.5)*binanglesize-np.pi
                    xlayer[frame][bin]=np.cos(anglef)*ri
                    ylayer[frame][bin]=np.sin(anglef)*ri
    xfront[frame,-1]=xfront[frame,0]
    xlayer[frame,-1]=xlayer[frame,0]
    yfront[frame,-1]=yfront[frame,0]
    ylayer[frame,-1]=ylayer[frame,0]

# plot cells
frame=frameplot
for cell in range (0,n[frame]):
    xx=x[frame][cell]
    yy=y[frame][cell]
    dd=d[frame][cell]/2.
    dep=depth[frame][cell]
    if(dep<layer and abs(xx)<25. and abs(yy-365.)<25.):
        c1=plt.Circle((xx,yy),dd,facecolor='gray',edgecolor='none')
        plt.gcf().gca().add_artist(c1)
    else:
        c1=plt.Circle((xx,yy),dd,facecolor='lightgray',edgecolor='none')
        plt.gcf().gca().add_artist(c1)

d=5.
cR0=min(1.,max(0.,-1./(0.012*d*radius[frameplot])))
cR1=min(1.,max(0.,-1./(0.012*d*radius[frameplot+1])))
#cR2=min(1.,max(0.,-1./(0.012*d*radius[frameplot+2])))
#cR3=min(1.,max(0.,-1./(0.012*d*radius[frameplot+3])))
cG0=min(1.,max(0.,1.-2.*abs(1./(0.012*d*radius[frameplot]))))
cG1=min(1.,max(0.,1.-2.*abs(1./(0.012*d*radius[frameplot+1]))))
#cG2=min(1.,max(0.,1.-2.*abs(1./(0.012*d*radius[frameplot+2]))))
#cG3=min(1.,max(0.,1.-2.*abs(1./(0.012*d*radius[frameplot+3]))))
cB0=min(1.,max(0.,1./(0.012*d*radius[frameplot])))
cB1=min(1.,max(0.,1./(0.012*d*radius[frameplot+1])))
#cB2=min(1.,max(0.,1./(0.012*d*radius[frameplot+2])))
#cB3=min(1.,max(0.,1./(0.012*d*radius[frameplot+3])))

print 1000./(5.*radius[frameplot])
print 1000./(5.*radius[frameplot+1])

#plt.plot(xfront[frameplot+3,:],yfront[frameplot+3,:],'-',linewidth=2,markersize=10,color=(cR3,cG3,cB3))
#plt.plot(xfront[frameplot+2,:],yfront[frameplot+2,:],'-',linewidth=2,markersize=10,color=(cR2,cG2,cB2))
plt.plot(xfront[frameplot+1,:],yfront[frameplot+1,:],'-',linewidth=2,markersize=10,color=(cR1,cG1,cB1))
plt.plot(xfront[frameplot,:],yfront[frameplot,:],'-',linewidth=2,markersize=10,color=(cR0,cG0,cB0))
#plt.plot(xlayer[frameplot,:],ylayer[frameplot,:],'-',linewidth=2)
plt.plot([-10.5,-10.5,10.5,10.5],[400.,0.,0.,400.],'w',linewidth=6)
plt.axes().set_aspect('equal')
plt.xlim([-25.,25.])
plt.ylim([360.,410.])
plt.xticks([])
plt.yticks([])
plt.axis('off')
plt.savefig('system_forward_slice_v3_flat_lambda8_lightoutside.pdf')
plt.savefig('system_forward_slice_v3_flat_lambda8_lightoutside.png')
plt.show() 

# print running time
t2=time.time()
print int(t2-t1), 'seconds'
