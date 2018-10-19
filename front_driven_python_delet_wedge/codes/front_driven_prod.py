import time
import math
import random
import numpy as np

t1=time.time()

################################################################
################### input from screen / bash ###################

# input lengthscales
T=input()        #     T = surface tension
L=input()        #     L = length of box
w0=input()       #    w0 = initial width of mutant clone
diam=input()     #  diam = diamter of cell (unit of length)
ndis=input()     #  ndis = number of lattice point per diam

# input dynamic params
diff=input()     # diffusion coeff of wild-type/mutant boundaries
steps=input()    # number of growth steps in simulation
skipb=input()    # number of steps between writing boundary data
skiph=input()    # number of steps between writing front data
dh0=input()      # distance flat front translates per growth step

# input evolution params
s=input()        # Selection coefficient: k_mut = k_wt*(1+s)

# input random seed
randseed=input() # random seed for boundary fluctuations
trials=input()   # number of indpendent runs

# input output file names
fileA=input()    # file to store front position
fileB=input()    # file to store bounary position

################################################################
###################### calc run quantities #####################

# initial random number generator
random.seed(randseed)

# calc run quantities
dx=diam/ndis
nx=int(L/dx)+1
xdiff=math.sqrt(3.*diff*dh0)

# consts for deriv calcs
dx2=2.*dx
dxsq=dx**2
dh=np.zeros(nx)

# for slanted wedge
dhdxinit=math.sqrt(-s*(2.+s))/(1.+s)
dhend=dx*dhdxinit

# array for saving data
fbA=open(fileA,'w')
fbB=open(fileB,'w')

################################################################
#################### loop over growth steps ####################
for trial in range(0,trials):
    print "starting trial ", trial
    
    # initialize front
    x=np.linspace(-L/2,L/2,nx)
    b=[-w0/2,w0/2]

    h=np.zeros(nx)
    for i in range(0,nx):
        xabs=math.fabs(x[i])
        if(xabs>w0/2.):
            h[i]=(xabs-w0/2.)*dhdxinit

    for step in range(0,steps):
        # store bound indices
        #  - X=0/1: to left/right of wild-type/mutant boundary
        #  - Y=A/B: to lattice site to left/right of each boundary 
        #  - biXY: index corresponding to boundary 
        #  - corrX: proximity to lattice site B (=0 at A, =1 at B)
        num0=(L/2.+b[0])/dx
        num1=(L/2.+b[1])/dx
        bi0A=int(num0)%nx
        bi1A=int(num1)%nx
        bi0B=(bi0A+1)%nx
        bi1B=(bi1A+1)%nx
        corr0=num0-int(num0)
        corr1=num1-int(num1)
        
        # loop over front positions
        for i in range(1,nx):
            # front +/- 1 lattice point
            h1=h[i-1]
            h2=h[i] 
            if(i+1<len(h)): h3=h[i+1]   
            else: h3=h[1]
        
            # calc 1st & 2nd derivatives
            hderiv1=(h3-h1)/dx2
            hderiv2=(h3-2.*h2+h1)/dxsq
            hderiv1sq=math.sqrt(1.+hderiv1**2)
            curv=hderiv2/hderiv1sq**3
            
            # calc flat front vels
            if(x[i]-b[0]<0. or x[i]-b[1]>0.): vrat0=1.
            else: vrat0=1.+s

            # modulate vels based on local curvature
            vrat=vrat0*(1.+T*curv)        
            dh[i]=vrat*dh0*hderiv1sq
            
            # update boundaries 
            if(i==bi0A): b[0]+=-(1.-corr0)*dh[i]*hderiv1/hderiv1sq
            elif(i==bi0B): b[0]+=-corr0*dh[i]*hderiv1/hderiv1sq
            if(i==bi1A): b[1]+=-(1.-corr1)*dh[i]*hderiv1/hderiv1sq
            elif(i==bi1B): b[1]+=-corr1*dh[i]*hderiv1/hderiv1sq

        # diffusion of boundaries w/ extinction
        if(b[1]>b[0]): 
            b[0]+=random.uniform(-1.,1.)*xdiff
            b[1]+=random.uniform(-1.,1.)*xdiff
        else: 
            b[0]=np.mean(b)
            b[1]=np.mean(b)
                
        # update front
        dh[0]=dh[nx-1]
        h+=dh
            
        # set end points to be consistent with wedge slope
        h[0]=h[1]+dhend
        h[nx-1]=h[nx-2]+dhend
        
        # save bounds
        if((step+1)%skipb==0): 
            h0=(1.-corr0)*h[bi0A]+corr0*h[bi0A]
            h1=(1.-corr1)*h[bi1A]+corr1*h[bi1A]
            print>>fbA,b[0],h0,b[1],h1
            fbA.flush()
        if((step+1)%skiph==0): 
            for i in range(0,nx): print>>fbB,x[i],h[i]
            fbB.flush()

# print running time
t2=time.time()
print int(t2-t1), 'seconds'
