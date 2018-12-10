      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!
      !!  Completely overdamped MD for dimers with growth
      !!
      !!
      !!  Cells (1) grow in 'growth layer', (2) are pushed in a
      !!    'propagation layer', (3) are held fixed in a 'boundary 
      !!    layer', & (4) are removed beyond the 'boundary layer'
      !!
      !!  Depths defined as distance to closest cell in front
      !!  
      !!  Front is defined as cells furthest from colony center-of-mass
      !!
      !!  Options - Restart: T/F = use restart file/start from scratch
      !!              Movie: T/F = do/do not output movie
      !!
      !!  F = b*m*dr/dt (m=1 implicit)   
      !!  T = b*I*dth/dt (I=inertia, th=orientation angle)   
      !!
      !!  Carl Schreck
      !!  11/30/2018
      !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      program main

      implicit none
      integer Ntot
      parameter(Ntot=2**16)
      double precision pi
      parameter(pi=3.1415926535897932d0)
      double precision x(Ntot),vx(Ntot),ax(Ntot),bx(Ntot),fx(Ntot),xa(2)
      double precision y(Ntot),vy(Ntot),ay(Ntot),by(Ntot),fy(Ntot),ya(2)
      double precision th(Ntot),vth(Ntot),ath(Ntot),bth(Ntot),fth(Ntot)
      double precision xp(Ntot),yp(Ntot),D(Ntot),alpha(Ntot),rate(Ntot)
      double precision inert(Ntot),depth(Ntot),rate0(Ntot),tdiv,corr
      double precision desync,kinetic,KE,V,cc,ss,alpha0,alphamax,width
      double precision maxdis,dt,ran2,layerwidth,layerdepth,propdepth
      double precision bounddepth,propdist,bounddist,radius,dd,ddsq
      double precision b,rateWT,dr(2),dk(2)
      integer N,seed,steps,i,j,k,countn,nl(12*Ntot,2),kstart,seedstart
      integer restartexist,dataskip,prodskip,div,layerskip,restskip
      integer nrem,forcelist(Ntot),proplist(Ntot),nprop,nsum
      character file1*199,file2*199,file3*199
      logical restart,movie
      common /f1com/ alpha
      common /f2com/ nl,countn 
      common /f3com/ proplist
      common /f5com/ alphamax,alpha0

      ! read geometric parameters
      read(*,*) alpha0
      read(*,*) alphamax

      ! read rates
      read(*,*) rateWT
      read(*,*) b

      ! read steps
      read(*,*) steps
      read(*,*) layerskip
      read(*,*) dataskip
      read(*,*) prodskip
      read(*,*) restskip
      read(*,*) dt

      ! read growth layer parameters
      read(*,*) layerwidth      
      read(*,*) layerdepth

      ! read layer parameters for force calc
      read(*,*) propdepth
      read(*,*) bounddepth

      ! read run parameters
      read(*,*) desync
      read(*,*) seed

      ! read output files
      read(*,*) file1
      read(*,*) file2
      read(*,*) file3

      ! read options
      read(*,*) movie
      read(*,*) restart

      write(*,*) "file1=", file1

      ! distances of propagation/boundary layer from front
      propdist=layerdepth+propdepth
      bounddist=layerdepth+propdepth+bounddepth

      width=0.2d0 ! width of neighborlist 
      tdiv=dlog10(2d0)/dlog10(1d0+dt*rateWT)*dt ! # steps in 1 generation
 
      ! initialize system from scratch or restart file
      call initialize(file1,file2,file3,restart,movie,rateWT,
     +     rate0,desync,kstart,seed,b,N,rate,depth,inert,d,x,y,
     +     th,vx,vy,vth,ax,ay,ath,bx,by,bth,xp,yp,width,nsum)

      ! loop over time-steps      
      do k=kstart+1,steps
         ! grow/divide cells
         call grow(dt,N,nsum,depth,layerdepth,rate,rate0,rateWT,width,
     +        D,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,xp,yp,seed,desync)

         ! calc propagation list
         call calc_proplist(N,nprop,depth,proplist,
     +        vx,vy,vth,ax,ay,ath,bx,by,bth,propdist)

         ! remove cells & make neighbor list
         call checklist(N,x,y,xp,yp,maxdis)
         if(maxdis.gt.width*d(1)) then
            call remove(N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,d,
     +           alpha,depth,rate,rate0,inert,proplist,bounddist)
            call makelist(N,x,y,d,xp,yp,width)
         endif

         ! calculate inertia of each cell
         call calc_inert(N,inert,D)
       
         ! Gear precictor-corrector
         call predict(dt,N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth)
         call force(N,x,y,th,d,V,fx,fy,fth)           
         call correct(dt,N,x,y,th,vx,vy,vth,ax,ay,ath,
     +        bx,by,bth,fx,fy,fth,inert,b)
         KE=kinetic(N,vx,vy,vth,inert)

         ! calc distance to front     
         if(mod(k,layerskip).eq.0) then
            call calcdepth_radial(N,x,y,d,layerwidth,depth)
         endif         
         
         ! output data to screen
         if(mod(k,dataskip).eq.0) then
            write(*,'(ES20.12,3I,3ES20.12)') dble(k)*dt/tdiv,
     +           nsum,N,nprop,V/dble(nprop),KE/dble(nprop)
         endif         

         ! save cell positions for movie
         if(movie.and.mod(k,prodskip).eq.0) then
            write(1,*) 2*N
            do i=1,N
               cc=dcos(th(i))
               ss=dsin(th(i))
               dd=alpha(i)-1d0
               dr(1)=(1d0+dd)/(1d0+dd**2)*dd**2*D(i)/2d0
               dr(2)=-(1d0+dd)/(1d0+dd**2)*D(i)/2d0
               do j=1,2
                  xa(j)=x(i)+dr(j)*cc
                  ya(j)=y(i)+dr(j)*ss
               enddo
               dk(1)=D(i)
               dk(2)=dd*D(i)
               write(1,'(4F)') xa(1),ya(1),dk(1),depth(i)
               write(1,'(4F)') xa(2),ya(2),dk(2),depth(i)
            enddo
            flush(1)
         endif         

         ! save restart file
         if(restart.and.mod(k,restskip).eq.0) then
            open(unit=2,file=TRIM(file2))
            write(2,'(4I)') k, N, nsum, seed
            do i=1,N
               write(2,'(17E26.18)') x(i),y(i),th(i),vx(i),vy(i),
     +              vth(i),ax(i),ay(i),ath(i),bx(i),by(i),bth(i),
     +              d(i),alpha(i),depth(i),rate(i),rate0(i)
            enddo
            flush(2)
            close(2)
         endif   

         ! output radius
         if(mod(k,dataskip).eq.0) then
            call calc_radius(N,x,y,depth,radius)
            write(3,'(E,2F)') dble(k)*dt/tdiv,radius            
            flush(3)      
         endif
      enddo
      
      end ! end main

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!  initialize cell position & momenta  !!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine initialize(file1,file2,file3,restart,movie,rateWT,
     +     rate0,desync,kstart,seed,b,N,rate,depth,inert,d,x,y,
     +     th,vx,vy,vth,ax,ay,ath,bx,by,bth,xp,yp,width,nsum)

      integer Ntot
      parameter(Ntot=2**16)
      double precision pi
      parameter(pi=3.1415926535897932d0)
      double precision x(Ntot),y(Ntot),th(Ntot),vx(Ntot),vy(Ntot)
      double precision vth(Ntot),ax(Ntot),ay(Ntot),ath(Ntot),bx(Ntot)
      double precision by(Ntot),bth(Ntot),fx(Ntot),fy(Ntot),fth(Ntot)
      double precision inert(Ntot),depth(Ntot),d(Ntot),alpha(Ntot),V
      double precision rate(Ntot),b,dd,ddsq,tmp,rateWT,ran2,rate0(Ntot)
      double precision alpha0,alphamax,desync,xp(Ntot),yp(Ntot),width
      integer N,kstart,seed,seedstart,i,proplist(Ntot),nsum
      logical restartexist,restart,movie
      character file1*199,file2*199,file3*199
      common /f1com/ alpha
      common /f3com/ proplist
      common /f5com/ alphamax,alpha0

      ! check if restart file exists
      inquire(file=file2,exist=restartexist)
      if(restart.and.restartexist) then  
        ! open files
         if(movie) open(unit=1,file=TRIM(file1),ACCESS="APPEND")
         if(restart) open(unit=2,file=TRIM(file2))
         open(unit=3,file=TRIM(file3),ACCESS="APPEND")
         
         ! read restart file
         read(2,*) kstart, N, nsum, seedstart
         do i=1,N
            read(2,*) x(i),y(i),th(i),vx(i),vy(i),vth(i),
     +           ax(i),ay(i),ath(i),bx(i),by(i),bth(i),
     +           d(i),alpha(i),depth(i),rate(i),rate0(i)
         enddo
         close(2)

         ! calculate inertia of each cell
         call calc_inert(N,inert,D)
 
         ! burn seeds
         do while (seed.ne.seedstart)
            tmp=ran2(seed)
         enddo
      else ! no restart file exists
         kstart=0

         ! open files
         open(unit=1,file=TRIM(file1))
         open(unit=3,file=TRIM(file3))

         ! random initial config
         N=2
         nsum=N
         do i=1,N
            d(i)=1d0
            x(i)=dble(i)-dble(N+1)/2d0
            y(i)=0d0
            th(i)=(ran2(seed)-0.5d0)*2d0*pi
            depth(i)=0d0
            proplist(i)=1
         enddo

         ! assign initial rates
         do i=1,N     
            rate0(i)=rateWT               
            rate(i)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0(i)
         enddo

         ! assign initial aspect ratios
         do i=1,N     
            alpha(i)=alpha0*(1d0+ran2(seed))
         enddo
         
         ! calculate inertia of each cell
         call calc_inert(N,inert,D)
         call makelist(N,x,y,d,xp,yp,width)
         call force(N,x,y,th,d,V,fx,fy,fth)  
         do i=1,N
            vx(i)=b*fx(i)
            vy(i)=b*fy(i)
            vth(i)=b*fth(i)/inert(i)
            ax(i)=0d0
            ay(i)=0d0
            ath(i)=0d0
            bx(i)=0d0
            by(i)=0d0
            bth(i)=0d0         
         enddo

      endif

      return
      end ! end initialize


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!    grow & divide cells    !!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine grow(dt,N,nsum,depth,layerdepth,rate,rate0,rateWT,
     +   width,D,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,xp,yp,seed,desync) 

      integer Ntot
      parameter(Ntot=2**16)
      double precision pi
      parameter(pi=3.1415926535897932d0)
      double precision x(Ntot),y(Ntot),th(Ntot),vx(Ntot),vy(Ntot)
      double precision vth(Ntot),ax(Ntot),ay(Ntot),ath(Ntot),bx(Ntot)
      double precision by(Ntot),bth(Ntot),depth(Ntot),d(Ntot),rate(Ntot)
      double precision alpha(Ntot),rate0(Ntot),xp(Ntot),yp(Ntot)
      double precision dt,ran2,alpha0,alphamax,desync,width,radius
      double precision rateWT,layerdepth,corr
      integer N,nsum,seed,i
      common /f1com/ alpha
      common /f5com/ alphamax,alpha0

      do i=1,N
         ! grow cell i
         if(depth(i).lt.layerdepth) then
            corr=(1d0+(alpha(i)-1d0)**2)/2d0/(alpha(i)-1d0)
            alpha(i)=alpha(i)+corr*dt*rate(i)            
         endif

         ! divide cell i
         if(alpha(i).gt.alphamax) then
            ! divide into 2 - N=current cels, nsum=total 
            N=N+1
            nsum=nsum+1

            ! set cell growth rates
            rate0(i)=rateWT
            rate0(N)=rateWT

            ! divide into 2 - 1st assigned index N+1
            D(N)=D(i)
            x(N)=x(i)+alpha0/2d0*dcos(th(i))
            y(N)=y(i)+alpha0/2d0*dsin(th(i))
            rate(N)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0(N)
            alpha(N)=alpha0
            vx(N)=vx(i)
            vy(N)=vy(i)
            vth(N)=vth(i)               
            ax(N)=ax(i)
            ay(N)=ay(i)
            ath(N)=ath(i)               
            bx(N)=bx(i)
            by(N)=by(i)
            bth(N)=bth(i)

            ! divide into 2 - 2nd assigned index i
            x(i)=x(i)-alpha0/2d0*dcos(th(i))
            y(i)=y(i)-alpha0/2d0*dsin(th(i))
            rate(i)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0(i) 
            alpha(i)=alpha0

            ! axial division 
            th(N)=th(i)
            th(i)=th(i)+pi
            
            ! update depth
            depth(N)=depth(i)+alpha0/2d0*dsin(th(i))
            depth(i)=depth(i)-alpha0/2d0*dsin(th(i))

            ! update neighbor list
            call makelistind(N,N,x,y,d,xp,yp,width)
         endif
      enddo

      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!  remove cells & association information  !!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine remove(N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,d,
     +     alpha,depth,rate,rate0,inert,proplist,bounddist)
      
      integer Ntot
      parameter(Ntot=2**16)
      double precision x(Ntot),y(Ntot),th(Ntot),vx(Ntot),vy(Ntot)
      double precision vth(Ntot),ax(Ntot),ay(Ntot),ath(Ntot),bx(Ntot)
      double precision by(Ntot),bth(Ntot),depth(Ntot),d(Ntot),rate(Ntot)
      double precision alpha(Ntot),inert(Ntot),bounddist,rate0(Ntot)
      integer N,nrem,i,j,proplist(Ntot)

      nrem=0
      do i=N,1,-1
         if(depth(i).gt.bounddist) then
            nrem=nrem+1
            do j=i+1,N
               x(j-1)=x(j)
               y(j-1)=y(j)
               th(j-1)=th(j)
               vx(j-1)=vx(j)
               vy(j-1)=vy(j)
               vth(j-1)=vth(j)
               ax(j-1)=ax(j)
               ay(j-1)=ay(j)
               ath(j-1)=ath(j)
               bx(j-1)=bx(j)
               by(j-1)=by(j)
               bth(j-1)=bth(j)
               d(j-1)=d(j)
               alpha(j-1)=alpha(j)
               depth(j-1)=depth(j)
               rate(j-1)=rate(j)
               rate0(j-1)=rate0(j)
               inert(j-1)=inert(j)
               proplist(j-1)=proplist(j)
            enddo
         endif
      enddo
      N=N-nrem
         
      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!  calc propagation list  !!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine calc_proplist(N,nprop,depth,proplist,
     +        vx,vy,vth,ax,ay,ath,bx,by,bth,propdist)
 
      integer Ntot
      parameter(Ntot=2**16)
      double precision vx(Ntot),vy(Ntot),vth(Ntot),ax(Ntot)
      double precision ay(Ntot),ath(Ntot),bx(Ntot),by(Ntot)
      double precision bth(Ntot),depth(Ntot),propdist
      integer N,nprop,i,proplist(Ntot)

      nprop=0
      do i=1,N
         if(depth(i).lt.propdist) then
            proplist(i)=1
            nprop=nprop+1
         else
            proplist(i)=0
            vx(i)=0d0
            vy(i)=0d0
            vth(i)=0d0
            ax(i)=0d0
            ay(i)=0d0
            ath(i)=0d0
            bx(i)=0d0
            by(i)=0d0
            bth(i)=0d0
         endif
      enddo

      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!  calc moment of interia of each cell  !!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine calc_inert(N,inert,D)
      integer Ntot
      parameter(Ntot=2**16)
      double precision inert(Ntot),D(Ntot),alpha(Ntot),dd,ddsq
      integer i,N
      common /f1com/ alpha
      
      do i=1,N
         dd=alpha(i)-1d0
         ddsq=dd*dd
         inert(i)=((1d0+ddsq**2)/(1d0+ddsq)+2d0*ddsq*
     +        (1d0+dd)**2/(1d0+ddsq)**2)*d(i)**2/8d0
      enddo

      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!  check max displacement to update list  !!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine checklist(N,x,y,xp,yp,maxdis)
      integer Ntot
      parameter(Ntot=2**16)
      double precision maxdis,x(Ntot),y(Ntot),xp(Ntot),yp(Ntot),df
      integer N

      df=2d0

      maxdis=0d0
      do i=1,N
	maxdis=max(dabs(x(i)-xp(i)),maxdis)
	maxdis=max(dabs(y(i)-yp(i)),maxdis)
      enddo
      maxdis=2d0*dsqrt(df*maxdis*maxdis)

      return
      end ! end checklist

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!   make neighbor list   !!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine makelist(N,x,y,d,xp,yp,width)
      integer Ntot
      parameter(Ntot=2**16)
      double precision x(Ntot),y(Ntot),xp(Ntot),yp(Ntot),d(Ntot)
      double precision rij,dij,rijsq,width,di_up(Ntot),alphamax
      double precision alpha0,xij,yij,dijlist
      integer countn,nl(12*Ntot,2),N,i,j
      common /f2com/ nl,countn
      common /f5com/ alphamax,alpha0

      countn=0      
      do i=2,N
         do j=1,i-1
            xij=x(i)-x(j)
            dij=alphamax*d(i) ! max distance - aspect ratio = 2
            dijlist=dij+width*d(1)
            if(dabs(xij).lt.dijlist) then
               yij=y(i)-y(j)
               rijsq=xij*xij+yij*yij
               if(rijsq.lt.dijlist**2) then
                  countn=countn+1
                  nl(countn,1)=i
                  nl(countn,2)=j
               endif
            endif
         enddo
      enddo
      
      do i=1,N
         xp(i)=x(i)
         yp(i)=y(i)
      enddo      
      
      return
      end ! end makelist
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!   make neighbor list only for cell i   !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine makelistind(i,N,x,y,d,xp,yp,width)
      integer Ntot
      parameter(Ntot=2**16)
      double precision x(Ntot),y(Ntot),xp(Ntot),yp(Ntot),d(Ntot)
      double precision rij,dij,rijsq,width,di_up(Ntot)
      double precision xij,yij,dijlist,alphamax,alpha0
      integer countn,nl(12*Ntot,2),N,i,j
      common /f2com/ nl,countn
      common /f5com/ alphamax,alpha0

      do j=1,i-1         
         xij=x(i)-x(j)        
         dij=alphamax*d(i) ! max distance = aspect ratio
         dijlist=dij+width*d(1)
         if(dabs(xij).lt.dijlist) then
            yij=y(i)-y(j)
            rijsq=xij*xij+yij*yij
            if(rijsq.lt.dijlist**2) then
               countn=countn+1
               nl(countn,1)=i
               nl(countn,2)=j
            end if
         endif
      enddo
      
      xp(i)=x(i)
      yp(i)=y(i)
      
      return
      end ! end makelistind
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!    predicts new positions and velocities    !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine predict(dt,N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth)     
      integer Ntot
      parameter(Ntot=2**16)
      integer N,i,proplist(Ntot)
      double precision x(Ntot),y(Ntot),vx(Ntot),vy(Ntot),ax(Ntot)
      double precision ay(Ntot),bx(Ntot),by(Ntot),th(Ntot),vth(Ntot)
      double precision ath(Ntot),bth(Ntot),dt,c1,c2,c3
      common /f3com/ proplist

      c1 = dt
      c2 = c1*dt/2d0
      c3 = c2*dt/3d0

      do i=1,N
         if(proplist(i).eq.1) then 
            x(i) = x(i) + c1*vx(i) + c2*ax(i) + c3*bx(i)
            y(i) = y(i) + c1*vy(i) + c2*ay(i) + c3*by(i)
            th(i) = th(i) + c1*vth(i) + c2*ath(i) + c3*bth(i)         
            vx(i) = vx(i) + c1*ax(i) + c2*bx(i)
            vy(i) = vy(i) + c1*ay(i) + c2*by(i)     
            vth(i) = vth(i) + c1*ath(i) + c2*bth(i)     
            ax(i) = ax(i) + c1*bx(i)
            ay(i) = ay(i) + c1*by(i)
            ath(i) = ath(i) + c1*bth(i)
         endif
      enddo

      end ! end prediction step


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!   corrects prediction   !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine correct(dt,N,x,y,th,vx,vy,vth,ax,ay,ath,
     +     bx,by,bth,fx,fy,fth,inert,b)
      integer Ntot
      parameter(Ntot=2**16)
      integer i,N,proplist(Ntot)
      double precision b,dt,x(Ntot),y(Ntot),vx(Ntot),vy(Ntot),ax(Ntot)
      double precision ay(Ntot),bx(Ntot),by(Ntot),th(Ntot),vth(Ntot)
      double precision ath(Ntot),bth(Ntot),fx(Ntot),fy(Ntot),fth(Ntot)
      double precision inert(Ntot),c1,c2,c3,cg0,cg2,cg3
      double precision gear0,gear2,gear3,corrx,corry,corrth
      common /f3com/ proplist

      gear0 = 3d0/8d0
      gear2 = 3d0/4d0
      gear3 = 1d0/6d0

      c1 = dt
      c2 = c1*dt/2d0
      c3 = c2*dt/2d0

      cg0 = gear0*c1
      cg2 = gear2*c1/c2
      cg3 = gear3*c1/c3

      do i=1,N
         if(proplist(i).eq.1) then 
            vxi = b*fx(i)
            vyi = b*fy(i)
            vthi = b*fth(i)/inert(i)
            corrx = vxi - vx(i)
            corry = vyi - vy(i)
            corrth = vthi - vth(i)
            x(i) = x(i) + cg0*corrx
            y(i) = y(i) + cg0*corry
            th(i) = th(i) + cg0*corrth        
            vx(i) = vxi
            vy(i) = vyi
            vth(i) = vthi
            ax(i) = ax(i) + cg2*corrx
            ay(i) = ay(i) + cg2*corry
            ath(i) = ath(i) + cg2*corrth
            bx(i) = bx(i) + cg3*corrx
            by(i) = by(i) + cg3*corry
            bth(i) = bth(i) + cg3*corrth
         endif
      enddo

      end ! end correction step

              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!    calc kinetic energy    !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      function kinetic(N,vx,vy,vth,inert)
      integer Ntot
      parameter(Ntot=2**16)
      integer i,N,proplist(Ntot)
      double precision vx(Ntot),vy(Ntot),vth(Ntot),inert(Ntot),kinetic
      common /f3com/ proplist

      kinetic=0d0
      do i=1,N
         if(proplist(i).eq.1) then         
            kinetic=kinetic+vx(i)**2+vy(i)**2+inert(i)*vth(i)**2   
         endif
      enddo   
      kinetic=kinetic/2d0

      end ! end kinetic energy calc

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!       calc force         !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine force(N,x,y,th,d,V,fx,fy,fth) ! dimer force
      integer Ntot
      parameter(Ntot=2**16)
      double precision x(Ntot),y(Ntot),th(Ntot),alpha(Ntot),D(Ntot)
      double precision radi_up(Ntot),fx(Ntot),fy(Ntot),fth(Ntot)
      double precision xa(Ntot,2),ya(Ntot,2),dk(Ntot,2),dr(Ntot,2)
      double precision c(Ntot),s(Ntot),f_x,f_y,fc,fr,LJ,dij,rij,xij,yij
      double precision dij_up,d1,fact,fthi,fthj,rijsq,dd,dd2,V,Vij
      integer countn,nl(12*Ntot,2),N,ki,kj,jj,up,down,proplist(Ntot)
      common /f1com/ alpha
      common /f2com/ nl,countn
      common /f3com/ proplist
      
      ! set forces & torques to zero
      do i=1,N
         fx(i)=0d0
         fy(i)=0d0
         fth(i)=0d0
      enddo
      V=0d0
      
      ! convert to from molecules to atoms
      do i=1,N
         c(i)=dcos(th(i))
         s(i)=dsin(th(i))
         dd=alpha(i)-1d0
         dd2=(1d0+dd)/(1d0+dd**2)*D(i)/2d0
         dr(i,1)=dd2*dd**2
         dr(i,2)=-dd2
         do k=1,2
            xa(i,k)=x(i)+dr(i,k)*c(i)
            ya(i,k)=y(i)+dr(i,k)*s(i)
         enddo
         dk(i,1)=D(i)
         dk(i,2)=dd*D(i)
         radi_up(i)=(dk(i,2)-2d0*dr(i,2))/2d0
      enddo

      ! inter-particle interactions   
      d1=D(1)
      do k=1,countn
         i=nl(k,1)
         j=nl(k,2)         
         if(proplist(i).eq.1.or.proplist(j).eq.1) then 
            dij_up=radi_up(i)+radi_up(j)
            xij=x(i)-x(j)
            if(dabs(xij).lt.dij_up) then 
               yij=y(i)-y(j)        
               rijsq=xij**2+yij**2
               if(rijsq.lt.dij_up*dij_up) then
                  do ki=1,2
                     do kj=1,2
                        dij=(dk(i,ki)+dk(j,kj))/2d0
                        xij=xa(i,ki)-xa(j,kj)
                        yij=ya(i,ki)-ya(j,kj)
                        rijsq=xij**2+yij**2
                        if(rijsq.lt.dij**2) then
                           rij=dsqrt(rijsq)

                           ! forces from Schreck et al, Soft Matter 2010
                           fc=(1d0-rij/dij)/dij     
                           Vij=(1d0-rij/dij)**2/2d0

                           ! modify so that force ~ overlap
                           fact=(dij/d1)**2
                           fr=fc/rij*fact
                           V=V+Vij*fact                     

                           ! calc forces & torques
                           f_x=fr*xij
                           f_y=fr*yij
                           if(proplist(i).eq.1) then
                              fx(i)=fx(i)+f_x
                              fy(i)=fy(i)+f_y
                              fth(i)=fth(i)+dr(i,ki)*(c(i)*f_y-s(i)*f_x)
                           endif
                           if(proplist(j).eq.1) then
                              fx(j)=fx(j)-f_x
                              fy(j)=fy(j)-f_y
                              fth(j)=fth(j)-dr(j,kj)*(c(j)*f_y-s(j)*f_x)
                           endif               
                        endif
                     enddo
                  enddo
               endif
            endif
         endif
      enddo
       
      return							
      end ! end force calc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!  check max displacement to update list  !!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine calcdepth_radial(N,x,y,d,layerwidth,depth)
      double precision pi
      parameter(pi=3.1415926535897932d0)
      integer Ntot
      parameter(Ntot=2**16)
      double precision x(Ntot),y(Ntot),d(Ntot),layerwidth,xR,xL,yR,yL
      double precision ymin,ymax,depth(Ntot),dx,rsq(Ntot),xcm,ycm
      double precision circum,binanglewidth,rsqmax(2**16),R,angle,drsq
      integer N,numbins,i,j,bin,dbin,nj,numfrontU,numfrontD,bini(Ntot)
      integer front(Ntot),binangle(Ntot)

      ! calc com of colony
      xcm=0d0
      ycm=0d0
      do i=1,N
         xcm=xcm+x(i)
         ycm=ycm+y(i)
      enddo
      xcm=xcm/dble(N)
      ycm=ycm/dble(N)
      
      ! calc angular bin size
      xR=0d0
      xL=0d0
      yR=0d0
      yL=0d0
      do i=1,N
         if(x(i)-xcm.gt.xR) then
            xR=x(i)-xcm
         elseif(x(i)-xcm.lt.xL) then
            xL=x(i)-xcm
         endif
         if(y(i)-ycm.gt.yR) then
            yR=y(i)-ycm
         elseif(y(i)-ycm.lt.yL) then
            yL=y(i)-ycm
         endif
      enddo
      R=(xR-xL+yR-yL)/4d0+1d0
      circum=2d0*pi*R
      numbins=int(circum/layerwidth)+1
      binanglewidth=2d0*pi/dble(numbins)

      ! calc radial distance to front
      do j=1,numbins
         rsqmax(j)=0d0
      enddo
      do i=1,N
         angle=datan((y(i)-ycm)/(x(i)-xcm))
         if(x(i)-xcm.lt.0d0.and.y(i)-ycm.gt.0d0) then
            angle=angle+pi
         else if(x(i)-xcm.lt.0d0.and.y(i)-ycm.le.0d0) then
            angle=angle-pi
         endif
         angle=angle+pi
         binangle(i)=int(angle/binanglewidth)+1
         rsq(i)=(x(i)-xcm)**2+(y(i)-ycm)**2         
         if(rsq(i).gt.rsqmax(binangle(i))) then
            rsqmax(binangle(i))=rsq(i)
         endif
      enddo      
      do i=1,N            
         depth(i)=dsqrt(rsqmax(binangle(i)))-dsqrt(rsq(i))
      enddo

      ! assign cells near front to be at front
      numfront=0
      do i=1,N
         if(depth(i).lt.D(i)) then 
            numfront=numfront+1
            front(numfront)=i
            depth(i)=0d0
         endif
      enddo
      
      ! calc distance to nearest cell at front
      do i=1,N
         do jj=1,numfront
            j=front(jj)
            dx=x(i)-x(j)
            if(dabs(dx).lt.depth(i)) then
               dy=y(i)-y(j)
               drsq=dx*dx+dy*dy
               if(drsq.lt.depth(i)**2) then
                  depth(i)=dsqrt(drsq)
               endif
            endif
         enddo
      enddo
 
      return
      end ! end depth calc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!      calc msd radius of front      !!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine calc_radius(N,x,y,depth,radius)
      
      integer Ntot
      parameter(Ntot=2**16)
      double precision x(Ntot),y(Ntot),depth(Ntot),radius
      double precision xave,yave,xsqave,ysqave
      integer count,i,N

      xave=0d0
      yave=0d0
      xsqave=0d0
      ysqave=0d0
      count=0
      do i=1,N
         if(depth(i).lt.1d0) then
            count=count+1
            xave=xave+x(i)
            yave=yave+y(i)
            xsqave=xsqave+x(i)*x(i)
            ysqave=ysqave+y(i)*y(i)
         endif
      enddo
      xave=xave/dble(count)
      yave=yave/dble(count)
      xsqave=xsqave/dble(count)
      ysqave=ysqave/dble(count)
      radius=dsqrt(xsqave-xave**2+ysqave-yave**2)+0.5d0

      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!    random number generator    !!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      double precision ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END ! end ran2
