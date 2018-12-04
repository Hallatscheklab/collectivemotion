      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!
      !!  Completely overdamped MD for dimers with growth
      !!
      !!
      !!  BC's = PBC in x
      !!  Cells (1) grow in growth layer, (2) are pushed in 
      !!    propagation layer, (3) are held fixed in boundary 
      !!    layer, & (4) are removed beyond boundary layer
      !!
      !!  Depths defined as distance to closest cell in front
      !!  
      !!  Options - Restart: T/F = use restart file/start from scratch
      !!              Movie: T/F = do/do not output movie
      !!             Bottom: T/F = keep/discard bottom
      !!
      !!  F = b*m*dr/dt (m=1 implicit)   
      !!  T = b*I*dth/dt (I=inertia, th=orientation angle)   
      !!
      !!  Carl Schreck
      !!  4/7/2016
      !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      program main

      implicit none
      integer Ntot,Ndivtot
      parameter(Ntot=2**21,Ndivtot=2**24)
      double precision pi
      parameter(pi=3.1415926535897932d0)
      double precision x(Ntot),vx(Ntot),ax(Ntot),bx(Ntot),fx(Ntot),xa(2)
      double precision y(Ntot),vy(Ntot),ay(Ntot),by(Ntot),fy(Ntot),ya(2)
      double precision th(Ntot),vth(Ntot),ath(Ntot),bth(Ntot),fth(Ntot)
      double precision xp(Ntot),yp(Ntot),D(Ntot),alpha(Ntot),rate(Ntot)
      double precision inert(Ntot),depth(Ntot),rate0(Ntot),tdiv,alpha0
      double precision b,exp,desync,kinetic,KE,V,ran2,cc,ss,corr,ddsq
      double precision layerwidth,layerdepth,dr(2),dk(2),maxdis,alphamax
      double precision width,propdepth,bounddepth,propdist,bounddist,dd
      double precision dt,att,rateWT,w0,s,xR,xL,yR,yL,R,xcm,ycm
      double precision angle,angle0,burnR
      integer N,seed,steps,i,j,k,countn,nl(12*Ntot,2),kstart,seedstart
      integer restartexist,dataskip,prodskip,div,layerskip,restskip
      integer nrem,forcelist(Ntot),proplist(Ntot),nprop,nsum,Nc,Ncb
      integer Np,Npb,Nf,Nfb,Nu,Nub,Nc2,Ncb2,celltype(Ntot),mutate
c      integer tree(Ntot,2),divindex(Ntot),ndiv,treeskip
c      integer tree(Ndivtot,2),divindex(Ndivtot),ndiv,treeskip
      integer tree(Ndivtot,2),ndiv,treeskip,ii
      integer divindex(Ntot),divtocellindex(Ndivtot)
      double precision divrad(Ndivtot),divdepth(Ndivtot),radius
      character file1*99,file2*99,file3*99,file4*99,file5*99
      logical restart,movie,bottom
      common /f1com/ exp,alpha
      common /f2com/ nl,countn 
      common /f3com/ proplist
      common /f4com/ bottom
      common /f5com/ alphamax,alpha0

      ! read geometric parameters
      read(*,*) alpha0
      read(*,*) alphamax
      read(*,*) att

      ! read rates
      read(*,*) rateWT
      read(*,*) b

      ! read steps
      read(*,*) steps
      read(*,*) layerskip
      read(*,*) dataskip
      read(*,*) prodskip
      read(*,*) restskip
      read(*,*) treeskip
      read(*,*) dt

      ! read radii where mutants introduced
      read(*,*) burnR

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
c     read(*,*) file4
c      read(*,*) file5

      ! read options
      read(*,*) movie
      read(*,*) restart
      read(*,*) bottom

      ! rescue parameterss
      read(*,*) s
      read(*,*) w0

      ! parameters
      exp=2d0     ! 2 = LS, 2.5 = Hertzian, >2.9 = RLJ
      width=0.2d0 ! width of neighborlist 

      ! calculate # steps until division
      tdiv=dlog10(2d0)/dlog10(1d0+dt*rateWT)*dt

      ! distances of propagation/boundary layer from front
      propdist=layerdepth+propdepth
      bounddist=layerdepth+propdepth+bounddepth
 
      ! initialize system from scratch or restart file
      call initialize(file1,file2,file3,file4,file5,restart,
     +     movie,rateWT,rate0,desync,kstart,seed,b,att,N,rate,
     +     depth,inert,d,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,
     +     0d0,celltype,0d0,mutate,tree,divindex,divtocellindex,
     +     divrad,divdepth,ndiv,xp,yp,width,nsum)
 
c      ! total # of cells (including those removed)
c      nsum=N



c      write(*,*) 111111
c      do i=1,N
c         write(1010101,*) sqrt(x(i)**2+y(i)**2), depth(i)
c      enddo
c      call calcdepth_radial(N,x,y,d,layerwidth,depth)
c      do i=1,N
c         write(2020202,*) sqrt(x(i)**2+y(i)**2), depth(i)
c      enddo
c      write(*,*) 222222


      ! loop over time-steps      
      do k=kstart+1,steps

c         ! mutate
c         if(mutate.eq.0) then
c            ! calc com of colony
c            xcm=0d0
c            ycm=0d0
c            do i=1,N
c               xcm=xcm+x(i)
c               ycm=ycm+y(i)
c            enddo
c            xcm=xcm/dble(N)
c            ycm=ycm/dble(N)
c      
c            ! calc colony radius
c            xR=0d0
c            xL=0d0
c            yR=0d0
c            yL=0d0
c            do i=1,N
c               if(x(i)-xcm.gt.xR) then
c                  xR=x(i)-xcm
c               elseif(x(i)-xcm.lt.xL) then
c                  xL=x(i)-xcm
c               endif
c               if(y(i)-ycm.gt.yR) then
c                  yR=y(i)-ycm
c               elseif(y(i)-ycm.lt.yL) then
c                  yL=y(i)-ycm
c               endif
c            enddo
c            R=(xR-xL+yR-yL)/4d0+1d0
c
c            ! mutate at radius burnR
c            if(R.gt.burnR) then   
c               mutate=1
c
c               angle0=w0/R
c               do i=1,N
c                  angle=datan((y(i)-ycm)/(x(i)-xcm))
c                  if(x(i)-xcm.lt.0d0.and.y(i)-ycm.gt.0d0) then
c                     angle=angle+pi
c                  else if(x(i)-xcm.lt.0d0.and.y(i)-ycm.le.0d0) then
c                     angle=angle-pi
c                  endif
c                  if(dabs(angle-pi/2d0).gt.angle0/2d0) then
c                     celltype(i)=0
c                     rate0(i)=rateWT               
c                  else
c                     celltype(i)=1
c                     rate0(i)=(1d0+s)*rateWT
c                  endif
c               enddo
c            endif
c         endif

         ! grow/divide cells
         call grow(dt,N,nsum,depth,layerdepth,rate,rate0,rateWT,width,
     +     att,D,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,xp,yp,seed,desync,
     +     celltype,s,tree,divindex,divtocellindex,divrad,divdepth,ndiv)

         ! calc propagation list
         call calc_proplist(N,nprop,depth,proplist,
     +        vx,vy,vth,ax,ay,ath,bx,by,bth,propdist)

         ! remove cells & make neighbor list
         call checklist(N,x,y,xp,yp,maxdis)
         if(maxdis.gt.width*d(1)) then
            call remove(N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,
     +           d,alpha,depth,rate,rate0,celltype,inert,
     +           proplist,bounddist,divindex,divtocellindex)
            call makelist(N,x,y,d,xp,yp,width,att)
         endif

         ! calculate inertia of each cell
         call calc_inert(N,inert,D)
       
         ! Gear precictor-corrector
         call predict(dt,N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth)
         call force(N,x,y,th,d,V,fx,fy,fth,att)           
         call correct(dt,N,x,y,th,vx,vy,vth,ax,ay,ath,
     +        bx,by,bth,fx,fy,fth,inert,b)
         KE=kinetic(N,vx,vy,vth,inert)

         ! calc distance to front     
         if(mod(k,layerskip).eq.0) then
            call calcdepth_radial(N,x,y,d,layerwidth,depth)
         endif         
         
         ! output data to screen
         if(mod(k,dataskip).eq.0) then
c            write(*,'(ES20.12,3I,2ES20.12)') dble(k)*dt/tdiv,
c     +           nsum, N, nprop, V/dble(nprop), KE/dble(nprop)
            write(*,'(ES20.12,3I,2ES20.12)') dble(k)*dt/tdiv,
     +           nsum, N, nprop, V/dble(nprop), KE/dble(nprop)            
         endif         
 
         ! save config
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
c               write(1,'(4F,I)') xa(1),ya(1),dk(1),depth(i),celltype(i)
c               write(1,'(4F,I)') xa(2),ya(2),dk(2),depth(i),celltype(i)
               write(1,'(4F,I)') xa(1),ya(1),dk(1),depth(i),celltype(i)
               write(1,'(4F,I)') xa(2),ya(2),dk(2),depth(i),celltype(i)
            enddo
            flush(1)
         endif         
 
         ! save restart file
         if(restart.and.mod(k,restskip).eq.0) then
            open(unit=2,file=TRIM(file2))
c            write(2,'(4I)') k, N, seed, mutate, ndiv
            write(2,'(6I)') k, N, nsum, seed, mutate, ndiv
            do i=1,N
c               write(2,'(17E26.18,I)') x(i),y(i),th(i),vx(i),
c     +              vy(i),vth(i),ax(i),ay(i),ath(i),bx(i),
c     +              by(i),bth(i),d(i),alpha(i),depth(i),
c     +              rate(i),rate0(i),celltype(i),divindex(i)
               write(2,'(17E26.18,2I)') x(i),y(i),th(i),vx(i),
     +              vy(i),vth(i),ax(i),ay(i),ath(i),bx(i),
     +              by(i),bth(i),d(i),alpha(i),depth(i),
     +              rate(i),rate0(i),celltype(i),divindex(i)               
            enddo
            write(2,*) ndiv
            do i=1,ndiv
               write(2,'(3I,2F)') tree(i,1),tree(i,2),divtocellindex(i),
     +              divrad(i),divdepth(i)
            enddo            
            flush(2)
            close(2)
         endif   

         ! save time-series data
         if(mod(k,dataskip).eq.0) then
c            call contacts(N,x,y,th,d,depth,layerdepth,
c     +           Np,Npb,Nc,Ncb,Nc2,Ncb2,Nf,Nfb,Nu,Nub) 
cc            write(3,'(E,12I)') dble(k)*dt/tdiv,Np,Nc,Nc2,Nf,Nu,
cc     +           Npb,Ncb,Ncb2,Nfb,Nub
c            write(3,'(E,12I)') dble(k)*dt/tdiv,Np,Nc,Nc2,Nf,Nu,
c     +           Npb,Ncb,Ncb2,Nfb,Nub            
c            flush(3)

c            write(4,*) nprop, V/dble(nprop) ! replace w growth layer
c            flush(4)

            ! calc radius
            call calc_radius(N,x,y,depth,radius)
            write(3,'(E,2F)') dble(k)*dt/tdiv,radius         
            flush(3)
         endif

c         ! save tree file
c         if(mod(k,treeskip).eq.0) then
c            call calc_radius(N,x,y,depth,radius)
c            write(5,*) ndiv,radius
c            do div=1,ndiv
c               if(tree(div,2)==0) then
c                  write(5,'(2I,2ES16.8,3I)')tree(div,1),tree(div,2),
c     +                 divrad(div),divdepth(div),0,0,0
c               else
c                  ii=divtocellindex(div) 
c                  if(ii.eq.-1) then
c                     write(5,'(2I,2ES16.8,3I)')tree(div,1),tree(div,2),
c     +                    divrad(div),divdepth(div),-1,-1,-1
c                  else
c                     write(5,'(2I,5ES16.8)')tree(div,1),tree(div,2),
c     +                   divrad(div),divdepth(div),x(ii),y(ii),depth(ii)
c                  endif
c               endif
c            enddo
c            flush(5)            
c         endif
      enddo
      
      end ! end main


      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!  initialize cell position & momenta  !!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine initialize(file1,file2,file3,file4,file5,restart,movie,
     +    rateWT,rate0,desync,kstart,seed,b,att,N,rate,depth,inert,d,x,
     +    y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,w0,celltype,s,mutate,tree,
     +    divindex,divtocellindex,divrad,divdepth,ndiv,xp,yp,width,nsum)

      integer Ntot,Ndivtot
      parameter(Ntot=2**21,Ndivtot=2**24)
c      parameter(Ntot=2**21)
      double precision pi
      parameter(pi=3.1415926535897932d0)
      double precision x(Ntot),y(Ntot),th(Ntot),vx(Ntot),vy(Ntot)
      double precision vth(Ntot),ax(Ntot),ay(Ntot),ath(Ntot),bx(Ntot)
      double precision by(Ntot),bth(Ntot),fx(Ntot),fy(Ntot),fth(Ntot)
      double precision inert(Ntot),depth(Ntot),d(Ntot),alpha(Ntot),V
      double precision rate(Ntot),exp,att,b,dd,ddsq,tmp,w0,rateWT
      double precision ran2,rate0(Ntot),alpha0,alphamax,desync,s
      double precision xp(Ntot),yp(Ntot),width
c      integer N,kstart,seed,restartexist,seedstart,i,proplist(Ntot)

      integer N,kstart,seed,seedstart,i,proplist(Ntot),nsum
      logical restartexist

c      integer celltype(Ntot),mutate,tree(Ntot,2),divindex(Ntot),ndiv         
      integer celltype(Ntot),mutate,ndiv,tree(Ndivtot,2)
c      integer divindex(Ndivtot)

      integer divindex(Ntot),divtocellindex(Ndivtot)
      double precision divrad(Ndivtot),divdepth(Ndivtot)
      character file1*99,file2*99,file3*99,file4*99,file5*99
      logical restart,movie
      common /f1com/ exp,alpha
      common /f3com/ proplist
      common /f4com/ bottom
      common /f5com/ alphamax,alpha0

      ! check if restart file exists
      inquire(file=file2,exist=restartexist)
      if(restart.and.restartexist) then  
         ! open files
         if(movie) open(unit=1,file=TRIM(file1),ACCESS="APPEND")
         if(restart) open(unit=2,file=TRIM(file2))
         open(unit=3,file=TRIM(file3),ACCESS="APPEND")
c         open(unit=4,file=TRIM(file4),ACCESS="APPEND")
c         open(unit=5,file=TRIM(file5),ACCESS="APPEND")
         
         ! read restart file
         read(2,*) kstart, N, nsum, seedstart, mutate
         do i=1,N
            read(2,*) x(i),y(i),th(i),vx(i),vy(i),vth(i),ax(i),ay(i),
     +           ath(i),bx(i),by(i),bth(i),d(i),alpha(i),depth(i),
     +           rate(i),rate0(i),celltype(i),divindex(i)
         enddo
            
         call calc_inert(N,inert,D)
         call makelist(N,x,y,d,xp,yp,width,att)
         call force(N,x,y,th,d,V,fx,fy,fth,att)  
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
         read(2,*) ndiv
         do i=1,ndiv
            read(2,*) tree(i,1),tree(i,2),
     +           divtocellindex(i),divrad(i),divdepth(i)
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
         mutate=0

         ! open files
         open(unit=1,file=TRIM(file1))
         open(unit=3,file=TRIM(file3))
c         open(unit=4,file=TRIM(file4))
c         open(unit=5,file=TRIM(file5))

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

         ! initial tree information
         do i=1,N
            tree(i,1)=0
            tree(i,2)=1
            divindex(i)=i         
            divtocellindex(i)=i
            divrad(i)=1d0
            divdepth(i)=0d0
         enddo
         ndiv=2
            
         ! initial growth rates
         do i=1,N
            if(dabs(x(i)).gt.w0/2d0) then
               celltype(i)=0
               rate0(i)=rateWT               
            else
               celltype(i)=1
               rate0(i)=(1d0+s)*rateWT
            endif
         enddo

         ! assign initial aspect ratios & rates
         do i=1,N     
            alpha(i)=alpha0*(1d0+dble(i-1)/2d0)
            rate(i)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0(i)
         enddo
         
         ! assign initial aspect ratios & rates
         do i=1,N     
            alpha(i)=alpha0*(1d0+ran2(seed))
         enddo         
         
         ! calculate inertia of each cell
         call calc_inert(N,inert,D)
         call makelist(N,x,y,d,xp,yp,width,att)
         call force(N,x,y,th,d,V,fx,fy,fth,att)  
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
     +     width,att,D,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,
     +     xp,yp,seed,desync,celltype,s,tree,divindex,
     +     divtocellindex,divrad,divdepth,ndiv)
 
      integer Ntot,Ndivtot
      parameter(Ntot=2**21,Ndivtot=2**24)
c      parameter(Ntot=2**21)
      double precision pi
      parameter(pi=3.1415926535897932d0)
      double precision exp,x(Ntot),y(Ntot),th(Ntot),vx(Ntot),vy(Ntot)
      double precision vth(Ntot),ax(Ntot),ay(Ntot),ath(Ntot),bx(Ntot),dt
      double precision by(Ntot),bth(Ntot),depth(Ntot),d(Ntot),rate(Ntot)
      double precision alpha(Ntot),rate0(Ntot),att,ran2,alpha0,alphamax
      double precision layerdepth,corr,desync,s,rateWT,xp(Ntot),yp(Ntot)
      double precision width,radius
      integer N,nsum,seed,i,celltype(Ntot)
c      integer tree(Ntot,2),divindex(Ntot),ndiv
c      integer tree(Ndivtot,2),divindex(Ndivtot),ndiv

      integer tree(Ndivtot,2),ndiv,divindex(Ntot)
      integer divtocellindex(Ndivtot) 
      double precision divrad(Ndivtot),divdepth(Ndivtot)

      common /f1com/ exp,alpha
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

            ! mutate
            if(celltype(i).eq.0) then         
               celltype(i)=0           
               celltype(N)=0           
               rate0(i)=rateWT
               rate0(N)=rateWT
            else if(celltype(i).eq.1) then
               celltype(i)=1
               celltype(N)=1
               rate0(i)=(1+s)*rateWT
               rate0(N)=(1+s)*rateWT
            endif

c            write(*,*) "AA"

            ! divide into 2 - 1st assigned index N+1
            D(N)=D(i)
            x(N)=x(i)+alpha0/2d0*dcos(th(i))
            y(N)=y(i)+alpha0/2d0*dsin(th(i))
            !th(N)=th(i)
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

c           write(*,*) "BB"

            ! divide into 2 - 2nd assigned index i
            x(i)=x(i)-alpha0/2d0*dcos(th(i))
            y(i)=y(i)-alpha0/2d0*dsin(th(i))
            rate(i)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0(i) 
            alpha(i)=alpha0

c           write(*,*) "CC"

            ! axial division 
            th(N)=th(i)
            th(i)=th(i)+pi
            
            ! update depth
            depth(N)=depth(i)+alpha0/2d0*dsin(th(i))
            depth(i)=depth(i)-alpha0/2d0*dsin(th(i))

            ! tree information
            tree(divindex(i),2)=0
            divtocellindex(divindex(i))=0
            do k=1,2
               tree(ndiv+k,1)=divindex(i)
               tree(ndiv+k,2)=1
            enddo
            divindex(i)=ndiv+1
            divindex(N)=ndiv+2
            divtocellindex(ndiv+1)=i
            divtocellindex(ndiv+2)=N

c           write(*,*) "DD"

            ! first div wrong after restart !

            ! spatial info for tree
            call calc_radius(N,x,y,depth,radius)
            divrad(ndiv+1)=radius
            divrad(ndiv+2)=radius
            divdepth(ndiv+1)=depth(i)
            divdepth(ndiv+2)=depth(N)

c            if(tree(ndiv+1,1).eq.4109) then
c               write(*,*) N, radius
c               do i=1,N
c                  write(*,*) x(i), y(i), depth(i)
c               enddo
c
c               write(*,*) 
c         write(*,'(2I,2F)')tree(ndiv+1,1),tree(ndiv+1,2),radius,depth(i)
c         write(*,'(2I,2F)')tree(ndiv+2,1),tree(ndiv+2,2),radius,depth(N)
c            endif

            ndiv=ndiv+2


c           write(*,*) "EE"

            ! update neighbor list
            call makelistind(N,N,x,y,d,xp,yp,width,att)

c           write(*,*) "FF"

c           write(*,*) "GG"

         endif
      enddo

      return
      end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!  initialize cell position & momenta  !!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine remove(N,x,y,th,vx,vy,vth,ax,ay,ath,
     +     bx,by,bth,d,alpha,depth,rate,rate0,celltype,
     +     inert,proplist,bounddist,divindex,divtocellindex)
      
      integer Ntot,Ndivtot
      parameter(Ntot=2**21,Ndivtot=2**24)
c      parameter(Ntot=2**21)
      double precision x(Ntot),y(Ntot),th(Ntot),vx(Ntot),vy(Ntot)
      double precision vth(Ntot),ax(Ntot),ay(Ntot),ath(Ntot),bx(Ntot)
      double precision by(Ntot),bth(Ntot),depth(Ntot),d(Ntot),rate(Ntot)
      double precision alpha(Ntot),inert(Ntot),bounddist,rate0(Ntot)
      integer N,nrem,i,j,proplist(Ntot),celltype(Ntot)
c      integer tree(Ntot,2),divindex(Ntot)
      integer tree(Ndivtot,2),divindex(Ntot),divtocellindex(Ndivtot)

      ! use correct array sizes for index arrays

c     double precision divrad(Ndivtot),divdepth(Ndivtot)

      nrem=0
      do i=N,1,-1
         if(depth(i).gt.bounddist) then
            nrem=nrem+1

            divtocellindex(divindex(i))=-1
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
               celltype(j-1)=celltype(j)
               inert(j-1)=inert(j)
               proplist(j-1)=proplist(j)

               !write(*,*) j, divindex(j), divtocellindex(divindex(j))

               divtocellindex(divindex(j))=j-1
               divindex(j-1)=divindex(j) 
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
      parameter(Ntot=2**21)
      double precision vx(Ntot),vy(Ntot),vth(Ntot),ax(Ntot),ay(Ntot)
      double precision ath(Ntot),bx(Ntot),by(Ntot),bth(Ntot),depth(Ntot)
      double precision propdist
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
!!!!!!!!!!!!!!  check max displacement to update list  !!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine calc_inert(N,inert,D)
      integer Ntot
      parameter(Ntot=2**21)
      double precision inert(Ntot),D(Ntot),alpha(Ntot),dd,ddsq,exp
      integer i,N
      common /f1com/ exp,alpha
      
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
      parameter(Ntot=2**21)
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

      subroutine makelist(N,x,y,d,xp,yp,width,att)
      integer Ntot
      parameter(Ntot=2**21)
      double precision x(Ntot),y(Ntot),xp(Ntot),yp(Ntot),d(Ntot)
      double precision rij,dij,rijsq,width,di_up(Ntot),alphamax
      double precision alpha0,att,xij,yij,dijlist
      integer countn,nl(12*Ntot,2),N,i,j
      common /f2com/ nl,countn
      common /f5com/ alphamax,alpha0

      countn=0      
      do i=2,N
         do j=1,i-1
            xij=x(i)-x(j)
            !xij=xij-dnint(xij/Lx)*Lx    
            dij=alphamax*d(i) ! max distance - aspect ratio = 2
            dijlist=dij+(width+att)*d(1)
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

      subroutine makelistind(i,N,x,y,d,xp,yp,width,att)
      integer Ntot
      parameter(Ntot=2**21)
      double precision x(Ntot),y(Ntot),xp(Ntot),yp(Ntot),d(Ntot)
      double precision rij,dij,rijsq,width,di_up(Ntot),att
      double precision xij,yij,dijlist,alphamax,alpha0
      integer countn,nl(12*Ntot,2),N,i,j
      common /f2com/ nl,countn
      common /f5com/ alphamax,alpha0

      do j=1,i-1         
         xij=x(i)-x(j)
         !xij=xij-dnint(xij/Lx)*Lx            
         dij=alphamax*d(i) ! max distance = aspect ratio
         dijlist=dij+(width+att)*d(1)
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
      parameter(Ntot=2**21)
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

      subroutine correct(dt,N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,
     +     fx,fy,fth,inert,b)
      integer Ntot
      parameter(Ntot=2**21)
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

            

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!           force            !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine force(N,x,y,th,d,V,fx,fy,fth,att) ! dimer force
      integer Ntot
      parameter(Ntot=2**21)
      double precision x(Ntot),y(Ntot),th(Ntot),alpha(Ntot),D(Ntot)
      double precision radi_up(Ntot),fx(Ntot),fy(Ntot),fth(Ntot),V,Vij
      double precision f_x,f_y,fc,fr,LJ,LJ0,exp,dij,rij,xij,yij,dij_up
      double precision fact,att,fthi,fthj,rijsq,c(Ntot),s(Ntot),dd,dd2
      double precision xa(Ntot,2),ya(Ntot,2),dk(Ntot,2),dr(Ntot,2),di1j1
      integer countn,nl(12*Ntot,2),N,ki,kj,jj,up,down,proplist(Ntot)
      common /f1com/ exp,alpha
      common /f2com/ nl,countn
      common /f3com/ proplist
      
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
      do k=1,countn
         i=nl(k,1)
         j=nl(k,2)         
         if(proplist(i).eq.1.or.proplist(j).eq.1) then 
            dij_up=radi_up(i)+radi_up(j)
            xij=x(i)-x(j)
            !xij=xij-dnint(xij/Lx)*Lx
            if(dabs(xij).lt.dij_up) then 
               yij=y(i)-y(j)        
               rijsq=xij**2+yij**2
               if(rijsq.lt.dij_up*dij_up) then
                  di1j1=(dk(i,1)+dk(j,1))/2d0
                  do ki=1,2
                     do kj=1,2
                        dij=(dk(i,ki)+dk(j,kj))/2d0
                        xij=xa(i,ki)-xa(j,kj)
                        !xij=xij-dnint(xij/Lx)*Lx
                        yij=ya(i,ki)-ya(j,kj)
                        rijsq=xij**2+yij**2
                        if(rijsq.lt.(dij+att)**2) then
                           rij=dsqrt(rijsq)
                           if(exp.eq.2d0) then
                              fc=(1d0-rij/dij)/dij     
                              Vij=(1d0-rij/dij)**2/exp
     +                             -(att/dij)**2/exp
                              fact=(dij/di1j1)**2
                           elseif(exp.lt.2.9) then
                              fc=(1d0-rij/dij)/dij     
                              Vij=(1d0-rij/dij)**2/exp
     +                             -(att/dij)**2/exp
                              fact=(dij/di1j1)**exp
                           else
                              LJ=(dij/rij)*(dij/rij)
                              LJ=LJ*LJ*LJ
                              LJ0=(dij/(dij+att))**6
                              fc=1d0/rij*LJ*(LJ-1d0)
                              Vij=(LJ-1d0)**2-(LJ0-1d0)**2
                              fact=(dij/di1j1)**2
                           endif                    
                           fr=fc/rij*fact
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
                           V=V+Vij*fact                     
                        endif
                     enddo
                  enddo
               endif
            endif
         endif
      enddo
 
      if(exp .gt. 2.9) then
         do i=1,N
            fx(i)=fx(i)/6d0
            fy(i)=fy(i)/6d0
            fth(i)=fth(i)/6d0 
         enddo         
         V=V/72d0
      endif
      
      return							
      end ! end force calc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!          contacts            !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine contacts(N,x,y,th,d,depth,layerdepth,
     +     Np,Npb,Nc,Ncb,Nc2,Ncb2,Nf,Nfb,Nu,Nub) 
      integer Ntot
      parameter(Ntot=2**21)
      double precision x(Ntot),y(Ntot),th(Ntot),alpha(Ntot),D(Ntot),exp
      double precision radi_up(Ntot),dij,rij,xij,yij,dij_up,dd,dd2
      double precision c(Ntot),s(Ntot),xa(Ntot,2),ya(Ntot,2),dk(Ntot,2)
      double precision dr(Ntot,2),di1j1,depth(Ntot),layerdepth
      integer countn,nl(12*Ntot,2),N,ki,kj,Np,Npb,Nc,Ncb
      integer Nf,Nfb,Nu,Nub,con(Ntot,2),Nc2,Ncb2,clist(Ntot,2)
      common /f1com/ exp,alpha
      common /f2com/ nl,countn
      
      ! count # cells in growth layer
      Np=0
      Npb=0
      do i=1,N   
         if(depth(i).lt.layerdepth) then
            Np=Np+1
            if(depth(i).gt.D(i)) then
               Npb=Npb+1
            endif 
         endif      
      enddo
      
      ! convert from dimers to monomers
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
      
      ! count contacts 
      Nc=0
      Ncb=0
      do i=1,N
         con(i,1)=0
         con(i,2)=0
      enddo
      do k=1,countn
         i=nl(k,1)
         j=nl(k,2)        
         if(depth(i).lt.layerdepth.or.depth(j).lt.layerdepth) then
            di1j1=(dk(i,1)+dk(j,1))/2d0
            do ki=1,2
               do kj=1,2
                  dij=(dk(i,ki)+dk(j,kj))/2d0
                  xij=xa(i,ki)-xa(j,kj)
                  !xij=xij-dnint(xij/Lx)*Lx
                  yij=ya(i,ki)-ya(j,kj)
                  rijsq=xij**2+yij**2
                  if(rijsq.lt.dij**2) then
                     Nc=Nc+1
                     if(depth(i).gt.D(i).or.depth(j).gt.D(j)) then
                        Ncb=Ncb+1                              
                     endif
                     con(i,ki)=con(i,ki)+1
                     con(j,kj)=con(j,kj)+1
                  endif                           
               enddo
            enddo
         endif
      enddo

      Nf=0
      Nu=0
      Nfb=0
      Nub=0
      do i=1,N
         if(depth(i).lt.layerdepth) then
            if(con(i,1).ge.2.and.con(i,2).ge.2) then
               clist(i,1)=1 
               clist(i,2)=1
            else
               if(con(i,1).le.2.and.con(i,2).le.2) then
                  Nf=Nf+1
                  if(depth(i).gt.D(i)) then
                     Nfb=Nfb+1
                  endif
                  clist(i,1)=0
                  clist(i,2)=0
               else 
                  Nu=Nu+1
                  if(depth(i).gt.D(i)) then
                     Nub=Nub+1
                  endif
                  if(con(i,1).le.1) then 
                     clist(i,1)=0
                     clist(i,2)=1
                  else
                     clist(i,1)=1
                     clist(i,2)=0
                  endif
               endif
            endif
         endif
      enddo

      ! count contacts 
      Nc2=0
      Ncb2=0
      do k=1,countn
         i=nl(k,1)
         j=nl(k,2)        
         if(depth(i).lt.layerdepth.or.depth(j).lt.layerdepth) then
            di1j1=(dk(i,1)+dk(j,1))/2d0
            do ki=1,2
               do kj=1,2
                  if(clist(i,ki).eq.1.and.clist(j,kj).eq.1)then
                     dij=(dk(i,ki)+dk(j,kj))/2d0
                     xij=xa(i,ki)-xa(j,kj)
                     yij=ya(i,ki)-ya(j,kj)
                     rijsq=xij**2+yij**2
                     if(rijsq.lt.dij**2) then
                        Nc2=Nc2+1
                        if(depth(i).gt.D(i).or.depth(j).gt.D(j)) then
                           Ncb2=Ncb2+1                              
                        endif
                        con(i,ki)=con(i,ki)+1
                        con(j,kj)=con(j,kj)+1
                     endif            
                  endif
               enddo
            enddo
         endif
      enddo

      return							
      end ! end contact calc     

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!    calc kinetic energy    !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      function kinetic(N,vx,vy,vth,inert)
      integer Ntot
      parameter(Ntot=2**21)
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


      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!  check max displacement to update list  !!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine calcdepth(N,x,y,d,Lx,layerwidth,depth)
      integer Ntot,LMAX,DEPTHMAX
      parameter(Ntot=2**21,LMAX=5000,DEPTHMAX=1000)
      double precision x(Ntot),y(Ntot),d(Ntot),Lx,layerwidth,xi,xij,yij
      double precision ymin,ymax,depth(Ntot),depthD(Ntot),depthU(Ntot)
      double precision offset,dx,drsq
      integer N,numbins,i,j,bin,dbin,nj,numfrontU,numfrontD
      integer bini(Ntot),npart(LMAX),ibin(LMAX,DEPTHMAX)
      integer frontU(Ntot),frontD(Ntot)
      logical bottom
      common /f4com/ bottom

      numbins=nint(Lx/layerwidth)
      
      ! calc vertical depths
      do bin=1,numbins
         npart(bin)=0
      enddo         
      offset=floor(Lx/2d0/layerwidth)+1
      do i=1,N
         xi=x(i)-dnint(x(i)/Lx)*Lx
         bin=floor(xi/layerwidth)+offset
         bini(i)=bin              ! x bin of cell i
         npart(bin)=npart(bin)+1  ! # cells in bin
         ibin(bin,npart(bin))=i   ! ibin = cell index in bin      
      enddo
      do i=1,N
         ymin=0d0
         ymax=0d0
         do dbin=-1,1
            bin=mod(numbins+bini(i)+dbin-1,numbins)+1
            do nj=1,npart(bin)                  
               j=ibin(bin,nj)
               xij=x(i)-x(j)
               xij=xij-dnint(xij/Lx)*Lx
               if(dabs(xij)<layerwidth) then
                  if(y(j).gt.ymax) then
                     ymax=y(j)
                  elseif(y(j).lt.ymin) then
                     ymin=y(j)
                  endif
               endif
            enddo
         enddo
         depthU(i)=ymax-y(i)
         depthD(i)=y(i)-ymin 
      enddo

      ! assign cells near front to be at front
      numfrontU=0
      numfrontD=0
      do i=1,N
         if(depthU(i).lt.D(i)) then 
            numfrontU=numfrontU+1
            frontU(numfrontU)=i
            depthU(i)=0d0
         endif
         if(depthD(i).lt.D(i)) then
            numfrontD=numfrontD+1
            frontD(numfrontD)=i
            depthD(i)=0d0
         endif
      enddo
      
      ! calc distance to nearest cell at front
      do i=1,N
         do jj=1,numfrontU
            j=frontU(jj)
            dx=x(i)-x(j)
            dx=dx-dnint(dx/Lx)
            if(dabs(dx).lt.depthU(i)) then
               dy=y(i)-y(j)
               drsq=dx*dx+dy*dy
               if(drsq.lt.depthU(i)**2) then
                  depthU(i)=dsqrt(drsq)
               endif
            endif
         enddo
         do jj=1,numfrontD
            j=frontD(jj)
            dx=x(i)-x(j)
            dx=dx-dnint(dx/Lx)
            if(dabs(dx).lt.depthD(i)) then
               dy=y(i)-y(j)
               drsq=dx*dx+dy*dy
               if(drsq.lt.depthD(i)**2) then
                  depthD(i)=dsqrt(drsq)
               endif
            endif
         enddo
         if(bottom) then
            depth(i)=min(depthU(i),depthD(i))   
         else
            depth(i)=depthU(i)
         endif
      enddo

      return
      end ! end depth calc


      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!  check max displacement to update list  !!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine calcdepth_radial(N,x,y,d,layerwidth,depth)
      double precision pi
      parameter(pi=3.1415926535897932d0)
      integer Ntot
      parameter(Ntot=2**21)
      double precision x(Ntot),y(Ntot),d(Ntot),layerwidth,xR,xL,yR,yL
      double precision ymin,ymax,depth(Ntot),dx,rsq(Ntot),xcm,ycm
      double precision circum,binanglewidth,rsqmax(2**21),R,angle,drsq
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
c      xR=0d0
c      xL=0d0
c      yR=0d0
c      yL=0d0
      xR=1d16
      xL=-1d16
      yR=1d16
      yL=-1d16
      do i=1,N
c         if(x(i)-xcm.gt.xR) then
c            xR=x(i)-xcm
c         elseif(x(i)-xcm.lt.xL) then
c            xL=x(i)-xcm
c         endif
c         if(y(i)-ycm.gt.yR) then
c            yR=y(i)-ycm
c         elseif(y(i)-ycm.lt.yL) then
c            yL=y(i)-ycm
c         endif
         if(x(i)-xcm.lt.xR.and.x(i).gt.0d0.and.dabs(y(i)).lt.2d0) then
            xR=x(i)-xcm
        elseif(x(i)-xcm.gt.xL.and.x(i).lt.0d0.and.dabs(y(i)).lt.2d0)then
            xL=x(i)-xcm
         endif
         if(y(i)-ycm.lt.yR.and.y(i).gt.0d0.and.dabs(x(i)).lt.2d0) then
            yR=y(i)-ycm
        elseif(y(i)-ycm.gt.yL.and.x(i).lt.0d0.and.dabs(y(i)).lt.2d0)then
            yL=y(i)-ycm
         endif
      enddo
      R=(xR-xL+yR-yL)/4d0+1d0

      circum=2d0*pi*R
      numbins=int(circum/layerwidth)+1
      binanglewidth=2d0*pi/dble(numbins)

      ! calc radial distance to front
      do j=1,numbins
c         rsqmax(j)=0d0
         rsqmax(j)=1d16
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
c         if(rsq(i).gt.rsqmax(binangle(i))) then
         if(rsq(i).lt.rsqmax(binangle(i))) then
            rsqmax(binangle(i))=rsq(i)
         endif
      enddo      
      do i=1,N            
c         depth(i)=sqrt(rsqmax(binangle(i)))-sqrt(rsq(i))
         depth(i)=sqrt(rsq(i))-sqrt(rsqmax(binangle(i)))
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
 
c      ! calc new depths from nearest cell at front
c      ! to replace previous
c      ! need ibin(,) & npart()
c      ! check bin
c      do i=1,N      
c         dbin=ceiling(depthU(i)/layerwidth)
c         do dbini=-dbin,dbin
c            bin=mod(numbins+bini(i)+dbini-1,numbins)+1
c            do jj=1,npart(bin)
c               j=ibin(bin,jj)
c               if(frontUind(j).eq.1) then
c                  dx=x(i)-x(j)
c                  dx=dx-dnint(dx/Lx)*Lx
c                  if(dabs(dx).lt.depthU(i)) then
c                     dy=y(i)-y(j)
c                     drsq=dx*dx+dy*dy
c                     if(drsq.lt.depthU(i)**2) then
c                        depthU(i)=dsqrt(drsq)
c                     endif
c                  endif
c               endif
c            enddo
c         enddo
c      enddo     

c      do i=1,N
c         depth(i)=0d0
c      enddo

      return
      end ! end depth calc


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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!      calc msd radius of front      !!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine calc_radius(N,x,y,depth,radius)
      
      integer Ntot
      parameter(Ntot=2**21)
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
