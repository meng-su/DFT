      program main
      implicit none
      real*8  beta,conv,pi
      parameter(pi=3.14159265358979323846264338327950d0)
      real*8  x,f,E,dx,rmax,rmin,xmax,xmin
      integer  nmax,lmax,n,node,ngrid,z,l,i,j,k
      integer  num,pwnum,atm,atn,val,ncore
      real*8,allocatable :: r(:),V(:),wf(:),Veff(:),Vxc(:),VH(:),
     & VHb(:),Eig(:,:),rho(:),wave(:,:),eigen(:),Vloc(:),
     & waveps(:,:),Vps(:),T(:,:)
      integer*4,allocatable :: atom(:,:),A(:,:) , conf(:,:)

      allocate(atom(4,4),A(4,4))
      z=26

c read the configuration of atom, you can change the parameters of the configuration files -- atom.txt
ccc
c Z: Nuclear charge number
c the file atom.txt:
c           s.    p.   d.   f 
c n = 1.    2
c n = 2.    2.    6
c n = 3.    2.    6.   6
c n = 4.    2
ccc
      open (1,file='atom.txt',action='read')
      read (1,*) A
      close(1)
      atom=transpose(A)
      print *, 'the occupation of electron in atom:'
      print '(4i5)',transpose( atom )
      nmax=size(atom,dim=1)
      lmax=size(atom,dim=2)
      allocate(Eig(nmax,lmax))

c prepare the radial grid
      ngrid=501
      allocate(r(ngrid),V(ngrid),wf(ngrid),Veff(ngrid),
     & Vxc(ngrid),VH(ngrid),VHb(ngrid),rho(ngrid)
     &  ,wave(nmax**2/2+nmax,ngrid),
     &   eigen(nmax**2/2+nmax) )
      rmin=10.0d0**(-6.0d0)/(z*1.0d0); rmax=30.0d0
      xmin=log(rmin); xmax=log(rmax)
      dx=(xmax-xmin)/(ngrid-1)*1.0d0
      do i=1,ngrid
         r(i)=dexp(dx*(i-1)+xmin)
      end do
      open(2,file='grid.txt',action='write')
      write (2,*) r
      close(2)

c prepare the initial Veff
      V=-z/r
      Veff=V
c prepare the convergence judgement condition
      VHb(:)=0
c conv=sum((VH-VHb)*r*dx)
      conv=1.0d0
c accelerate the convergence speed
      beta=0.5

c start the iteration of solving K-S Eq.
      do while (conv>10.0d0**(-5.0d0))
c      do k=1,5
c calculate each orbital of K-S Eq.  and calculate the density    
           rho=0.0d0
           k=1 
           do n=1,nmax,1
              do l=0,n-1,1
c                if (atom(n,l+1)/=0) then
                 node=n-l-1
                 call radial
     >                (r,Veff,z,l,node,ngrid
     <                                 ,E,wf)
                 Eig(n,l+1)=E
                 eigen(k)=E
                 wave(k,:)=wf
                 k=k+1
                 rho=1/(4*pi)*wf*wf/r*atom(n,l+1)+rho
c                 print '(a5,3f20.7)','E=',Eig(n,l+1) 
c                end if
              end do
           end do
c           print '(a5,3f20.7)','E=',Eig 
           do j=1,ngrid,1
               VH(j)=(4.0d0*pi)*(sum(rho(1:j)*r(1:j)**3.0d0
     &          *dx)/r(j)+sum(rho(j:ngrid)*r(j:ngrid)**2.0d0*dx))
           end do
           VH=beta*VHb+(1-beta)*VH
           conv=sum((VH-VHb)*r*dx)
           call LDA
     >         (rho,ngrid
     <               ,Vxc)
           Veff=V+Vxc+VH
           VHb=VH
c           print '(a5,3f20.7)','E=',Eig 
      end do
      print '(a5)','E='
      print '(4f20.7)',transpose(Eig)
      open (1,file='orbital.txt',action='write')
      k=1
      do n=1,nmax,1
           do l=0,n-1,1
                if (atom(n,l+1)/=0) then
                    write (1,'(a6,i4,a6,f15.7)') 'occ=',atom(n,l+1),
     &                    'E=',Eig(n,l+1)
                    write (1,*) wave(k,:)
                 if (n==3 .and. l==2) then
                   val=k
                   open(2,file='wavection.txt')
                   write(2,*)wave(k,:)*sqrt(r)
                end if
                    k=k+1
                end if
            end do
      end do
      print '(a10,i10)','val=',val
      close(1)
c      deallocate(A) 
      allocate(conf(3,5),Vloc(ngrid),waveps
     & (k,ngrid),Vps(ngrid),T(ngrid,ngrid)  )
      Vloc=Veff
c      print '(a10,f12.7)','Veff(1)=',Veff(1)
c      print '(a10,f12.7)','Vloc(1)=',Vloc(1)
      conf(1,:)=(/4,7,3,5,6 /)
      conf(2,:)=(/410,434,403,428,433/)
      conf(3,:)=(/0,0,1,1,2 /)

c      conf=transpose(A)
      pwnum=5;
      num=size(wave,dim=1)
      atm=size(atom,dim=1)
      atn=size(atom,dim=2)
      do i=1,ngrid-1
           if ((r(i)-2.0d0)*(r(i+1)-2.0d0)<=0.0d0) then
              ncore=i
           end if
      end do
      do i=1,num
           wave(i,:)=wave(i,:)*sqrt(r)
      end do
ccc calculate the pseudopotential
c      call ps
c     > (wave,conf,r,num,pwnum,ngrid,eigen,Vloc,atom,
c     >  atm,atn,val,ncore,
c     <  waveps,Vps,T )

      
      end program main

     


