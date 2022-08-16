c this sub solve the radial K-S Eq.      
      subroutine radial
     > (r,V,z,l,node,ngrid
     <              ,E,wf)
      integer*4 ngrid,node,z,l
      integer*4 ilc,ncross,i,j,k
      real*8  r(ngrid),V(ngrid)
      real*8  wf(ngrid)
      real*8  E
      real*8  dx,Emid,a,M,yilc,judge
      real*8  y(ngrid),f(ngrid),x(ngrid)
      x=dlog(r)
      dx=x(2)-x(1)
      Emax=0.0d0
      Emin=-100.0*z*z
c      Emin=-10.0d0**(5.0d0)
      Emid=(Emax+Emin)/2.0d0
c      print '(a7,f20.7)','Emid=',Emid
      ilc=ngrid
      y(1)=r(1)**(l+1.0d0) *(1.0d0 - 2.0d0*z*r(1)/(2.0d0*l+2.0d0))
     & / sqrt(r(1));
      y(2)=r(2)**(l+1.0d0) *(1.0d0 - 2.0d0*z*r(2)/(2.0d0*l+2.0d0)) 
     & / sqrt(r(2));
      a=137.0d0;
      do while (abs(Emax-Emin)>=10.0d0**(-4.0d0))
            ncross=0;
            M=1;
            f=(2*M*r*r*(V-Emid)+(l+0.5d0)**2.0d0)/12.0d0*dx*dx
            do i=1,ngrid-1,1
                if ((V(i)-Emid)*(V(i+1)-Emid)<=0) then
                      ilc=i
                end if
            end do
c            print '(a7,i7)','ncross=',ncross
            if (ilc<=2 .or. ilc>=ngrid-2) then
               Emin=-maxval(abs(V));
               Emid=(Emax+Emin)/2.0d0;
            end if
c            print '(a7,f20.7)','Emin=',Emin
c            print '(a7,f20.7)','Emax=',Emax
c            print '(a7,f20.7)','Emid=',Emid
            f=1.0d0-f;
            do i=3,ilc,1
                y(i)=((12-10*f(i-1))*y(i-1)-f(i-2)*y(i-2))/f(i);
                if (abs(y(i))>10.0d0**10.0d0) then
                   norm=abs(y(i));
                   do k=1,i,1
                      y(k)=y(k)/norm;
                   end do
                end if
            if (y(i-1)*y(i)<=0) then
               ncross=ncross+1
            end if
            end do
            yilc=y(ilc)
            if (ncross/=node) then
               if (ncross>node) then
                  Emax=Emid;
               end if
               if (ncross<node) then
                  Emin=Emid;
               end if
               Emid=(Emax+Emin)/2.0d0;
               cycle
            end if
            y(ngrid)=dx/Z;
            y(ngrid-1)= (12-10*f(ngrid))*y(ngrid)/f(ngrid-1);
            do i=ngrid-2,ilc,-1
                y(i)=((12-10*f(i+1))*y(i+1)-f(i+2)*y(i+2))/f(i)
                if (abs(y(i))>10.0d0**10.0d0) then
                   norm=abs(y(i));
                do k=ngrid,i,-1
                   y(k)=y(k)/norm;
                end do
                end if
            end do
            yilc=yilc/y(ilc);
            do i=ilc,ngrid,1
                y(i)=y(i)*yilc;
            end do
            judge=(y(ilc-1)+y(ilc+1)-(14-12*f(ilc))*y(ilc))
     &       /(dx*r(ilc))
            if (judge*y(ilc)>0) then
               Emax=Emid;
            end if
            if (judge*y(ilc)<0) then
               Emin=Emid;
            end if
            if (judge*y(ilc)==0) then
               E=Emid;
               stop
            end if
               Emid=(Emax+Emin)/2.0d0;
c               print '(a7,f20.7)','Emid=',Emid
c               print '(a7,f20.7)','Emin=',Emin
c               print '(a7,f20.7)','Emax=',Emax
c               print '(a10,f20.7)','a=',(-338.0000610-338.0000305)/2.0d0
      end do
            y=y/sqrt(sum(y**2.0d0*r*r*dx))
            E=Emid
            wf=y
      end subroutine radial
