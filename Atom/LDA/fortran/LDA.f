      subroutine LDA
     > (rho,ngrid,
     <       Vxc )
      implicit none
      integer ngrid
      integer i
      real*8 rho(ngrid),Vxc(ngrid)
      real*8,allocatable :: rs(:),epsx(:),epsc(:),ecp(:)
      real*8 cp,cof,rp,pi
      allocate(rs(ngrid),epsx(ngrid),epsc(ngrid),ecp(ngrid))
      pi=3.14159265358979323846264338327950d0
      do i=1,ngrid,1
         if (rho(i)==0.0d0) then
            rho(i)=10.0d0**(-10.0d0)
         end if   
      end do
      rs=(3/(4*pi*rho))**(1.0d0/3.0d0)
      cp=0.0450d0; rp=21.0d0;
      cof=3.0d0/2.0d0/pi*(9*pi/4)**(1.0d0/3.0d0)
      epsx=-cof/rs
      epsc=-cp*dlog(1+rp/rs);
      Vxc=4.0d0/3.0d0*epsx + epsc; 

      end 
