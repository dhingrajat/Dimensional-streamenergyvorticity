**********4-CHIP CASE***************************************************
*****    THIS PROGRAM SOLVES FOR CONJUGATE MIXED CONVECTIVE, LAMINAR
*****    AND FULLY DEVELOPED FLOW OVER A PROTRUDED HEATED BLOCK
*****    IN A VERTICAL CHANNEL (FLUID IS NEWTONIAN)
*****    RADIATION EFFECT IS CONSIDERED SOLVED USING ABSORPTION FACTOR
*****    AN IMPLICIT FINITE DIFFERENCE FALSE TRANSIENT METHOD
*****    HAS BEEN APPLIED USING STREAM FUNCTION-VORTICITY APPROACH.
*****    (GAUSS SEIDEL ITERATION METHOD WITH SUCCESSIVE OVER
*****    RELAXATION HAS BEEN USED FOR SOLVING THE DISCRETIZED EQUATIONS)
************************************************************************
*****
*****    Line 23 , 24 , 271 , 716 , 774 , 805 , 2220
*****    conv_Q , coeff_aV , solve_T , coeff_aT , output
*****
************************************************************************
*****  MAIN PROGRAM FOR COMPLETE SOLUTION
*icvt=1 upwinding and icvt=0 is central difference

      include 'psemiimp.in'

* Same Everywhere in Solid and Fluid zone
* N = Total Number of grids in x-direction
      double precision trv,axl
      trv = 1.2d0
      axl = 12.d0

      dx1=trv/dfloat(N-1)
      dx2=axl/dfloat(M-1)

      nc(1) = int((ycp1/dx2)+1)
      nc(2) = int(nc(1)+xcp/dx1)
      nc(3) = int(nc(2)+(ycp2-ycp1)/dx2)
      nc(4) = int(nc(3)+xcp/dx1)

      nc(5) = int(nc(4)+(ycp3-ycp2)/dx2)
      nc(6) = int(nc(5)+xcp/dx1)
      nc(7) = int(nc(6)+(ycp2-ycp1)/dx2)
      nc(8) = int(nc(7)+xcp/dx1)

      nc(9) = int(nc(8)+(ycp3-ycp2)/dx2)
      nc(10)= int(nc(9)+xcp/dx1)
      nc(11)= int(nc(10)+(ycp2-ycp1)/dx2)
      nc(12)= int(nc(11)+xcp/dx1)

      nc(13)= int(nc(12)+(ycp3-ycp2)/dx2)
      nc(14)= int(nc(13)+xcp/dx1)
      nc(15)= int(nc(14)+(ycp2-ycp1)/dx2)
      nc(16)= int(nc(15)+xcp/dx1)

      nc(17)= int(nc(16)+(axl-(3.d0*ycp3+ycp2-3.d0*ycp1))/dx2)
      nc(18)= int(nc(17)+(trv-2.d0*t_w)/dx1)
      nc(19)= int(nc(18)+axl/dx2)
      nc(20)= int(nc(19)+((trv-2.d0*t_w)/dx1)-1)

C      do k = 1,L
C      M1(k)= int((ycp1+(k-1)*(ycp2-ycp1)+(k-1)*(ycp3-ycp2))/dx2)+1
C      M2(k)= int((ycp2+(k-1)*(ycp2-ycp1)+(k-1)*(ycp3-ycp2))/dx2)+1
C      enddo

      M1(1)= int(ycp1/dx2)+1
      M2(1)= int(ycp2/dx2)+1
      M1(2)= int(ycp3/dx2)+1
	
      do k = 1,L
      M1(k)=M1(1)+(k-1)*(M1(2)-M1(1))
      M2(k)=M2(1)+(k-1)*(M1(2)-M1(1))
*      write(*,33)M1(k),M2(k)
*33    format(1x,i5,2x,i5)
      enddo
      
      N0=int(t_w/dx1)+1
      N1=int(xcp/dx1)+N0
      N2=N-N0+1

      beta=(dx2/dx1)**2.d0
      ae=beta
      aw=beta
      ap=-2.d0*(1.d0+beta)

      Tre=426.136d0
C      Tre=qv*0.3d0*0.6d0/t_fl
* Defining few constants to be used frequently
      a1=(1.d0/(dx1**2.d0))+(1.d0/(dx2**2.d0))

* Defining Emissivity and Grid Area dx, dy or (dx+dy)/2
* Area of dx1 surface
      do i=2,nc(1)-1
      da(i)=dx2
      em(i)=emb
      enddo
      do i=nc(2)+1,nc(3)-1
      da(i)=dx2
      em(i)=emc
      enddo
      do i=nc(4)+1,nc(5)-1
      da(i)=dx2
      em(i)=emb
      enddo
      do i=nc(6)+1,nc(7)-1
      da(i)=dx2
      em(i)=emc
      enddo
      do i=nc(8)+1,nc(9)-1
      da(i)=dx2
      em(i)=emb
      enddo
      do i=nc(10)+1,nc(11)-1
      da(i)=dx2
      em(i)=emc
      enddo
      do i=nc(12)+1,nc(13)-1
      da(i)=dx2
      em(i)=emb
      enddo
      do i=nc(14)+1,nc(15)-1
      da(i)=dx2
      em(i)=emc
      enddo
      do i=nc(16)+1,nc(17)-1
      da(i)=dx2
      em(i)=emb
      enddo
      do i=nc(18)+1,nc(19)-1
      da(i)=dx2
      em(i)=emb
      enddo
* Area of dx2 surface
      do i=nc(1)+1,nc(2)-1
      da(i)=dx1
      em(i)=emc
      enddo
      do i=nc(3)+1,nc(4)-1
      da(i)=dx1
      em(i)=emc
      enddo
      do i=nc(5)+1,nc(6)-1
      da(i)=dx1
      em(i)=emc
      enddo
      do i=nc(7)+1,nc(8)-1
      da(i)=dx1
      em(i)=emc
      enddo
      do i=nc(9)+1,nc(10)-1
      da(i)=dx1
      em(i)=emc
      enddo
      do i=nc(11)+1,nc(12)-1
      da(i)=dx1
      em(i)=emc
      enddo
      do i=nc(13)+1,nc(14)-1
      da(i)=dx1
      em(i)=emc
      enddo
      do i=nc(15)+1,nc(16)-1
      da(i)=dx1
      em(i)=emc
      enddo
      do i=nc(17)+1,nc(18)-1
      da(i)=dx1
      em(i)=eme
      enddo
      do i=nc(19)+1,nc(20)
      da(i)=dx1
      em(i)=eme
      enddo
* Area of (dx1+dx2)/2 surface
* Emissivity at corners are different
      i=1
      da(i)=(dx1+dx2)/2.d0
      em(i)=(emb*dx2+eme*dx1)/(dx1+dx2)
      do i=1,19
      i1=nc(i)
      da(i1)=(dx1+dx2)/2.d0
      if ((i.eq.1).or.(i.eq.4).or.(i.eq.5).or.(i.eq.8))then
      i1=nc(i)
      em(i1)=(emb*dx2+emc*dx1)/(dx1+dx2)
      else if ((i.eq.9).or.(i.eq.12).or.(i.eq.13).or.(i.eq.16)) then
      i1=nc(i)
      em(i1)=(emb*dx2+emc*dx1)/(dx1+dx2)
      else if ((i.eq.2).or.(i.eq.3).or.(i.eq.6).or.(i.eq.7))then
      i1=nc(i)
      em(i1)=emc
      else if ((i.eq.10).or.(i.eq.11).or.(i.eq.14).or.(i.eq.15))then
      i1=nc(i)
      em(i1)=emc
      else if ((i.eq.17).or.(i.eq.18).or.(i.eq.19))then
      i1=nc(i)
      em(i1)=(emb*dx2+eme*dx1)/(dx1+dx2)
      endif
      enddo

      open(1,file='a_f.dat')
      do i=1,nc(20)
      do j=1,nc(20)
      read(1,10)a_f(i,j)
10    format(d16.8)
      enddo
      enddo
      close(1)
	
      call init
C      call restart
      call vel_bc
      call si_bc
      call vort_bc
      call temp_bc
      iteration=1
* iteration=1 denotes the 1st iteration
      rms=1.d0
      dt=1.d-3
* Main Loop Starts
      do while(rms.gt.1e-11)

      if(iteration.eq.1)then
      call tstep
      endif

      do i=1,N
      do j=1,M
      vort_old(i,j)=vort(i,j)
      T_old(i,j)=T(i,j)
      enddo
      enddo
	
      do i=1,N
      do j=1,M
      rhs(i,j)=T_old(i,j)
      enddo
      enddo
* Calculating the coefficient of temperature and vorticity
      call coeff_aT
      call solve_T
      call coeff_aV
      call solve_vort
      call solve_si
      call velocity
* Check RMS values
      sum1=0.d0
      sum2=0.d0
      do i=2,N-1
      do j=2,M-1
      sum1=sum1+(vort(i,j)-vort_old(i,j))**2.d0
      sum2=sum2+(T(i,j)-T_old(i,j))**2.d0
      enddo
      enddo
      rms1=dsqrt(sum1/dfloat(N*M))
      rms2=dsqrt(sum2/dfloat(N*M))
      rms=dmax1(rms1,rms2)
      rms=rms2

      write(*,*)rms1,rms2
      write(*,*)
      iteration=iteration+1
	
      enddo
* Main loop finished
      call output

      stop
      end
      
C============================================
C     Subroutine starts here      ===========
C     first subroutine call last in main code
C============================================
C============================================
      subroutine output
      include 'psemiimp.in'
	  
      double precision pro,pro1,pro2,pro3,ads

      open(3,file='leftfluid.dat')

      do i=1,nc(17)
      write(3,30)(i-1)*da(i),ph(i)
30    format(1x,f14.10,2x,d22.10)
      enddo

      close(3)

      open(4,file='rightemp.dat')

      do j=1,M
      write(4,40)(j-1)*dx2,T(N2,j)
40    format(1x,f14.10,2x,d22.10)
      enddo

      close(4)
* Write the Velocities (u & v) for all grid points
      open(5,file='vel_out.dat')
      write(5,50)
50    format(1x,'DIMENSIONLESS VELOCITY U (HORIZONTAL-X DIRECTION) AND V
     &(VERTICAL Y-DIRECTION) AT ALL CORRESPONDING GRID POINTS')
      write(5,*)
      write(5,60)
60    format(7x,'X',10x,'Y',15x,'U',8x,'V')
      write(5,*)
      do i=N0,N2
      do j=1,M
      bad=0

      if(i.le.N1)then
      do k=1,L
      if((j.gt.M1(k)).and.(j.lt.M2(k)))then
      bad=1
      endif
      enddo
      endif

      if (bad.eq.0)then
      write(5,70)(i-1)*dx1,(j-1)*dx2,u(1,i,j),u(2,i,j)
70    format (1x,f10.6,2x,f10.6,2x,f22.10,1x,f22.10)
      endif
      enddo
      write(5,*)
      enddo
      close(5)

* Write the Vorticity, Stream function
* For all Grid Points

      open(6,file='vort_si.dat')
      write(6,80)
80    format(1x,'DIMENSIONLESS VORTICITY,STREAM FUNCTION')
      write(6,*)
      write(6,90)
90    format(11x,'X',21X,'Y',25x,'VORTICITY',11x,'STREAM FUNCTION')
      write(6,*)

      do i=N0,N2
      do j=1,M
      bad=0

      if(i.lt.N1)then
      do k=1,L
      if((j.gt.M1(k)).and.(j.lt.M2(k)))then
      bad=1
      endif
      enddo
      endif

      if (bad.eq.0) then
      write(6,100)(i-1)*dx1,(j-1)*dx2,vort(i,j),si(i,j)
100   format(4x,f14.10,8x,f14.10,16x,f22.10,8x,f22.10)
      endif
      enddo
      write(6,*)
      enddo
      close(6)

* Write the Temperature for all grid points
      open(7,file='temperature.dat')
      write(7,110)
110   format(1x,'Dimensionless Temperature at grid points')
      write(7,*)
      write(7,120)
120   format(7x,'X',10x,'Y',15x,'T')
      write(7,*)
      do i=1,N
      do j=1,M
      write(7,130)(i-1)*dx1,(j-1)*dx2,T(i,j)
130   format (1x,f10.6,2x,f10.6,2x,f22.10)
      enddo
      write(7,*)
      enddo
      close(7)

* Write the Nusselt number over chip surface
      open(16,file='localnu.dat')
      write(16,111)
111   format(1x,'Local Nu number at chip surface')
      write(16,*)
      write(16,121)
121   format(7x,'X',10x,'Y',15x,'T')
      write(16,*)

      ads = 0.d0
      do k = 1,L

             pro  = 0.d0
             pro1 = 0.d0
             pro2 = 0.d0
             pro3 = 0.d0
             if (k.eq.1) then
               k1 = 1
             else 
               k1 = 2
             endif
             ads = ads + (k1-1)*(ycp2-ycp1)

      do i = N0,N1
         j = M1(k)

             avN = (T(i,j)-T(i,j-1))*dx1/(2.d0*dx2*T(i,j))
             pro1 = pro1 + avN
             if (i.eq.N0) then
               i2 = 0
             else
               i2 = 1
             endif
             ads = ads + i2*dx1
             write(16,131)ads,avN
131          format (1x,f10.6,2x,f22.10)

      enddo
      write(16,*)
             pro1 = pro1/(N1-N0+1.d0)

      do j = M1(k)+1,M2(k)
         i = N1

             avN = (T(i,j)-T(i+1,j))*dx2/(2.d0*dx1*T(i,j))
             pro2 = pro2 + avN
             ads = ads + dx2
             write(16,132)ads,avN
132          format (1x,f10.6,2x,f22.10)

      enddo
      write(16,*)
             pro2 = pro2/(M2(k)-M1(k))

      do i = N0,N1-1
         j = M2(k)

             i1=N1-i
             avN = (T(i1,j)-T(i1,j+1))*dx1/(2.d0*dx2*T(i1,j))
             pro3 = pro3 + avN
             ads = ads + dx1
             write(16,133)ads,avN
133          format (1x,f10.6,2x,f22.10)

      enddo
      write(16,*)
             pro3 = pro3/(N1-N0)
             pro = (pro1 + pro2 + pro3) / 3.d0

             write(16,134)pro
134          format(14x,f22.10)
             write(16,135)
135          format(1x,'Chip change to next chip')

      enddo

      close(16)

      open(11,file='ve.dat')
      do i=N0,N2
      do j=1,M
      bad=0

      if(i.lt.N1)then
      do k=1,L
      if((j.gt.M1(k)).and.(j.lt.M2(k)))then
      bad=1
      endif
      enddo
      endif

      if (bad.eq.0)then
      write(11,370)u(1,i,j),u(2,i,j)
370   format (2x,f22.10,1x,f22.10)
      endif
      enddo
      enddo
      close(11)

* Write the Vorticity, Stream function
* For all Grid Points

      open(12,file='vo.dat')
      do i=N0,N2
      do j=1,M
      bad=0

      if(i.lt.N1)then
      do k=1,L
      if((j.gt.M1(k)).and.(j.lt.M2(k)))then
      bad=1
      endif
      enddo
      endif

      if (bad.eq.0)then
      write(12,380)vort(i,j),si(i,j)
380   format(16x,f22.10,8x,f22.10)
      endif
      enddo
      enddo
      close(12)

* Write the Temperature for all grid points
      open(13,file='temp.dat')

      do i=1,N
      do j=1,M
      write(13,330)T(i,j)
330   format (2x,f22.10)
      enddo
      enddo
      close(13)

* Write the Temperature for all grid points
      open(14,file='leftwall.dat')

      do j=1,M
      write(14,400)(j-1)*dx2,T(N0,j)
400   format(1x,f14.10,2x,d22.10)
      enddo

      close(14)

* Write the velocity for all grid points
      open(15,file='velexit.dat')

      do i=N0,N2
      write(15,410)(i-N0)*dx1,u(2,i,M)
410   format(1x,f14.10,2x,d22.10)
      enddo

      close(15)

      return
      end
c========================================
      subroutine velocity
      include 'psemiimp.in'
* Calulating velocities (U&V) Using Stream function definition
      do i=N0+1,N2-1
      do j=2,M-1
      bad=0

      if(i.le.N1)then
      do k=1,L
      if((j.ge.M1(k)).and.(j.le.M2(k)))then
      bad=1
      endif
      enddo
      endif

      if (bad.eq.0)then
      u(1,i,j)=(si(i,j+1)-si(i,j-1))/(2.d0*dx2)
      u(2,i,j)=(si(i-1,j)-si(i+1,j))/(2.d0*dx1)
      endif 
      enddo
      enddo

      call vel_bc

      return
      end 
c========================================
      subroutine solve_si
      include 'psemiimp.in'
* Solve Poission's Equation for Stream Function for getting 
* Stream-function values (si) using Gauss-Seidel iteration 
* with Successive over relaxation

      write(*,140)
140   format(1x,'si')
************************
      rms=1.d0
      iter_s=1
      do while(rms.gt.1.e-8)
*****  EVALUATE SI VALUES AT ALL GRID POINTS EXCEPT BONDARY
      do i=N0+1,N2-1
      do j=2,M-1
      bad=0

      if(i.le.N1)then

      if((j.ge.M1(1)).and.(j.le.M2(1)))then
      bad=1
      elseif ((j.ge.M1(2)).and.(j.le.M2(2))) then
      bad=1
      elseif ((j.ge.M1(3)).and.(j.le.M2(3))) then
      bad=1
      elseif ((j.ge.M1(4)).and.(j.le.M2(4))) then
      bad=1
      endif

      endif

      if (bad.eq.0)then
      resid=-(dx2**2.d0)*vort(i,j)-ae*si(i+1,j)-aw*si(i-1,j)-
     &si(i,j+1)-si(i,j-1)-ap*si(i,j)

      si(i,j)=si(i,j)+1.d0*resid/ap
      endif
      enddo
      enddo

*****  NOW INCLUDING BOUNDARY POINTS
      call si_bc
*****  CHECK RMS VALUE

      sum=0.d0
      do i=N0+1,N2-1
      do j=2,M-1
      bad=0

      if(i.le.N1)then

      if((j.ge.M1(1)).and.(j.le.M2(1)))then
      bad=1
      elseif ((j.ge.M1(2)).and.(j.le.M2(2))) then
      bad=1
      elseif ((j.ge.M1(3)).and.(j.le.M2(3))) then
      bad=1
      elseif ((j.ge.M1(4)).and.(j.le.M2(4))) then
      bad=1
      endif

      endif

      if (bad.eq.0)then
      resid=-(dx2**2.d0)*vort(i,j)-ae*si(i+1,j)-aw*si(i-1,j)-
     &si(i,j+1)-si(i,j-1)-ap*si(i,j)
      sum=sum+resid**2.d0
      endif

      enddo
      enddo

      rms=dsqrt(sum/dfloat(N*M))
      iter_s=iter_s+1
      enddo

      write(*,150)iter_s
150   format(1x,i5)
      
      return
      end
c========================================
      subroutine solve_vort
      include 'psemiimp.in'
* Solving Vorticity-Transport Equation to get Vorticity using
* Gauss Seidel Iteration with Successive Over Relaxation
      write(*,160)
160   format(1x,'vorticity')
**************************
      do i=N0+1,N2-1
      do j=2,M-1
      bad=0

      do k=1,L
      if ((j.ge.M1(k)).and.(j.le.M2(k)).and.(i.le.N1))then
      bad=1
      endif
      enddo

      if (bad.eq.0)then
      call source(i,j,T,Sr)
      rhs(i,j)=dt*Sr*Gr/(Re**2.d0)+vort_old(i,j)
      endif
      enddo
      enddo
**************************
      rms=1.d0
      iter_v=1
      do while(rms.gt.1.e-8)
*****  EVALUATE VORTICITY AT ALL GRID PTS EXCEPT BOUNDARY
      do i=N0+1,N2-1
      do j=2,M-1
      bad=0

      do k=1,L
      if((j.ge.M1(k)).and.(j.le.M2(k)).and.(i.le.N1))then
      bad=1
      endif
      enddo

      if (bad.eq.0)then
      call conv_Q(i,j,vort,Fc)
      call diff_U(i,j,vort,Fu)
      res=rhs(i,j)+dt*(-Fc+Fu/Re)-vort(i,j)
      vort(i,j)=vort(i,j)+1.d0*res/apV(i,j)
      endif
      enddo
      enddo
*****  NOW INCLUDING BOUNDARY
      call vort_bc
*****  CHECK RMS VALUE
      sum=0.d0
      do i=N0+1,N2-1
      do j=2,M-1
      bad=0

      if(i.le.N1)then
      do k=1,L
      if((j.ge.M1(k)).and.(j.le.M2(k)))then
      bad=1
      endif
      enddo
      endif

      if (bad.eq.0)then
      call conv_Q(i,j,vort,Fc)
      call diff_U(i,j,vort,Fu)
      res=rhs(i,j)+dt*(-Fc+Fu/Re)-vort(i,j)
      sum=sum+res**2.d0
      endif
      enddo
      enddo
      rms=dsqrt(sum/dfloat(N*M))
      iter_v=iter_v+1
      enddo
      write(*,170)iter_v
170   format(1x,i5)

      return
      end
c========================================
      subroutine source(i,j,phi,Sr)
      include 'psemiimp.in'
      dimension phi(N,M)
* Subroutine to calculate Buoyancy Term
* Need to be solved in Vorticity-Transport Equation

      phix=(phi(i+1,j)-phi(i-1,j))/(2.d0*dx1)
      Sr=phix

      return
      end
c========================================
      subroutine conv_Q(i,j,phi,Fc)
      include 'psemiimp.in'
      dimension phi(N,M)

      icvt = 1
	  
      if ((i.gt.N0).and.(i.lt.N2))then
      chip=0
      do k=1,L
      if ((j.ge.M1(k)).and.(j.le.M2(k)).and.(i.le.N1))then
      chip=1
      endif
      enddo
      if (chip.eq.0)then

      if (icvt.eq.0)then
      Fc=u(1,i,j)*(phi(i+1,j)-phi(i-1,j))/(2.d0*dx1)
      Fc=Fc+(u(2,i,j)*(phi(i,j+1)-phi(i,j-1))/(2.d0*dx2))
      else if (icvt.eq.1)then

      if (u(1,i,j).ge.0.d0) then
      cux= u(1,i,j)*(phi(i,j)-phi(i-1,j))/dx1
      else
      cux= u(1,i,j)*(phi(i+1,j)-phi(i,j))/dx1
      endif

      if (u(2,i,j).ge.0.d0) then
      cvy= u(2,i,j)*(phi(i,j)-phi(i,j-1))/dx2
      else
      cvy= u(2,i,j)*(phi(i,j+1)-phi(i,j))/dx2
      endif

      Fc=cux+cvy
      endif

      else
      Fc=0.d0
      endif
      else
      Fc=0.d0
      endif

      return
      end
c========================================
      subroutine diff_U(i,j,phi,Fu)
      include 'psemiimp.in'
      dimension phi(N,M)
* Subrotine for Calculating diffusion term in
* Vorticity-Transport Equations

      Fu=(phi(i+1,j)-2.d0*phi(i,j)+phi(i-1,j))/
     &(dx1**2.d0)+
     &(phi(i,j+1)-2.d0*phi(i,j)+phi(i,j-1))/
     &(dx2**2.d0)

      return
      end
c========================================
      subroutine coeff_aV
      include 'psemiimp.in'

      icvt = 1
       
      do i=N0+1,N2-1
      do j=2,M-1

      chip=0

      do k=1,L
        if ((j.ge.M1(k)).and.(j.le.M2(k)).and.(i.le.N1))then
          chip=1
        endif
      enddo

      if (chip.eq.0)then
         if (icvt.eq.1) then
      a=(1.d0/dx1)*(dmax1(u(1,i,j),0.d0)+dmax1(-u(1,i,j),0.d0))+
     & (1.d0/dx2)*(dmax1(u(2,i,j),0.d0)+dmax1(-u(2,i,j),0.d0))
      apV(i,j)=1.d0+dt*(a+2.d0*a1/Re)
         elseif (icvt.eq.0)then
           apV(i,j)=1.d0+(2.d0*a1*dt/Re)
         endif
      endif

      enddo
      enddo

      return
      end
c========================================
      subroutine solve_T
      include 'psemiimp.in'

        double precision p1,p2,p3
          icvt = 1

      write(*,180)
180   format(1x,'temp')

* For radiation surface temp in 1-D
* rs at n level
* ph at n+1 level

      do j1=1,nc(1)
        ph(j1)=T(N0,j1)
      enddo
      do j1=nc(1)+1,nc(2)
        ph(j1)=T(N0+j1-nc(1),M1(1))
      enddo
      do j1=nc(2)+1,nc(3)
        ph(j1)=T(N1,M1(1)+j1-nc(2))
      enddo
      do j1=N0,N1-1
        ph(nc(4)+N0-j1)=T(j1,M2(1))
      enddo
      do j1=nc(4)+1,nc(5)
        ph(j1)=T(N0,M2(1)+j1-nc(4))
      enddo
      do j1=nc(5)+1,nc(6)
        ph(j1)=T(j1-nc(5)+N0,M1(2))
      enddo
      do j1=nc(6)+1,nc(7)
        ph(j1)=T(N1,M1(2)+j1-nc(6))
      enddo
      do j1=N0,N1-1
        ph(nc(8)+N0-j1)=T(j1,M2(2))
      enddo
      do j1=nc(8)+1,nc(9)
        ph(j1)=T(N0,M2(2)+j1-nc(8))
      enddo
      do j1=nc(9)+1,nc(10)
        ph(j1)=T(N0+j1-nc(9),M1(3))
      enddo
      do j1=nc(10)+1,nc(11)
        ph(j1)=T(N1,M1(3)+j1-nc(10))
      enddo
      do j1=N0,N1-1
        ph(nc(12)+N0-j1)=T(j1,M2(3))  
      enddo
      do j1=nc(12)+1,nc(13)
        ph(j1)=T(N0,M2(3)+j1-nc(12))
      enddo
      do j1=nc(13)+1,nc(14)
        ph(j1)=T(N0+j1-nc(13),M1(4))
      enddo
      do j1=nc(14)+1,nc(15)
        ph(j1)=T(N1,M1(4)+j1-nc(14))
      enddo
      do j1=N0,N1-1
        ph(nc(16)+N0-j1)=T(j1,M2(4))  
      enddo
      do j1=nc(16)+1,nc(17)
        ph(j1)=T(N0,M2(4)+j1-nc(16))
      enddo
      do j1=nc(17)+1,nc(18)
        ph(j1)=T(N0+j1-nc(17),M)
      enddo
      do j1=1,M-1
        ph(nc(19)+1-j1)=T(N2,j1)
      enddo
      do j1=N0+1,N2-1
        ph(nc(20)+N0+1-j1)=T(j1,1)
      enddo

      rms=1.d0
      iter_t = 0
      do while(rms.gt.1.e-11)
*****  EVALUATE TEMP AT ALL GRID PTS EXCEPT BOUNDARY
      do j=2,M-1
      do i=2,N-1

* Calculating in fluid and solids (wall and chip)

* Inside Wall
      if((i.lt.N0).or.(i.gt.N2))then
        p1 =(T(i+1,j)-2.d0*T(i,j)+T(i-1,j))/(dx1**2.d0)
        p2 =(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3 = p1+p2
        Fd = p3*t_wa*rcf/(t_fl*Re*Pr*rcw)
      endif

* Inside Chip
      if ((i.gt.N0).and.(i.lt.N1))then
        do k=1,L
         if ((j.gt.M1(k)).and.(j.lt.M2(k)))then
           p1 =(T(i+1,j)-2.d0*T(i,j)+T(i-1,j))/(dx1**2.d0)
           p2 =(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
           p3 =p1+p2
* Channel width is always 1 so not mention
           Fd = p3*t_ch*rcf/(t_fl*Re*Pr*rcc)
           Fd = Fd+(rcf/(rcc*Re*Pr))*(1.d0/(xcp*(ycp2-ycp1)))
         endif
        enddo
      endif

* Inside Fluid
      if ((i.gt.N0).and.(i.lt.N2))then
       chip=0
        do k=1,L
         if ((j.ge.M1(k)).and.(j.le.M2(k)).and.(i.le.N1))then
           chip=1
         endif
        enddo

       if (chip.eq.0)then
        p1 =(T(i+1,j)-2.d0*T(i,j)+T(i-1,j))/(dx1**2.d0)
        p2 =(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3 =p1+p2
        Fd = p3/(Re*Pr)

         if (icvt.eq.0)then
           p1 =u(1,i,j)*(T(i+1,j)-T(i-1,j))/(2.d0*dx1)
           p2 =u(2,i,j)*(T(i,j+1)-T(i,j-1))/(2.d0*dx2)
           Fc =p1+p2
         else if (icvt.eq.1)then
           if (u(1,i,j).ge.0.d0) then
             Fc1= u(1,i,j)*(T(i,j)-T(i-1,j))/dx1
           else
             Fc1= u(1,i,j)*(T(i+1,j)-T(i,j))/dx1
           endif

           if (u(2,i,j).ge.0.d0) then
             Fc2= u(2,i,j)*(T(i,j)-T(i,j-1))/dx2
           else
             Fc2= u(2,i,j)*(T(i,j+1)-T(i,j))/dx2
           endif

         Fc=Fc1+Fc2
         endif
        Fd=Fd-Fc
       endif
      endif

* Chip and Left Wall Interface
      if (i.eq.N0)then
       if ((j.ge.M1(1)).and.(j.le.M2(1)))then
        p1 =(t_ch*T(i+1,j)-(t_ch+t_wa)*T(i,j)+t_wa*T(i-1,j))/(dx1**2.d0)
        p2 =(t_ch+t_wa)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3 =((2.d0*p1)+p2+t_fl/(xcp*(ycp2-ycp1)))/(Re*Pr*t_fl)
        Fd = p3/((rcw/rcf)+(rcc/rcf))
       elseif ((j.ge.M1(2)).and.(j.le.M2(2)))then
        p1 =(t_ch*T(i+1,j)-(t_ch+t_wa)*T(i,j)+t_wa*T(i-1,j))/(dx1**2.d0)
        p2 =(t_ch+t_wa)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3 =((2.d0*p1)+p2+t_fl/(xcp*(ycp2-ycp1)))/(Re*Pr*t_fl)
        Fd = p3/((rcw/rcf)+(rcc/rcf))
       elseif ((j.ge.M1(3)).and.(j.le.M2(3)))then
        p1 =(t_ch*T(i+1,j)-(t_ch+t_wa)*T(i,j)+t_wa*T(i-1,j))/(dx1**2.d0)
        p2 =(t_ch+t_wa)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3 =((2.d0*p1)+p2+t_fl/(xcp*(ycp2-ycp1)))/(Re*Pr*t_fl)
        Fd = p3/((rcw/rcf)+(rcc/rcf))
       elseif ((j.ge.M1(4)).and.(j.le.M2(4)))then
        p1 =(t_ch*T(i+1,j)-(t_ch+t_wa)*T(i,j)+t_wa*T(i-1,j))/(dx1**2.d0)
        p2 =(t_ch+t_wa)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3 =((2.d0*p1)+p2+t_fl/(xcp*(ycp2-ycp1)))/(Re*Pr*t_fl)
        Fd = p3/((rcw/rcf)+(rcc/rcf))
* Left wall side emitting radiation
       elseif (j.lt.M1(1))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1= (rad1*T(i,j)+(Tin**4.d0))*sig*emb*dx2
        j1=j
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2= (rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad= rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1 =(t_fl*T(i+1,j)-(t_fl+t_wa)*T(i,j)+t_wa*T(i-1,j))/(dx1**2.d0)
        p2 =(t_fl+t_wa)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3 =(2.d0*p1)+p2
        Fd =(p3-rad)/(t_fl*Re*Pr*(1.d0+(rcw/rcf)))

       elseif ((j.gt.M2(1)).and.(j.lt.M1(2)))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1 =(rad1*T(i,j)+(Tin**4.d0))*sig*emb*dx2
        j1=nc(4)+j-M2(1)
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2= (rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad= rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1 =(t_fl*T(i+1,j)-(t_fl+t_wa)*T(i,j)+t_wa*T(i-1,j))/(dx1**2.d0)
        p2 =(t_fl+t_wa)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3 =(2.d0*p1)+p2
        Fd =(p3-rad)/(t_fl*Re*Pr*(1.d0+(rcw/rcf)))

       elseif ((j.gt.M2(2)).and.(j.lt.M1(3)))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1 =(rad1*T(i,j)+(Tin**4.d0))*sig*emb*dx2
        j1=nc(8)+j-M2(2)
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2= (rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad= rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1 =(t_fl*T(i+1,j)-(t_fl+t_wa)*T(i,j)+t_wa*T(i-1,j))/(dx1**2.d0)
        p2 =(t_fl+t_wa)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3 =(2.d0*p1)+p2
        Fd =(p3-rad)/(t_fl*Re*Pr*(1.d0+(rcw/rcf)))

       elseif ((j.gt.M2(3)).and.(j.lt.M1(4)))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1=(rad1*T(i,j)+(Tin**4.d0))*sig*emb*dx2
        j1=nc(12)+j-M2(3)
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2= (rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad= rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1 =(t_fl*T(i+1,j)-(t_fl+t_wa)*T(i,j)+t_wa*T(i-1,j))/(dx1**2.d0)
        p2 =(t_fl+t_wa)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3 =(2.d0*p1)+p2
        Fd =(p3-rad)/(t_fl*Re*Pr*(1.d0+(rcw/rcf)))

       elseif (j.gt.M2(4))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1=(rad1*T(i,j)+(Tin**4.d0))*sig*emb*dx2
        j1=nc(16)+j-M2(4)
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2= (rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad= rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1 =(t_fl*T(i+1,j)-(t_fl+t_wa)*T(i,j)+t_wa*T(i-1,j))/(dx1**2.d0)
        p2 =(t_fl+t_wa)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3 =(2.d0*p1)+p2
        Fd =(p3-rad)/(t_fl*Re*Pr*(1.d0+(rcw/rcf)))
       endif
*ph at n+1 level i,j need to define i1
*rs at n level j1 starts from 1 to nc(12)
      endif
* Right Wall
      if (i.eq.N2) then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1 =(rad1*T(i,j)+(Tin**4.d0))*sig*emb*dx2
        j1=nc(19)-j+1
        rad=0.d0
         do i1=1,nc(20)
           rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &          (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &          (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &          (4.d0*Tre*(Tin**3.d0))
           rad2= (rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
           rad= rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1 =(t_fl*T(i-1,j)-(t_fl+t_wa)*T(i,j)+t_wa*T(i+1,j))/(dx1**2.d0)
        p2 =(t_fl+t_wa)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3 =(2.d0*p1)+p2
        Fd= (p3-rad)/(t_fl*Re*Pr*(1.d0+(rcw/rcf)))
      endif
* Top face of chip
      if (j.eq.M2(1)) then
       if ((i.gt.N0).and.(i.le.N1))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1 =(rad1*T(i,j)+(Tin**4.d0))*sig*emc*dx1
        j1=nc(4)-i+N0
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2 =(rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad =rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1=(t_fl*T(i,j+1)-(t_fl+t_ch)*T(i,j)+t_ch*T(i,j-1))/(dx2**2.d0)
        p2=(t_fl+t_ch)*(T(i+1,j)-2.d0*T(i,j)+T(i-1,j))/(dx1**2.d0)
        p3=(2.d0*p1)+p2-rad+t_fl/(xcp*(ycp2-ycp1))
        Fd =p3/(t_fl*Re*Pr*(1.d0+(rcc/rcf)))
       endif
      endif

      if (j.eq.M2(2)) then
       if ((i.gt.N0).and.(i.le.N1))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1= (rad1*T(i,j)+(Tin**4.d0))*sig*emc*dx1
        j1=nc(8)-i+N0
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2 =(rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad =rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1=(t_fl*T(i,j+1)-(t_fl+t_ch)*T(i,j)+t_ch*T(i,j-1))/(dx2**2.d0)
        p2=(t_fl+t_ch)*(T(i+1,j)-2.d0*T(i,j)+T(i-1,j))/(dx1**2.d0)
        p3=(2.d0*p1)+p2-rad+t_fl/(xcp*(ycp2-ycp1))
        Fd =p3/(t_fl*Re*Pr*(1.d0+(rcc/rcf)))
       endif
      endif

      if (j.eq.M2(3)) then
       if ((i.gt.N0).and.(i.le.N1))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1= (rad1*T(i,j)+(Tin**4.d0))*sig*emc*dx1
        j1=nc(12)-i+N0
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2 =(rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad =rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1=(t_fl*T(i,j+1)-(t_fl+t_ch)*T(i,j)+t_ch*T(i,j-1))/(dx2**2.d0)
        p2=(t_fl+t_ch)*(T(i+1,j)-2.d0*T(i,j)+T(i-1,j))/(dx1**2.d0)
        p3=(2.d0*p1)+p2-rad+t_fl/(xcp*(ycp2-ycp1))
        Fd =p3/(t_fl*Re*Pr*(1.d0+(rcc/rcf)))
       endif
      endif

      if (j.eq.M2(4)) then
       if ((i.gt.N0).and.(i.le.N1))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1= (rad1*T(i,j)+(Tin**4.d0))*sig*emc*dx1
        j1=nc(16)-i+N0
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2 =(rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad =rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1=(t_fl*T(i,j+1)-(t_fl+t_ch)*T(i,j)+t_ch*T(i,j-1))/(dx2**2.d0)
        p2=(t_fl+t_ch)*(T(i+1,j)-2.d0*T(i,j)+T(i-1,j))/(dx1**2.d0)
        p3=(2.d0*p1)+p2-rad+t_fl/(xcp*(ycp2-ycp1))
        Fd =p3/(t_fl*Re*Pr*(1.d0+(rcc/rcf))) 
       endif
      endif

* Bottom face of chip
      if (j.eq.M1(1)) then
       if ((i.gt.N0).and.(i.le.N1))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1 =(rad1*T(i,j)+(Tin**4.d0))*sig*emc*dx1
        j1=nc(1)+i-N0
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2 =(rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad =rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1=(t_ch*T(i,j+1)-(t_fl+t_ch)*T(i,j)+t_fl*T(i,j-1))/(dx2**2.d0)
        p2=(t_fl+t_ch)*(T(i+1,j)-2.d0*T(i,j)+T(i-1,j))/(dx1**2.d0)
        p3=(2.d0*p1)+p2-rad+t_fl/(xcp*(ycp2-ycp1))
        Fd =p3/(t_fl*Re*Pr*(1.d0+(rcc/rcf)))
       endif
      endif

      if (j.eq.M1(2)) then
       if ((i.gt.N0).and.(i.le.N1))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1 =(rad1*T(i,j)+(Tin**4.d0))*sig*emc*dx1
        j1=nc(5)+i-N0
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2 =(rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad =rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1=(t_ch*T(i,j+1)-(t_fl+t_ch)*T(i,j)+t_fl*T(i,j-1))/(dx2**2.d0)
        p2=(t_fl+t_ch)*(T(i+1,j)-2.d0*T(i,j)+T(i-1,j))/(dx1**2.d0)
        p3=(2.d0*p1)+p2-rad+t_fl/(xcp*(ycp2-ycp1))
        Fd =p3/(t_fl*Re*Pr*(1.d0+(rcc/rcf)))
       endif
      endif

      if (j.eq.M1(3)) then
       if ((i.gt.N0).and.(i.le.N1))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1 =(rad1*T(i,j)+(Tin**4.d0))*sig*emc*dx1
        j1=nc(9)+i-N0
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2 =(rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad =rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1=(t_ch*T(i,j+1)-(t_fl+t_ch)*T(i,j)+t_fl*T(i,j-1))/(dx2**2.d0)
        p2=(t_fl+t_ch)*(T(i+1,j)-2.d0*T(i,j)+T(i-1,j))/(dx1**2.d0)
        p3=(2.d0*p1)+p2-rad+t_fl/(xcp*(ycp2-ycp1))
        Fd =p3/(t_fl*Re*Pr*(1.d0+(rcc/rcf)))
       endif
      endif

      if (j.eq.M1(4)) then
       if ((i.gt.N0).and.(i.le.N1))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1 =(rad1*T(i,j)+(Tin**4.d0))*sig*emc*dx1
        j1=nc(13)+i-N0
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2 =(rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad =rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1=(t_ch*T(i,j+1)-(t_fl+t_ch)*T(i,j)+t_fl*T(i,j-1))/(dx2**2.d0)
        p2=(t_fl+t_ch)*(T(i+1,j)-2.d0*T(i,j)+T(i-1,j))/(dx1**2.d0)
        p3=(2.d0*p1)+p2-rad+t_fl/(xcp*(ycp2-ycp1))
        Fd =p3/(t_fl*Re*Pr*(1.d0+(rcc/rcf)))
       endif
      endif
* Front face of chip
      if (i.eq.N1) then
       if ((j.gt.M1(1)).and.(j.lt.M2(1)))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1 =(rad1*T(i,j)+(Tin**4.d0))*sig*emc*dx2
        j1=nc(2)+j-M1(1)
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2 =(rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad =rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1=(t_ch*T(i-1,j)-(t_fl+t_ch)*T(i,j)+t_fl*T(i+1,j))/(dx1**2.d0)
        p2=(t_fl+t_ch)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3=(2.d0*p1)+p2-rad+t_fl/(xcp*(ycp2-ycp1))
        Fd =p3/(t_fl*Re*Pr*(1.d0+(rcc/rcf)))
       endif
      endif

      if (i.eq.N1) then
       if ((j.gt.M1(2)).and.(j.lt.M2(2)))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1 =(rad1*T(i,j)+(Tin**4.d0))*sig*emc*dx2
        j1=nc(6)+j-M1(2)
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2 =(rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad =rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1=(t_ch*T(i-1,j)-(t_fl+t_ch)*T(i,j)+t_fl*T(i+1,j))/(dx1**2.d0)
        p2=(t_fl+t_ch)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3=(2.d0*p1)+p2-rad+t_fl/(xcp*(ycp2-ycp1))
        Fd =p3/(t_fl*Re*Pr*(1.d0+(rcc/rcf)))
       endif
      endif

      if (i.eq.N1) then
       if ((j.gt.M1(3)).and.(j.lt.M2(3)))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1 =(rad1*T(i,j)+(Tin**4.d0))*sig*emc*dx2
        j1=nc(10)+j-M1(3)
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2 =(rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad =rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1=(t_ch*T(i-1,j)-(t_fl+t_ch)*T(i,j)+t_fl*T(i+1,j))/(dx1**2.d0)
        p2=(t_fl+t_ch)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3=(2.d0*p1)+p2-rad+t_fl/(xcp*(ycp2-ycp1))
        Fd =p3/(t_fl*Re*Pr*(1.d0+(rcc/rcf)))
       endif
      endif

      if (i.eq.N1) then
       if ((j.gt.M1(4)).and.(j.lt.M2(4)))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1 =(rad1*T(i,j)+(Tin**4.d0))*sig*emc*dx2
        j1=nc(14)+j-M1(4)
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2 =(rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad =rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1=(t_ch*T(i-1,j)-(t_fl+t_ch)*T(i,j)+t_fl*T(i+1,j))/(dx1**2.d0)
        p2=(t_fl+t_ch)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3=(2.d0*p1)+p2-rad+t_fl/(xcp*(ycp2-ycp1))
        Fd =p3/(t_fl*Re*Pr*(1.d0+(rcc/rcf)))
       endif
      endif

      res=rhs(i,j)+dt*Fd-T(i,j)
      T(i,j)=T(i,j)+1.d0*res/apT(i,j)

* For radiation surface temp in 1-D
* rs at n level
* ph at n+1 level
      if (i.eq.N0)then
       if (j.le.M1(1))then
        ph(j)=T(i,j)
       else if ((j.ge.M2(1)).and.(j.le.M1(2)))then
        j1=nc(4)+j-M2(1)
        ph(j1)=T(i,j)
       else if ((j.ge.M2(2)).and.(j.le.M1(3)))then
        j1=nc(8)+j-M2(2)
        ph(j1)=T(i,j)
       else if ((j.ge.M2(3)).and.(j.le.M1(4)))then
        j1=nc(12)+j-M2(3)
        ph(j1)=T(i,j)
       else if (j.ge.M2(4)) then
        j1=nc(16)+j-M2(4)
        ph(j1)=T(i,j)
       endif
      endif

      if (j.eq.M1(1)) then
       if ((i.gt.N0).and.(i.le.N1)) then
        j1=nc(1)+i-N0
        ph(j1)=T(i,j)
       endif
      endif
      if (j.eq.M1(2)) then
       if ((i.gt.N0).and.(i.le.N1)) then
        j1=nc(5)+i-N0
        ph(j1)=T(i,j)
       endif
      endif
      if (j.eq.M1(3)) then
       if ((i.gt.N0).and.(i.le.N1)) then
        j1=nc(9)+i-N0
        ph(j1)=T(i,j)
       endif
      endif
      if (j.eq.M1(4)) then
       if ((i.gt.N0).and.(i.le.N1)) then
        j1=nc(13)+i-N0
        ph(j1)=T(i,j)
       endif
      endif

      if (j.eq.M2(1)) then
       if ((i.gt.N0).and.(i.le.N1)) then
        j1=nc(4)-i+N0
        ph(j1)=T(i,j)
       endif
      endif
      if (j.eq.M2(2)) then
       if ((i.gt.N0).and.(i.le.N1)) then
        j1=nc(8)-i+N0
        ph(j1)=T(i,j)
       endif
      endif
      if (j.eq.M2(3)) then
       if ((i.gt.N0).and.(i.le.N1)) then
        j1=nc(12)-i+N0
        ph(j1)=T(i,j)
       endif
      endif
      if (j.eq.M2(4)) then
       if ((i.gt.N0).and.(i.le.N1)) then
        j1=nc(16)-i+N0
        ph(j1)=T(i,j)
       endif
      endif

      if (i.eq.N1) then
       if ((j.gt.M1(1)).and.(j.lt.M2(1))) then
        j1=nc(2)+j-M1(1)
        ph(j1)=T(i,j)
       endif
      endif
      if (i.eq.N1) then
       if ((j.gt.M1(2)).and.(j.lt.M2(2))) then
        j1=nc(6)+j-M1(2)
        ph(j1)=T(i,j)
       endif
      endif
      if (i.eq.N1) then
       if ((j.gt.M1(3)).and.(j.lt.M2(3))) then
        j1=nc(10)+j-M1(3)
        ph(j1)=T(i,j)
       endif
      endif
      if (i.eq.N1) then
       if ((j.gt.M1(4)).and.(j.lt.M2(4))) then
        j1=nc(14)+j-M1(4)
        ph(j1)=T(i,j)
       endif
      endif

      if (i.eq.N2)then
       j1=nc(19)-j+1
       ph(j1)=T(i,j)
      endif

      enddo
      enddo
* Substituting BC's
      call temp_bc
        do i=N0,N2
          j1=nc(17)+i-N0
          ph(j1)=T(i,M)
        enddo
        do i=N0+1,N2
          j1=nc(20)+N0+1-i
          ph(j1)=T(i,1)
        enddo
      ph(1)=T(N0,1)

* RMS check will now
      sum=0.d0
      do j=2,M-1
        do i=2,N-1

* Calculating in fluid and solids (wall and chip)
* Inside Wall
        if((i.lt.N0).or.(i.gt.N2))then
          p1 =(T(i+1,j)-2.d0*T(i,j)+T(i-1,j))/(dx1**2.d0)
          p2 =(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
          p3 =p1+p2
          Fd = p3*t_wa*rcf/(t_fl*Re*Pr*rcw)
        endif

* Inside Chip
        if ((i.gt.N0).and.(i.lt.N1))then
         do k=1,L
          if ((j.gt.M1(k)).and.(j.lt.M2(k)))then
           p1=(T(i+1,j)-2.d0*T(i,j)+T(i-1,j))/(dx1**2.d0)
           p2=(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
           p3=p1+p2
* Channel width is always 1 so not mention
           Fd = p3*t_ch*rcf/(t_fl*Re*Pr*rcc)
           Fd = Fd+(rcf/(rcc*Re*Pr))*(1.d0/(xcp*(ycp2-ycp1)))
          endif
         enddo
        endif

* Inside Fluid
      if ((i.gt.N0).and.(i.lt.N2))then
         chip=0
         do k=1,L
          if ((j.ge.M1(k)).and.(j.le.M2(k)).and.(i.le.N1))then
           chip=1
          endif
         enddo

       if (chip.eq.0)then
        p1 =(T(i+1,j)-2.d0*T(i,j)+T(i-1,j))/(dx1**2.d0)
        p2 =(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3 =p1+p2
        Fd = p3/(Re*Pr)

         if (icvt.eq.0)then
           p1 =u(1,i,j)*(T(i+1,j)-T(i-1,j))/(2.d0*dx1)
           p2 =u(2,i,j)*(T(i,j+1)-T(i,j-1))/(2.d0*dx2)
           Fc =p1+p2
         else if (icvt.eq.1)then
           if (u(1,i,j).ge.0.d0) then
             Fc1= u(1,i,j)*(T(i,j)-T(i-1,j))/dx1
           else
             Fc1= u(1,i,j)*(T(i+1,j)-T(i,j))/dx1
           endif

           if (u(2,i,j).ge.0.d0) then
             Fc2= u(2,i,j)*(T(i,j)-T(i,j-1))/dx2
           else
             Fc2= u(2,i,j)*(T(i,j+1)-T(i,j))/dx2
           endif

         Fc=Fc1+Fc2
         endif
        Fd=Fd-Fc
       endif
      endif

* Chip and Left Wall Interface
      if (i.eq.N0)then
       if ((j.ge.M1(1)).and.(j.le.M2(1)))then
        p1 =(t_ch*T(i+1,j)-(t_ch+t_wa)*T(i,j)+t_wa*T(i-1,j))/(dx1**2.d0)
        p2 =(t_ch+t_wa)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3 =((2.d0*p1)+p2+t_fl/(xcp*(ycp2-ycp1)))/(Re*Pr*t_fl)
        Fd = p3/((rcw/rcf)+(rcc/rcf))
       elseif ((j.ge.M1(2)).and.(j.le.M2(2)))then
        p1 =(t_ch*T(i+1,j)-(t_ch+t_wa)*T(i,j)+t_wa*T(i-1,j))/(dx1**2.d0)
        p2 =(t_ch+t_wa)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3 =((2.d0*p1)+p2+t_fl/(xcp*(ycp2-ycp1)))/(Re*Pr*t_fl)
        Fd = p3/((rcw/rcf)+(rcc/rcf))
       elseif ((j.ge.M1(3)).and.(j.le.M2(3)))then
        p1 =(t_ch*T(i+1,j)-(t_ch+t_wa)*T(i,j)+t_wa*T(i-1,j))/(dx1**2.d0)
        p2 =(t_ch+t_wa)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3 =((2.d0*p1)+p2+t_fl/(xcp*(ycp2-ycp1)))/(Re*Pr*t_fl)
        Fd = p3/((rcw/rcf)+(rcc/rcf))
       elseif ((j.ge.M1(4)).and.(j.le.M2(4)))then
        p1 =(t_ch*T(i+1,j)-(t_ch+t_wa)*T(i,j)+t_wa*T(i-1,j))/(dx1**2.d0)
        p2 =(t_ch+t_wa)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3 =((2.d0*p1)+p2+t_fl/(xcp*(ycp2-ycp1)))/(Re*Pr*t_fl)
        Fd = p3/((rcw/rcf)+(rcc/rcf))
* Left wall side emitting radiation 
        elseif (j.lt.M1(1))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1= (rad1*T(i,j)+(Tin**4.d0))*sig*emb*dx2
        j1=j
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2= (rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad= rad+rad2
          enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1 =(t_fl*T(i+1,j)-(t_fl+t_wa)*T(i,j)+t_wa*T(i-1,j))/(dx1**2.d0)
        p2 =(t_fl+t_wa)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3 =(2.d0*p1)+p2
        Fd =(p3-rad)/(t_fl*Re*Pr*(1.d0+(rcw/rcf)))

       elseif ((j.gt.M2(1)).and.(j.lt.M1(2)))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1 =(rad1*T(i,j)+(Tin**4.d0))*sig*emb*dx2
        j1=nc(4)+j-M2(1)
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2= (rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad= rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1 =(t_fl*T(i+1,j)-(t_fl+t_wa)*T(i,j)+t_wa*T(i-1,j))/(dx1**2.d0)
        p2 =(t_fl+t_wa)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3 =(2.d0*p1)+p2
        Fd =(p3-rad)/(t_fl*Re*Pr*(1.d0+(rcw/rcf)))

       elseif ((j.gt.M2(2)).and.(j.lt.M1(3)))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1 =(rad1*T(i,j)+(Tin**4.d0))*sig*emb*dx2
        j1=nc(8)+j-M2(2)
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2= (rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad= rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1 =(t_fl*T(i+1,j)-(t_fl+t_wa)*T(i,j)+t_wa*T(i-1,j))/(dx1**2.d0)
        p2 =(t_fl+t_wa)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3 =(2.d0*p1)+p2
        Fd =(p3-rad)/(t_fl*Re*Pr*(1.d0+(rcw/rcf)))

       elseif ((j.gt.M2(3)).and.(j.lt.M1(4)))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1=(rad1*T(i,j)+(Tin**4.d0))*sig*emb*dx2
        j1=nc(12)+j-M2(3)
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2= (rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad= rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1 =(t_fl*T(i+1,j)-(t_fl+t_wa)*T(i,j)+t_wa*T(i-1,j))/(dx1**2.d0)
        p2 =(t_fl+t_wa)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3 =(2.d0*p1)+p2
        Fd =(p3-rad)/(t_fl*Re*Pr*(1.d0+(rcw/rcf)))

       elseif (j.gt.M2(4))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1=(rad1*T(i,j)+(Tin**4.d0))*sig*emb*dx2
        j1=nc(16)+j-M2(4)
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2= (rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad= rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1 =(t_fl*T(i+1,j)-(t_fl+t_wa)*T(i,j)+t_wa*T(i-1,j))/(dx1**2.d0)
        p2 =(t_fl+t_wa)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3 =(2.d0*p1)+p2
        Fd =(p3-rad)/(t_fl*Re*Pr*(1.d0+(rcw/rcf)))
       endif
*ph at n+1 level i,j need to define i1
*rs at n level j1 starts from 1 to nc(12)
      endif
* Right Wall
      if (i.eq.N2) then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1 =(rad1*T(i,j)+(Tin**4.d0))*sig*emb*dx2
        j1=nc(19)-j+1
        rad=0.d0
         do i1=1,nc(20)
           rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &          (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &          (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &          (4.d0*Tre*(Tin**3.d0))
           rad2= (rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
           rad= rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1 =(t_fl*T(i-1,j)-(t_fl+t_wa)*T(i,j)+t_wa*T(i+1,j))/(dx1**2.d0)
        p2 =(t_fl+t_wa)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3 =(2.d0*p1)+p2
        Fd= (p3-rad)/(t_fl*Re*Pr*(1.d0+(rcw/rcf)))
      endif
* Top face of chip
      if (j.eq.M2(1)) then
       if ((i.gt.N0).and.(i.le.N1))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1 =(rad1*T(i,j)+(Tin**4.d0))*sig*emc*dx1
        j1=nc(4)-i+N0
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2 =(rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad =rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1=(t_fl*T(i,j+1)-(t_fl+t_ch)*T(i,j)+t_ch*T(i,j-1))/(dx2**2.d0)
        p2=(t_fl+t_ch)*(T(i+1,j)-2.d0*T(i,j)+T(i-1,j))/(dx1**2.d0)
        p3=(2.d0*p1)+p2-rad+t_fl/(xcp*(ycp2-ycp1))
        Fd =p3/(t_fl*Re*Pr*(1.d0+(rcc/rcf)))
       endif
      endif

      if (j.eq.M2(2)) then
       if ((i.gt.N0).and.(i.le.N1))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1= (rad1*T(i,j)+(Tin**4.d0))*sig*emc*dx1
        j1=nc(8)-i+N0
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2 =(rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad =rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1=(t_fl*T(i,j+1)-(t_fl+t_ch)*T(i,j)+t_ch*T(i,j-1))/(dx2**2.d0)
        p2=(t_fl+t_ch)*(T(i+1,j)-2.d0*T(i,j)+T(i-1,j))/(dx1**2.d0)
        p3=(2.d0*p1)+p2-rad+t_fl/(xcp*(ycp2-ycp1))
        Fd =p3/(t_fl*Re*Pr*(1.d0+(rcc/rcf)))
       endif
      endif

      if (j.eq.M2(3)) then
       if ((i.gt.N0).and.(i.le.N1))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1= (rad1*T(i,j)+(Tin**4.d0))*sig*emc*dx1
        j1=nc(12)-i+N0
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2 =(rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad =rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1=(t_fl*T(i,j+1)-(t_fl+t_ch)*T(i,j)+t_ch*T(i,j-1))/(dx2**2.d0)
        p2=(t_fl+t_ch)*(T(i+1,j)-2.d0*T(i,j)+T(i-1,j))/(dx1**2.d0)
        p3=(2.d0*p1)+p2-rad+t_fl/(xcp*(ycp2-ycp1))
        Fd =p3/(t_fl*Re*Pr*(1.d0+(rcc/rcf)))
       endif
      endif

      if (j.eq.M2(4)) then
       if ((i.gt.N0).and.(i.le.N1))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1= (rad1*T(i,j)+(Tin**4.d0))*sig*emc*dx1
        j1=nc(16)-i+N0
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2 =(rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad =rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1=(t_fl*T(i,j+1)-(t_fl+t_ch)*T(i,j)+t_ch*T(i,j-1))/(dx2**2.d0)
        p2=(t_fl+t_ch)*(T(i+1,j)-2.d0*T(i,j)+T(i-1,j))/(dx1**2.d0)
        p3=(2.d0*p1)+p2-rad+t_fl/(xcp*(ycp2-ycp1))
        Fd =p3/(t_fl*Re*Pr*(1.d0+(rcc/rcf))) 
       endif
      endif

* Bottom face of chip
      if (j.eq.M1(1)) then
       if ((i.gt.N0).and.(i.le.N1))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1 =(rad1*T(i,j)+(Tin**4.d0))*sig*emc*dx1
        j1=nc(1)+i-N0
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2 =(rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad =rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1=(t_ch*T(i,j+1)-(t_fl+t_ch)*T(i,j)+t_fl*T(i,j-1))/(dx2**2.d0)
        p2=(t_fl+t_ch)*(T(i+1,j)-2.d0*T(i,j)+T(i-1,j))/(dx1**2.d0)
        p3=(2.d0*p1)+p2-rad+t_fl/(xcp*(ycp2-ycp1))
        Fd =p3/(t_fl*Re*Pr*(1.d0+(rcc/rcf)))
       endif
      endif

      if (j.eq.M1(2)) then
       if ((i.gt.N0).and.(i.le.N1))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1 =(rad1*T(i,j)+(Tin**4.d0))*sig*emc*dx1
        j1=nc(5)+i-N0
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2 =(rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad =rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1=(t_ch*T(i,j+1)-(t_fl+t_ch)*T(i,j)+t_fl*T(i,j-1))/(dx2**2.d0)
        p2=(t_fl+t_ch)*(T(i+1,j)-2.d0*T(i,j)+T(i-1,j))/(dx1**2.d0)
        p3=(2.d0*p1)+p2-rad+t_fl/(xcp*(ycp2-ycp1))
        Fd =p3/(t_fl*Re*Pr*(1.d0+(rcc/rcf)))
       endif
      endif

      if (j.eq.M1(3)) then
       if ((i.gt.N0).and.(i.le.N1))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1 =(rad1*T(i,j)+(Tin**4.d0))*sig*emc*dx1
        j1=nc(9)+i-N0
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2 =(rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad =rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1=(t_ch*T(i,j+1)-(t_fl+t_ch)*T(i,j)+t_fl*T(i,j-1))/(dx2**2.d0)
        p2=(t_fl+t_ch)*(T(i+1,j)-2.d0*T(i,j)+T(i-1,j))/(dx1**2.d0)
        p3=(2.d0*p1)+p2-rad+t_fl/(xcp*(ycp2-ycp1))
        Fd =p3/(t_fl*Re*Pr*(1.d0+(rcc/rcf)))
       endif
      endif
	  
      if (j.eq.M1(4)) then
       if ((i.gt.N0).and.(i.le.N1))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1 =(rad1*T(i,j)+(Tin**4.d0))*sig*emc*dx1
        j1=nc(13)+i-N0
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2 =(rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad =rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1=(t_ch*T(i,j+1)-(t_fl+t_ch)*T(i,j)+t_fl*T(i,j-1))/(dx2**2.d0)
        p2=(t_fl+t_ch)*(T(i+1,j)-2.d0*T(i,j)+T(i-1,j))/(dx1**2.d0)
        p3=(2.d0*p1)+p2-rad+t_fl/(xcp*(ycp2-ycp1))
        Fd =p3/(t_fl*Re*Pr*(1.d0+(rcc/rcf)))
       endif
      endif
* Front face of chip
      if (i.eq.N1) then
       if ((j.gt.M1(1)).and.(j.lt.M2(1)))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1 =(rad1*T(i,j)+(Tin**4.d0))*sig*emc*dx2
        j1=nc(2)+j-M1(1)
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2 =(rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad =rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1=(t_ch*T(i-1,j)-(t_fl+t_ch)*T(i,j)+t_fl*T(i+1,j))/(dx1**2.d0)
        p2=(t_fl+t_ch)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3=(2.d0*p1)+p2-rad+t_fl/(xcp*(ycp2-ycp1))
        Fd =p3/(t_fl*Re*Pr*(1.d0+(rcc/rcf)))
       endif
      endif

      if (i.eq.N1) then
       if ((j.gt.M1(2)).and.(j.lt.M2(2)))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1 =(rad1*T(i,j)+(Tin**4.d0))*sig*emc*dx2
        j1=nc(6)+j-M1(2)
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2 =(rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad =rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1=(t_ch*T(i-1,j)-(t_fl+t_ch)*T(i,j)+t_fl*T(i+1,j))/(dx1**2.d0)
        p2=(t_fl+t_ch)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3=(2.d0*p1)+p2-rad+t_fl/(xcp*(ycp2-ycp1))
        Fd =p3/(t_fl*Re*Pr*(1.d0+(rcc/rcf)))
       endif
      endif

      if (i.eq.N1) then
       if ((j.gt.M1(3)).and.(j.lt.M2(3)))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1 =(rad1*T(i,j)+(Tin**4.d0))*sig*emc*dx2
        j1=nc(10)+j-M1(3)
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2 =(rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad =rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1=(t_ch*T(i-1,j)-(t_fl+t_ch)*T(i,j)+t_fl*T(i+1,j))/(dx1**2.d0)
        p2=(t_fl+t_ch)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3=(2.d0*p1)+p2-rad+t_fl/(xcp*(ycp2-ycp1))
        Fd =p3/(t_fl*Re*Pr*(1.d0+(rcc/rcf)))
       endif
      endif

      if (i.eq.N1) then
       if ((j.gt.M1(4)).and.(j.lt.M2(4)))then
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &       (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     &       (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     &       (4.d0*Tre*(Tin**3.d0))
        rad1 =(rad1*T(i,j)+(Tin**4.d0))*sig*emc*dx2
        j1=nc(14)+j-M1(4)
        rad=0.d0
         do i1=1,nc(20)
          rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &         (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &         (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &         (4.d0*Tre*(Tin**3.d0))
          rad2 =(rad2*ph(i1)+(Tin**4.d0))*sig*em(i1)*da(i1)*a_f(i1,j1)
          rad =rad+rad2
         enddo
        rad = 2.d0*(rad1-rad)/(dx1*dx2*Tre)
        p1=(t_ch*T(i-1,j)-(t_fl+t_ch)*T(i,j)+t_fl*T(i+1,j))/(dx1**2.d0)
        p2=(t_fl+t_ch)*(T(i,j+1)-2.d0*T(i,j)+T(i,j-1))/(dx2**2.d0)
        p3=(2.d0*p1)+p2-rad+t_fl/(xcp*(ycp2-ycp1))
        Fd =p3/(t_fl*Re*Pr*(1.d0+(rcc/rcf)))
       endif
      endif

      res=rhs(i,j)+dt*Fd-T(i,j)
      sum=sum+res**2.d0

      enddo
      enddo
      rms=dsqrt(sum/dfloat(N*M))
      iter_t=iter_t+1
      enddo

      write(*,190)iter_t
190   format(1x,i5)

      return
      end
c========================================
      subroutine coeff_aT
      include 'psemiimp.in'
	  
      icvt = 1

* For calculating the coefficient apT
* For energy equation i.e. solve_T
      do j=1,nc(1)
      rs(j)=rhs(N0,j)
      enddo
      do j=nc(1)+1,nc(2)
      rs(j)=rhs(N0+j-nc(1),M1(1)) 
      enddo
      do j=nc(2)+1,nc(3)
      rs(j)=rhs(N1,M1(1)+j-nc(2))
      enddo
      do j=N0,N1-1
      rs(nc(4)+N0-j)=rhs(j,M2(1))  
      enddo
	  
      do j=nc(4)+1,nc(5)
      rs(j)=rhs(N0,M2(1)+j-nc(4))
      enddo
      do j=nc(5)+1,nc(6)
      rs(j)=rhs(j-nc(5)+N0,M1(2))  
      enddo
      do j=nc(6)+1,nc(7)
      rs(j)=rhs(N1,M1(2)+j-nc(6))
      enddo
      do j=N0,N1-1
      rs(nc(8)+N0-j)=rhs(j,M2(2))  
      enddo

      do j=nc(8)+1,nc(9)
      rs(j)=rhs(N0,M2(2)+j-nc(8))
      enddo
      do j=nc(9)+1,nc(10)
      rs(j)=rhs(N0+j-nc(9),M1(3))
      enddo
      do j=nc(10)+1,nc(11)
      rs(j)=rhs(N1,M1(3)+j-nc(10))
      enddo
      do j=N0,N1-1
      rs(nc(12)+N0-j)=rhs(j,M2(3))  
      enddo

      do j=nc(12)+1,nc(13)
      rs(j)=rhs(N0,M2(3)+j-nc(12))
      enddo
      do j=nc(13)+1,nc(14)
      rs(j)=rhs(N0+j-nc(13),M1(4))
      enddo
      do j=nc(14)+1,nc(15)
      rs(j)=rhs(N1,M1(4)+j-nc(14))
      enddo
      do j=N0,N1-1
      rs(nc(16)+N0-j)=rhs(j,M2(4))  
      enddo

      do j=nc(16)+1,nc(17)
      rs(j)=rhs(N0,M2(4)+j-nc(16))
      enddo
      do j=nc(17)+1,nc(18)
      rs(j)=rhs(N0+j-nc(17),M)
      enddo
      do j=1,M-1
      rs(nc(19)+1-j)=rhs(N2,j)
      enddo
      do j=N0+1,N2-1
      rs(nc(20)+N0+1-j)=rhs(j,1)
      enddo

* apT for Complete Geometry
      do i=2,N-1
      do j=2,M-1

* Calculating in fluid and solids (wall and chip)
* Inside Wall
      if((i.lt.N0).or.(i.gt.N2))then
        apT(i,j)=1.d0+(2.d0*dt*a1*t_wa*rcf/(t_fl*Re*Pr*rcw))
      endif
* Inside Chip
      if ((i.gt.N0).and.(i.lt.N1))then
        do k=1,4
          if ((j.gt.M1(k)).and.(j.lt.M2(k)))then
            apT(i,j)=1.d0+(2.d0*dt*a1*t_ch*rcf/(t_fl*rcc*Re*Pr))
          endif
        enddo
      endif
* Inside Fluid
      if ((i.gt.N0).and.(i.lt.N2))then
      chip=0
        do k=1,4
          if ((j.ge.M1(k)).and.(j.le.M2(k)).and.(i.le.N1))then
            chip=1
          endif
        enddo
        if (chip.eq.0)then
          if (icvt.eq.1) then
          a=(1.d0/dx1)*(dmax1(u(1,i,j),0.d0)+dmax1(-u(1,i,j),0.d0))+
     &      (1.d0/dx2)*(dmax1(u(2,i,j),0.d0)+dmax1(-u(2,i,j),0.d0))
          apT(i,j)=1.d0+dt*(a+2.d0*a1/(Re*Pr))
          else if(icvt.eq.0) then
          apT(i,j)=1.d0+2.d0*dt*a1/(Re*Pr)
          endif
        endif
      endif
* Radiating Surfaces and Chip Wall surface
* Left wall
      if (i.eq.N0)then
       if ((j.ge.M1(1)).and.(j.le.M2(1)))then
        apT(i,j)=2.d0*dt*a1*(t_wa+t_ch)/(t_fl*Re*Pr)
        apT(i,j)=1.d0+apT(i,j)/((rcw/rcf)+(rcc/rcf))
       elseif ((j.ge.M1(2)).and.(j.le.M2(2)))then
        apT(i,j)=2.d0*dt*a1*(t_wa+t_ch)/(t_fl*Re*Pr)
        apT(i,j)=1.d0+apT(i,j)/((rcw/rcf)+(rcc/rcf))
       elseif ((j.ge.M1(3)).and.(j.le.M2(3)))then
        apT(i,j)=2.d0*dt*a1*(t_wa+t_ch)/(t_fl*Re*Pr)
        apT(i,j)=1.d0+apT(i,j)/((rcw/rcf)+(rcc/rcf))
       elseif ((j.ge.M1(4)).and.(j.le.M2(4)))then
        apT(i,j)=2.d0*dt*a1*(t_wa+t_ch)/(t_fl*Re*Pr)
        apT(i,j)=1.d0+apT(i,j)/((rcw/rcf)+(rcc/rcf))

       elseif (j.lt.M1(1))then
        j1=j
        i1=j1
        rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &  (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &  (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &  (4.d0*Tre*(Tin**3.d0))
        rad2= rad2*sig*em(i1)*da(i1)*a_f(i1,j1)

        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &  (4.d0*rhs(i,j)*rhs(i,j)*Tre*Tre*Tre*Tin)+
     &  (6.d0*rhs(i,j)*Tre*Tre*Tin*Tin)+
     &  (4.d0*Tre*Tin*Tin*Tin)
        rad1= ((rad1*sig*emb*dx2)-rad2)/(dx1*dx2*Tre)
        apT(i,j)=2.d0*(rad1+(a1*(t_wa+t_fl)))
        apT(i,j)=1.d0+dt*apT(i,j)/((1.d0+(rcw/rcf))*t_fl*Re*Pr)

       elseif (j.gt.M2(4))then
        j1=nc(16)+j-M2(4)
        i1=j1
        rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &  (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &  (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &  (4.d0*Tre*(Tin**3.d0))
        rad2= rad2*sig*em(i1)*da(i1)*a_f(i1,j1)

        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &  (4.d0*rhs(i,j)*rhs(i,j)*Tre*Tre*Tre*Tin)+
     &  (6.d0*rhs(i,j)*Tre*Tre*Tin*Tin)+
     &  (4.d0*Tre*Tin*Tin*Tin)
        rad1=((rad1*sig*emb*dx2)-rad2)/(dx1*dx2*Tre)
        apT(i,j)=2.d0*(rad1+(a1*(t_wa+t_fl)))
        apT(i,j)=1.d0+dt*apT(i,j)/((1.d0+(rcw/rcf))*t_fl*Re*Pr)

       elseif ((j.gt.M2(1)).and.(j.lt.M1(2)))then
        j1=nc(4)+j-M2(1)
        i1=j1
        rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &  (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &  (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &  (4.d0*Tre*(Tin**3.d0))
        rad2= rad2*sig*em(i1)*da(i1)*a_f(i1,j1)

        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &  (4.d0*rhs(i,j)*rhs(i,j)*Tre*Tre*Tre*Tin)+
     &  (6.d0*rhs(i,j)*Tre*Tre*Tin*Tin)+
     &  (4.d0*Tre*Tin*Tin*Tin)
        rad1=((rad1*sig*emb*dx2)-rad2)/(dx1*dx2*Tre)
        apT(i,j)=2.d0*(rad1+(a1*(t_wa+t_fl)))
        apT(i,j)=1.d0+dt*apT(i,j)/((1.d0+(rcw/rcf))*t_fl*Re*Pr)

       elseif ((j.gt.M2(2)).and.(j.lt.M1(3)))then
        j1=nc(8)+j-M2(2)
        i1=j1
        rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &  (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &  (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &  (4.d0*Tre*(Tin**3.d0))
        rad2= rad2*sig*em(i1)*da(i1)*a_f(i1,j1)
	  
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &  (4.d0*rhs(i,j)*rhs(i,j)*Tre*Tre*Tre*Tin)+
     &  (6.d0*rhs(i,j)*Tre*Tre*Tin*Tin)+
     &  (4.d0*Tre*Tin*Tin*Tin)
        rad1=((rad1*sig*emb*dx2)-rad2)/(dx1*dx2*Tre)
        apT(i,j)=2.d0*(rad1+(a1*(t_wa+t_fl)))
        apT(i,j)=1.d0+dt*apT(i,j)/((1.d0+(rcw/rcf))*t_fl*Re*Pr)
	  
       elseif ((j.gt.M2(3)).and.(j.lt.M1(4)))then
        j1=nc(12)+j-M2(3)
        i1=j1
        rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     &  (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     &  (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     &  (4.d0*Tre*(Tin**3.d0))
        rad2= rad2*sig*em(i1)*da(i1)*a_f(i1,j1)
	  
        rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     &  (4.d0*rhs(i,j)*rhs(i,j)*Tre*Tre*Tre*Tin)+
     &  (6.d0*rhs(i,j)*Tre*Tre*Tin*Tin)+
     &  (4.d0*Tre*Tin*Tin*Tin)
        rad1=((rad1*sig*emb*dx2)-rad2)/(dx1*dx2*Tre)
        apT(i,j)=2.d0*(rad1+(a1*(t_wa+t_fl)))
        apT(i,j)=1.d0+dt*apT(i,j)/((1.d0+(rcw/rcf))*t_fl*Re*Pr)

       endif
      endif
* Right Wall
      if (i.eq.N2)then

      j1=nc(19)-j+1	  
      i1=j1
      rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     & (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     & (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad2= rad2*sig*em(i1)*da(i1)*a_f(i1,j1)

      rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     & (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     & (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad1=((rad1*sig*emb*dx2)-rad2)/(dx1*dx2*Tre)

      apT(i,j)=2.d0*((a1*(t_wa+t_fl))-rad1)
      apT(i,j)=1.d0+dt*apT(i,j)/((1.d0+(rcw/rcf))*t_fl*Re*Pr)

      endif
* Top face of chips
      if (j.eq.M2(1))then
      if ((i.gt.N0).and.(i.le.N1))then

      j1=nc(4)-i+N0
      i1=j1
      rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     & (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     & (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad2= rad2*sig*em(i1)*da(i1)*a_f(i1,j1)

      rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     & (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     & (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad1=((rad1*sig*emc*dx1)-rad2)/(dx1*dx2*Tre)

      apT(i,j)=((t_fl+t_ch)*a1+rad1)/(1.d0+(rcc/rcf))
      apT(i,j)=1.d0+(2.d0*dt*apT(i,j))/(t_fl*Re*Pr)
      endif
      endif

      if (j.eq.M2(2))then
      if ((i.gt.N0).and.(i.le.N1))then

      j1=nc(8)-i+N0
      i1=j1
      rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     & (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     & (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad2= rad2*sig*em(i1)*da(i1)*a_f(i1,j1)

      rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     & (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     & (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad1=((rad1*sig*emc*dx1)-rad2)/(dx1*dx2*Tre)

      apT(i,j)=((t_fl+t_ch)*a1+rad1)/(1.d0+(rcc/rcf))
      apT(i,j)=1.d0+(2.d0*dt*apT(i,j))/(t_fl*Re*Pr)
      endif
      endif


      if (j.eq.M2(3))then
      if ((i.gt.N0).and.(i.le.N1))then

      j1=nc(12)-i+N0
      i1=j1
      rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     & (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     & (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad2= rad2*sig*em(i1)*da(i1)*a_f(i1,j1)

      rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     & (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     & (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad1=((rad1*sig*emc*dx1)-rad2)/(dx1*dx2*Tre)

      apT(i,j)=((t_fl+t_ch)*a1+rad1)/(1.d0+(rcc/rcf))
      apT(i,j)=1.d0+(2.d0*dt*apT(i,j))/(t_fl*Re*Pr)
      endif
      endif
	  
      if (j.eq.M2(4))then
      if ((i.gt.N0).and.(i.le.N1))then

      j1=nc(16)-i+N0
      i1=j1
      rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     & (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     & (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad2= rad2*sig*em(i1)*da(i1)*a_f(i1,j1)

      rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     & (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     & (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad1=((rad1*sig*emc*dx1)-rad2)/(dx1*dx2*Tre)

      apT(i,j)=((t_fl+t_ch)*a1+rad1)/(1.d0+(rcc/rcf))
      apT(i,j)=1.d0+(2.d0*dt*apT(i,j))/(t_fl*Re*Pr)
      endif
      endif
* Bottom face of chips
      if (j.eq.M1(1))then
      if ((i.gt.N0).and.(i.le.N1))then

      j1=nc(1)+i-N0
      i1=j1
      rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     & (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     & (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad2= rad2*sig*em(i1)*da(i1)*a_f(i1,j1)
	  
      rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     & (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     & (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad1=((rad1*sig*emc*dx1)-rad2)/(dx1*dx2*Tre)

      apT(i,j)=((t_fl+t_ch)*a1-rad1)/(1.d0+(rcc/rcf))
      apT(i,j)=1.d0+(2.d0*dt*apT(i,j))/(t_fl*Re*Pr)
      endif
      endif

      if (j.eq.M1(2))then
      if ((i.gt.N0).and.(i.le.N1))then

      j1=nc(5)+i-N0
      i1=j1
      rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     & (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     & (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad2= rad2*sig*em(i1)*da(i1)*a_f(i1,j1)
	  
      rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     & (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     & (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad1=((rad1*sig*emc*dx1)-rad2)/(dx1*dx2*Tre)

      apT(i,j)=((t_fl+t_ch)*a1-rad1)/(1.d0+(rcc/rcf))
      apT(i,j)=1.d0+(2.d0*dt*apT(i,j))/(t_fl*Re*Pr)
      endif
      endif

      if (j.eq.M1(3))then
      if ((i.gt.N0).and.(i.le.N1))then

      j1=nc(9)+i-N0
      i1=j1
      rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     & (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     & (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad2= rad2*sig*em(i1)*da(i1)*a_f(i1,j1)

      rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     & (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     & (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad1=((rad1*sig*emc*dx1)-rad2)/(dx1*dx2*Tre)

      apT(i,j)=((t_fl+t_ch)*a1-rad1)/(1.d0+(rcc/rcf))
      apT(i,j)=1.d0+(2.d0*dt*apT(i,j))/(t_fl*Re*Pr)
      endif
      endif
	  
      if (j.eq.M1(4))then
      if ((i.gt.N0).and.(i.le.N1))then

      j1=nc(13)+i-N0
      i1=j1
      rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     & (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     & (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad2= rad2*sig*em(i1)*da(i1)*a_f(i1,j1)

      rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     & (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     & (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad1=((rad1*sig*emc*dx1)-rad2)/(dx1*dx2*Tre)

      apT(i,j)=((t_fl+t_ch)*a1-rad1)/(1.d0+(rcc/rcf))
      apT(i,j)=1.d0+(2.d0*dt*apT(i,j))/(t_fl*Re*Pr)
      endif
      endif
* Front face of chip
      if (i.eq.N1) then
      if ((j.gt.M1(1)).and.(j.lt.M2(1)))then

      j1=nc(2)+j-M1(1)
      i1=j1
      rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     & (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     & (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad2= rad2*sig*em(i1)*da(i1)*a_f(i1,j1)

      rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     & (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     & (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad1=((rad1*sig*emc*dx2)-rad2)/(dx1*dx2*Tre)

      apT(i,j)=((t_fl+t_ch)*a1+rad1)/(1.d0+(rcc/rcf))
      apT(i,j)=1.d0+(2.d0*dt*apT(i,j))/(t_fl*Re*Pr)
      endif

      if ((j.gt.M1(2)).and.(j.lt.M2(2)))then

      j1=nc(6)+j-M1(2)
      i1=j1
      rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     & (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     & (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad2= rad2*sig*em(i1)*da(i1)*a_f(i1,j1)

      rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     & (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     & (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad1=((rad1*sig*emc*dx2)-rad2)/(dx1*dx2*Tre)
      
      apT(i,j)=((t_fl+t_ch)*a1+rad1)/(1.d0+(rcc/rcf))
      apT(i,j)=1.d0+(2.d0*dt*apT(i,j))/(t_fl*Re*Pr)
      endif

      if ((j.gt.M1(3)).and.(j.lt.M2(3)))then

      j1=nc(10)+j-M1(2)
      i1=j1
      rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     & (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     & (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad2= rad2*sig*em(i1)*da(i1)*a_f(i1,j1)

      rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     & (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     & (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad1=((rad1*sig*emc*dx2)-rad2)/(dx1*dx2*Tre)

      apT(i,j)=((t_fl+t_ch)*a1+rad1)/(1.d0+(rcc/rcf))
      apT(i,j)=1.d0+(2.d0*dt*apT(i,j))/(t_fl*Re*Pr)
      endif

      if ((j.gt.M1(4)).and.(j.lt.M2(4)))then

      j1=nc(14)+j-M1(2)
      i1=j1
      rad2=((rs(i1)**3.d0)*(Tre**4.d0))+
     & (4.d0*rs(i1)*rs(i1)*(Tre**3.d0)*Tin)+
     & (6.d0*rs(i1)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad2= rad2*sig*em(i1)*da(i1)*a_f(i1,j1)

      rad1=((rhs(i,j)**3.d0)*(Tre**4.d0))+
     & (4.d0*rhs(i,j)*rhs(i,j)*(Tre**3.d0)*Tin)+
     & (6.d0*rhs(i,j)*(Tre**2.d0)*Tin*Tin)+
     & (4.d0*Tre*(Tin**3.d0))
      rad1=((rad1*sig*emc*dx2)-rad2)/(dx1*dx2*Tre)
      
      apT(i,j)=((t_fl+t_ch)*a1+rad1)/(1.d0+(rcc/rcf))
      apT(i,j)=1.d0+(2.d0*dt*apT(i,j))/(t_fl*Re*Pr)
      endif

      endif

      enddo
      enddo
	  
      return
      end
c========================================
      subroutine tstep
      include 'psemiimp.in'

      cux=100.d0
      cuy=100.d0
      do i=N0+1,N2-1
      do j=2,M-1
      bad=0
      if(i.le.N1)then
      do k=1,L
      if((j.ge.M1(k)).and.(j.le.M2(k)))then
      bad=1
      endif
      enddo
      endif
* Courant Fredricks Levy Criterion
      if (bad.eq.0)then
      cux=dmin1(cux,dabs(dx1/u(1,i,j)))
      cuy=dmin1(cuy,dabs(dx2/u(2,i,j)))
      endif
      enddo
      enddo
      cumax=dmin1(cux,cuy)
      dtc=cumax
* Fourier Stability Criterion
      dtf=(0.25d0*Re)/((1.d0/(dx1**2.d0))+(1.d0/(dx2**2.d0)))
* Minimum of Both
      dt=dmin1(dtc,dtf)

      return
      end
C========================================
      subroutine temp_bc
      include 'psemiimp.in'

* At inlet side 2 corners 
* Inflow and Outflow conditions
* H is Always = 1 so not taken as non-dimensional parameter*
      do i=1,N0-1
      T(i,1)=T(i,2)/(1+(ht*dx2/t_wa))
      T(i,M)=T(i,M-1)/(1+(ht*dx2/t_wa))
      enddo
      do i=N0,N2
      T(i,1)=0.d0
      T(i,M)=T(i,M-1)
      enddo
* H is Always = 1 so not taken as non-dimensional parameter*
      do i=N2+1,N
      T(i,1)=T(i,2)/(1+(ht*dx2/t_wa))
      T(i,M)=T(i,M-1)/(1+(ht*dx2/t_wa))
      enddo
* At right and left side of walls
* H is Always = 1 so not taken as non-dimensional parameter*
      do j=2,M-1
      T(1,j)=T(2,j)/(1+(ht*dx1/t_wa))
      T(N,j)=T(N-1,j)/(1+(ht*dx1/t_wa))
      enddo

      return
      end
C========================================
      subroutine vort_bc
      include 'psemiimp.in'

      do i=N0,N2
* Inflow and Outflow conditions
*      vort(i,1)=6.d0-12.d0*x(i)
      vort(i,1)=0.d0
      vort(i,M)=vort(i,M-1)
* At Chip's Lower and Upper Surface
      if(i.le.N1)then
      do k=1,L
      vort(i,M1(k))=2.d0*(si(i,M1(k))-si(i,M1(k)-1))/(dx2**2.d0)
      vort(i,M2(k))=2.d0*(si(i,M2(k))-si(i,M2(k)+1))/(dx2**2.d0)
      enddo
      endif
* Condition started from N0 only
      enddo

      do j=2,M-1
* Started from 2 as 1 is covered above same for M
      bad=0
      do k=1,L
      if ((j.gt.M1(k)).and.(j.lt.M2(k)))then
      bad=1
      endif
      enddo
      if (bad.eq.0) then
* At Left Wall i.e. at N0
      vort(N0,j)=2.d0*(si(N0,j)-si(N0+1,j))/(dx1**2.d0)
      else
* At Chip's Vertical Surface i.e. at N1
      vort(N1,j)=2.d0*(si(N1,j)-si(N1+1,j))/(dx1**2.d0)
      endif
* At Right Wall i.e. at N2
      vort(N2,j)=2.d0*(si(N2,j)-si(N2-1,j))/(dx1**2.d0)
      enddo

      return
      end
C========================================
      subroutine si_bc
      include 'psemiimp.in'

      do i=N0,N2
* Inflow and Outflow conditions
*      si(i,1)=2.d0*x(i)**3.d0-3.d0*x(i)**2.d0
      si(i,1)=-1.d0*x(i)
      si(i,M)=si(i,M-1)
* At Chip's Lower and Upper Surface
      if(i.le.N1)then
      do k=1,L
      si(i,M1(k))=0.d0
      si(i,M2(k))=0.d0
      enddo
      endif
* Condition started from N0 only
      enddo

      do j=2,M-1
* Started from 2 as 1 is covered above same for M
      bad=0
      do k=1,L
      if ((j.gt.M1(k)).and.(j.lt.M2(k)))then
      bad=1
      endif
      enddo
      if (bad.eq.0)then
* At Left Wall i.e. at N0
      si(N0,j)=0.d0
      else
* At Chip's Vertical Surface i.e. at N1
      si(N1,j)=0.d0
      endif
* At Right Wall i.e. at N2
      si(N2,j)=si(N2,1)
      enddo

      return
      end
C========================================
      subroutine vel_bc
      include 'psemiimp.in'

* For Parabolic velocity profile at Inlet (Channel)
* x(i) = Grid Point distance in x-direction x(N0)=0 and x(N2)=1
      do i=N0,N2
      x(i)=0.d0+dx1*dfloat(i-N0)
      enddo
      
      do i=N0,N2
* Inflow and Outflow conditions
      u(1,i,1)=0.d0
      u(1,i,M)=u(1,i,M-1)
*      u(2,i,1)=6.d0*x(i)-6.d0*x(i)**2.d0
      u(2,i,1)=1.d0
      u(2,i,M)=u(2,i,M-1)
* At Chip's Lower and Upper Surface
      if(i.le.N1)then
      do k=1,L
      u(1,i,M1(k))=0.d0
      u(2,i,M1(k))=0.d0
      u(1,i,M2(k))=0.d0
      u(2,i,M2(k))=0.d0
      enddo
      endif
* Condition started from N0 only
      enddo

      do j=2,M-1
* Started from 2 as 1 is covered above same for M
      bad=0
      do k=1,L
      if ((j.gt.M1(k)).and.(j.lt.M2(k)))then
      bad=1
      endif
      enddo
      if (bad.eq.0)then
* At Left Wall i.e. at N0
      u(1,N0,j)=0.d0
      u(2,N0,j)=0.d0
      else
* At Chip's Vertical Surface i.e. at N1
      u(1,N1,j)=0.d0
      u(2,N1,j)=0.d0
      endif
* At Right Wall i.e. at N2
      u(1,N2,j)=0.d0
      u(2,N2,j)=0.d0
      enddo

      return
      end

C========================================
      subroutine restart
      include 'psemiimp.in'

      open(8,file='vo.dat')
      do i=N0,N2
      do j=1,M
      bad=0

      if(i.lt.N1)then
      do k=1,L
      if((j.gt.M1(k)).and.(j.lt.M2(k)))then
      bad=1
      endif
      enddo
      endif

      if (bad.eq.0)then
      read(8,220)vort(i,j),si(i,j)
220   format(16x,f22.10,8x,f22.10)
      endif
      enddo
      enddo
      close(8)

      open(9,file='ve.dat')
      do i=N0,N2
      do j=1,M
      bad=0

      if(i.lt.N1)then
      do k=1,L
      if((j.gt.M1(k)).and.(j.lt.M2(k)))then
      bad=1
      endif
      enddo
      endif

      if (bad.eq.0)then
      read(9,230)u(1,i,j),u(2,i,j)
230   format(2x,f22.10,1x,f22.10)
      endif
      enddo
      enddo
      close(9)

      open(10,file='temp.dat')
      do i=1,N
      do j=1,M
      read(10,300)T(i,j)
300   format(2x,f22.10)
      enddo
      enddo
      close(10)

      return
      end

C========================================
      subroutine init
      include 'psemiimp.in'

* Defining 0.d0 everywhere for vorticity, stream function
* temperature, velocity (U, V)
      do i=1,N
      do j=1,M
      vort(i,j)=0.d0
      si(i,j)=0.d0
      T(i,j)=0.d0
      u(1,i,j)=0.d0	  
      u(2,i,j)=0.d0
      enddo
      enddo

      return
      end
