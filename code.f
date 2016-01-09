        include 'data.in'

        double precision trv,axl
        trv = 1.2d0
        axl = 12.d0

        dx=1.2/(N-1)
        dy=axl/dfloat(M-1)
        M1(1)= int(ycp1/dy)+1
        M2(1)= int(ycp2/dy)+1
        M1(2)= int(ycp3/dy)+1
        do k = 1,L
          M1(k)=M1(1)+(k-1)*(M1(2)-M1(1))
          M2(k)=M2(1)+(k-1)*(M1(2)-M1(1))
*         write(*,33)M1(k),M2(k)
*33       format(1x,i5,2x,i5)
        enddo
      
        N0=int(t_w/dx)+1
        N1=int(xcp/dx)+N0
        N2=N-N0+1
        a1=(1.0/dx**2.d0)+(1.0/dx2**2.d0)

* Initilaization
        do i=1,N
          do j=1,M
            W(i,j)=0.0  !Vorticity
            H(i,j)=0.0  !Stream function
            T(i,j)=0.0  !Energy
            u(i,j)=0.0  !u-velocity	  
            v(i,j)=0.0  !v-velocity
          enddo
        enddo

C If Required call restart

* Boundary Velocity
        do i=N0,N2
          u(i,1)=0.0
          u(i,M)=u(i,M-1)
            if(i.le.N1)then
              do k=1,L
                u(i,M1(k))=0.0
                u(i,M2(k))=0.0
                v(i,M1(k))=0.0
                v(i,M2(k))=0.0
              enddo
            endif
          v(i,1)=vinlet
          v(i,M)=v(i,M-1)
        enddo
        do j=2,M-1
          ib=0
          do k=1,L
            if ((j.gt.M1(k)).and.(j.lt.M2(k)))then
              ib=1
            endif
          enddo
          if (ib.eq.0)then
            u(N0,j)=0.0
            v(N0,j)=0.0
          else
            u(N1,j)=0.0
            v(N1,j)=0.0
          endif
          u(N2,j)=0.0
          v(N2,j)=0.0
        enddo

* Boundary Stream-Function || x(N0)=0 and x(N2)=1
        do i=N0,N2
          x(i)=0.0+dx*dfloat(i-N0)
        enddo
        do i=N0,N2
          H(i,1)=-1.0*x(i)
          H(i,M)=H(i,M-1)
          if(i.le.N1)then
            do k=1,L
              H(i,M1(k))=0.0
              H(i,M2(k))=0.0
            enddo
          endif
        enddo
        do j=2,M-1
          ib=0
          do k=1,L
            if ((j.gt.M1(k)).and.(j.lt.M2(k)))then
              ib=1
            endif
          enddo
          if (ib.eq.0)then
            H(N0,j)=0.0
          else
            H(N1,j)=0.0
          endif
          H(N2,j)=H(N2,1)
        enddo

* Boundary Vorticity
        do i=N0,N2
          W(i,1)=0.0
          W(i,M)=W(i,M-1)
          if (i.le.N1)then
            do k=1,L
              W(i,M1(k))=2.0*(H(i,M1(k))-H(i,M1(k)-1))/dy**2
              W(i,M2(k))=2.0*(H(i,M2(k))-H(i,M2(k)+1))/dy**2
            enddo
          endif
        enddo
        do j=2,M-1
          ib=0
          do k=1,L
            if ((j.gt.M1(k)).and.(j.lt.M2(k)))then
              ib=1
            endif
          enddo
          if (ib.eq.0)then
            W(N0,j)=2.0*(H(N0,j)-H(N0+1,j))/dx**2
          else
            W(N1,j)=2.0*(H(N1,j)-H(N1+1,j))/dx**2
          endif
          W(N2,j)=2.0*(H(N2,j)-H(N2-1,j))/dx**2
        enddo

* Boundary Temperature
        do i=1,N0-1
          T(i,1)=T(i,2)
          T(i,M)=T(i,M-1)
        enddo
        do i=N0,N2
          T(i,1)=300.0
          T(i,M)=T(i,M-1)
        enddo
        do i=N2+1,N
          T(i,1)=T(i,2)
          T(i,M)=T(i,M-1)
        enddo
        do j=2,M-1
          T(1,j)=T(2,j)
          T(N,j)=T(N-1,j)
        enddo

* Main Iterative
        rms=1.0
        dt=0.001
        do while(rms.gt.1e-5)
          rms=0.0
          do i=1,N
            do j=1,M
              W1(i,j)=W(i,j)
              T1(i,j)=T(i,j)
            enddo
          enddo

* Energy loop
          rmt=1.0
          do while(rmt.gt.1e-4)
            rmt=0.0
            do i=1,N
              do j=1,M
                RS(i,j)=T(i,j)
              enddo
            enddo
            do i=1,N
              do j=1,M
                if (i.lt.N0)then
                  if ((i.eq.1).and.(j.eq.1))then
                    T(i,j)=(T1(i,j)/dt+t1/(r1*c1)*(2.0*T(i+1,j)
     &                /dx**2+2.0*T(i,j+1)/dy**2))
     &                /(1/dt+2.0*b1/(r1*c1)*a1)                  
                  elseif ((j.eq.1).and.(i.gt.1))then
                    T(i,j)=(T1(i,j)/dt+b1/(r1*c1)*((T(i+1,j)+T(i-1,j))
     &                /dx**2+2.0*T(i,j+1)/dy**2))
     &                /(1/dt+2.0*b1/(r1*c1)*a1)
                  elseif ((i.eq.1).and.(j.gt.1).and.(j.lt.M))then
                    T(i,j)=(T1(i,j)/dt+b1/(r1*c1)*(2.0*T(i+1,j)
     &                /dx**2+(T(i,j+1)+T(i,j-1))/dy**2))
     &                /(1/dt+2.0*b1/(r1*c1)*a1)
                  elseif ((i.eq.1).and.(j.eq.M))then
                    T(i,j)=(T1(i,j)/dt+b1/(r1*c1)*(2.0*T(i+1,j)
     &                /dx**2+2.0*T(i,j-1)/dy**2))
     &                /(1/dt+2.0*b1/(r1*c1)*a1)
                  elseif ((j.eq.M).and.(i.gt.1))then
                    T(i,j)=(T1(i,j)/dt+b1/(r1*c1)*((T(i+1,j)+T(i-1,j))
     &                /dx**2+2.0*T(i,j-1)/dy**2))
     &                /(1/dt+2.0*b1/(r1*c1)*a1)
                  else
                    T(i,j)=(T1(i,j)/dt+b1/(r1*c1)*((T(i+1,j)+T(i-1,j))
     &                /dx**2+(T(i,j+1)+T(i,j-1))/dy**2))
     &                /(1/dt+2.0*b1/(r1*c1)*a1)
                  endif
                elseif (i.gt.N2)then
                  if ((i.eq.N).and.(j.eq.1))then
                    T(i,j)=(T1(i,j)/dt+b1/(r1*c1)*(2.0*T(i-1,j)
     &                /dx**2+2.0*T(i,j+1)/dy**2))
     &                /(1/dt+2.0*b1/(r1*c1)*a1)  
                  elseif ((j.eq.1).and.(i.lt.N))then
                    T(i,j)=(T1(i,j)/dt+b1/(r1*c1)*((T(i+1,j)+T(i-1,j))
     &                /dx**2+2.0*T(i,j+1)/dy**2))
     &                /(1/dt+2.0*b1/(r1*c1)*a1)
                  elseif ((i.eq.N).and.(j.gt.1).and.(j.lt.M))then
                    T(i,j)=(T1(i,j)/dt+b1/(r1*c1)*(2.0*T(i-1,j)
     &                /dx**2+(T(i,j+1)+T(i,j-1))/dy**2))
     &                /(1/dt+2.0*b1/(r1*c1)*a1)
                  elseif ((i.eq.N).and.(j.eq.M))then
                    T(i,j)=(T1(i,j)/dt+b1/(r1*c1)*(2.0*T(i-1,j)
     &                /dx**2+2.0*T(i,j-1)/dy**2))
     &                /(1/dt+2.0*b1/(r1*c1)*a1)
                  elseif ((j.eq.M).and.(i.lt.N))then
                    T(i,j)=(T1(i,j)/dt+b1/(r1*c1)*((T(i+1,j)+T(i-1,j))
     &                /dx**2+2.0*T(i,j-1)/dy**2))
     &                /(1/dt+2.0*b1/(r1*c1)*a1)
                  else
                    T(i,j)=(T1(i,j)/dt+b1/(r1*c1)*((T(i+1,j)+T(i-1,j))/dx**2+
     &                (T(i,j+1)+T(i,j-1))/dy**2))/(1/dt+2.0*b1/(r1*c1)*a1)
                  endif
                elseif (i.eq.N0)then
                  if (j.eq.1)then
                    T(i,j)=300.0
                  elseif (j.eq.M)then
                    T(i,j)=(T1(i,j)*(r1*c1+r*c)/(2.0*dt)+((b1*T(i-1,j)+
     &                b*T(i+1,j))/dx**2+(b1+b)*T(i,j-1)/dy**2))
     &                /((r1*c1+r*c)/(2.0*dt)+(b1+b)*a1)
                  else
                    ib=0
                    do k=1,L
                      if ((j.gt.M1(k)).and.(j.lt.M2(k)))then
                        ib=1
                      endif
                    enddo
                    if (ib.eq.0)then
                      T(i,j)=(T1(i,j)*(r1*c1+r*c)/(2.0*dt)+((b1*T(i-1,j)+
     &                  b*T(i+1,j))/dx**2+(b1+b)/2.0*(T(i,j+1)+T(i,j-1))/
     &                  dy**2))/((r1*c1+r*c)/(2.0*dt)+(b1+b)*a1)
                    else
                      T(i,j)=(T1(i,j)*(r1*c1+r2*c2)/(2.0*dt)+((b1*T(i-1,j)+
     &                  b2*T(i+1,j))/dx**2+(b1+b2)/2.0*(T(i,j+1)+T(i,j-1))/
     &                  dy**2)+q/(2.0*r2*c2))/((r1*c1+r2*c2)/(2.0*dt)+
     &                  (b1+b2)*a1)
                    endif
                  endif
                elseif (i.eq.N2)then
                  if (j.eq.1)then
                    T(i,j)=300.0
                  elseif (j.eq.M)then
                    T(i,j)=(T1(i,j)*(r1*c1+r*c)/(2.0*dt)+((b*T(i-1,j)+
     &                b1*T(i+1,j))/dx**2+(b1+b)*T(i,j-1)/dy**2))
     &                /((r1*c1+r*c)/(2.0*dt)+(b1+b)*a1)
                  else
                    T(i,j)=(T1(i,j)*(r1*c1+r*c)/(2.0*dt)+((b*T(i-1,j)+
     &                b1*T(i+1,j))/dx**2+(b1+b)/2.0*(T(i,j+1)+T(i,j-1))/
     &                dy**2))/((r1*c1+r*c)/(2.0*dt)+(b1+b)*a1)
                  endif
                else
                  if (j.eq.1)then
                    T(i,j)=300.0
                  elseif (j.eq.M)then
                    T(i,j)=(T1(i,j)/dt+b/(r*c)*((T(i+1,j)+T(i-1,j))/
     &                dx**2+2.0*T(i,j-1)/dy**2)-u(i,j)*(T(i+1,j)-
     &                T(i-1,j))/(2.0*dx))/(1/dt+2.0*b/(r*c)*a1)
                  else
                    ib=0
                    do k=1,L
                      if ((j.ge.M1(k)).and.(j.le.M2(k)).and.(i.le.N1))
     &                  then
                        ib=1  
                      endif
                    enddo
                    if (ib.eq.0)then
                      T(i,j)=(T1(i,j)/dt+b/(r*c)*((T(i+1,j)+T(i-1,j))/
     &                  dx**2+(T(i,j+1)+T(i,j-1))/dy**2)-u(i,j)*
     &                  (T(i+1,j)-T(i-1,j))/(2.0*dx)-v(i,j)*(T(i,j+1)
     &                  -T(i,j-1))/(2.0*dy))/(1/dt+2.0*b/(r*c)*a1)
                    else
                      ib=0
                      do k=1,L
                        if (j.eq.M1(k))then
                          ib=1
                        elseif (j.eq.M2(k))then
                          ib=2
                        elseif (i.eq.N1)then
                          ib=3
                        endif
                      enddo
                      if (ib.eq.1)then
                        T(i,j)=(T1(i,j)*(r2*c2+r*c)/(2.0*dt)+((b2+b)
     &                    /2.0*(T(i+1,j)+T(i-1,j))/dx**2+(b*T(i,j-1)
     &                    +b2*T(i,j+1))/dy**2)+q/(2.0*r2*c2))/((r2*
     &                    c2+r*c)/(2.0*dt)+(b2+b)*a1)
                      elseif (ib.eq.2)then
                        T(i,j)=(T1(i,j)*(r2*c2+r*c)/(2.0*dt)+((b2+b)
     &                    /2.0*(T(i+1,j)+T(i-1,j))/dx**2+(b2*T(i,j-1)
     &                    +b*T(i,j+1))/dy**2)+q/(2.0*r2*c2))/((r2*c2+
     &                    r*c)/(2.0*dt)+(b2+b)*a1)
                      elseif (ib.eq.3)then
                        T(i,j)=(T1(i,j)*(r2*c2+r*c)/(2.0*dt)+((b2*
     &                    T(i-1,j)+b*T(i+1,j))/dx**2+(b2+b)/2.0*
     &                    (T(i,j+1)+T(i,j-1))/dy**2)+q/(2.0*r2*c2))/
     &                    ((r2*c2+r*c)/(2.0*dt)+(b2+b)*a1)
                      else
                        T(i,j)=(q/(r2*c2)+T1(i,j)/dt+b1/(r1*c1)*
     &                    ((T(i+1,j)+T(i-1,j))/dx**2+(T(i,j+1)+
     &                    T(i,j-1))/dy**2) )/(1/dt+2.0*(b1/r1*c1)*a1)
                      endif
                    endif
                  endif
                endif
              enddo
            enddo
            do i=1,N
              do j=1,M
                rmt=dmax1(rmt,abs(T(i,j)-RS(i,j)))
              enddo
            enddo
          enddo
                
        
main loop--        enddo



