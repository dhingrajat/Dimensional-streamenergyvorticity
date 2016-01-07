      implicit double precision(a-h,o-z)
      integer bad,chip,icvt
* icvt defines the convective scheme = 1/0 (upwind/central difference)
* N is number of grid points in x direction for trv=1.2

      parameter(N=61,M=601,L=4,Re=500.d0,Pr=0.71d0,Gr=8.65d5,ht=0.d0)
* Thermal Conductivities   ------------   Wall Thickness
      parameter(t_ch=120.d0,t_fl=2.64d-2,t_wa=79.8d0,t_w=0.1d0)
* Chip related dimensions  ------------   Heat generation
      parameter(xcp=0.3d0,ycp1=2.d0,ycp2=2.6d0,ycp3=3.d0,qv=1.d5)
* Emissivities board - chip - exit & entry
      parameter(emb=0.55d0,emc=0.55d0,eme=1.d0,sig=5.67d-8)
* density*thermal capacity of chip - fluid - wall
      parameter(Tin=300,rcc=1638.d0,rcf=1.182d0,rcw=1575.66d0)

      common dt,rms,N0,N1,N2,a1,Tre,iter_s,res,ae,aw,ap,beta,a
      common iteration,iter_t,iter_v,M1(L),M2(L),nc(20),x(N),cvy
      common dx1,dx2,bad,chip,da(1462),em(1462),rad,rad1,rad2,Fc1,Fc2
      common apT(N,M),apV(N,M),rs(1462),ph(1462),rhs(N,M),si(N,M)
      common u(2,N,M),a_f(1462,1462),vort_old(N,M),vort(N,M),T(N,M),
     & T_old(N,M)
