!      HotCross Distribution - Version 1.0

       subroutine density                                               

!      tag: subroutine ../vents2/lavelle/NewCross/density_3d.f - version 6.4.97

!                                                                       
! ----- perturbation density - rho prime
!                                                                       
       use tracer_module
       use grid_module

       include 'array_size.inc'

       real rho_rho, xisw
       integer i,j,k,isw

       isw = 0
       xisw = real (isw)
!$OMP PARALLEL DO

       do  k = 1, nz 
       do  j = 1, ny   
       do  i = 1, nx                            

          rhop(i,j,k) = rho_rho(s(i,j,k), t(i,j,k), 0.) - rho_bkg(k) 
          rrho(i,j,k)=rhop(i,j,k)+rho_bkg(k)

       enddo;enddo;enddo

!$OMP END PARALLEL DO

       do k = 1, nz
       do j = 1, ny

        i = 0
        rhop(i,j,k) = rho_rho(s(i,j,k),t(i,j,k),0.) - rho_bkg(k)
        rrho(i,j,k) = rhop(i,j,k) + rho_bkg(k)
        i = nx+1
        rhop(i,j,k) = rho_rho(s(i,j,k),t(i,j,k),0.) - rho_bkg(k)
        rrho(i,j,k) = rhop(i,j,k) + rho_bkg(k)

       enddo;enddo

       do k = 1, nz
       do i = 1, nx

        j = 0
        rhop(i,j,k) = rho_rho(s(i,j,k),t(i,j,k),0.) - rho_bkg(k)
        rrho(i,j,k) = rhop(i,j,k) + rho_bkg(k)
        j = ny+1
        rhop(i,j,k) = rho_rho(s(i,j,k),t(i,j,k),0.) - rho_bkg(k)
        rrho(i,j,k) = rhop(i,j,k) + rho_bkg(k)

       enddo;enddo

       do  j = 1, ny
       do  i = 1, nx

        k = 0
        rhop(i,j,k) = rho_rho(s(i,j,k),t(i,j,k),0.) - rho_bkg(k)
        rrho(i,j,k) = rhop(i,j,k) + rho_bkg(k)
        k = nz+1
        rhop(i,j,k) = rho_rho(s(i,j,k),t(i,j,k),0.) - rho_bkg(k)
        rrho(i,j,k) = rhop(i,j,k) + rho_bkg(k)

       enddo;enddo
       
!       call view_3d(rhop,0,nx+1,0,ny+1,0,nz+1); write(6,*) 'rhop';stop
       
       end subroutine density

! --------------------------------------------------------------------
! ------------------ Equation of state of seawater -------------------
! --------------------------------------------------------------------

       real function rho_rho(s,t,p)

! ---  seawater density for 0<T<40 C and 0<S<42 0/00 and 0<p<1000 bars
! ---  from Gill, pp. 599-600
! ---  JW Lavelle, 1986,  checked against UNESCO and PMEL
! ---  equal to rho, named changed to rho_rho to free variable in code

       real  kw, kst0, kstp
       real    p, s, t
       real  rho0, rho1

!       if( (t.gt.40).or.(s.gt.42).or.(t.lt.0.0).or.(s.lt.0.0).or.          !dummed to speed calculation
!     1      (p.lt.0.0).or.(p.gt.1000.0) ) then
!             write(6,*) 'density parameters out of range'
!             stop
!       endif

       rho0 = 999.842594d0 + 6.793952d-2*t - 9.095290d-3*t*t 
     1  +1.001685d-4*t**3 - 1.120083d-6*t**4 + 6.536332d-9*t**5

       rho1 = rho0 + s*(0.824493d0 - 4.0899d-3*t + 7.6438d-5*t*t 
     1  - 8.2467d-7*t**3 + 5.3875d-9*t**4) + s**1.5*( - 5.72466d-3 
     2  + 1.0227d-4*t - 1.6546d-6*t*t) + 4.8314d-4*s*s

       rho_rho = rho1

! if (p.ewq. 0) skip the following

!       kw = 19652.21d0 + 148.4206d0*t - 2.327105d0*t*t 
!     1  + 1.360477d-2*t**3 - 5.155288d-5*t**4

!       kst0 = kw + s*( 54.6746d0 - 0.603459d0*t + 1.09987d-2*t*t 
!     1  - 6.1670d-5*t**3) + s**1.5*( 7.944d-2 + 1.6483d-2*t 
!     2  - 5.3009d-4*t*t)

!          kstp = kst0 + p*( 3.239908d0 + 1.43713d-3*t + 1.16092d-4*t*t 
!     1    - 5.77905d-7*t**3) + p*s*( 2.2838d-3 - 1.0981d-5*t 
!     2    - 1.6078d-6*t*t) + 1.91075d-4*p*s**1.5 + p*p*( 8.50935d-5 
!     3    - 6.12293d-6*t + 5.2787d-8*t*t) + p*p*s*( - 9.9348d-7 
!     4    + 2.0816d-8*t + 9.1697d-10*t*t)

!       rho_rho = rho1 / ( 1.0d0 - p / kstp )

       return

       end

