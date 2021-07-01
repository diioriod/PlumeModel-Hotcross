!      HotCross Distribution - Version 1.0
!      tag: subroutine ../vents2/lavelle/NewCross/viscosity3d.f - version 6.4.97
!           version that kills vertical smag contribution

       subroutine viscosity3d

       use viscosity_module    !deprecated, ss => ai
       use grid_module
       use velocity_module
       use tracer_module
       use constants_module
       use signals_module
       include 'array_size.inc'

       real coef, Rf, eps30
       real s11,s12,s22,s33,s13,s23
       integer i,j,k
            
!           write(6,*) cs_h, cs_v; stop
            ls =   sqrt(dx * dy)
            coef = ismag_h * (cs_h * ls)**2 
            eps30 = 1.0e-30

            ss = 0

!$OMP PARALLEL  DO

       do   k =  1 , nz 
       do   j =  1 , ny 
       do   i =  1 , nx 

!2014 revision

              s11 =  ( um(i,j,k) - um(i-1,j,k) ) / dx 
              s22 =  ( vm(i,j,k) - vm(i,j-1,k) ) / dy
! this is du/dy + dv/dx calculated over 2 grid cells and averaged over 2 columns/rows
!              s12 =  0.25 * ( um(i,j+1,k) + um(i-1,j+1,k) - um(i,j-1,k) - um(i-1,j-1,k)  ) / dy
!     &             + 0.25 * ( vm(i+1,j,k) + vm(i+1,j-1,k) - vm(i-1,j,k) - vm(i-1,j-1,k)  ) / dx 
              s33 = ( wm(i,j,k) - wm(i,j, k-1) ) / dz 

! this is (du/dy + dv/dx)/2 calculated over 2 grid cells and averaged over 2 columns/rows
              s12 =  0.25 * ( ( um(i,j+1,k) + um(i-1,j+1,k) - um(i,j-1,k) - um(i-1,j-1,k)  ) / dy
     &                      + ( vm(i+1,j,k) + vm(i+1,j-1,k) - vm(i-1,j,k) - vm(i-1,j-1,k)  ) / dx ) / 2.0d0

              ss(i,j,k) = 2.0d0*(s11**2) + 2.0d0*(s22**2) + 4.0d0*(s12**2) + 2.0d0*(s33**2) 

!             ss(i,j,k) = ( ( um(i,j,k) - um(i-1,j,k) ) / dx ) **2  + ( ( vm(i,j,k) - vm(i,j-1,k) ) / dy ) **2
!     &                + ( ( um(i,j+1,k) + um(i-1,j+1,k) - um(i,j-1,k) - um(i-1,j-1,k)  ) / dy / 4.0d0  
!     &                +   ( vm(i+1,j,k) + vm(i+1,j-1,k) - vm(i-1,j,k) - vm(i-1,j-1,k)  ) / dx / 4.0d0 ) **2 / 2.0d0  
!     &                +  ( ( wm(i,j+1,k) + wm(i,j+1,k-1) - wm(i,j-1,k) - wm(i,j-1,k-1)  ) / dy / 4.0d0 ) **2 / 2.0d0  
!     &                +  ( ( wm(i+1,j,k) + wm(i+1,j,k-1) - wm(i-1,j,k) - wm(i-1,j,k-1)  ) / dx / 4.0d0) **2 / 2.0d0  

!old             epsh(i,j,k) =  (2.0d0*ss(i,j,k) )
!old              ss(i,j,k ) =  coef * sqrt(2.0d0* ss(i,j,k) + eps30 )
!                epsh(i,j,k) = (epsh(i,j,k) * (ss(i,j,k) + axmin)  !epsh now = AH*S_H^2

              epsh(i,j,k) = ss(i,j,k)  !store S_H^2
              ss(i,j,k ) =  coef * sqrt( ss(i,j,k) + eps30)   !ss now =  mixing coefficient AH without axmin
              epsh(i,j,k) = epsh(i,j,k) * (ss(i,j,k) + axmin)  ! epsh now = S_H^2*A_H

!2014 Revision

       enddo;enddo;enddo

!$OMP END PARALLEL DO

            !call stat_3d(epsh,0, nx+1,0,ny+1,0,nz+1)

            ls =   dz
          coef = ismag_v * (cs_v * ls)**2

!$OMP PARALLEL  DO

      do k =  1 , nz
      do j =  1 , ny
      do i =  1 , nx

! 2014 revision

! this is (du/dz + dw/dx)/2 calculated over 2 grid cells and averaged over 2 columns/rows
         S13 = 0.25 * ( ( um(i,j,k-1) + um(i-1,j,k-1) - um(i,j,k+1) - um(i-1,j,k+1) ) / dz
     &             +    ( wm(i+1,j,k) + wm(i+1,j,k-1) - wm(i-1,j,k) - wm(i-1,j,k-1) ) / dx )/2.0d0 
 
! this is (dv/dz + dw/dy)/2 calculated over 2 grid cells and averaged over 2 columns/rows
         S23 = 0.25 * ( ( vm(i,j,k-1) + vm(i,j-1,k-1) - vm(i,j,k+1) - vm(i,j-1,k+1) ) / dz
     &            +     ( wm(i,j+1,k) + wm(i,j+1,k-1) - wm(i,j-1,k) - wm(i,j-1,k-1) ) / dy )/2.0d0

         ss_v(i,j,k) =  4.0d0*(s13**2) + 4.0d0*(s23**2)

!     &         + ( ( wm(i,j,k) - wm(i,j, k-1) ) / dz ) **2

!     &         +  (( (vm(i,j,k-1) + vm(i,j-1,k-1) + vm(i,j,k) + vm(i,j-1,k))
!     &             - (vm(i,j,k) + vm(i,j-1,k) + vm(i,j,k+1) + vm(i,j-1,k+1) ) ) / dz / 4.0
!     &            ) **2 / 2.0

!     &         +  (( (um(i,j,k-1) + um(i-1,j,k-1)+ um(i,j,k) + um(i-1,j,k))
!     &              - (um(i,j,k) + um(i-1,j,k) + um(i,j,k+1) + um(i-1,j,k+1) )  ) / dz / 4.0
!     &            ) **2 / 2.0

! this the density gradient calculated over 2 grid cells
         Ri(i,j,k) = - g * ( ( rhop(i,j,k-1) + rho_bkg(k-1) ) - ( rhop(i,j,k+1) + rho_bkg(k+1) ) ) / 2.0d0
     &            / dz / rho0 / (ss_v(i,j,k) + eps30)

         Rf = 3  ! put as Prandtl number for now

        epsv(i,j,k)  =  ss_v(i,j,k) 
!        ss_v(i,j,k) =  coef * sqrt( ss_v(i,j,k) + eps30 )* sqrt ( 1.0 + eps - 0.167 )

!  !Di Iorio Eq. 13:
!            epsv(i,j,k)  =   ss_v(i,j,k)    !here epsv temporarily represents S_V^2 
!  !Di Iorio Eq. 8:
            ss_v(i,j,k) =  coef * sqrt( ss_v(i,j,k)+eps30 )*csqrt( 1.0 + eps - Ri(i,j,k)/Rf ) !now ss_v is vertical mixing coef without azmin
!  !Di Iorio Eq. 13: dissipation = S_V^2 * A_V
            epsv(i,j,k)  =  epsv(i,j,k) * (ss_v(i,j,k)+azmin)  !Di Iorio Eq. 13

!           epsv(i,j,k)  =  (coef * sqrt ( 2.0d0 * ss_v(i,j,k) ) * sqrt ( 1.0 + eps - Rf ) + azmin ) * ss_v(i,j,k) 

!2104 revision

      enddo;enddo;enddo

!$OMP END PARALLEL DO

                do k = 1, nz 
                do j = 1, ny

                   ss(0,j,k)   = ss(1,j,k) * bcondx + ( 1 - bcondx ) * ss(nx,j,k)
                   ss(nx+1,j,k)= ss(nx,j,k) * bcondx + ( 1 - bcondx ) * ss(1,j,k)
 
                enddo;enddo 

                do   k = 1, nz
                do   i = 1 , nx
       
                   ss(i,0,k) = ss(i,1,k) * bcondy + ( 1 - bcondy ) * ss(i,ny,k)
                   ss(i,ny+1,k) = ss(i,ny,k)  * bcondy + ( 1 - bcondy ) * ss(i,1,k)

                enddo;enddo

                do  j = 1, ny 
                do  i = 1, nx

                   ss_v(i,j,0) = ss_v(i,j,1)
                   ss_v(i,j,nz+1) = ss_v(i,j,nz)                                            

                enddo;enddo

                return

                end
