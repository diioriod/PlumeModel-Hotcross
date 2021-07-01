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

              ss(i,j,k) = ( ( um(i,j,k) - um(i-1,j,k) ) / dx ) **2  + ( ( vm(i,j,k) - vm(i,j-1,k) ) / dy ) **2

     &                + ( ( um(i,j+1,k) + um(i-1,j+1,k) - um(i,j-1,k) - um(i-1,j-1,k)  ) / dy / 4.0d0  

     &                +   ( vm(i+1,j,k) + vm(i+1,j-1,k) - vm(i-1,j,k) - vm(i-1,j-1,k)  ) / dx / 4.0d0 ) **2 / 2.0d0  

     &                +  ( ( wm(i,j+1,k) + wm(i,j+1,k-1) - wm(i,j-1,k) - wm(i,j-1,k-1)  ) / dy / 4.0d0 ) **2 / 2.0d0  

     &                +  ( ( wm(i+1,j,k) + wm(i+1,j,k-1) - wm(i-1,j,k) - wm(i-1,j,k-1)  ) / dx / 4.0d0) **2 / 2.0d0  

            epsh(i,j,k) =  (2.0 * ss(i,j,k) )
             ss(i,j,k ) = coef * sqrt( 2.0d0 * ss(i,j,k) + eps30 )
            epsh(i,j,k) =  epsh(i,j,k) * ss(i,j,k )

       enddo;enddo;enddo

!$OMP END PARALLEL DO

            !call stat_3d(epsh,0, nx+1,0,ny+1,0,nz+1)

            ls =   dz
          coef = ismag_v * (cs_v * ls)**2

!$OMP PARALLEL  DO

      do k =  1 , nz
      do j =  1 , ny
      do i =  1 , nx

            ss_v(i,j,k) =

     &         + ( ( wm(i,j,k) - wm(i,j, k-1) ) / dz ) **2

     &         +  (( (vm(i,j,k-1) + vm(i,j-1,k-1) + vm(i,j,k) + vm(i,j-1,k))
     &             - (vm(i,j,k) + vm(i,j-1,k) + vm(i,j,k+1) + vm(i,j-1,k+1) ) ) / dz / 4.0
     &            ) **2 / 2.0

     &         +  (( (um(i,j,k-1) + um(i-1,j,k-1)+ um(i,j,k) + um(i-1,j,k))
     &              - (um(i,j,k) + um(i-1,j,k) + um(i,j,k+1) + um(i-1,j,k+1) )  ) / dz / 4.0
     &            ) **2 / 2.0

                Ri(i,j,k) = - g * ( ( rhop(i,j,k-1) + rho_bkg(k-1) ) - ( rhop(i,j,k+1) + rho_bkg(k+1) ) ) / 2.0d0
     &            /  dz / rho0 / ( 2.0 * (ss_v(i,j,k)  + eps30 ) )

                Rf = .16666

            epsv(i,j,k)  =  (2.0 * ss_v(i,j,k) )
            ss_v(i,j,k ) =  coef * sqrt( 2.0d0 * ss_v(i,j,k) + eps30 )* sqrt ( 1.0 + eps - Rf )
            epsv(i,j,k)  =  epsv(i,j,k) * ss_v(i,j,k )

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
