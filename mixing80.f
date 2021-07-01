       subroutine mixing  !actuall diffusivity

       use geometry_module!, only: mp_, np_, hzp
       use grid_module!, only: mtime,dx,dy,dz
       use signals_module!, only:  bit5k, bit2k
       use viscosity_module!, only:  axmin, aymin, prandtl, azzmin, ss_v, ss, chi_v, chi_h

       use tracer_module!  , only:   kxt,kyt,kzt,t,t_bkg
       include 'array_size.inc'
       real tprime(0:nx+1,0:ny+1,0:nz+1)

          integer i,j,k

          kxt = 0.0d0
          kyt = 0.0d0
          kzt = 0.0d0

!         write(6,*) 'ss'
          
!$OMP PARALLEL DO

          do k = 1, nz    !-1   ! set to nz-1 so that artifically large horizontal dispersion won't occur
                           ! at seafloor just be cause of large vertical shear there; alternative is to split smag formula
                           ! into horizontal and vertical production terms
                           ! or change bottom boundary so partial slip and not so much vertical shear
                           ! or modulate vertical shear term in a way that boudnary layer doesn't make large bottom
                           ! horiztonal mixing
          do j = 0, ny+1
          do i = 0, nx
                kxt(i,j,k) = bit5k * (( .5d0 * ( ss(i,j,k) + ss(i+1,j,k) ) ) / prandtl + axmin / prandtl)
          enddo;enddo;enddo

!$OMP END PARALLEL DO

!$OMP PARALLEL DO

          do k = 1, nz!-1
          do j = 0, ny
          do i = 0, nx+1
                kyt(i,j,k) = bit5k * (( .5d0 * ( ss(i,j,k) + ss(i,j+1,k) ) ) / prandtl + aymin / prandtl)
          enddo;enddo;enddo

!$OMP END PARALLEL DO

                kzt = 0.0
!$OMP PARALLEL DO

          do k = 1, nz-1 
          do j = 0, ny+1
          do i = 0, nx+1
                kzt(i,j,k) = bit2k * ((0.5d0 * ( ss_v(i,j,k) + ss_v(i,j,k+1) ) ) / prandtl
     &                  +  azzmin(k) / prandtl)
           enddo;enddo;enddo

!$OMP END PARALLEL DO


!$OMP PARALLEL DO

            do k = 0, nz+1
            do j = 0, ny+1
            do i = 0, nx+1
               tprime(i,j,k) = t(i,j,k)- t_bkg(k)
            enddo;enddo;enddo

!$OMP END PARALLEL DO

               chi_v = 0; chi_h = 0

!$OMP PARALLEL DO
! Aug 2014 corrected the calculation of chi with /2dz /2dx /2dy since difference is over 2 grid cells
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
               chi_v(i,j,k) = kzt(i,j,k)* ( (tprime(i,j,k-1) - tprime(i,j,k+1))/ (hzp(k)*2.0d0* dz))**2 
            enddo;enddo;enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
               chi_h(i,j,k) = kxt(i,j,k)* ( (tprime(i+1,j,k) - tprime(i-1,j,k))* mp_(i,j)/ (2.0d0*dx) )**2
     &                      + kyt(i,j,k)* ( (tprime(i,j+1,k) - tprime(i,j-1,k))* np_(i,j)/ (2.0d0*dy) )**2 
            enddo;enddo;enddo

!$OMP END PARALLEL DO


      end subroutine mixing
