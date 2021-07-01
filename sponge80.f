!on 3-4-08, changed minimum ncut from 12 to 10.

       subroutine sponge_maker

!        descendant of sponge_maker_xy_revised.f (work velocity replaces word transport)
     
         use signals_module,    only: bcondx, bcondy, tau, sponge_factor, spongex, spongey, spongee
         use velocity_module,   only: alphae, alphau, alphav, alphaw, alphaee        
         use grid_module, only: dt
         
         real missing_value, beta
         include 'array_size.inc'

         integer ncut, i, j,k
         integer, save:: ict
         !integer spongex,spongey,spongez,spongee
         integer spongez
         data ict/0/

         missing_value = -1.0d+34; beta = 0.67d0 
         beta=0.95

         ! write(6,*) "sponge_factor-beta", beta

          spongex=1; spongey=1; spongez=1; spongee = 1

          !if (bcondx .eq.0) spongex=1
          !if (bcondy .eq.0) spongey=0

          ncut = min(8, nint (.10 * real(nx)))

          alphau = 0

        do j = 0, ny +1
            do i = -1, ncut
              alphau(i,j) = beta*real(ncut-i)/real(ncut+1)
            enddo
            do i = nx + 1 - ncut -1,  nx+1
              alphau(i,j) = beta * real( i - ( nx+1 - ncut -1) ) / real (ncut +1)
            enddo 
        enddo

       if (bcondx .eq. 0) then
          do j = 0, ny +1
              alphau(-1,j) = alphau(nx-1,j)
              alphau(nx+1,j) = alphau(1,j)
          enddo
       endif

         alphau = alphau**2 

!         write(6,*) 'alphau with bcondx = ', bcondx
!         call view_2d( alphau, -1, nx +1, 0, ny+1)

         alphav = 0
                                                                                                                            
         ncut = min(8, nint (.10*real(ny)))
         do i =  0, nx +1
           do j = -1, ncut
              alphav(i,j) = beta * real(ncut-j) / real(ncut+1)
           enddo

           do j = ny + 1 - ncut -1,  ny+1
              alphav(i,j) = beta * real( j - ( ny+1 - ncut -1) ) / real (ncut +1)
           enddo
         enddo

       if (bcondy .eq. 0) then

          do i = 0, nx +1
             alphav(i,-1) = alphav(i, ny - 1)
             alphav(i, ny + 1) = alphav(i,1)
          enddo

!       elseif (bcondy .eq. 1) then

!            do i = 0, nx+1
!              alphav(i,0)    = alphav(i,1)
!              alphav(i,-1)   = alphav(i,1)
!              alphav(i,ny)   = alphav(i,ny-1)
!              alphav(i,ny+1) = alphav(i,ny-1)
!            enddo

!       else
!           write(6,*) 'bcondy must be either 0 or 1'
!           stop 'bcondy violation'

        endif

         alphav = alphav**2

         ncut = min(10,nint (.20 * real(nz)))
         alphaw = 0
            do k = -1, ncut
              alphaw(k) = beta *real(ncut-k)/real(ncut)
           enddo
         alphaw = spongez*alphaw**2   !this is lid sponge

         alphau = sponge_factor / dt * alphau       !+ rayleigh
         alphav = sponge_factor / dt * alphav       !+ rayleigh
         alphaw = sponge_factor / dt * alphaw

         alphae = missing_value
         alphau = spongex * alphau
         alphav = spongey * alphav

         do j = 0,  ny+1 
            do i = 0,  nx+1 
!                alphae(i,j) = (alphau(i,j) + alphau(i-1,j) + alphav(i,j) + alphav(i,j-1) )/2.0d0 !+ rfriction
                 alphae(i,j) = max((alphau(i,j) + alphau(i-1,j)),  (alphav(i,j) + alphav(i,j-1)) )/2.0d0
            enddo
         enddo
!but lid on alphae 
              
           do j = 0, ny+1
           do i = 0, nx+1
           do k = 0, nz+1
                alphaee(i,j,k) =  max(alphae(i,j), (alphaw(k)+ alphaw(k-1))/2.0)
           enddo;enddo;enddo
             !call view_3d(alphaee,0,nx+1,0,ny+1,0,nz+1)
             !stop 88888

             alphaee = spongee * alphaee    !this is a 3-d sponge

             !alphaw_perimeter = spongez * alphae
             alphae = spongee * alphae
             
       if (bcondx .eq. 0) then
          do j = 0, ny +1
              alphae(0,j)      = alphae(nx,j)
              alphae(nx+1,j)   = alphae(1,j)
          enddo;endif

       if (bcondy .eq. 0) then
          do i = 0, nx +1
             alphae(i,0)       = alphae(i,ny )
             alphae(i,ny+1)    = alphae(i,1)
          enddo;endif


         return
         end 
