!copy of mpdata71-hotcross_spongee.f
!CONVECTIVE ADJUSTMENT, if need, should be in gaia, not here.
!uc, etc reasssignments just before smol loop endif make iord_3 possible
!fully omp, without excess writes. 3-12-12
 
!mpdata5ababyyy.f uses fully implict at switch3

!hzw and differencing in x and y for bit5
  
!2-4-08 changed bd conditions so gradient ot tracer is 0.
!daughter of tracer_smolarkiewicz33gg.f

                                                                                                                               
!QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ

        subroutine tracer31 ( c, q_c, c_bkg, cmin ,cmax, c_total0, c_end, ucc, vcc, wcc, tau3, output_unit_number )
                                                  
!       constraints:

                  !  3 - dimensional upstream corrected transport algorithm
                  !  requires mass conserving input velocity field
                  !  if Q_W i.e. water mass is not zero must add term to correction velocities

!       conventions:

                  !  c is the tracer field 
                  !  q_c is the source intensity array
                  !  cmin is the minimum value allowed for the field
                  !  u, v, w are the advection velocities
                  !  kxt, kyt, kzt are the diffusivity arrays 
                  !  tau3 is the time step (in seconds) of integration
!       tests:
!                 !conservation without advection or diffusion
!                 !advection only - uniform c
!                 !advection only spot c      
                                                                                                                        
       use tracer_smol_module
       use constants_module
       use signals_module
       use geometry_module
       use grid_module,  only: dx, dy, ds, modulo_x, modulo_y, mtime 
       use tracer_module, only: kxt,kyt,kzt 
       use velocity_module, only: alphaee, alphae, alphau, alphav, alphaw, u, v, w, um, vm, wm
       use viscosity_module

       include 'array_size.inc'                       
                  
       real  c(0:nx+1,0:ny+1,0:nz+1)
       real, intent(in)   :: q_c(0:nx+1,0:ny+1,0:nz+1)
       real, intent(in)   :: ucc(-1:nx+1,0:ny+1,0:nz+1), vcc(0:nx+1,-1:ny+1,0:nz+1), wcc(0:nx+1,0:ny+1,-1:nz+1)
       real, intent(in)   :: c_bkg(0:nz+1), cmin ,cmax, tau3
       real, intent(inout):: c_total0
       integer  c_end(0:nx+1,0:ny+1), output_unit_number
 
       real worku(-1:nx+1,0:ny+1,0:nz+1), workv(0:nx+1,-1:ny+1,0:nz+1), workw(0:nx+1,0:ny+1,-1:nz+1)

       real block_a, block_b, block_c
       real ce(0:nx+1,0:ny+1,0:nz+1)

       real xfac, yfac, zfac, ca, cb 

       real sum, sum_sponge, sum_source, c_total

       real epsa
       integer switch1, switch2, switch3, nox, noy, iord          

       integer i, j, k, iord_number
       integer, save::  icta

!       allocate ( uc(-1:nx+1,0:ny+1,0:nz+1), vc(0:nx+1,-1:ny+1,0:nz+1), wc(0:nx+1,0:ny+1,-1:nz+1))
!       allocate ( cp(0:nx+1,0:ny+1,0:nz+1) , cc(0:nx+1,0:ny+1,0:nz+1))                        
!       allocate ( ud(-1:nx+1,0:ny+1,0:nz+1), vd(0:nx+1,-1:ny+1,0:nz+1), wd(0:nx+1,0:ny+1,-1:nz+1))

!QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
       data icta/0/
 
        uc = 0;vc=0;wc=0;cp=0;cc=0;ud=0;vd=0;wd=0

        epsa = 1.0d-30

        iord_number = 3     !3 gives two correction steps- moved to *case

        switch1 = 1  !not used -advection
        switch2 = 1  !diffusion
        switch3 = 1  !sponge

        ismol_v =  1
        ismol_hx = 1
        ismol_hy = 1
        nox = 1; noy =1

!        uc = missing_value; vc = missing_value; wc = missing_value
 
!$OMP PARALLEL WORKSHARE

        uc(0:nx,0:ny+1,0:nz+1) = hznx(0:nx,0:ny+1,0:nz+1) * ucc(0:nx,0:ny+1,0:nz+1) * tau3 / dx * 
     &     switch1*bit_z_hx(0:nx,0:ny+1,0:nz+1)        
        vc(0:nx+1,0:ny,0:nz+1) = hzmy(0:nx+1,0:ny,0:nz+1) * vcc(0:nx+1,0:ny,0:nz+1) * tau3 / dy * 
     &     switch1*bit_z_hy(0:nx+1,0:ny,0:nz+1)     
        wc(0:nx+1,0:ny+1,0:nz) = -wcc(0:nx+1,0:ny+1,0:nz) / mp_k(0:nx+1,0:ny+1,0:nz) 
     &                                                         / np_k(0:nx+1,0:ny+1,0:nz) * tau3 / ds * 
     &     switch1 * bit_z_hz(0:nx+1,0:ny+1,0:nz) 

!$OMP END PARALLEL WORKSHARE

         cp = c 

!END LOOP NEW_VEL

! )))))))))))))))))))))))))))

!               massdiv = 0.0

!!$OMP PARALLEL DO

!          do  k = 1, nz
!          do  j = 1, ny
!          do  i = 1, nx

!             massdiv(i,j,k) =
                       !out minus in
!     &                 ( hznx(i,j,k) * ucc(i,j,k) * bit_z_hx(i,j,k) - hznx(i-1,j,k) * ucc(i-1,j,k) * bit_z_hx(i-1,j,k) ) 
!     &               + ( hzmy(i,j,k) * vcc(i,j,k) * bit_z_hy(i,j,k) - hzmy(i,j-1,k) * vcc(i,j-1,k) * bit_z_hy(i,j-1,k) ) 
!     &               + ( wcc(i,j,k-1) * bit_z_hz(i,j,k-1) - wcc(i,j,k) * bit_z_hz(i,j,k) ) / mp_(i,j) / np_(i,j)

! for test case  bit_z_hz = 1 everywhere

!          enddo;enddo;enddo

!!$OMP END PARALLEL DO

!            idum = setvbuf3f(804,2,100)
!            write(804,*) maxval( massdiv), minval(massdiv)
  
        xfac = nox * bit5k  / dx / dx * tau3 * switch2   
        yfac = noy * bit5k  / dy / dy * tau3 * switch2         
        zfac = bit2k  / ds / ds * tau3 * switch2

        !CALL MIXING

!$OMP PARALLEL DO 

         do  k = 1, nz                                 
         do  j = 1, ny                              
         do  i = 1, nx

             c(i,j,k) = c(i,j,k) + (

     &          + xfac * ( hznx(i,j,k)  *  kxt(i,j,k)   * bit_z_hx(i,j,k) * mpx_(i,j)   * ( cp(i+1,j,k) - cp(i,j,k) )
     &                 -  hznx(i-1,j,k) *  kxt(i-1,j,k) * bit_z_hx(i-1,j,k) * mpx_(i-1,j) * ( cp(i,j,k) - cp(i-1,j,k) )  )

     &          + yfac * ( hzmy(i,j,k)   * kyt(i,j,k) * bit_z_hy(i,j,k)  * npy_(i,j)   * ( cp(i,j+1,k) - cp(i,j,k) )
     &                   - hzmy(i,j-1,k) * kyt(i,j-1,k) * bit_z_hy(i,j-1,k) *  npy_(i,j-1) * ( cp(i,j,k) - cp(i,j-1,k) )  )

     &          + zfac * ( kzt(i,j,k-1) / hzw(k-1)* bit_z_hz(i,j,k-1) *
!    &                        ( (cp(i,j,k-1) -  c_bkg(k-1)) - (cp(i,j,k)- c_bkg(k))   )
     &                        ( (cp(i,j,k-1) +  c_bkg(k)) - (c_bkg(k-1) + cp(i,j,k))   )
     &                     -  kzt(i,j,k) / hzw(k)* bit_z_hz(i,j,k) *
!    &                        ( (cp(i,j,k) -  c_bkg(k)) - (cp(i,j,k+1)-  c_bkg(k+1))  )
     &                        ( (cp(i,j,k) +  c_bkg(k+1)) - (c_bkg(k) + cp(i,j,k+1))  )

     &                    ) / mp_(i,j) / np_(i,j)
                                                             !or look to Lavelle, 1997, to beat the slow change 
!in background profile

     &           ) / hzmn(i,j,k) * bit_z_hp(i,j,k) 


         ! CAUTION : a non-zero kzt at z = 0, z = p_end and non zero gradient of c_bkg/c at p_end could lead to 
         ! unphysical growth/decay in C load

          enddo;enddo;enddo

!$OMP END PARALLEL DO

         iord = 1 

!        ce is the change in c caused by advection

         call mpdata (ce, cp,  uc,  vc,  wc,  hzmn, iord, cmin, cmax, tau3)   

!$OMP PARALLEL WORKSHARE

          c = (c + ce) * bit_z_hp 

!$OMP END PARALLEL WORKSHARE


         c = c *  bit_z_hp   !NEW!NEW 

         CALL BOUNDARY_CONDITIONS (c,c_end,c_bkg)


! --- START SMOLARKIEWICZ -----------------------------------

       if (ismol .eq. 1) then

          do iord = 2, iord_number
          

!$OMP PARALLEL WORKSHARE

          cc = c * bit_z_hp           !NEW on 2-15-2012 


       ! --- ANTI DIFFUSION VELOCITY IN THE X DIRECTION----------------

!!ck          ud = missing_value
!!ck          vd = missing_value
!!ck          wd = missing_value

!$OMP END PARALLEL WORKSHARE

!$OMP PARALLEL PRIVATE (ca,cb,block_a,block_b) 

!$OMP DO

          do k = 1, nz 
          do j = 1, ny 
          do i = 0, nx           

               ca =  cc(i+1,j,k)
               cb =  cc(i,j,k)

!              block_c =  bit_z_hy(i,j,k)* bit_z_hy(i+1,j,k)*bit_z_hy(i,j-1,k)* bit_z_hy(i+1,j-1,k)
!     &                    * bit_z_hz(i,j,k)* bit_z_hz(i+1,j,k) * bit_z_hx(i,j,k+1)

               ud(i,j,k) = ( abs ( uc(i,j,k) )  -  uc(i,j,k) * uc(i,j,k)  / hzmnx(i,j,k) ) 
     &                     * (ca - cb) / (ca + cb + epsa) !* block_c  

               ca =  (cc(i+1,j+1,k) + cc(i,j+1,k) ) 
               cb =  (cc(i+1,j-1,k) + cc(i,j-1,k) )

               block_a = bit_z_hy(i+1,j,k) * bit_z_hy(i,j,k) * bit_z_hy(i+1,j-1,k) * bit_z_hy(i,j-1,k) 

               ud(i,j,k) =  ud(i,j,k) - noy * 0.5 * uc(i,j,k) / hzmnx(i,j,k) * (ca - cb) / ( ca + cb + epsa) 
     &            * 0.25 * ( vc(i+1,j,k) + vc(i,j,k) + vc(i+1,j-1,k) + vc(i,j-1,k) ) * block_a

               ca =  ( cc(i+1,j,k-1) + cc(i,j,k-1))
               cb =  ( cc(i+1,j,k+1) + cc(i,j,k+1))

               block_b =  bit_z_hz(i+1,j,k) * bit_z_hz(i,j,k) * bit_z_hz(i+1,j,k-1) * bit_z_hz(i,j,k-1)  

               ud(i,j,k) = ud(i,j,k) - 0.5 * uc(i,j,k) / hzmnx(i,j,k) * (ca - cb) / (ca + cb + epsa) 
     &            * 0.25 * ( wc(i+1,j,k) + wc(i,j,k) + wc(i+1,j,k-1) + wc(i,j,k-1) ) * block_b

               ud(i,j,k) = nox * ismol_hx * ud(i,j,k) * bit_z_hx(i,j,k) 

  10      continue

          enddo; enddo; enddo

!$OMP END DO 
!!NOWAIT

! --- ANTI DIFFUSION VELOCITY IN THE Y DIRECTION----------------

!$OMP DO

          do  k = 1, nz 
          do  j = 0, ny           
          do  i = 1, nx 

                ca = cc(i,j+1,k)
                cb = cc(i,j,k)
 
                vd(i,j,k) =  (abs ( vc(i,j,k) )  -  vc(i,j,k) * vc(i,j,k) / hzmny(i,j,k) )  
     &                          * noy * (ca - cb) / ( ca + cb + epsa)


                ca =  (cc(i+1,j+1,k) + cc(i+1,j,k))
                cb =  (cc(i-1,j+1,k) + cc(i-1,j,k))

                block_a = bit_z_hx(i,j,k) * bit_z_hx(i,j+1,k) * bit_z_hx(i-1,j,k) *  bit_z_hx(i-1,j+1,k) 

                vd(i,j,k) =  vd(i,j,k) - 0.5d0 * vc(i,j,k) / hzmny(i,j,k) * (ca - cb) / (ca + cb + epsa) 
     &             * 0.25d0 * ( uc(i,j+1,k) + uc(i,j,k) + uc(i-1,j+1,k) + uc(i-1,j,k) ) * block_a
               
                ca = ( cc(i,j+1,k-1) + cc(i,j,k-1))
                cb = ( cc(i,j+1,k+1) + cc(i,j,k+1))

                block_b =  bit_z_hz(i,j,k) *  bit_z_hz(i,j+1,k) *  bit_z_hz(i,j,k-1) * bit_z_hz(i,j+1,k-1)

                vd(i,j,k) =  vd(i,j,k) - 0.5d0 * vc(i,j,k) / hzmny(i,j,k) * (ca - cb) / (ca + cb + epsa)  
     &              * 0.25d0* ( wc(i,j,k) + wc(i,j+1,k) + wc(i,j,k-1) + wc(i,j+1,k-1) ) * block_b

                vd(i,j,k) = noy * ismol_hy * vd(i,j,k) * bit_z_hy(i,j,k) 

          enddo; enddo; enddo

!$OMP END DO 
!!NOWAIT

! -----ANTI DIFFUSION VELOCITY IN THE Z DIRECTION-----------------------

!$OMP DO

          do  k = 0, nz             
          do  j = 1, ny  
          do  i = 1, nx 

               ca = cc(i,j,k )
               cb = cc(i,j,k+1) 
               
               wd(i,j,k) = ( abs ( wc(i,j,k) ) -  wc(i,j,k) * wc(i,j,k) /  hzmn(i,j,k) ) *  (ca - cb) / 
     &                   ( ca + cb + epsa) 

               ca = ( cc(i+1,j,k) + cc(i+1,j,k+1) ) 
               cb = ( cc(i-1,j,k) + cc(i-1,j,k+1) )

               block_a = bit_z_hx(i,j,k) * bit_z_hx(i,j,k+1) *  bit_z_hx(i-1,j,k) * bit_z_hx(i-1,j,k+1)

               wd(i,j,k) =  wd(i,j,k) - nox * 0.5 * wc(i,j,k)  /  hzmn(i,j,k) * (ca - cb) / ( ca + cb + epsa)  
     &               *  0.25 *   ( uc(i,j,k) + uc(i,j,k+1) +  uc(i-1,j,k+1) + uc(i-1,j,k) ) * block_a

               ca =  ( cc(i,j+1,k+1) + cc(i,j+1,k) )
               cb =  ( cc(i,j-1,k+1) + cc(i,j-1,k) )

               block_b =  bit_z_hy(i,j,k) *  bit_z_hy(i,j,k+1) *   bit_z_hy(i,j-1,k) * bit_z_hy(i,j-1,k+1)          

               wd(i,j,k) =  wd(i,j,k) - noy * 0.5 * wc(i,j,k) /  hzmn(i,j,k) * (ca - cb) / ( ca + cb + epsa) 
     &              * 0.25  *  ( vc(i,j,k) + vc(i,j,k+1) +  vc(i,j-1,k+1) + vc(i,j-1,k) ) * block_b

               wd(i,j,k) =       ismol_v  * wd(i,j,k) * bit_z_hz(i,j,k) 

          enddo; enddo; enddo

!$OMP END DO
!$OMP END PARALLEL

! ce is the change caused by advection

              call mpdata (ce, cc, ud, vd, wd, hzmn, iord, cmin, cmax, tau3 )

              c = (c + ce) * bit_z_hp
    
              CALL BOUNDARY_CONDITIONS(c, c_end,c_bkg)

              !write(6,*) 'in corrector, iord =', iord

!$OMP PARALLEL WORKSHARE
              uc = ud
              vc = vd
              wc = wd
!$OMP END PARALLEL WORKSHARE
       enddo     !ends the iord loop



       endif
              

       ! -------  END SMOLARKIEWICZ ------------------------------------

       ! )))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
!$OMP PARALLEL  DO 

         do  k = 1, nz
         do  j = 1, ny
         do  i = 1, nx

         ! sponge the values back to background on the perimeter
         ! semi-implict with sponge back in surface layer 

!                     c(i,j,k) =  (c(i,j,k) - switch3 * alphae(i,j) * tau3 * (0* cp(i,j,k)  -  c_bkg(k) ) 
!1OLD    &                                    -  switch3 * alphaz(k)   * tau3 * ( cp(i,j,k) - c_bkg(k) ) 
!     &                          ) / (1.0d0 +  switch3 * alphae(i,j) * tau3 )  !fully implicit

                      c(i,j,k) =  c(i,j,k) - switch3 * 0.5d0 * alphaee(i,j,k) * tau3 * ( cp(i,j,k) - c_bkg(k) )


!!!  DOES SMOLARKIEWICZ TERM NEED CORRECTION BECAUSE OF THIS SPONGE??  Would only be corrected in sponge region where doesnt 
!!!  matter.
 
                 c(i,j,k) = c(i,j,k) 
     &                      + tau3  *  q_c(i,j,k) / (dx * dy * ds ) / hzmn(i,j,k) * bit_z_hp(i,j,k)
!    &                      + tau3  *  cp(i,j,k) * massdiv(i,j,k) / (dx * dy * ds ) / hzmn(i,j,k) * bit_z_hp(i,j,k)

         enddo;enddo;enddo

!$OMP END PARALLEL DO


            if(mod(mtime,100).eq. 0) then                                                                                        

            sum = 0.0d0
                                                                                                           
          do  k = 1, nz
          do  j = 1, ny
          do  i = 1, nx
                                                                                                                               
             sum = sum +  ( c(i,j,k) - c_bkg(k) ) * hzmn(i,j,k) * bit_z_hp(i,j,k)
                                                                                                                         
          enddo;enddo;enddo

             sum = sum * dx * dy * ds

             c_total = sum

             sum = 0.0d0
                                                                                                                              
! change in C over time period
                                                                                                                    
          do  k = 1, nz
          do  j = 1, ny
          do  i = 1, nx
                                                                                                                      
             sum = sum - switch3 * alphae(i,j) * ( c(i,j,k) - c_bkg(k) ) * hzmn(i,j,k) * bit_z_hp(i,j,k)
                                                                                                                             
          enddo;enddo;enddo

             sum_sponge = sum  * dx * dy * ds * tau3 

             sum = 0.0d0
                                                    
          do  k = 1, nz
          do  j = 1, ny
          do  i = 1, nx 

               sum = sum + tau3 * q_c(i,j,k) * bit_z_hp(i,j,k)
         
          enddo;enddo;enddo

             sum_source = sum 

!          write(output_unit_number, '(4d20.10,i8)')   c_total - c_total0 - sum_sponge - sum_source,  
!     &                c_total, sum_sponge, sum_source, mtime
          endif  
                                                                                                              
          c_total0 = c_total
           

       end subroutine tracer31
      
       subroutine boundary_conditions (c, c_end,c_bkg) 

       use signals_module, only:  bcondx,bcondy,bdbit
       use grid_module  ,  only: dx,dy,ds,  modulo_x, modulo_y
       !use tracer_module,  only:  c_bkg
        include 'array_size.inc'
        integer i,j,k, c_end(0:nx+1,0:ny+1)
       real c(0:nx+1,0:ny+1,0:nz+1), c_bkg(0:nz+1)

          do  j = 1, ny
          do  i = 1, nx           !BOX TOP

              k = 0
              c(i,j,k)   = c(i,j,k+1)                    !no gradient
              k = c_end(i,j)
              c(i,j,k+1) = c(i,j,k)   

          enddo
          enddo

           do k = 0, nz+1
           do j = 0, ny+1                                  !BOX END
                                                                                                                
               c(0,j,k) =  ( 1 - bcondx ) * c(nx,j,k)
     &                + bcondx *  ( 1 - bdbit(1) ) *  c(1,j,k)
     &                + bcondx *  bdbit(1) * c(1,j,k)   !( 2.0 * c_bkg(k) - c(1,j,k))  !gradient of tracer is same as gradient of background
                                                                                                                        
           c(nx+1,j,k) = ( 1 - bcondx ) * c(1,j,k)
     &                + bcondx *  ( 1 - bdbit(2) ) * c(nx,j,k)
     &                + bcondx *  bdbit(2) * c(nx,j,k)  !( 2.0 * c_bkg(k) - c(nx,j,k))
                                                                                                                   
          enddo;enddo

          do  k = 0, nz+1
          do  i = 0, nx+1                                 !BOX FRONT
                                                                                                               
           c(i,0,k) =    ( 1 - bcondy ) * c(i,ny,k)
     &                + bcondy *  ( 1 - bdbit(3) ) * c(i,1,k)
     &                + bcondy * bdbit(3) * c(i,1,k)         !( 2.0 *c_bkg(k) -  c(i,1,k))
           c(i,ny+1,k) = ( 1 - bcondy ) * c(i,1,k)
     &                + bcondy *  ( 1 - bdbit(4) ) * c(i,ny,k)
     &                + bcondy * bdbit(4) * c(i,ny,k)        !( 2.0 *c_bkg(k) -  c(i,ny,k))

          enddo;enddo

             !c(0,   0,   1:nz)    = c_bkg(1:nz) + (1-bcondx) * (1-bcondy) * c(nx,   ny, 1:nz)
             !c(nx+1,ny+1,1:nz)    = c_bkg(1:nz) + (1-bcondx) * (1-bcondy) * c(1,     1, 1:nz)
             !c(0,   ny+1,1:nz)    = c_bkg(1:nz) + (1-bcondx) * (1-bcondy) * c(nx,    1, 1:nz)
             !c(nx+1,0,   1:nz)    = c_bkg(1:nz) + (1-bcondx) * (1-bcondy) * c(1,    ny, 1:nz)

             if ( (modulo_x .eq. 1 .and. modulo_y .eq. 0) .or. (modulo_x .eq. 0 .and. modulo_y .eq. 1)) then

              c(0,   0,   1:nz)    =  modulo_x * c(nx,   0, 1:nz) + modulo_y * c(0,    ny, 1:nz)
              c(0,   ny+1,1:nz)    =  modulo_x * c(nx,ny+1, 1:nz) + modulo_y * c(0,     1, 1:nz)
              c(nx+1,0,   1:nz)    =  modulo_x * c(1,    0, 1:nz) + modulo_y * c(nx+1, ny, 1:nz)
              c(nx+1,ny+1,1:nz)    =  modulo_x * c(1, ny+1, 1:nz) + modulo_y * c(nx+1,  1, 1:nz)

             endif

      end subroutine boundary_conditions


      subroutine mpdata (cb, ca, uca, vca, wca, capG, iord, cmin, cmax, tau3 )

! Here uc,vc, wc = capU, capV, and capW
! ca is k-1 psi value and cb is k psi value, capG=(hzmn) is the 

            use constants_module,   only: missing_value
            include 'array_size.inc'

            real, intent(out)  :: cb(0:nx+1,0:ny+1,0:nz+1)
            real, intent(in)   :: ca(0:nx+1,0:ny+1,0:nz+1), capG(0:nx+1,0:ny+1,0:nz+1), cmin, cmax, tau3
            real, intent(in)   :: uca(-1:nx+1,0:ny+1,0:nz+1), vca(0:nx+1,-1:ny+1,0:nz+1), wca(0:nx+1,0:ny+1,-1:nz+1)

            real worku(-1:nx+1,0:ny+1,0:nz+1), workv(0:nx+1,-1:ny+1,0:nz+1),  workw(0:nx+1,0:ny+1,-1:nz+1), 
     &           massdiv(0:nx+1, 0:ny+1, 0:nz+1)

            real upi(-1:nx+1,  0:ny+1,  0:nz+1), umi(-1:nx+1, 0:ny+1, 0:nz+1),
     &           vpj( 0:nx+1, -1:ny+1,  0:nz+1), vmj( 0:nx+1,-1:ny+1, 0:nz+1),
     &           wpk( 0:nx+1,  0:ny+1, -1:nz+1), wmk( 0:nx+1, 0:ny+1,-1:nz+1)                        

            integer, intent(in) :: iord
            integer i,j,k

!            upi = missing_value
!            umi = missing_value
!            vpj = missing_value
!            vmj = missing_value
!            wpk = missing_value
!            wmk = missing_value


!$OMP PARALLEL DO

          do  k = 1, nz
          do  j = 1, ny
          do  i = 0, nx

            worku(i, j, k) = abs( uca(i, j, k))
            upi(i,j,k) =  0.5d0 * ( uca(i,j,k) + worku(i,j,k) )
            umi(i,j,k) =  0.5d0 * ( uca(i,j,k) - worku(i,j,k) )

            !the implicit form is 11.5 % slower

!           worku(0:nx, 1:ny, 1:nz) = abs( uca(0:nx, 1:ny, 1:nz))
!           upi(0:nx, 1:ny, 1:nz) =  0.5d0 * ( uca(0:nx, 1:ny, 1:nz) + worku(0:nx, 1:ny, 1:nz) )
!           umi(0:nx, 1:ny, 1:nz) =  0.5d0 * ( uca(0:nx, 1:ny, 1:nz) - worku(0:nx, 1:ny, 1:nz) )  

         enddo;enddo;enddo

!$OMP END PARALLEL DO

!$OMP PARALLEL workshare

!          do  k = 1, nz
!          do  j = 0, ny
!          do  i = 1, nx

!            workv(i,j,k) =  abs( vca(i,j,k))
!            vpj(i,j,k)   =  0.5d0 * ( vca(i,j,k) + workv(i,j,k) )
!            vmj(i,j,k)   =  0.5d0 * ( vca(i,j,k) - workv(i,j,k) )

           workv(1:nx, 0:ny, 1:nz) = abs( vca(1:nx, 0:ny, 1:nz))
           vpj(1:nx, 0:ny, 1:nz) =  0.5d0 * ( vca(1:nx, 0:ny, 1:nz) + workv(1:nx, 0:ny, 1:nz) )
           vmj(1:nx, 0:ny, 1:nz) =  0.5d0 * ( vca(1:nx, 0:ny, 1:nz) - workv(1:nx, 0:ny, 1:nz) )

!          enddo;enddo;enddo

!!$OMP END PARALLEL  workshare

!!$OMP workshare

!          do  k = 0, nz
!          do  j = 1, ny
!          do  i = 1, nx

!            workw(i,j,k) =  abs( wca(i,j,k) )
!            wpk(i,j,k)   =  0.5d0 * ( wca(i,j,k) + workw(i,j,k) )
!            wmk(i,j,k)   =  0.5d0 * ( wca(i,j,k) - workw(i,j,k) )

              workw(1:nx, 1:ny, 0:nz) = abs( wca(1:nx, 1:ny, 0:nz) )
              wpk(1:nx, 1:ny, 0:nz) =  0.5d0 * ( wca(1:nx, 1:ny, 0:nz) + workw(1:nx, 1:ny, 0:nz) )
              wmk(1:nx, 1:ny, 0:nz) =  0.5d0 * ( wca(1:nx, 1:ny, 0:nz) - workw(1:nx, 1:ny, 0:nz) )

!          enddo;enddo;enddo

!$OMP END PARALLEL WORKSHARE


!     &                + ( upi(0:nx-1,1:ny,1:nz) * ca(0:nx-1,1:ny,1:nz) - umi(1:nx,1:ny,1:nz)   * ca(2:nx+1,1:ny,1:nz) )
!     &                + ( vpj(1:nx,0:ny-1,1:nz) * ca(1:nx,0:ny-1,1:nz) - vmj(1:nx,1:ny,1:nz)   * ca(1:nx,2:ny+1,1:nz) )
!     &                + ( wpk(1:nx,1:ny,1:nz)   * ca(1:nx,1:ny,2:nz+1) - wmk(1:nx,1:ny,0:nz-1) * ca(1:nx,1:ny,0:nz-1) )

!     &              - ( ( upi(1:nx,1:ny,1:nz)   * ca(1:nx,1:ny,1:nz)   - umi(0:nx-1,1:ny,1:nz) * ca(1:nx,1:ny,1:nz)     )
!     &                + ( vpj(1:nx,1:ny,1:nz)   * ca(1:nx,1:ny,1:nz)   - vmj(1:nx,0:ny-1,1:nz) * ca(1:nx,1:ny,1:nz)   )
!     &                + ( wpk(1:nx,1:ny,0:nz-1) * ca(1:nx,1:ny,1:nz)   - wmk(1:nx,1:ny,1:nz)   * ca(1:nx,1:ny,1:nz)   ) )

!      if ( iord .eq. 1) then


            cb = missing_value   !setting the value(missing or 0) within this routine absolutely required 

!$OMP  PARALLEL WORKSHARE

            cb(1:nx,1:ny,1:nz) = 0*ca(1:nx,1:ny,1:nz) + (

     &                + ( upi(0:nx-1,1:ny,1:nz) * ca(0:nx-1,1:ny,1:nz) + umi(0:nx-1,1:ny,1:nz) * ca(1:nx,1:ny,1:nz)   )
     &                + ( vpj(1:nx,0:ny-1,1:nz) * ca(1:nx,0:ny-1,1:nz) + vmj(1:nx,0:ny-1,1:nz) * ca(1:nx,1:ny,1:nz)   )
     &                + ( wmk(1:nx,1:ny,1:nz)   * ca(1:nx,1:ny,1:nz)   + wpk(1:nx,1:ny,1:nz)   * ca(1:nx,1:ny,2:nz+1) )
     &              - ( ( upi(1:nx,1:ny,1:nz)   * ca(1:nx,1:ny,1:nz)   + umi(1:nx,1:ny,1:nz)   * ca(2:nx+1,1:ny,1:nz) )
     &                + ( vpj(1:nx,1:ny,1:nz)   * ca(1:nx,1:ny,1:nz)   + vmj(1:nx,1:ny,1:nz)   * ca(1:nx,2:ny+1,1:nz) )
     &                + ( wmk(1:nx,1:ny,0:nz-1) * ca(1:nx,1:ny,0:nz-1) + wpk(1:nx,1:ny,0:nz-1) * ca(1:nx,1:ny,1:nz)   ) )

     &                               ) / capG(1:nx,1:ny,1:nz)

!$OMP END  PARALLEL WORKSHARE


!      endif

!           massdiv = 0

!          do  k = 1, nz
!          do  j = 1, ny
!          do  i = 1, nx

!             massdiv(i,j,k) =

!     &                  ( (upi(i,j,k) + umi(i,j,k)) -  (upi(i-1,j,k) + umi(i-1,j,k)) )
!     &                + ( (vpj(i,j,k) + vmj(i,j,k)) -  (vpj(i,j-1,k) + vmj(i,j-1,k)) )
!     &                + ( -(wpk(i,j,k) + wmk(i,j,k)) +  (wpk(i,j,k-1) + wmk(i,j,k-1)) )

!          enddo;enddo;enddo

             !write(6,*) 'massdiv'

             !call view_3d(massdiv(1:nx,1:ny,1:nz),1,nx,1,ny,1,nz)

             if ( iord .gt. 1 ) then

!make cell interfaces impermeable to correction velocities under certain conditions 

!                  where ( cb(1:nx,1:ny,1:nz) .le. cmin )

!                           umi(0:nx-1,1:ny,1:nz) = 0.0d0
!                           upi(1:nx,1:ny,1:nz)   = 0.0d0
!                           vmj(1:nx,0:ny-1,1:nz) = 0.0d0
!                           vpj(1:nx,1:ny,1:nz)   = 0.0d0
!                           wpk(1:nx,1:ny,0:nz-1) = 0.0d0
!                           wmk(1:nx,1:ny,1:nz)   = 0.0d0

!                 endwhere

!                 where ( cb(1:nx,1:ny,1:nz) .ge. cmax )                

!                           upi(0:nx-1,1:ny,1:nz) = 0.0d0
!                           umi(1:nx,1:ny,1:nz)   = 0.0d0
!                           vpj(1:nx,0:ny-1,1:nz) = 0.0d0
!                           vmj(1:nx,1:ny,1:nz)   = 0.0d0
!                           wpk(1:nx,1:ny,1:nz)   = 0.0d0
!                           wmk(1:nx,1:ny,0:nz-1) = 0.0d0

!                 endwhere                                                                                              


!$OMP PARALLEL DO

                do k = 1, nz
                do j = 1, ny
                do i = 1, nx

                    if ( (ca(i,j,k) + cb(i,j,k)) .le. cmin ) then

                           umi(i-1,j,k) = 0.0d0
                           upi(i,j,k)   = 0.0d0
                           vmj(i,j-1,k) = 0.0d0
                           vpj(i,j,k)   = 0.0d0
                           wpk(i,j,k-1) = 0.0d0
                           wmk(i,j,k)   = 0.0d0

                    elseif ( (ca(i,j,k) + cb(i,j,k)) .ge. cmax ) then

                           upi(i-1,j,k) = 0.0d0
                           umi(i,j,k)   = 0.0d0
                           vpj(i,j-1,k) = 0.0d0
                           vmj(i,j,k)   = 0.0d0
                           wpk(i,j,k)   = 0.0d0
                           wmk(i,j,k-1) = 0.0d0

                    endif

                enddo;enddo

               enddo

!$OMP END PARALLEL DO

!$OMP PARALLEL DO 

!             cb(1:nx,1:ny,1:nz) = ca(1:nx,1:ny,1:nz) + (

!     &                + ( upi(0:nx-1,1:ny,1:nz) * ca(0:nx-1,1:ny,1:nz) + umi(0:nx-1,1:ny,1:nz) * ca(1:nx,1:ny,1:nz)   )
!     &                + ( vpj(1:nx,0:ny-1,1:nz) * ca(1:nx,0:ny-1,1:nz) + vmj(1:nx,0:ny-1,1:nz) * ca(1:nx,1:ny,1:nz)   )
!     &                + ( wmk(1:nx,1:ny,1:nz)   * ca(1:nx,1:ny,1:nz)   + wpk(1:nx,1:ny,1:nz)   * ca(1:nx,1:ny,2:nz+1) )
!     &              - ( ( upi(1:nx,1:ny,1:nz)   * ca(1:nx,1:ny,1:nz)   + umi(1:nx,1:ny,1:nz)   * ca(2:nx+1,1:ny,1:nz) )
!     &                + ( vpj(1:nx,1:ny,1:nz)   * ca(1:nx,1:ny,1:nz)   + vmj(1:nx,1:ny,1:nz)   * ca(1:nx,2:ny+1,1:nz) )
!     &                + ( wmk(1:nx,1:ny,0:nz-1) * ca(1:nx,1:ny,0:nz-1) + wpk(1:nx,1:ny,0:nz-1) * ca(1:nx,1:ny,1:nz)   ) )

!     &                               ) / capG(1:nx,1:ny,1:nz) 

         do k = 1, nz
         do j = 1, ny
         do i = 1, nx

             cb(i,j,k) = 0*ca(i,j,k) + (

     &                + ( upi(i-1,j,k) * ca(i-1,j,k) + umi(i-1,j,k) * ca(i,j,k)   )
     &                + ( vpj(i,j-1,k) * ca(i,j-1,k) + vmj(i,j-1,k) * ca(i,j,k)   )
     &                + ( wmk(i,j,k)   * ca(i,j,k)   + wpk(i,j,k)   * ca(i,j,k+1) )
     &              - ( ( upi(i,j,k)   * ca(i,j,k)   + umi(i,j,k)   * ca(i+1,j,k) )
     &                + ( vpj(i,j,k)   * ca(i,j,k)   + vmj(i,j,k)   * ca(i,j+1,k) )
     &                + ( wmk(i,j,k-1) * ca(i,j,k-1) + wpk(i,j,k-1) * ca(i,j,k)   ) )

     &                               ) / capG(i,j,k)

         enddo;enddo;enddo

!$OMP END PARALLEL DO

         endif

      end subroutine mpdata


