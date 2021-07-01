
! changed to admit the abs(frequencies) higher than (1/10.0) hrs - from the Thomson record

! YOU ARE USING SITE_FORCING:  Thomson's 2001 data = but just mean, semi-diurnal, and < 10hr periods

! then iterate on the size of spectral coefficient needed
! dt = time step in code

        integer np
        parameter (np = 6600)  ! nh = 10*24)       !nh arrives from case file       

        integer, save:: ns, ntt, m, time0
        real, save::   ramp_hrs

        real, save ::  om(1:np), ur(1:np), ui(1:np), vr(1:np), vi(1:np)
        real, save ::  sr(1:np), si(1:np), freq(1:np), period(1:np), weight(1:np), dummy

        real utest,vtest,time_,scaling,factor_c 
        real, allocatable, save:: ut(:),vt(:)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%SETUP%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (mtime .eq. 0) then

            missing_value = -1.0e+34

            time0 = 0.0

            !time0 is defined in case

            ramp_hrs = 0.25 * 24.0    
            ramp_hrs = 2.0

            write(6,*) 'ramp_hrs = ', ramp_hrs

            ns = np/2 - 1
                                                                                               
            if ( mod(np,2) .ne. 0) stop 'np in site_forcing.inc incorrect'
                                                                             

            ntt = nt + nint(-time0/dt)

            write(6,*) 'number of time steps including start up;', ntt
 
            allocate ( ut(-1:ntt+1),vt(-1:ntt+1))

            ut =  0; vt = 0

!construct a forcing time series at the dt increment

                om(1) = 2.0*pi/12.42/3600.0

               !om(1) = 2.0*pi/12.5/3600.0

           do l = 0, ntt+1
                                                                                                                          
              time_ = dt * l + time0
!Commenting out ut and vt to test with no background flow. <TEMPORARY>
!              ut(l)=0
!	      vt(l)=0
              ut(l) = 0.032   + 0.036  * sin(om(1) * time_) 
              vt(l) =           0.0024 * cos(om(1) * time_)

               factor_c = min (1.0,    (l  * dt) / (ramp_hrs * 3600.0)  )
                                                                                                                          
              write(7001,'(4e14.6)') ut(l),vt(l), factor_c * ut(l),factor_c * vt(l)

              if ( mod(nint(time_),3600) .eq. 0) write(7002,'(5e14.6)') ut(l),vt(l), factor_c * ut(l),factor_c * vt(l),time_/3600.0
              if ( mod(nint(time_),300) .eq. 0)  write(7003,'(4e14.6)') ut(l),vt(l),time_,time_/3600.0

              ut(l) =  factor_c * ut(l)
              vt(l) =  factor_c * vt(l)

           enddo

              rewind(7001);rewind(7002)

        endif
!&&&&&&&&&&&& END OF setup &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

         tau = leap * dt                   
 
         time_ = time + time0

!        write(6,*) 'mtime, delta, leap, tau',  mtime, delta, leap, tau

         l = nint(time/dt)

              ubkg_m  = ut(l-1)           !ubkg(-1:nx+1, 0:ny+1),  ubkg_p(-1:nx+1, 0:ny+1),  ubkg_m(-1:nx+1, 0:ny+1)
              ubkg    = ut(l)   
              ubkg_p  = ut(l+1) 

              vbkg_m  = vt(l-1)
              vbkg    = vt(l)
              vbkg_p  = vt(l+1)                                                                                   

              ebkg =0.0d0
              ebkg_p = 0.0d0
              ebkg_m = 0.0d0
                
              !write(701,'(7e14.6,i10)')  ut(l+1), ut(l), ut(l-1), vt(l-1), vt(l), vt(l+1), time_, l

         Fbx = missing_value
         Fby = missing_value
         Fbz = 0.0

!NOTE- THESE FBX AND FBY SHOULD BE THE CREATED AS SOLUTIONS TO TRUNCATED EQUATIONS WITH THE HMN FACTORS
! IN OTHER WORDS, SERVE bODY FORCE TO MOMENTUM.F WITH ITS GEOMETRIC FACTORS IN PLACE
         
           !call stat_2d(alphau*ubkg_m,-1,nx+1,0,ny+1)
           !write(6,*) 'fc=',fc

        do k = 1, nz
        do j = 1, ny
        do i = 0, nx 

           Fbx(i,j,k) =

!     &             alphau(i,j) *  ubkg_m(i,j)

     &          + (ubkg_p(i,j) - (1-delta) * ubkg(i,j) - delta * ubkg_m(i,j)) / tau
                                                                                                                        
     &          - bit7 * fc * ( vbkg(i+1,j) + vbkg(i+1,j-1)  +  vbkg(i,j) + vbkg(i,j-1) ) / 4.0d0


!WARNING - fix terms when non-orthogonal grid

!     &    + bit9 * (
!     &           dndxp_(i+1,j) * ( vbit (i+1,j,k) + vbit(i+1,j-1,k) ) * hzp_(i+1,j,k) * ( vbit (i+1,j,k) + vbit(i+1,j-1,k) )
!     &         +  dndxp_(i,j) * ( vbit (i,j,k) + vbit(i,j-1,k) ) * hzp_(i,j,k) * ( vbit (i,j,k) + vbit(i,j-1,k) )
!     &             ) / 8.0d0
                                                                                                                   
!     &    - bit9 * (
!     &            dmdyp_(i+1,j) * ( ubit (i+1,j,k) + ubit(i,j,k) ) *  hzp_(i+1,j,k) * ( vbit (i+1,j,k) + vbit(i+1,j-1,k) )
!     &        +   dmdyp_(i,j) * ( ubit (i,j,k) + ubit(i-1,j,k) ) *  hzp_(i,j,k) *   ( vbit (i,j,k) + vbit(i,j-1,k)     )
!     &             ) / 8.0d0


        enddo;enddo

      do j = 0, ny
      do i = 1, nx

          Fby(i,j,k) =

!     &         alphav(i,j) * vbkg_m(i,j)

     &     +   (vbkg_p(i,j) - (1-delta) * vbkg(i,j) - delta * vbkg_m(i,j) ) / tau
                                                                                                                                 
     &     + bit7 *  fc * ( ubkg(i,j+1)  + ubkg(i-1,j+1) +  ubkg(i,j)   + ubkg(i-1,j) )  / 4.0d0

!WARNING - FIX THESE TERMS FOR NON-ORTHOGONAL GRID

!     &    - bit9 * (
!     &          dndxp_(i,j) * ( vbit (i,j,k) + vbit(i,j-1,k) ) * hzp_(i,j,k) * ( ubit (i,j,k) + ubit(i-1,j,k) )
!     &         + dndxp_(i,j+1) * ( vbit (i,j+1,k) + vbit(i,j,k) ) * hzp_(i,j+1,k) * ( ubit (i,j+1,k) + ubit(i-1,j+1,k) )
!     &             ) / 8.0d0
                                                                                                                                 
!     &    + bit9 * (
!     &          dmdyp_(i,j) *   ( ubit (i,j,k) + ubit(i-1,j,k) ) *       hzp_(i,j,k) * ( ubit (i,j,k) + ubit(i-1,j,k) )
!     &         + dmdyp_(i,j+1) * ( ubit (i,j+1,k) + ubit(i-1,j+1,k) ) * hzp_(i,j+1,k) * ( ubit (i,j+1,k) + ubit(i-1,j+1,k) )
!     &             ) / 8.0d0


      enddo;enddo
      enddo

           bit1x=1;bit1y=1

           Fbx(0:nx,1:ny,0:nz) =  bit1x *  Fbx(0:nx,1:ny,0:nz)  
           Fby(1:nx,0:ny,0:nz) =  bit1y *  Fby(1:nx,0:ny,0:nz) 

!          Fbz(1:nx,1:ny,nz-1) =  bit1z *  isource(1:nx,1:ny) * ( wp(1:nx,1:ny,nz-1) - wsource ) / 3600.0
!             make  bit1z = 1 in case, kill wp(1:nx,1:ny,nz:nz+1)= 0, check sign on this term
!             make wsource = 0.2

           !write(188,'(2e15.7)') Fbx(nx/2,ny/2,3), Fby(nx/2,ny/2,3)

         !write(6,*) 'gonna show Fbx', nx,ny,nz

         !if (mtime. eq. 10) then

             !call view_3d(Fbx,-1,nx+1,0,ny+1,0,nz+1)
             !call view_3d(Fby, 0,nx+1,-1,ny+1,0,nz+1)
             !stop 8888

        !endif

           fbz = 0


!           k = 0                                                                                                                         

!           fbz(nx/2-k:nx/2+k+1, ny/2-k:ny/2+k+1, 10:20) =  0* 1.0d-6  *
!     &          sin( .5d0 * 3.14159260 * time /(2.0*3600.0d0) )



       return
