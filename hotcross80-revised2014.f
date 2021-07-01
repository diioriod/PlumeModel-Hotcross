! has w sponge on sides now, also on top
 
!  HotCross Distribution - Version 8.0

! --------------------------------------------------------------------------------------------------
! -------- A 3-D CONVECTION PLUME MODEL WITH CROSS FLOW ----------------------------
! --------------------------------------------------------------------------------------------------
! ---------J.W. Lavelle,June, 2013 ----------------------------------------
! -------- NOAA/Pacific Marine Environmental Lab---------------------------------
! -------- Seattle, Washington 98115 --------------------------------------------
! -------- j.william.lavelle@lavelle-associates.org, 206-527-1062 -----------------------------------
! --------------------------------------------------------------------------------------------------
! -------- See Lavelle, J.W. (1997):  Buoyancy-driven plumes in rotating, stratified cross flows:---
! -------- Plume dependence on rotation, turbulent mixing, and cross-flow strength. ----------------
! -------- J. Geophys. Res., 102(C2), 3405-3420, for model physics. --------------------------------
! --------------------------------------------------------------------------------------------------
! -------- and 
!--------- Lavelle, J. W., Daniela Di Iorio, and Peter Rona(2013), A turbulent convection model with--
!--------- an observational context for a deep-sea hydrothermal plume in a time-variable cross-flow.--
!--------- J. Geophysical Res.- Oceans, submitted
! --------------------------------------------------------------------------------------------------\  
!compile:
!pgf95  -mp -r8 -O3  -Mextend -fast -Kieee -Minline -Mdclchk  -o hotcross  modules80.f hotcross80.f 
!valgrind --leak-check=full -v --track-origins=yes --read-var-info=yes hotcross
!source80.f smagorinsky_split_mixing.f potential_density80.f initialize_profiles80.f
!hs3crt_composite80.o mixing80.f body_force80.f sponge80.f mpdata-hotcross80.f write_snapshota.f
!view_3d.f view_2d.f  stat_2d.f stat_3d.f  -L/usr/local/lib -L/usr/bin -lnetcdf  -lcurl


!#put_block.o put_block1.o put_block2.o

!run
!hotcross > hotcross.out &

! Requires: 
!   array_size.inc
!   case.inc
!   source.inc 
!   stratification.inc 
!   site_forcing.inc
!   write.inc 

      	program hotcross

      	use velocity_module
       	use tracer_module
       	use viscosity_module
       	use constants_module
        use signals_module
        use grid_module
        use force_module
        use source_module
        use background_module
        use geometry_module
        use tracer_smol_module

        include 'array_size.inc'                                 

        real sum

        real worku(-1:nx+1,0:ny+1,0:nz+1), workv(0:nx+1,-1:ny+1,0:nz+1),  workw(0:nx+1,0:ny+1,-1:nz+1)
               
        real lat, q0_t, sumq
        integer i, j, k, l, m, numa
        integer(8) tvalues(8)
      	integer  OMP_get_thread_num, OMP_get_num_threads
        integer idum, setvbuf3f
        integer getcwd
        character *80 path_name, program_name
        character(len=8)   wdate(2),wdate_t
        character(len=10)  wtime(2),wtime_t
        character(len=5)    zone(2),zone_t

               character * 90 output_file, title, history0
               common/start_info1/output_file,title, history0
               real megawatts, vent_temperature_anomaly, vent_salinity_anomaly, vent_tracer_anomaly
               real settling_vel 
               integer min_or_max_t, min_or_max_s, min_or_max_c 

               tvalues = 0
               wtime_t = '000000.000'
               zone_t = '+0800'
               wdate_t = '21030101'

               call date_and_time(wdate_t,wtime_t,zone_t,tvalues)
                     wdate(1) = wdate_t; wtime(1)= wtime_t

        if (nx .eq. 192) call OMP_set_num_threads(8)

        call OMP_set_num_threads(8)

!$omp PARALLEL

         numa = OMP_get_num_threads()

         if ( OMP_get_thread_num() .eq. 0) then

           write(6,*) 'number of threads being used: ', numa

        endif

!$omp END PARALLEL


         CALL initialize_constants

!        call initialize_dummy_variables

                 bit_z_hp =1;  bit_z_hx =1; bit_z_hy=1; bit_z_hz=1
                 mp_= 1.0; np_ = 1.0; mp_k = 1.0;  np_k= 1.0; mpx_ = 1.0 ;  npy_ = 1.0
                 hzmn= 1.0; hzmnx=1.0d0; hzmny=1.0d0; hznx = 1.0d0;  hzmy = 1.0d0
                 hzp = 1.0d0;  hzw= 1.0d0

        call system ('head -200 makefile')

	! -------- the computational environment --------------------    

        include 'case.inc'

              allocate (p(0:nx+1,0:ny+1,0:nz+1))
              allocate (u(-1:nx+1,0:ny+1,0:nz+1), up(-1:nx+1,0:ny+1,0:nz+1), um(-1:nx+1,0:ny+1,0:nz+1))
              allocate (v(0:nx+1,-1:ny+1,0:nz+1), vp(0:nx+1,-1:ny+1,0:nz+1), vm(0:nx+1,-1:ny+1,0:nz+1))
              allocate (w(0:nx+1,0:ny+1,-1:nz+1), wm(0:nx+1,0:ny+1,-1:nz+1), wp(0:nx+1,0:ny+1,-1:nz+1))
              p=0.;w=0.;up=0.;u=0.;um=0.;vp=0.;v=0.;vm=0.;wp=0.;w=0;wm=0.

              allocate (kxt(-1:nx+1,0:ny+1,0:nz+1), kyt(0:nx+1,-1:ny+1,0:nz+1), kzt(0:nx+1,0:ny+1,-1:nz+1))
              allocate (s(0:nx+1,0:ny+1,0:nz+1), Q_s(0:nx+1,0:ny+1,0:nz+1))
              allocate (t(0:nx+1,0:ny+1,0:nz+1), Q_t(0:nx+1,0:ny+1,0:nz+1))
              allocate (c(0:nx+1,0:ny+1,0:nz+1), Q_c(0:nx+1,0:ny+1,0:nz+1))
              kxt=0.;kyt=0.;kzt=0.;s=0.;q_s=0.;t=0.;q_t=0.;c=0.;Q_c=0.

              allocate ( uc(-1:nx+1,0:ny+1,0:nz+1), vc(0:nx+1,-1:ny+1,0:nz+1), wc(0:nx+1,0:ny+1,-1:nz+1))
              allocate ( cp(0:nx+1,0:ny+1,0:nz+1) , cc(0:nx+1,0:ny+1,0:nz+1))
              allocate ( ud(-1:nx+1,0:ny+1,0:nz+1), vd(0:nx+1,-1:ny+1,0:nz+1), wd(0:nx+1,0:ny+1,-1:nz+1))
              uc = 0.;vc=0.;wc=0.;cp=0.;cc=0.;ud=0.;vd=0.;wd=0.

        xstart = - real(nx) / 3.0 * dx
        ystart = - real(ny) / 2.0 * dy

        ds = dz 

        bdbit(1:4) = 1       !bdbit= 0 if closed, 1 if open
        bdbit(3:4) = 1

        open_box = 0                                       ! default value; means box is closed

        open_box = 1
                                                   ! otherwise bottom_slip = +1 will require stilling call to EKMAN
        write(6,'(1x,a80)') title
        !call date_time (history0)

           !file_number=2000
           !inquire (unit= file_number,exist=unit_new, opened=unit_opened)
           !    if (unit_new .and. .not. unit_opened) then
           !        open ( unit = file_number , status = 'unknown')
           !        idum = setvbuf3f(file_number,2,5)
           !        if (idum .ne. 0) stop 'setvbuf messed up'
           !    else
           !        stop 'unit \'file_number\' is already open'
           !    endif

           open(unit=2000)  !status = "UNKNOWN" is default
               !idum = setvbuf3f(2000,2,10)
                idum = setvbuf3f(6,2,10)
           open( unit = 65)! , status = "UNKNOWN")
d!              idum = setvbuf3f(65,2,10)



           !file_number=6006
           !inquire (unit= file_number,exist=unit_new, opened=unit_opened)
           !    if (unit_new .and. .not. unit_opened) then
           !        open ( unit = file_number , status = 'unknown')
           !        idum = setvbuf3f(file_number,2,5)
           !        if (idum .ne. 0) stop 'setvbuf messed up'
           !    else
           !        stop 'unit \'file_number\' is already open'
           !    endif


        idum = getcwd(path_name)
        call getarg(0,program_name)

             write(6,'(a11,a90)') "path name: ", path_name
             write(6,*) "program name: ", program_name

        if (bcondx .gt. 1 .or. bcondx .lt. 0 ) stop 'bcondx must be either 1(non-clyclic) or 0(cyclic) for x-direction'
        if (bcondy .gt. 1 .or. bcondy .lt. 0 ) stop 'bcondy must be either 1(non-clyclic) or 0(cyclic) for y-direction' 
        if (bottom_slip .ne. 1 .and. bottom_slip .ne. -1) 
     &            stop 'bottom_slip must be either -1(zero vel at bottom) or 1(free slip at bottom)'
        if (side_slip .ne. 1 .and. side_slip .ne. -1) 
     &            stop 'side_slip must be either -1(no tangential vel along sides) or 1 (free slip along sides)'

        fc =   2.0 * omega * sin (pi * lat / 180.0)
        fc_y = real(bit9) * 2.0 * omega * cos (pi * lat / 180.0)                                                     

        pgrady_con = 0.0       
        pgradx_con =  0.0
          
       open( unit = 81 , status = 'old', file = 'cfl.out')
            idum = setvbuf3f(81,2,10)
       open(unit=411)!,  status = 'unknown')
d           idum= setvbuf3f(411,2,10) 
       open( unit = 401 , status = 'old', file = 'temp.diagnostics')
d           idum = setvbuf3f(401,2,10)
       open( unit = 402 , status = 'old', file = 'salt.diagnostics')
d           !idum=  setvbuf3f(402,2,10)
       open( unit = 403 , status = 'old', file = 'tracer.diagnostics')
d           !idum=  setvbuf3f(403,2,10)                                               
        CALL INITIALIZE_PROFILES !(back)
!           write(6,*) t_bkg,s_bkg;stop

        CALL DENSITY

!       rho0 = rho_rho (s_bkg(nz),t_bkg(nz),0.0d0)

        write(6,*) 'rho0 =', rho0

        min_or_max_t  = -1
        min_or_max_c  = -1
        min_or_max_s  = -1

        CALL INITIALIZE

        if ( vent_salinity_anomaly .ge. 0 ) then
              s_ref = minval ( s_bkg( 1:nz) ) * (1- 1.0d-8)    
              min_or_max_s  = -1
              write (6,*)  'min_or_max_s=',  min_or_max_s, '   indicates non-negative salinity anomaly at source'
        else
              s_ref = maxval ( s_bkg( 1:nz) ) * (1.0 + 1.0d-8)
              min_or_max_s  = +1
              write (6,*)  'min_or_max_s=',  min_or_max_s , '   indicates negative salinity anomaly at source'
        endif

        write(6,*) 't_ref,s_ref,c_ref', t_ref,s_ref,c_ref
        write(6,'(a6,d20.7,a30) ') ' q0_t = ', q0_t,'  J / s = W  total heat discharge'

        write(6,*) 'vent_temperature_anomaly: ', vent_temperature_anomaly
        write(6,*) 'vent_salinity_anomaly: ', vent_salinity_anomaly
        write(6,*) 'vent_tracer_anomaly: ', vent_tracer_anomaly
        write(6,*) 'bcondx,bcondy', bcondx,bcondy

        CALL GUTENBERG

        write(6,*) 'gutenberg'

        m = -99                              ! this is a signal to initialize the netcdf blocks
        
        time = -dt

        CALL SPONGE_MAKER

        alphae_complement = 1.0d0

        where (alphae .gt. 0.0) alphae_complement = 0.0d0         

        include 'write.inc'

        t_total=0.;s_total=0.;c_total=0.


      	! ----------------- begin time stepping -----------------------         

      	do 12 m = 0, nt

           time = time + dt
           mtime = m

           leap = 2
           delta = 1
           
           if ( m .eq. 0 ) then 
                 leap = 1     
                 delta =0
            endif

        CALL BODY_FORCE                 !UPDATE EXTERNAL FORCES

        CALL FORCES                     !UPDATE INTERNAL FORCES (cORIOLIS, NON-LINEAR TERMS, VISCOUS TERMS, DENSITY FORCE) 

        CALL PRESSURE                   !CALCULATE SEA SURFACE ELEVATION(ETA), OR PRESSURE ON RIGID LID(ETA)

        CALL MOMENTUM_U                 !CALCULATE NEW U VELOCITY CAUSED BY FORCES AND ETA

        CALL MOMENTUM_V                 !CALCULATE NEW V VELOCITY CAUSED BY FORCES AND ETA

        CALL MOMENTUM_W                 !CALCULATE NEW W VELOCITY CAUSED BY FORCES AND ETA

!       CALL JOULES_3D

        CALL VISCOSITY3D                !CALCULATE VISCOSITY COEFFICIENTS USING VELOCITY SHEARS 

        CALL SOURCE_terms               !CALCULATE NEW INSTANTANEOUS EFFLUENT FLUXES

        CALL MIXING                     !CALCULATE DIFFUSIVITY COEFFICIENTS

        ws_c = 0.0d0

             p_end = nz

        if ( toggle(6) .eq. 1 ) then    !CALCULATE POTENTIAL TEMPERATURE LOOP

             worku =  0.5 * ( up + u )        
             workv =  0.5 * ( vp + v )        
             workw =  0.5 * ( wp + w )

             tmin = minval(t_bkg); tmax = 1.0e+10        !maxval(t_bkg)

             CALL TRACER31 (t, Q_t, t_bkg, tmin, tmax, t_total, p_end,  worku, workv,  workw, dt, 401 )

             do j = 1,ny;do i=1,nx

                 t(i,j,nz+1) = (tsource + t_bkg(nz)) * isource_map(i,j) + real ( 1 - isource_map(i,j))  * t(i,j,nz+1)

             enddo;enddo
               
        endif                                             
  
        if ( toggle(7) .eq. 1 ) then     !CALCULATE SALINITY LOOP

            smin = minval(s_bkg)
            smax = maxval(s_bkg)
            smin = 29.2              !salinity at source
            smin = 0

            CALL TRACER31 (s, Q_s, s_bkg, smin, smax, s_total, p_end,  worku, workv,  workw, dt, 402 )

            do j = 1, ny;do i=1,nx

                  s(i,j,nz+1) = (ssource + s_bkg(nz)) * isource_map(i,j) + real ( 1 - isource_map(i,j)) * s(i,j,nz+1)

            enddo;enddo

        endif

        if ( bit3 .eq. 1) CALL DENSITY             !CALCULATE NEW DENSITY FIELD BASED ON UPDATED TEMP AND SALINITY 

        settling_vel = 0.0d0

        cmin = 0.0; cmax = 1.0e+10; c_bkg= 0.0

       	if (toggle(8) .eq.1 ) then                 !CALCULATE PASSIVE TRACER TRANPORT/DISPERION

            CALL TRACER31 (c, Q_c, c_bkg, cmin, cmax, c_total, p_end,  worku, workv,  workw, dt, 403 )
        endif


        CALL REFRESH             !SHIFT VALUES ONE TIME STEP INTO PAST

        if ( mod(mtime,10) .eq. 0) CALL CFL_VALUES        !TRACK THE STABILITY CRITERION

!-----------output control in write.inc

            include 'write.inc'   !IF SIGNAL TO WRITE IS YES, THEN WRITE OUTPUT NETCDF FILES
        write(2000,'(7e15.6)')  time, u(nx/2,ny/2,nz-2),v(nx/2,ny/2,nz-2),w(nx/2,ny/2,nz-2),t(nx/2,ny/2,nz-2), 
     &               s(nx/2,ny/2,nz-2), rhop(nx/2,ny/2,nz-2)

 12     continue                  !END TIME INTEGRATION LOOP
 
       m = -998                   ! this is a signal to close the netcdf writes
       include 'write.inc'

       write(6,*) 'hotcross run for', nt , 'iterations'
                write(6,*) 'started  hotcross:   ',  wdate(1),'   ', wtime(1)
       call date_and_time(wdate_t,wtime_t,zone_t,tvalues)
                write(6,*) 'finished hotcross:   ',  wdate_t,'   ', wtime_t

       end program hotcross      

! --------------------------------------------------------------------------------
! -------------------- SUBROUTINE GRACE ------------------------------------------
! --------------------------------------------------------------------------------

      subroutine cfl_values

        use grid_module
        use signals_module
        use velocity_module
        use viscosity_module
        include 'array_size.inc'
        real maxw, maxv, maxu, minu, minv, minw
        integer istop,i,j,k,l,m

!!        data istop/0/

!!        m = -998                   !this is a signal to close netcdf write blocks
!!        include 'write.inc'

!!        write(6,*) 'termination occurred on times step: ' , mtime
!!        write(6,*) 'TERMINATION ABNORMAL BUT GRACEFUL'
!!        istop = 1        
                 
!Q      ENTRY  CFL_VALUES

                minu = minval ( u (0:nx,1:ny,1:nz) ) * dt/dx
                maxu = maxval ( u (0:nx,1:ny,1:nz) ) * dt/dx
                                                                                                                          
                minv = minval ( v (1:nx,0:ny,1:nz) ) * dt/dy
                maxv = maxval ( v (1:nx,0:ny,1:nz) ) * dt/dy
                                                                                                                            
                minw = minval ( w (1:nx,1:ny,0:nz) ) * dt/dz
                maxw = maxval ( w (1:nx,1:ny,0:nz) ) * dt/dz

                write(81,'(6e15.7,i10)') maxu, minu, maxv,minv, maxw, minw, mtime
                                                                                                                 
          if (max(maxu,maxv,-minu,-minv,maxw,-minw) .gt. 0.75) then

                write(6,*) 'EMERGENCY BAILOUT AT TIME STEP = ', MTIME
                CALL put_snapshot_open (0,nx+1,0,ny+1,0,nz+1,output_file_2)
                CALL put_snapshot      (0,nx+1,0,ny+1,0,nz+1,output_file_2)
                CALL put_snapshot_close(0,nx+1,0,ny+1,0,nz+1,output_file_2)
                rewind (81)
                stop 'stop on impending cfl violation'

          endif

!                 do k = 1,nz;do j=1,ny;do i=1,nx
 
!                  if ( ss_v(i,j,k)/ ss_v(i,j,k) .ne. 1.0d0) then
                   
!                  write(6,*) 'EMERGENCY BAILOUT with NaN AT TIME STEP = ', MTIME
!!

!                  CALL put_snapshot_open (0,nx+1,0,ny+1,0,nz+1,output_file_2)
!                  CALL put_snapshot      (0,nx+1,0,ny+1,0,nz+1,output_file_2)
!                  CALL put_snapshot_close(0,nx+1,0,ny+1,0,nz+1,output_file_2)
!                  rewind (81)
!                  stop 'stop on NaN cfl violation'

!                 endif

!                 enddo;enddo;enddo
                                                  
!               if( istop .eq. 1) stop

!Q         return

       	end subroutine cfl_values

! ------------------------------------------------------------------------------------------------
! ----------------- SUBROUTINE INITIALIZE --------------------------------------------------------
! ------------------------------------------------------------------------------------------------

       subroutine initialize

       use constants_module
       use grid_module
       use velocity_module
       use tracer_module
       use viscosity_module
       use signals_module
       include 'array_size.inc'

       real c_zero,c_top, c_fluxh, q_vertical(0:nz+1)
       real t_bkg_end(0:nz+1), s_bkg_end(0:nz+1), c_bkg_end(0:nz+1)
       real harvest(-1:nx+1,-1:ny+2,-1:nz+2 )
       integer i,j,k,m       

       ! INITIALIZE_VELOCITY
 
       random_amp1=0     
       call random_number (harvest)
       u  =  random_amp1 * ( 2. * ( harvest (-1:nx+1,0:ny+1,0:nz+1) -.5 ) )
       call random_number (harvest)
       up = random_amp1 * ( 2. * ( harvest (-1:nx+1,0:ny+1,0:nz+1) -.5 ) )
       call random_number (harvest)
       um = random_amp1 * ( 2. * ( harvest (-1:nx+1,0:ny+1,0:nz+1) -.5 ) )
       call random_number (harvest)
       v  = random_amp1 * ( 2. * ( harvest (0:nx+1,-1:ny+1,0:nz+1) -.5 ) )
       call random_number (harvest)
       vp = random_amp1 * ( 2. * ( harvest (0:nx+1,-1:ny+1,0:nz+1) -.5 ) )
       call random_number (harvest)
       vm = random_amp1 * ( 2. * ( harvest (0:nx+1,-1:ny+1,0:nz+1) -.5 ) )
       call random_number (harvest)
       w  = random_amp1 * ( 2. * ( harvest (0:nx+1,0:ny+1,-1:nz+1) -.5 ) )
       call random_number (harvest)
       wp = random_amp1 * ( 2. * ( harvest (0:nx+1,0:ny+1,-1:nz+1) -.5 ) )
       call random_number (harvest)
       wm = random_amp1 * ( 2. * ( harvest (0:nx+1,0:ny+1,-1:nz+1) -.5 ) )
       call random_number (harvest)
       p  = random_amp1 * ( 2. * ( harvest (0:nx+1,0:ny+1,0:nz+1) -.5 ) )

!      call random_number (harvest)
!      ut = random_amp1 * ( 2 * ( harvest (-1:nx+1,0:ny+1,0:nz+1) -.5 ) )
!      call random_number (harvest)
!      vt = random_amp1 * ( 2 * ( harvest (0:nx+1,-1:ny+1,0:nz+1) -.5 ) )
!      call random_number (harvest)
!      wt = random_amp1 * ( 2 * ( harvest (0:nx+1,0:ny+1,-1:nz+1) -.5 ) )

! ----------  VELOCITY ----------------------------------------------------------------

        do 700 k = 0 , nz +1                     !  "boundary layer" made with z-dependent viscosity

             pgradx(k) = bit8 * pgradx_con
             pgrady(k) = bit8 * pgrady_con

  700   continue

        !CALL EKMAN                               ! RETURNS U_BKG AND V_BKG

        if (bit7 .eq. 0) then

            do k = 0, nz + 1  		         ! added to keep non-rotating case symmetric 6/7/95
                v_bkg(k) = 0.0
                u_bkg(k) = 0.0
            enddo

            write(6,*) "WARNING - since rotation is zero, pressure gradient profile has been tailored to
     &      provide a u_bkg profile identical to that with rotation at same latitude and a v_bkg of zero"

            do k = 1, nz 
                pgrady(k) = 0.0
                pgradx(k) = bit8 * (  azzmin(k) * (u_bkg(k+1) - u_bkg(k) ) -  azzmin(k-1) * (u_bkg(k) - u_bkg(k-1) ) ) / dz / dz
            enddo
                                    
        endif

!       call random_number (harvest)
!       wt = random_amp1 * ( 2 * ( harvest (0:nx+1,0:ny+1,-1:nz+1) -.5 ) )

        do  k = 0, nz+1
        do  j = 0, ny+1
        do  i = -1, nx +1

                um(i,j,k) = u_bkg(k) + random_amp1 * ( 2 * ( harvest (i,j,k) -.5 ) )

        enddo;enddo;enddo

        call random_number (harvest)

        do  k = 0, nz+1
        do  j = 0, ny+1
        do  i = -1, nx +1

                u(i,j,k) = u_bkg(k) + random_amp1 * ( 2 * ( harvest (i,j,k) -.5 ) )

        enddo;enddo;enddo

        call random_number (harvest)
        do k = 0, nz+1
        do j = -1, ny+1
        do i = 0, nx +1

                vm(i,j,k) = v_bkg(k) + random_amp1 * ( 2 * ( harvest (i,j,k) -.5 ) )

        enddo;enddo;enddo

        call random_number (harvest)
        do k = 0, nz+1
        do j = -1, ny+1
        do i = 0, nx +1

                v(i,j,k) = v_bkg(k) + random_amp1 * ( 2 * ( harvest (i,j,k) -.5 ) )

        enddo;enddo;enddo

        call fill_value(u,-1,nx+1,0,ny+1,0,nz+1)
        call fill_value(v,0,nx+1,-1,ny+1,0,nz+1)
        call fill_value(w,0,nx+1,0,ny+1,-1,nz+1)
        call fill_value(up,-1,nx+1,0,ny+1,0,nz+1)
        call fill_value(vp,0,nx+1,-1,ny+1,0,nz+1)
        call fill_value(wp,0,nx+1,0,ny+1,-1,nz+1)
        call fill_value(p,0,nx+1,0,ny+1,0,nz+1)

         t_ref = minval(t_bkg(1:nz)) * (1.0 - 1.0d-8)             !corrector now cannot make lower temp
         s_ref = minval(s_bkg(1:nz)) * (1.0 - 1.0d-8)  
         c_ref = minval(c_bkg(1:nz)) * (1.0 - 1.0d-8)   


	! INITIALIZE  TRACERS

         if( toggle(9)  .eq. 1) then
 
		allocate ( cc1 (0:nx+1, 0: ny+1, 0:nz+1) )
		allocate (q_cc1(0:nx+1, 0: ny+1, 0:nz+1) )
                q_cc1 = 0.0;cc1 =0.0
                dqcc1dt = 0*2.9 e-3 * dQTdt * rho0 * cpp
        	cc1_lambda = 0.0
        	cc1_ws =  2.0d-7          !1.0d-5/5.0       !( 1.0d-5) / 10
        	cc1_bkg = 1.0d-30
                c_fluxh = 3.40d-5
                c_zero = 45.0*1000
                c_top =  33.0*1000      !now units of nmol/m3

                q_vertical = 0
                cc1_bkg = cc1_bkg 

                do 10 k = 0, nz+1
                do 10 j = 0, ny+1
                do 10 i = 0, nx+1 
                    cc1(i,j,k) = cc1_bkg(k)
  10            continue

		cc1_ref = minval( cc1_bkg(1:nz) )

        endif

        if( toggle(10)  .eq. 1) then

                allocate ( cc2 (0:nx+1, 0: ny+1, 0:nz+1) )
                allocate (q_cc2(0:nx+1, 0: ny+1, 0:nz+1) )

                q_cc2 = 0.0;cc2=0
                dqcc2dt = 1.0e+6
                dqcc2dt = 0*2.9 e-3 * dQTdt * rho0 * cpp      !dqTdt is already dividied by ncells
                                                            ! dqTdt had previously been dividded by  cpp*rho 
							    ! 2.9 e-3 is Del_aluminum of (4.5e+6)/cpp/rho/350 oC
                cc2_lambda = 0.0
                cc2_ws =  1.0d-5     ! 1.0e-4
                cc2_bkg = 1.0e-30
                c_top =0
                c_zero= 2*16.0*1000          !now units of nmol/m3

!                call profile_maker ( cc2_bkg, nz, c_zero, c_top, q_vertical, 2, dx, dy, dz, cc2_ws, azzmin)
                do 20 k =0, nz+1
                do 20 j= 0, ny+1
                do 20 i = 0, nx+1
   20                cc2(i,j,k) = cc2_bkg(k)
                cc2_ref = minval( cc2_bkg(1:nz) )

        endif

        if( toggle(11)  .eq. 1) then

                allocate ( cc3 (0:nx+1, 0: ny+1, 0:nz+1) )
                allocate (q_cc3(0:nx+1, 0: ny+1, 0:nz+1) )

                q_cc3 = 0.0;cc3= 0.0
                dqcc3dt = 0*1.0e+6
                cc3_lambda = 0.0
                cc3_ws = 0.0
                cc3_bkg = 1.0e-30
                do 30 k = 0, nz+1
                do 30 j = 0, ny+1
                do 30 i = 0, nx+1
                    cc3(i,j,k) = cc3_bkg(k)
                    if ( (zstart + k *dz ) .gt. (zend - 30))  cc3(i,j,k) = 10*1000
                    if ( (zstart + k *dz ) .gt. (zend - 30))  cc3_bkg(k) = 10*1000
   30           continue
                cc3_ref = minval( cc3_bkg(1:nz) )

        endif

        if( toggle(12)  .eq. 1) then

                allocate ( cc4 (0:nx+1, 0: ny+1, 0:nz+1) )
                allocate (q_cc4(0:nx+1, 0: ny+1, 0:nz+1) )

                q_cc4 = 0.0;cc4=0.0
                dqcc4dt = 1.0e+6
                cc4_lambda = 0.0
                cc4_ws = 0.0
                cc4_bkg = 1.0e-30
                do  k = 0, nz+1
                do  j = 0, ny+1
                do  i = 0, nx+1
                   cc4(i,j,k) = cc4_bkg(k)
                enddo;enddo;enddo
                cc4_ref = minval( cc4_bkg(1:nz) )

        endif

       return
       end
!
! **********************************************************************************************************
! ----------------------------------------------------------------------------------------------
! ----------------- SUBROUTINE MOMENTUM_U  ------------------------------------------------------
! -----------------------------------------------------------------------------------------------

       subroutine momentum_u                       

       use force_module, only: fx, fbx
       use constants_module
       use grid_module
       use velocity_module
       use viscosity_module
       use signals_module
       include 'array_size.inc'
       integer i,j,k,m

             up = missing_value

!$OMP PARALLEL DO

       do  k = 1, nz
       do  j = 1, ny
       do  i = 0, nx

             up(i,j,k) = (1 - delta) * u(i,j,k)  + delta * um(i,j,k) 
     &                  + dt * leap * ( -  (p(i+1,j,k) - p(i,j,k) ) / dx  - pgradx(k) + fx(i,j,k) )

      enddo;enddo;enddo

!$OMP END PARALLEL DO


! --- boundary conditions --------------                                    

      do  j = 1, ny
      do  i = 0, nx

          up(i,j,0)    = up(i,j,1)                  !BOX top
          up(i,j,nz+1) = bottom_slip * up(i,j,nz)                      

      enddo;enddo

      do k = 1, nz
      do i = 0, nx                                            !box front and back

           up(i, 0, k)  = side_slip * up(i,1,k)  * bcondy + ( 1 - bcondy ) * up (i,ny,k)
           up(i,ny+1,k) = side_slip * up(i,ny,k) * bcondy + ( 1 - bcondy ) * up (i,1,k)
!NEW
            up(i, 0, k)  =  up(i, 1, k)
            up(i,ny+1,k) =  up(i,ny,k)

        enddo; enddo

        do k = 1, nz                   !box ends
        do j = 1, ny

          up(-1,j,k) = up(0,j,k)
          up(nx+1,j,k) = up(nx,j,k)   ! zero derivatives

        enddo;enddo

        return                                                           
        end

! ------------------------------------------------------------------------------------------------
! ------------------ SUBROUTINE MOMENTUM_V -------------------------------------------------------
! ------------------------------------------------------------------------------------------------

       subroutine momentum_v                       

       use force_module, only: fy, fby
       use constants_module
       use grid_module
       use velocity_module
       use viscosity_module
       use signals_module
       include 'array_size.inc'
       integer i,j,k,m


           vp = missing_value

!$OMP PARALLEL DO

      do  k = 1, nz
      do  j = 0, ny          !1,ny-1
      do  i = 1, nx

             vp(i,j,k) = (1 - delta) * v(i,j,k) + delta * vm(i,j,k) 
     &                   + dt * leap * ( -  (p(i,j+1,k) - p(i,j,k))/dy  -  pgrady(k) + fy(i,j,k) )

      enddo;enddo;enddo

!$OMP END PARALLEL DO

! --- boundary conditions ---                                    

      do  j = 0, ny
      do  i = 1, nx

          vp(i,j,0) =  vp(i,j,1)                  !BOX top
          vp(i,j,nz+1) =  bottom_slip * vp(i,j,nz)              

      enddo;enddo
                                                         
        do k = 1, nz                                          !BOX x-ends
        do j = 0, ny
!NEW
          vp(0,j,k)    =   vp(1,j,k)                                !box ends
          vp(nx+1,j,k) =   vp(nx,j,k)
!NEW-NEW

        enddo; enddo

        do k = 1, nz
        do i = 1, nx                               !box y-ends 
!NEW
          vp(i,ny+1,k) =  vp(i,ny,k)    
          vp(i,-1,k)   =  vp(i,0,k)     !makes vp = vbkg(i,j) at boundary 

        enddo; enddo

      return                                                           
      end

! --------------------------------------------------------------------------
! ---------SUBROUTINE MOMENTUM_w -------------------------------------------
! ---------------------------------------------------------------------------

       subroutine momentum_w                       

       use force_module, only : fz, fbz
       use constants_module
       use grid_module
       use velocity_module
       use viscosity_module
       use signals_module
       use tracer_module, only: Q_t, isource_map
       include 'array_size.inc'
       integer i,j,k,m           
       real sumup, sumdown, wup, wramp_hrs

!$OMP PARALLEL DO
   
       do  k = 1, nz-1
       do  j = 1, ny
       do  i = 1, nx

             wp(i,j,k) = (1 - delta) * w(i,j,k) + delta *wm(i,j,k) 
     &                              + leap * dt * ( -  (p(i,j,k+1) - p(i,j,k)) / dz + fz(i,j,k) )

       enddo;enddo;enddo

!$OMP END PARALLEL DO

           sumup = 0.0; sumdown = 0.0

      do  j = 1, ny
      do  i = 1, nx

          if ( alphae(i,j) .eq. 0) then 
             if (isource_map(i,j) .eq. 1) then
                  sumup = sumup + 1.0
             else 
                  sumdown = sumdown + 1.0
             endif
          endif 
      enddo;enddo
 
         wramp_hrs = 6.0d0

         wup = wsource *  min (1.0, (mtime  * dt) / (wramp_hrs * 3600.0)  )

      do  j = 1, ny
      do  i = 1, nx                  !box top

         wp(i,j,-1)   = 0.0    ! - wp(i,j,1)
         wp(i,j,nz+1) = 0.0    ! - wp(i,j,nz-1)
         wp(i,j,0)    = 0.0
         wp(i,j,nz)   = 0.0      

         if( alphae(i,j) .eq. 0) then
             if (isource_map(i,j) .eq. 1) then
                  wp(i,j,nz)   = -wup
                  wp(i,j,nz+1) = -wup
             else
                  wp(i,j,nz)   = wup * sumup / sumdown
                  wp(i,j,nz+1) = wup * sumup / sumdown
            endif
         endif

      enddo;enddo


        
        do k = 0, nz
        do i = 1, nx                                            !box y-ends

           wp(i, 0, k)  = side_slip * wp(i,1,k)  * bcondy + ( 1 - bcondy ) * wp(i,ny,k)
           wp(i,ny+1,k) = side_slip * wp(i,ny,k) * bcondy + ( 1 - bcondy ) * wp(i,1,k)

!NEW

           wp(i, 0, k)  =  0* wp(i, 1, k)
           wp(i,ny+1,k) =  0* wp(i,ny,k)   


         enddo; enddo

        do k = 1, nz                                          !BOX x-ends
        do j = 0, ny

          wp(0,j,k) = 0                                 !wp(1,j,k)                             !  0.0
          wp(nx+1,j,k) = 2.0 * wp(nx,j,k) -  wp(nx-1,j,k)
!NEW
          wp(0,j,k)    = 0 * wp(1,j,k) 
          wp(nx+1,j,k) = 0 * wp(nx,j,k) 

        enddo; enddo

      return                                                           
      end
! ________________________________________________________________________________________________

      subroutine forces

       use force_module
       use constants_module
       use grid_module
       use velocity_module
       use viscosity_module   !deprecated, ss => ai
       use signals_module
       use tracer_module, only: rhop

       include 'array_size.inc'
       integer i,j,k
       real coe_2zz, coe_5xx, coe_5yy, coe_4x, coe_4y, coe_4z, coe_7, coe_9, coe_10 
!       real s11,s12,s13,s21,s22,s23,s31,s32,s33
       real avxp,avxm,avyp,avym,ahp,ahm
!       real s13m,s13p,s23p,s23m

       real dissipation_x, dissipation_y, dissipation_z

       coe_2zz = bit2 / dz / dz
       coe_5xx = bit5 / dx / dx      
       coe_5yy = bit5 / dy / dy

       coe_4x = bit4 / 4.0 / dx
       coe_4y = bit4 / 4.0 / dy
       coe_4z = bit4 / 4.0 / dz

       coe_7 =  bit7 * fc / 4.0d0
       coe_9 =  bit9 * fc_y / 4.0d0
       coe_10 = bit10 


!$OMP PARALLEL WORKSHARE

          axu(0:nx+1,0:ny+1,0:nz+1) =  ss(0:nx+1,0:ny+1,0:nz+1) + axmin
          ayu(0:nx,0:ny,0:nz+1)     = .25d0 * ( ss(0:nx,0:ny,0:nz+1) + ss(1:nx+1,0:ny,0:nz+1) + ss(0:nx,1:ny+1,0:nz+1) 
     &                                     + ss(1:nx+1,1:ny+1,0:nz+1) ) + aymin 
          axv(0:nx,0:ny,0:nz+1)     = .25d0 * ( ss(0:nx,0:ny,0:nz+1) + ss(1:nx+1,0:ny,0:nz+1) + ss(0:nx,1:ny+1,0:nz+1)  
     &                                     + ss(1:nx+1,1:ny+1,0:nz+1) ) + axmin
          ayv(0:nx+1,0:ny+1,0:nz+1) =  ss(0:nx+1,0:ny+1,0:nz+1) + aymin
          axw(0:nx,0:ny+1,0:nz)     = .25d0 * ( ss(0:nx,0:ny+1,0:nz) + ss(1:nx+1,0:ny+1,0:nz) + ss(0:nx,0:ny+1,1:nz+1)  
     &                                     + ss(1:nx+1,0:ny+1,1:nz+1) ) + axmin  
          ayw(0:nx+1,0:ny,0:nz)     = .25d0 * ( ss(0:nx+1,0:ny,0:nz) + ss(0:nx+1,1:ny+1,0:nz) + ss(0:nx+1,0:ny,1:nz+1)  
     &                                     + ss(0:nx+1,1:ny+1,1:nz+1) ) +aymin

!$OMP END PARALLEL WORKSHARE

          do k = 0,nz

              azu(1:nx-1,1:ny,k) = .25 * ( ss_v(1:nx-1,1:ny,k) + ss_v(2:nx,1:ny,k) + ss_v(1:nx-1,1:ny,k+1) 
     &              + ss_v(2:nx,1:ny, k+1) ) + azzmin(k)
              azv(1:nx,1:ny-1,k) = .25 * ( ss_v(1:nx,1:ny-1,k) + ss_v(1:nx,2:ny,k) + ss_v(1:nx,1:ny-1,k+1) 
     &              + ss_v(1:nx,2:ny, k+1) ) + azzmin(k)

          enddo
 
!          do k = 1,nz+1
!2014 REVISION
!               azw(0:nx+1,0:ny+1,k) =  ss_v(0:nx+1,0:ny+1,k) + azzmin(k)
!                azw(0:nx+1,0:ny+1,k) =  ss(0:nx+1,0:ny+1,k) + axmin           
!2014 REVISION
!          enddo

          do i = 0, 7
              azu(i,1:ny,1:nz+1) = azu(i,1:ny,1:nz+1) + real (7 - i ) *.01 / 7.0    !intensifies vertical mixing at entrance
              azv(i,0:ny,1:nz+1) = azv(i,0:ny,1:nz+1) + real (7 - i ) *.01 / 7.0
          enddo

          fx = missing_value
          fy = missing_value

!$OMP PARALLEL DO

       do  k = 1, nz            
       do  j = 1, ny
       do  i = 0, nx

         !dissipation_x = d/dx[2 AH du/dx] + d/dy[ AH (du/dy+dv/dx)] + d/dz[ Av (dw/dx + du/dz)]

          dissipation_x =      

     &     2.0 *  bit5 * (axu(i+1,j,k) * ( um(i+1,j,k) - um(i,j,k) ) - axu(i,j,k)   * ( um(i,j,k) - um(i-1,j,k) ) ) /dx/dx
     &          + bit5 * (ayu(i,j,k) * ( um(i,j+1,k) - um(i,j,k) ) - ayu(i,j-1,k) * ( um(i,j,k) - um(i,j-1,k) ) )  /dy/dy
     &          + bit5 * (ayu(i,j,k) * ( vm(i+1,j,k) - vm(i,j,k) ) - ayu(i,j-1,k) * ( vm(i+1,j-1,k) - vm(i,j-1,k) ))/dx/dy
     &          + bit2 * (azu(i,j,k) * ( um(i,j,k+1) - um(i,j,k) ) - azu(i,j,k-1) * ( um(i,j,k) - um(i,j,k-1) )   ) /dz/dz
!     &          + bit2 * (azu(i,j,k-1) * ( wm(i+1,j,k-1) - wm(i,j,k-1) ) - azu(i,j,k) * ( wm(i+1,j,k) - wm(i,j,k) ) ) /dx/dz
     &          + bit2 * (azu(i,j,k) * ( wm(i+1,j,k) - wm(i,j,k) ) - azu(i,j,k-1) * ( wm(i+1,j,k-1) - wm(i,j,k-1) ) ) /dx/dz

          fx(i,j,k) =  dissipation_x 

!2014 rev -dissipation code below replaced

!    &                  coe_2zz * ( azu(i,j,k) * (um(i,j,k+1)-um(i,j,k))  - azu(i,j,k-1) * (um(i,j,k)-um(i,j,k-1)) )

!2014 rev

     &                 + bit3 * 0.0d0

     &                - 0.5 *alphau(i,j) * (um(i,j,k)- ubkg_m(i,j))

     &                 - coe_4x * ( (u(i+1,j,k)+u(i,j,k))**2 - (u(i,j,k)+u(i-1,j,k))**2 ) 

     &                 - coe_4z * ( (w(i,j,k)+w(i+1,j,k)) * (u(i,j,k+1)+u(i,j,k)) 
     &                            - (w(i,j,k-1)+w(i+1,j,k-1)) * (u(i,j,k)+u(i,j,k-1)) )  

     &                 - coe_4y * ( (v(i,j,k)+v(i+1,j,k)) * (u(i,j+1,k)+u(i,j,k))
     &                            - (v(i,j-1,k)+v(i+1,j-1,k)) * (u(i,j,k)+u(i,j-1,k)) ) 

!2014 rev -dissipation code below replaced

!     &                 + coe_5xx * (axu(i+1,j,k) * ( um(i+1,j,k) - um(i,j,k)   ) - axu(i,j,k) * ( um(i,j,k) - um(i-1,j,k) ) )
!     &                 + coe_5yy * (  ayu(i,j,k) * ( um(i,j+1,k) - um(i,j,k) ) - ayu(i,j-1,k) * ( um(i,j,k) - um(i,j-1,k) ) )

!2014 rev
             
     &                 + coe_7 * ( v(i+1,j,k) + v(i,j,k) + v(i+1,j-1,k) + v(i,j-1,k) ) 
     &                 - coe_9 * ( w(i+1,j,k) + w(i+1,j,k-1) + w(i,j,k) + w(i,j,k-1) )     

     &                 + coe_10 *  fbx(i,j,k) 

       enddo;enddo;enddo   

!$OMP END PARALLEL DO

!$OMP PARALLEL DO


                !call view_3d(fx,-1,nx+1,0,ny+1,0,nz+1)

       do 20 k = 1, nz            
       do 20 j = 0, ny      !1,ny-1
       do 20 i = 1, nx

        !dissipation_y = d/dx[ AH (du/dy+dv/dx)] + d/dy[2 AH dv/dy] + d/dz[ Av (dw/dy + dv/dz)]
      
         dissipation_y =

     &    + bit5 * (  axv(i,j,k) * ( vm(i+1,j,k) - vm(i,j,k) ) - axv(i-1,j,k) * ( vm(i,j,k) - vm(i-1,j,k) ) )  / dx/dx
     &    + bit5 * (  axv(i,j,k) * ( um(i,j+1,k) - um(i,j,k) ) - axv(i-1,j,k) * ( um(i-1,j+1,k) - um(i-1,j,k) ))/dx/dy
     &    + 2.0 * bit5 * (ayv(i,j+1,k) * ( vm(i,j+1,k) - vm(i,j,k) ) - ayv(i,j,k) * ( vm(i,j,k) - vm(i,j-1,k) ) )  / dy/dy
     &    + bit2 * (  azv(i,j,k) * ( vm(i,j,k+1 ) - vm(i,j,k) ) - azv(i,j,k-1) * ( vm(i,j,k) - vm(i,j,k-1) )    ) /dz/dz
!     &    + bit2 * (  azv(i,j,k-1) * ( wm(i,j+1,k-1) - wm(i,j,k-1) ) - azv(i,j,k) * ( wm(i,j+1,k) - wm(i,j,k) )) /dy/dz
     &    + bit2 * (  azv(i,j,k) * ( wm(i,j+1,k) - wm(i,j,k) ) - azv(i,j,k-1) * ( wm(i,j+1,k-1) - wm(i,j,k-1) )) /dy/dz

         fy(i,j,k) = dissipation_y

!2014 rev -- disspation code below replaced
!    &                 coe_2zz * (azv(i,j,k) * (vm(i,j,k+1)-vm(i,j,k)) - azv(i,j,k-1) * (vm(i,j,k)-vm(i,j,k-1)) )
!2014 rev

     &                 + bit3 * 0.0d0

     &                 - 0.5*alphav(i,j) * (vm(i,j,k)-vbkg_m(i,j))

     &                 - coe_4y * ( (v(i,j+1,k)+v(i,j,k)) **2 - (v(i,j,k)+v(i,j-1,k)) **2 ) 

     &                 - coe_4z * ( (v(i,j,k+1)+v(i,j,k)) * (w(i,j+1,k)+w(i,j,k)) -
     &                              (v(i,j,k)+v(i,j,k-1)) * (w(i,j+1,k-1)+w(i,j,k-1)) ) 

     &                 - coe_4x * ( (v(i+1,j,k)+v(i,j,k)) * (u(i,j,k)+u(i,j+1,k)) -
     &                              (v(i,j,k)+v(i-1,j,k)) * (u(i-1,j,k)+u(i-1,j+1,k)) )

!2014 rev -dissipation code below replaced

!     &                 + coe_5xx * ( axv(i,j,k) * ( vm(i+1,j,k) - vm(i,j,k) ) -  axv(i-1,j,k) * ( vm(i,j,k) - vm(i-1,j,k) ) )
!     &                 + coe_5yy * ( ayv(i,j+1,k) * ( vm(i,j+1,k) - vm(i,j,k) ) -  ayv(i,j,k) * ( vm(i,j,k) - vm(i,j-1,k) ) )

!2014 rev


     &                 - coe_7 * ( u(i,j,k)  + u(i-1,j,k) + u(i,j+1,k) + u(i-1,j+1,k) )


     &                 + coe_10 *  fby(i,j,k)

  20   continue   

!$OMP END PARALLEL DO

!$OMP PARALLEL DO
       do 30  k = 0, nz       !revision from    0, nz   
       do 30  j = 1, ny
       do 30  i = 1, nx


         ahp =  ss(i,j,k+1) + axmin 
         ahm =  ss(i,j,k) + axmin

         avxp =  0.25* (ss_v(i,j,k+1) + ss_v(i,j,k) +ss_v(i+1,j,k+1) + ss_v(i+1,j,k)) + azzmin(k)
         avxm =  0.25* (ss_v(i,j,k+1) + ss_v(i,j,k) +ss_v(i-1,j,k+1) + ss_v(i-1,j,k)) + azzmin(k)
         avyp =  0.25* (ss_v(i,j,k+1) + ss_v(i,j,k) +ss_v(i,j+1,k+1) + ss_v(i,j+1,k)) + azzmin(k)
         avym =  0.25* (ss_v(i,j,k+1) + ss_v(i,j,k) +ss_v(i,j-1,k+1) + ss_v(i,j-1,k)) + azzmin(k)

         !s13p =  (wm(i+1,j,k) - wm(i,j,k)) /dx + (um(i,j,k+1)-um(i,j,k))/dz 
         !s13m =  (wm(i,j,k) - wm(i-1,j,k)) /dx + (um(i-1,j,k+1)-um(i-1,j,k))/dz
         !s23p =  (wm(i,j+1,k) - wm(i,j,k)) /dy + (vm(i,j,k+1)-vm(i,j,k))/dz
         !s23m =  (wm(i,j,k) - wm(i,j-1,k)) /dy + (vm(i,j-1,k+1)-vm(i,j-1,k))/dz

        !dissipation_z =  d/dx[ Av (dw/dx + du/dz)] + d/dy [ Av (dw/dy + dv/dz)} + d/dz [2 AH dw/dz} 

         dissipation_z =  

     &       2.0 * bit5 * ( ahp * (wm(i,j,k+1) - wm(i,j,k) ) - ahm * (wm(i,j,k)-wm(i,j,k-1)) )/dz/dz
     &      + bit2 *  ( avxp * (wm(i+1,j,k)- wm(i,j,k) ) - avxm * (wm(i,j,k)- wm(i-1,j,k) ) ) / dx /dx
     &      + bit2 *  ( avxp * (um(i,j,k)- um(i,j,k+1) ) - avxm * (um(i-1,j,k)- um(i-1,j,k+1) ) ) / dz /dx
     &      + bit2 *  ( avyp * (wm(i,j+1,k)- wm(i,j,k) ) - avym * (wm(i,j,k)- wm(i,j-1,k) ) ) / dy /dy
     &      + bit2 *  ( avyp * (vm(i,j,k)- vm(i,j,k+1) ) - avym * (vm(i,j-1,k)- vm(i,j-1,k+1) ) ) / dz /dy

          fz(i,j,k) =  dissipation_z


!2014 rev -substituted
!   &                coe_2zz * ( azw(i,j,k+1) * (wm(i,j,k+1) - wm(i,j,k) ) - azw(i,j,k) * (wm(i,j,k)-wm(i,j,k-1)) )

!2014 rev

     &                 + bit3 * g * ( rhop(i,j,k) + rhop(i,j,k+1) ) / rho0 / 2.0

     &                 - 0.0 * alphaw(k)  *  wm(i,j,k)    !- wbkg_m(i,j))

     &                 - 0.5 * alphaee(i,j,k) * wm(i,j,k)

     &                 - coe_4z * ( (w(i,j,k+1)+w(i,j,k))**2 - (w(i,j,k)+w(i,j,k-1))**2 ) 

     &                 - coe_4x * ( (w(i+1,j,k)+w(i,j,k)) * (u(i,j,k)+u(i,j,k+1)) -
     &                            (w(i,j,k) + w(i-1,j,k)) * (u(i-1,j,k)+u(i-1,j,k+1))  ) 

     &                 - coe_4y * ( (w(i,j+1,k)+w(i,j,k)) * (v(i,j,k)+v(i,j,k+1)) -
     &                            (w(i,j,k)+w(i,j-1,k)) * (v(i,j-1,k)+v(i,j-1,k+1))  ) 

!2014 rev-substituted
!     &                 + coe_5xx * (axw(i,j,k) * ( wm(i+1,j,k) - wm(i,j,k) ) - axw(i-1,j,k) * ( wm(i,j,k) - wm(i-1,j,k) ) )
!     &                 + coe_5yy * (ayw(i,j,k) * ( wm(i,j+1,k) - wm(i,j,k) ) - ayw(i,j-1,k) * ( wm(i,j,k) - wm(i,j-1,k) ) )


!2014 rev
     



     &                 + coe_9 * ( u(i,j,k+1) + u(i-1,j,k+1) +  u(i,j,k) + u(i-1,j,k) ) 
     &                 + coe_10 *  fbz(i,j,k) 


  30   continue

!$OMP END PARALLEL DO


      end subroutine forces


! -------------------------------------------------------------------------------------------------
! ------------------ SUBROUTINE PRESSURE ----------------------------------------------------------
! -------------------------------------------------------------------------------------------------

       subroutine pressure                   

       use force_module
       use constants_module
       use grid_module
       use velocity_module
       use viscosity_module
       use signals_module
       use tracer_module

       include 'array_size.inc'
       integer i,j,k

       real bdxs(1:ny,1:nz),bdxf(1:ny,1:nz),bdys(1:nx,1:nz),bdyf(1:nx,1:nz),
     &     bdzs(1:nx,1:ny),bdzf(1:nx,1:ny), r(1:nx,1:ny,1:nz), qcheck(1:nx,1:ny,1:nz)
       real, save:: work(2*nx*ny*nz+6*nx+5*(ny+nz)+46)          !work at least 2*l*m*n +6*l+5*(m+n) +46
       real, save:: pterb,xl,yl,zl,xs,ys,zs,lambda,sumr, ave
       integer, save:: ierror

!       real dil(0:nx+1,0:ny+1,0:nz+1)

!        div_old = 0.0d0      
!        sumd = 0.0d0                                    

!        do 10 k = 1, nz
!        do 10 j = 1, ny
!        do 10 i = 1, nx

!           dil(i,j,k) = (um(i,j,k)-um(i-1,j,k)) / dx  + (wm(i,j,k) -wm(i,j,k-1)) / dz
!     1          +(vm(i,j,k) - vm(i,j-1,k)) / dy

!            sumd = sumd + dil(i,j,k) * dx * dy * dz

!   10  continue

!            div_p = maxval ( dil(1:nx,1:ny,1:nz) )
!            div_m = minval ( dil(1:nx,1:ny,1:nz) )

!               global_div = sumd                
!               write(400,4001) mtime,div_p,div_m,global_div
!4001           format(1x,i10,3e15.7)


! ------------ Establish the boundary conditions
          
         do 20 k = 1,nz
         do 20 j = 1,ny

            i = 0
!           bdxs(j,k) =  - dudt(0,j,k) - pgradx(k) + fx(i,j,k)
            bdxs(j,k) =  - ( u(i,j,k)  - um(i,j,k) ) / dt   - pgradx(k) + fx(i,j,k) 
            bdxs(j,k) =  - ( ubkg(i,j)  - ubkg_m(i,j) ) / dt   - pgradx(k) + fx(i,j,k)


            i = nx
!           bdxf(j,k) = - dudt(nx,j,k)  - pgradx(k) + fx(i,j,k)
            bdxf(j,k) =  - ( u(i,j,k) - um(i,j,k) ) / dt   -  pgradx(k) + fx(i,j,k)
            bdxf(j,k) =  - ( ubkg(i,j) - ubkg_m(i,j) ) / dt   -  pgradx(k) + fx(i,j,k)

   20    continue

         do 30 k = 1, nz
         do 30 i = 1, nx
          
            j = 0
!           bdys(i,k) =  - dvdt(i,0,k) -  pgrady(k) + fy(i,j,k)
            bdys(i,k) =  - ( v(i,j,k) - vm(i,j,k) ) / dt   -  pgrady(k) + fy(i,j,k)
            bdys(i,k) =  - ( vbkg(i,j) - vbkg_m(i,j) ) / dt   -  pgrady(k) + fy(i,j,k)
            j = ny
!           bdyf(i,k) = - dvdt(i,ny,k) -  pgrady(k) + fy(i,j,k)
            bdyf(i,k) = - ( v(i,j,k) - vm(i,j,k) ) / dt   - pgrady(k) + fy(i,j,k)
            bdyf(i,k) = - ( vbkg(i,j) - vbkg_m(i,j) ) / dt   - pgrady(k) + fy(i,j,k)

   30   continue


        do 35 j = 1, ny
        do 35 i = 1, nx
 
           k = 0
!          bdzs(i,j) = - dwdt(i,j,0) + fz(i,j,k)
           bdzs(i,j) = - ( w(i,j,k) - wm(i,j,k) ) / dt   + fz(i,j,k)

           k = nz
!          bdzf(i,j) = - dwdt(i,j,nz) + fz(i,j,k)
           bdzf(i,j) = - ( w(i,j,k) - wm(i,j,k) ) / dt   + fz(i,j,k)

   35   continue 

            !bdxs = 0
            !bdxf = 0
            !bdys = 0
            !bdyf = 0
            !bdzs = 0
            !bdzf = 0
            r = 0.0d0 

!$OMP PARALLEL DO

        do 40 k = 1,nz
        do 40 j = 1,ny
        do 40 i = 1,nx

           r(i,j,k) = 
    
     &     ( ( um(i,j,k) - um(i-1,j,k) ) / dx + ( vm(i,j,k) - vm(i,j-1,k) ) / dy + ( wm(i,j,k) - wm(i,j,k-1)) / dz ) / dt / leap +

     &     ( ( fx(i,j,k) - fx(i-1,j,k) ) / dx + ( fy(i,j,k) - fy(i,j-1,k) ) / dy + ( fz(i,j,k) - fz(i,j,k-1) ) / dz )

  40    continue

!$OMP END PARALLEL DO


d        sumr=0.0d0
d        rmax = r(1,1,1)
d        rmin = r(1,1,1)

d        do 50 k = 1, nz
d        do 50 j = 1, ny
d        do 50 i = 1, nx

d           rmax = max(rmax,r(i,j,k))
d           rmin = min (rmin,r(i,j,k))
d  50       sumr = sumr + r(i,j,k)*dx*dz*dy

d      	sumb=0.0
d       subx=0.0
d       suby=0.0
d	subz=0.0

d        do 60 k = 1, nz
d        do 60 j = 1,ny

d            subx = subx - bdxs(j,k)*dy*dz + bdxf(j,k)*dy*dz
d  60        sumb = sumb - bdxs(j,k)*dy*dz + bdxf(j,k)*dy*dz

d        do 62 k = 1,nz
d        do 62 i = 1,nx

d            suby = suby - bdys(i,k)*dz*dx + bdyf(i,k)*dx*dz
d  62        sumb = sumb  - bdys(i,k)*dz*dx + bdyf(i,k)*dx*dz

d       do 64 j = 1,ny
d       do 64 i = 1,nx

d           subz = subz  - bdzs(i,j)*dx*dy + bdzf(i,j)*dx*dy 
d 64        sumb = sumb  - bdzs(i,j)*dx*dy + bdzf(i,j)*dx*dy 

d        diff = sumr - sumb           !cyclic in x and y only

d        write(65,650) mtime,diff, sumr,sumb,subx,suby,subz
d 650    format ( 1x, i10,6e15.7)

        xl = nx*dx
        yl = ny*dy
        zl = nz*dz
        xs = 0.0d0
        ys = 0.0d0
        zs = 0.0d0

        lambda = ( -1.0 d-66 )     
        
         if( mtime.eq.0) CALL HS3CRI(xs,xl,nx,3,ys,yl,ny,3,zs,zl,nz,3,lambda,nx,ny,ierror,work)
     
             if(ierror.ne.0) then 

                  write (6,*) ierror
                  stop

             endif 
             qcheck = r
        CALL HS3CRT(bdxs,bdxf,bdys,bdyf,bdzs,bdzf,nx,ny,r,pterb,work)


        !call check_HS3CRT

            !write(6,*)'max and  min',  maxval ( r(1: nx, 1:ny,1:nz )), minval ( r(1: nx, 1:ny,1:nz ))

!$OMP PARALLEL DO
        do  k = 1, nz
        do  j = 1, ny
        do  i = 1, nx

          p(i,j,k) = r(i,j,k) - r(1,1,1)
                           
        enddo;enddo;enddo

!$OMP END PARALLEL DO

        sumr = 0.0
        do k = 1, nz
        do j = 1, ny
        do i = 1, nx  
          sumr =  p(i,j,k) + sumr
        enddo;enddo;enddo

        ave = sumr /real(nx)/real(ny)/real(nz)

        p = p - ave


        do 80 j = 1, ny
        do 80 i = 1, nx

             p(i,j,0)    = p(i,j,1)  - bdzs(i,j) * dz
 80          p(i,j,nz+1) = p(i,j,nz) + bdzf(i,j) * dz

        if (open_box .eq. 1) then 

        do  k = 1, nz
        do  j = 1, ny

              p(0,j,k) =    ( p(1,j,k) - bdxs(j,k) * dx ) 
              p(nx+1,j,k) = ( p(nx,j,k) + bdxf(j,k)* dx ) 

        enddo; enddo

        else 

        do  k = 1, nz
        do  j = 1, ny

              p(0,j,k) =    ( p(1,j,k) - bdxs(j,k) * dx ) * bcondx + ( 1 - bcondx ) * p(nx,j,k)
              p(nx+1,j,k) = ( p(nx,j,k) + bdxf(j,k)* dx ) * bcondx + ( 1 - bcondx ) * p(1,j,k)

        enddo; enddo
        endif

        do  i = 1, nx
        do  k = 1, nz

             p(i,0,k) =    ( p(i,1,k)  - bdys(i,k) * dy ) * bcondy + ( 1 - bcondy ) * p(i,ny,k)
             p(i,ny+1,k) = ( p(i,ny,k) + bdyf(i,k) * dy ) * bcondy + ( 1 - bcondy ) * p(i,1,k)

        enddo; enddo

      contains

        subroutine check_hs3crt

          real error(1:nx,1:ny,1:nz)
          error = 0.0

        do k = 2, nz-1
        do j = 2, ny-1
        do i = 2, nx-1

          error(i,j,k) = ((r(i+1,j,k)- r(i,j,k)) - (r(i,j,k)- r(i-1,j,k)))/dx/dx +  
     &      ((r(i,j+1,k)- r(i,j,k)) - (r(i,j,k)- r(i,j-1,k)))/dy/dy + 
     &      ((r(i,j,k+1)- r(i,j,k)) - (r(i,j,k)- r(i,j,k-1)))/dz/dz - qcheck(i,j,k)

        enddo;enddo;enddo

          call stat_3d(error,1,nx,1,ny,1,nz); call stat_3d(qcheck,1,nx,1,ny,1,nz)

        end subroutine check_hs3crt

                                            
        end subroutine pressure

! ----------------------------------------------------------------------------------------------
! ------------------------- SUBROUTINE REFRESH -------------------------------------------------
! ----------------------------------------------------------------------------------------------

       subroutine refresh

       use grid_module
       use velocity_module
       use constants_module
       include 'array_size.inc'
       integer i,j,k

       if ( leap .eq. 1 ) then

       do k = 0,nz+1
       do j = 0,ny+1
       do i = -1,nx+1

            u(i,j,k) = up(i,j,k)
            um(i,j,k) = u(i,j,k)

      enddo;enddo;enddo

       do k = 0, nz+1
       do j = -1, ny+1
       do i = 0, nx+1

            v(i,j,k) = vp(i,j,k)
            vm(i,j,k) = v(i,j,k)

       enddo;enddo;enddo

       do k = -1, nz+1
       do j = 0 ,ny+1
       do i = 0, nx+1

            w(i,j,k) = wp(i,j,k)
            wm(i,j,k) = w(i,j,k)

       enddo;enddo;enddo

       else if (leap .eq. 2) then

!!$OMP PARALLEL DO

       do  k = 0,nz+1
       do  j = 0,ny+1
       do  i = -1,nx+1
 
            u(i,j,k)  = u(i,j,k) + assel * (  up(i,j,k) - 2 * u(i,j,k) + um(i,j,k) )        
            um(i,j,k) = u(i,j,k)
!           ut(i,j,k) = ( u(i,j,k) + up(i,j,k) ) / 2.0d0        !i=-1, and nx+1 not defined for ut
            u(i,j,k) = up(i,j,k)

       enddo;enddo;enddo

!!$OMP END PARALLEL DO
!!$OMP PARALLEL DO

       do  k = -1, nz+1
       do  j = 0 ,ny+1
       do  i = 0, nx+1

            w(i,j,k) = w(i,j,k) + assel * (  wp(i,j,k) - 2 * w(i,j,k) + wm(i,j,k) )             
            wm(i,j,k) = w(i,j,k)
!           wt(i,j,k) = ( w(i,j,k) + wp(i,j,k) ) / 2.0d0       !k=-1, and nz+1 not defined for wt
            w(i,j,k) = wp(i,j,k)

       enddo;enddo;enddo
!!$OMP END PARALLEL DO

!!$OMP PARALLEL DO
       do  k = 0, nz+1
       do  j = -1, ny+1
       do  i = 0, nx+1

            v(i,j,k) = v(i,j,k) + assel * (  vp(i,j,k) - 2 * v(i,j,k) + vm(i,j,k) )
            vm(i,j,k) = v(i,j,k)
!           vt(i,j,k) = ( v(i,j,k) + vp(i,j,k) ) / 2.0d0       !j=-1, and ny+1 not defined for vt
            v(i,j,k) = vp(i,j,k)
      
       enddo;enddo;enddo
!!$OMP END PARALLEL DO
       endif

       return 
       end

!  ------------------------------------------------------------------------------------------------
!  ----------------- SUBROUTINE FILLVALUE ---------------------------------------------------------
!  ------------------------------------------------------------------------------------------------

       subroutine fill_value(var,nxs,nxe,nys,nye,nzs,nze)

       integer nxs,nxe,nys,nye,nzs,nze,i,j,k
       real var(nxs:nxe,nys:nye,nzs:nze),FillValue
       data FillValue / -1.0e+34 /       

           do i = nxs,nxe
              var(i,nys,nzs) = FillValue
              var(i,nys,nze) = FillValue
              var(i,nye,nzs) = FillValue
              var(i,nye,nze) = FillValue
           enddo

           do j =nys,nye
              var(nxs,j,nzs) = FillValue
              var(nxs,j,nze) = FillValue
              var(nxe,j,nzs) = FillValue
              var(nxe,j,nze) = FillValue
           enddo

           do  k =nzs,nze
              var(nxs,nys,k) = FillValue
              var(nxs,nye,k) = FillValue
              var(nxe,nys,k) = FillValue
              var(nxe,nye,k) = FillValue
           enddo

           return
           end        

	! )))))))))))))))))))))))))))))))))))))))))))))))))))
 
        subroutine gutenberg

        use velocity_module
        use tracer_module
        use viscosity_module
        use constants_module
        use signals_module
        use grid_module
        use source_module
        include 'array_size.inc'

      	write(6,*) 'nx = ', nx
      	write(6,*) 'ny = ', ny
      	write(6,*) 'nz = ', nz
        write(6,*) 'nt = ', nt
      	write(6,*) 'dx = ', dx
      	write(6,*) 'dy = ', dy      
      	write(6,*) 'dz = ', dz
      	write(6,*) 'dt = ', dt

        write(6,*) 'open_box = ', open_box
        write(6,*) 'ls = ', ls
        write(6,*) 'Cs_h =', cs_h
        write(6,*) 'Cs_v =', cs_v
        write(6,*) 'Prandtl No = ', prandtl

        write(6,*) 'bcondx,bcondy ',  bcondx, bcondy
        write(6,*) 'bottom_slip ',  bottom_slip 
        write(6,*) 'side_slip ',  side_slip

        write(6,*) 'omega = ', omega

!       write(6,*) 'b0_s = ', b0_s
!       write(6,'(a50,e14.5)') 'buoyancy flux from heat in units of m2/s3', b0_t
!       write(6,'(a50,e14.5)') 'buoyancy flux from salt in units of m2/s3', b0_s

        write(6,*) 'source_width = ', source_width
        write(6,*) 'source_length = ', source_length
!       write(6,*) 'rotation_freq = ', rotation_freq
!       write(6,*) 'fluid_height = ', fluid_height
        write(6,*) 'random_amp = ', random_amp
        write(6,*) 'random_amp1 = ', random_amp1
        write(6,*) 'source_on_time = ', source_on_time
        write(6,*) 'sponge_factor = ', sponge_factor                                                            
      	write(6,*) 'dQtdt =', dqtdt
        write(6,*) 'dQsdt =', dqsdt
        write(6,*) 'dQcdt=', dQcdt
        write(6,*) 'dQcc1dt =', dqcc1dt
        write(6,*) 'dQcc2dt=', dQcc2dt
        write(6,*) 'dQcc3dt =', dqcc3dt
        write(6,*) 'dQcc4dt=', dQcc4dt
!        write(6,*) 
!       write(6,*) 'source_diameter=', source_diameter, '(m)'
!       write(6,*) 'B0             =', b0 ,  ' (m^2/s^3)'
!       write(6,*) 'random_amp     =', random_amp
!       write(6,*)
        write(6,*) 'fc =', fc
        write(6,*) 'fc_y =', fc_y
        write(6,*) 'axmin  =',axmin
        write(6,*) 'aymin  =',aymin
        write(6,*) 'azmin  =',azmin
        write(6,*) 'azmax  =',azmax, 'if azmax = azmin, then mixing_height_scale matters not'
        write(6,*) 'kxmin  =',kxmin
        write(6,*) 'kymin  =',kymin
        write(6,*) 'kzmin  =',kzmin
        write(6,*) 'kzmax  =',kzmax, 'if bit_kz = 0, values given kxmin and kzmax matters not'
        write(6,*) 'mixing_height_scale =', mixing_height_scale, 'height in meters for bottom taper of viscosity' 
        !write(6,*) 'mixing_delta  =', mixing_delta, 'mixing_delta = 0(1) means no(yes) to Smagorinsky term in vertical viscosity'
        write(6,*) 'ismag_h',ismag_h
        write(6,*) 'ismag_v',ismag_v,' 0(1) is off(on) for vertical shear turb. generation terms in smagorinky' 
        write(6,*) 'bit2k = ',bit2k, 'bit2k = 0(1) means no(yes) to vertical diffusivity'
        write(6,*) 'assel = ',assel

      	write(6,*) 'bit2  =',bit2
        write(6,*) 'bit2k = ',bit2k
     	write(6,*) 'bit5  =',bit5
        write(6,*) 'bit5k  =',bit5k
      	write(6,*) 'bit4  =',bit4
      	write(6,*) 'bit3  =',bit3       
      	write(6,*) 'bit7  =',bit7
        write(6,*) 'bit8  =',bit8
        write(6,*) 'bit9  =',bit9
        write(6,*) 'bit10 =',bit10
      	write(6,*) 'ismol =',ismol

        return
        end
