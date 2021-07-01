      !     HotCross Distribution - Version 1.0

      module tracer_smol_module

         private nx, ny, nz

         include 'array_size.inc'

         real, allocatable, public, save ::  uc(:,:,:), vc(:,:,:), wc(:,:,:)
         real, allocatable, public, save ::  cp(:,:,:) , cc(:,:,:)
         real, allocatable, public, save ::  ud(:,:,:), vd(:,:,:), wd(:,:,:)

      end module tracer_smol_module

      ! *******************************************************

      module geometry_module

         private nx,ny,nz
         include 'array_size.inc'

         real, public:: mp_k( 0:nx+1, 0:ny+1,-1:nz+1 ), np_k( 0:nx+1, 0:ny+1, -1:nz+1 )
         real, public:: mp_( 0:nx+1, 0:ny+1 ), np_( 0:nx+1, 0:ny+1 )
         real, public:: mpx_( -1:nx+1, 0:ny+1 ), npy_( 0:nx+1, -1:ny+1 )
         real, public:: hzmnx(-1:nx+1, 0:ny+1,0:nz+1), hznx(-1:nx+1, 0:ny+1,0:nz+1)
         real, public:: hzmny(0:nx+1, -1:ny+1,0:nz+1), hzmy(0:nx+1, -1:ny+1,0:nz+1) 
         real, public:: hzmn(0:nx+1,0:ny+1,0:nz+1)
         real, public:: hzp(0:nz+1), hzw(-1:nz+1)
         integer,public :: bit_z_hp(0:nx+1,0:ny+1,0:nz+1), bit_z_hx(-1:nx+1,0:ny+1,0:nz+1), bit_z_hy(0:nx+1,-1:ny+1,0:nz+1),
     &          bit_z_hz(0:nx+1,0:ny+1,-1:nz+1)

!        contains

!             subroutine initialize_dummy_variables 
!                 bit_z_hp =1;  bit_z_hx =1; bit_z_hy=1; bit_z_hz=1
!                 mp_= 1.0; np_ = 1.0; mp_k = 1.0;  np_k= 1.0; mpx_ = 1.0 ;  npy_ = 1.0
!                 hzmn= 1.0; hzmnx=1.0d0; hzmny=1.0d0; hznx = 1.0d0;  hzmy = 1.0d0
!                 hzp = 1.0d0;  hzw= 1.0d0
!             end subroutine initialize_dummy_variables

        end module geometry_module

        module constants_module

            real, public:: g, pi, omega, rho0, eps, Cpp, Cs, eps10, fc, fc_y, assel, missing_value, cs_h, cs_v

        contains

            subroutine initialize_constants

               pi = 4.0d0 * atan2 ( 1.0d0 , 1.0d0 )
               g = 9.81d0
               omega =  2.0d0 * pi / 24.0 / 3600.0d0
               eps   =  1.0d-10
               eps10   =  1.0d-10
               Cs = 0.2d0
               Cpp = 4200.0d0 
               rho0 = 1027.65d0
               assel = .01d0   !.15
               missing_value = -1.0e+34

            end subroutine initialize_constants

            real function rad (xx)
               real xx
               rad = pi / 180.0d0 * xx

            end function rad



        end module constants_module


        ! *********************************************************

        module grid_module

        private nx,ny,nz
        include 'array_size.inc'

        real, public:: dx, dy, dz, ds
        real, public::  xstart, ystart, zstart, zend
        real, public:: time,dt
        integer, public:: mtime, leap, nt, modulo_x, modulo_y
        real, public:: z_centers(0:nz+1)
        end module grid_module

        ! ))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))

        module signals_module

        real,    public, save:: sponge_factor, tau, c_d
        integer, public, save:: cdf_ids(100)
        integer, public, save:: bit1,bit2,bit3,bit4,bit5,bit6,bit7,bit8,bit9, bit10, bit11,bit12,bit13,bit14,bit15
        integer, public, save:: bit5k, bit2k
        integer, public, save:: delta, ismol, bcondx, bcondy, bottom_slip, side_slip, slipon, back, open_box 
        integer, public, save:: ismol_hx,ismol_hy,ismol_v
        integer, public, save:: toggle(13)
        integer, public, save:: number_of_output_blocks
        integer, public, save:: ismag_vertical, ismag_v, ismag_h
        integer, public, save:: spongex,spongey,spongee
        integer, public, save:: bdbit(1:4)
        integer,public, save:: nzls,nzle,nxls,nxle,nyls,nyle, ncapture1, ncapture2,ncapstart, 
     &                     ncapture3, ncapture4, ncapture5, ncapture6
        character * 90, public, save:: output_file_1, output_file_2, output_file_3,output_file_4


        end module signals_module

        ! ))))))))))))))))))))))))))))))))))))))))))))

        module force_module

        private
        include 'array_size.inc'

        real, public::  fx(-1:nx+1,0:ny+1,0:nz+1), fy(0:nx+1,-1:ny+1,0:nz+1), fz(0:nx+1,0:ny+1,-1:nz+1)
        real, public:: Fbx(-1:nx+1,0:ny+1,0:nz+1),Fby(0:nx+1,-1:ny+1,0:nz+1),Fbz(0:nx+1,0:ny+1,-1:nz+1)

        end module force_module

        ! ))))))))))))))))))))))))))))))))))))))))))))

        module velocity_module

        private
        include 'array_size.inc'

        !real, public:: u(-1:nx+1,0:ny+1,0:nz+1), up(-1:nx+1,0:ny+1,0:nz+1), um(-1:nx+1,0:ny+1,0:nz+1)
        !real, public:: v(0:nx+1,-1:ny+1,0:nz+1), vp(0:nx+1,-1:ny+1,0:nz+1), vm(0:nx+1,-1:ny+1,0:nz+1)
        !real, public:: w(0:nx+1,0:ny+1,-1:nz+1), wp(0:nx+1,0:ny+1,-1:nz+1), wm(0:nx+1,0:ny+1,-1:nz+1)

        real, allocatable, public,   save:: u(:,:,:), up(:,:,:), um(:,:,:)
        real, allocatable, public,   save:: v(:,:,:), vp(:,:,:), vm(:,:,:)
        real, allocatable, public,   save:: w(:,:,:), wp(:,:,:), wm(:,:,:)
        real, allocatable, public,   save:: p(:,:,:)

        real, public:: u_bkg(0:nz+1)    !,ut(-1:nx+1,0:ny+1,0:nz+1), dudt(-1:nx+1,0:ny+1,0:nz+1)
        real, public:: ubkg(-1:nx+1, 0:ny+1),  ubkg_p(-1:nx+1, 0:ny+1),  ubkg_m(-1:nx+1, 0:ny+1)
        real, public:: vbkg( 0:nx+1,-1:ny+1),  vbkg_p( 0:nx+1,-1:ny+1),  vbkg_m( 0:nx+1,-1:ny+1)
        real, public::  v_bkg(0:nz+1)    !,vt(0:nx+1,-1:ny+1,0:nz+1), dvdt(0:nx+1,-1:ny+1,0:nz+1)

        real, public:: wsource
                             !,wt(0:nx+1,0:ny+1,-1:nz+1) , dwdt(0:nx+1,0:ny+1,-1:nz+1), wsource
        real, public:: ebkg( 0:nx+1, 0:ny+1),  ebkg_p( 0:nx+1, 0:ny+1),  ebkg_m( 0:nx+1, 0:ny+1)
        real, public:: alphau(-1:nx+1,0:ny+1), alphav(0:nx+1,-1:ny+1), alphae(0:nx+1,0:ny+1), alphaw(-1:nz+1), 
     &                alphaee(0:nx+1,0:ny+1,0:nz+1), alphae_complement(0:nx+1,0:ny+1)!,alphaw_perimeter(0:nx+1,0:ny+1)

        real, public:: spongez ( 0:nx+1, 0:ny+1, -1:nz+1 )
        real, public:: pgradx(0:nz+1), pgrady(0:nz+1), pgradx_con, pgrady_con
        real, public:: random_amp1 

        contains

           subroutine velocity_module_allocate

              !allocate (p(0:nx+1,0:ny+1,0:nz+1))
              !allocate (u(-1:nx+1,0:ny+1,0:nz+1), up(-1:nx+1,0:ny+1,0:nz+1), um(-1:nx+1,0:ny+1,0:nz+1))
              !allocate (v(0:nx+1,-1:ny+1,0:nz+1), vp(0:nx+1,-1:ny+1,0:nz+1), vm(0:nx+1,-1:ny+1,0:nz+1))
              !allocate (w(0:nx+1,0:ny+1,-1:nz+1), wm(0:nx+1,0:ny+1,-1:nz+1), wp(0:nx+1,0:ny+1,-1:nz+1))
 
           end subroutine velocity_module_allocate

        end module velocity_module

        ! )))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))


        module source_module

        private
        include 'array_size.inc'

        real, public:: b0_s, b0_t, source_width, source_length, rotation_freq, fluid_height, random_amp, source_on_time 

        end module source_module
        ! )))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))

        module background_module

        private
        include 'array_size.inc'

        real, public:: alpha, beta, N2, s0, t0, c0                  

        end module background_module
     
        !)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))

        module tracer_module

        private
        include 'array_size.inc'

        real, public:: s_bkg(0:nz+1), t_bkg(0:nz+1), c_bkg(0:nz+1), s_ref, t_ref, c_ref
        
        real, allocatable, public, save:: kxt(:,:,:), kyt(:,:,:), kzt(:,:,:)
        real, allocatable, public, save:: t(:,:,:), Q_t(:,:,:)
        real, allocatable, public, save:: s(:,:,:), Q_s(:,:,:)
        real, allocatable, public, save:: c(:,:,:), Q_c(:,:,:)
        real, allocatable, public, save:: cc1(:,:,:), Q_cc1(:,:,:)
        real, allocatable, public, save:: cc2(:,:,:), Q_cc2(:,:,:)
        real, allocatable, public, save:: cc3(:,:,:), Q_cc3(:,:,:)
        real, allocatable, public, save:: cc4(:,:,:), Q_cc4(:,:,:)

        !real, public:: kxt(-1:nx+1,0:ny+1,0:nz+1), kyt(0:nx+1,-1:ny+1,0:nz+1), kzt(0:nx+1,0:ny+1,-1:nz+1)
        !real, public:: t(0:nx+1,0:ny+1,0:nz+1), Q_t(0:nx+1,0:ny+1,0:nz+1)
        !real, public:: s(0:nx+1,0:ny+1,0:nz+1), Q_s(0:nx+1,0:ny+1,0:nz+1)
        !real, public:: c(0:nx+1,0:ny+1,0:nz+1), Q_c(0:nx+1,0:ny+1,0:nz+1)
         real, public:: t_total, dQtdt, tsource,tmin,tmax
         real, public:: s_total, dQSdt, ssource,smin,smax
         real, public:: c_total, dQcdt, ws_c,   cmin,cmax
         real, public:: t_total0,s_total0,c_total0



!       real, public:: spongex ( 0:nx+1, 0:ny+1, 0:nz+1 ), spongey ( 0:nx+1, 0:ny+1, 0:nz+1 ), spongez ( 0:nx+1, 0:ny+1, 0:nz+1 )
        real, public:: rhop(0:nx+1,0:ny+1,0:nz+1), rho_bkg(0:nz+1), rrho(0:nx+1,0:ny+1,0:nz+1), ref_pressure
        real, public:: source_diameter, b0

!        real, public, dimension(:,:,:), allocatable :: cc1, cc2, cc3, cc4, Q_cc1, Q_cc2, Q_cc3, Q_cc4

!       real, public:: cc1(0:nx+1,0:ny+1,0:nz+1),q_cc1(0:nx+1,0:ny+1,0:nz+1)
!       real, public:: cc2(0:nx+1,0:ny+1,0:nz+1),q_cc2(0:nx+1,0:ny+1,0:nz+1)
!       real, public:: cc3(0:nx+1,0:ny+1,0:nz+1),q_cc3(0:nx+1,0:ny+1,0:nz+1)
!       real, public:: cc4(0:nx+1,0:ny+1,0:nz+1),q_cc4(0:nx+1,0:ny+1,0:nz+1)

        real, public:: cc1_bkg(0:nz+1), cc1_ref, cc1_total, dQcc1dt, cc1_ws, cc1_lambda 
        real, public:: cc2_bkg(0:nz+1), cc2_ref, cc2_total, dQcc2dt, cc2_ws, cc2_lambda
        real, public:: cc3_bkg(0:nz+1), cc3_ref, cc3_total, dQcc3dt, cc3_ws, cc3_lambda
        real, public:: cc4_bkg(0:nz+1), cc4_ref, cc4_total, dQcc4dt, cc4_ws, cc4_lambda

        integer, public:: isource_map(0:nx+1,0:ny+1), p_end(0:nx+1,0:ny+1)

       contains

           subroutine tracer_module_allocate

!              allocate (kxt(-1:nx+1,0:ny+1,0:nz+1), kyt(0:nx+1,-1:ny+1,0:nz+1), kzt(0:nx+1,0:ny+1,-1:nz+1))
!              allocate (s(0:nx+1,0:ny+1,0:nz+1), Q_s(0:nx+1,0:ny+1,0:nz+1))
!              allocate (t(0:nx+1,0:ny+1,0:nz+1), Q_t(0:nx+1,0:ny+1,0:nz+1))
!              allocate (c(0:nx+1,0:ny+1,0:nz+1), Q_c(0:nx+1,0:ny+1,0:nz+1))

              !if (toggle(9)  .eq. 1) allocate (cc1(0:nx+1,0:ny+1,0:nz+1), Q_cc1(0:nx+1,0:ny+1,0:nz+1))
              !if (toggle(10) .eq. 1) allocate (cc2(0:nx+1,0:ny+1,0:nz+1), Q_cc2(0:nx+1,0:ny+1,0:nz+1))
              !if (toggle(11) .eq. 1) allocate (cc3(0:nx+1,0:ny+1,0:nz+1), Q_cc3(0:nx+1,0:ny+1,0:nz+1))
              !if (toggle(12) .eq. 1) allocate (cc4(0:nx+1,0:ny+1,0:nz+1), Q_cc4(0:nx+1,0:ny+1,0:nz+1))

           end subroutine tracer_module_allocate

        end module tracer_module

        !))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))

        module viscosity_module

        private
        include 'array_size.inc'

        real, public::  axu(0:nx+1,0:ny+1,0:nz+1),ayu(0:nx,0:ny,0:nz+1),    azu(-1:nx+1,0:ny+1,-1:nz+1)
        real, public::  axv(0:nx,0:ny,0:nz+1),    ayv(0:nx+1,0:ny+1,0:nz+1),azv(0:nx+1,0:ny,-1:nz+1)
        real, public::  axw(0:nx,0:ny+1,0:nz),    ayw(0:nx+1,0:ny,0:nz),    azw(0:nx+1,0:ny+1,0:nz+1)
        real, public::  xcrossterms(0:nx+1,0:ny+1,0:nz+1)

        real, public::  Ri(0:nx+1,0:ny+1,0:nz+1), ai(0:nx+1,0:ny+1,0:nz+1), azzmin(0:nz+1)
!        real, public::  Ritheta(0:nx+1,0:ny+1,0:nz+1)
        real, public::  ss_v(0:nx+1,0:ny+1,0:nz+1), ss(0:nx+1,0:ny+1,0:nz+1)
        real, public::  epsv(0:nx+1,0:ny+1,0:nz+1), epsh(0:nx+1,0:ny+1,0:nz+1), chi_v(0:nx+1,0:ny+1,0:nz+1), 
     &                  chi_h(0:nx+1,0:ny+1,0:nz+1)
  	real, public::  axmin, axmax, aymin, aymax, azmin, azmax, prandtl, ls, kzmin, kzmax, kxmin, kymin
        real, public::  mixing_height_scale
        integer, public:: bit_kz, mixing_delta                      

      	end module viscosity_module

       !))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))

        subroutine paws
                                                                                                                         
            write(6,*) 'press enter to continue'
            read(5,*)
                                                                                                                   
        end subroutine paws














