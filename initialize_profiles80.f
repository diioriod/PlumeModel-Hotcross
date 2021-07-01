!      HotCross Distribution - Version 1.0

       subroutine initialize_profiles 

!      tag: subroutine ../vents2/lavelle/NewCross/initialize_profiles.f - version 6.4.97

       use viscosity_module
       use velocity_module
       use constants_module
       use grid_module
       use tracer_module
       use signals_module
       use background_module
       include 'array_size.inc'

       integer i,j,k,background_number,isw,ndata,kup
       real height_scale, flux_t, flux_s, flux_c
       real  zz,  rho_rho, inflexion_to_bottom , spread_fraction
       real zdata(1:3000), sdata(1:3000), tdata(1:3000)
       external rho_rho

!         go to 4000   !activate to read in real data 

!         open ( unit = 1010, status='unknown', file ='site_stratification_fromobs.inc')
!         open ( unit = 1010, status='unknown', file ='site_stratification_linear.inc')

          sdata= 0.0; tdata=0.0; ndata = 0

          do k= 1,5 
              read(1010,*) !skip header
          enddo
             
          do k = 1, 3000
            read (1010,'(3F15.7)', end=999, err = 998) zdata(k), tdata(k), sdata(k) 
            write(666,*) zdata(k), tdata(k), sdata(k) 
            ndata = ndata + 1
         enddo
 998       stop 'error reading site_stratification_data.inc'
 999     continue

4000      continue

          isw = 0


       do k = 0, nz 

          zz = ( zstart - dz/2 + k * dz )

          z_centers(k) = zz

          call locator(zdata, ndata, zz, kup)

          t_bkg(k) = tdata(kup) +
     &                 ( tdata(kup+1) - tdata(kup)) /( zdata(kup+1) - zdata(kup)) * ( zz - zdata(kup))

          s_bkg(k) = sdata(kup) +
     &                 ( sdata(kup+1) - sdata(kup)) /( zdata(kup+1) - zdata(kup)) * ( zz - zdata(kup))


          azzmin(k) = ( azmin + ( azmax - azmin ) * exp ( - ( ( zend -zstart ) - k * dz ) / mixing_height_scale ) )
          write(6001,'(2i10,4F15.7)') nz, kup, zz, t_bkg(k), s_bkg(k),azzmin(k)

       enddo

          !write(6,*)  azmin,azmax,zend,zstart,mixing_height_scale; stop

          t_bkg(0)=t_bkg(1);t_bkg(nz+1)=t_bkg(nz)
          s_bkg(0)=s_bkg(1);s_bkg(nz+1)=s_bkg(nz)

          z_centers(0) = zstart
          z_centers(nz) = zend                      !false z_center, but used only in rho_rho

          kzmin  = azmin / prandtl
          kzmax  = azmax / prandtl

          kxmin = axmin / prandtl
          kymin = aymin / prandtl
                                         
          bit_kz = 1                     ! to preserve T and S background profile shapes
                                         ! in absence of the following heroic efforts to preserve them
                                         ! e.g. when t_bkg(k) and s_bkg(k) are chosen without
                                         ! making them consistent with azzmin(k) shape 
                                         ! set kzt = 0 using bit_kz = 0

          go to 4002 
          flux_t = azzmin(0) / prandtl * ( t_bkg(1) - t_bkg(0) ) / dz    !only t_bkg(1), t_bkg(0) coming out of background
          flux_s = azzmin(0) / prandtl * ( s_bkg(1) - s_bkg(0) ) / dz
          flux_c = azzmin(0) / prandtl * ( c_bkg(1) - c_bkg(0) ) / dz

          write(6,*) 'flux_t,flux_s,flux_c', flux_t, flux_s, flux_c

          do k = 1, nz +1

!!             t_bkg(k) = t_bkg(0) +  k * dz * flux_t / kzmin  -  mixing_height_scale / kzmin * flux_t * log( (kzmin + 
!!     &                     ( kzmax - kzmin) * exp ( - ( ( zend - zstart) - ( k -.5 ) * dz ) / mixing_height_scale ) ) / 
!!     &                     ( kzmin + ( kzmax - kzmin )* exp ( - ( ( zend - zstart)  +.5 * dz ) / mixing_height_scale ) ) ) 

              t_bkg(k) = t_bkg(k-1) + flux_t * dz / azzmin(k-1) * prandtl

 !!            s_bkg(k) = s_bkg(0) +  k * dz * flux_s / kzmin  -  mixing_height_scale / kzmin * flux_s * log( (kzmin +
 !!    &                     ( kzmax - kzmin) * exp ( - ( ( zend - zstart) - ( k -.5 ) * dz ) / mixing_height_scale ) ) /
 !!   &                     ( kzmin + ( kzmax - kzmin )* exp ( - ( ( zend - zstart)  +.5 * dz ) / mixing_height_scale ) ) )

               s_bkg(k) = s_bkg(k-1) + flux_s * dz / azzmin(k-1) * prandtl

!!!               c_bkg(k) = c_bkg(0) +  k * dz * flux_c / kzmin  - mixing_height_scale / kzmin * flux_c * log( (kzmin +
!!!     &                     ( kzmax - kzmin) * exp ( - ( ( zend - zstart) - ( k -.5 ) * dz ) / mixing_height_scale ) ) /
!!!     &                     ( kzmin + ( kzmax - kzmin )* exp ( - ( ( zend - zstart)   +.5 * dz ) / mixing_height_scale ) ) )

              c_bkg(k) = c_bkg(k-1) + flux_c * dz / azzmin(k-1) * prandtl
                                                                                                                   
        enddo
4002    continue

!          bit_kz = 0

           s_bkg ( 0 ) = s_bkg ( 1)
           s_bkg ( nz+1) = s_bkg ( nz)
           t_bkg ( 0 ) = t_bkg ( 1)
           t_bkg ( nz+1) = t_bkg ( nz)
           c_bkg ( 0 ) = c_bkg ( 1)
           c_bkg ( nz+1) = c_bkg ( nz)

              !write(6,*) ' initial profiles', prandtl
              !do k = 0, nz+1
              !            write(6002,'(4e14.7,i10)')  t_bkg(k), s_bkg(k), c_bkg(k), azzmin(k), k
              !enddo

              !stop

!              write(6,*) ' K(dT/dz) profiles'
!               do k = 1, nz+1
!                          write(6,'(e14.7,i10)')  (t_bkg(k)-t_bkg(k-1))* azzmin(k-1),k
!              enddo

         do k = 0, nz+1
            write(59,'(4e14.5)') Z_CENTERS(K), s_bkg(k),t_bkg(k),c_bkg(k)            
         enddo

         do i = 0, nx+1
         do j = 0, ny+1
         do k = 0, nz+1

            t(i,j,k) = t_bkg(k)
            s(i,j,k) = s_bkg(k)
            c(i,j,k) = c_bkg(k)

!           tm(i,j,k) = t_bkg(k)
!           sm(i,j,k) = s_bkg(k)
!           cm(i,j,k) = c_bkg(k)

         enddo;enddo;enddo

         do k = 0, nz+1

            rho_bkg (k) = rho_rho (s_bkg(k), t_bkg(k), isw*z_centers(k)/10)

         enddo

      	return
	end

        subroutine locator(xedges, nx, xpoint, i )

! a modified locate.f code for xedges indices running 1:nx
! subroutine to locate cell index, from numerical recipes, p.111

       integer nx, i
       real xedges(1:nx), xpoint
       integer ilo, ihi, imid

       ilo = 1
       ihi = nx

  10   if( (ihi-ilo) .gt. 1) then

          imid = (ihi+ilo)/2
          if(xpoint .ge.  xedges(imid)) then
             ilo = imid
          else
             ihi = imid
          endif
        go to 10
        endif

        if (xpoint .eq. xedges(1)) then
            i = 1
        else if (xpoint .eq. xedges(nx)) then
            i = nx - 1
        else
            i = ilo
        endif

        return

        end subroutine locator
