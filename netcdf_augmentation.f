!       HotCross Distribution - Version 1.0
!       tag: subroutine ../vents2/lavelle/NewCross/netcdf_augmentation.f - version 6.4.97

        subroutine start_n_count ( start, s1, s2, s3, s4, count, c1, c2, c3, c4 )

         integer count(4), start(4), c1, c2, c3, c4, s1, s2, s3, s4

         start(1) = s1
         start(2) = s2
         start(3) = s3
         start(4) = s4

         count(1) = c1
         count(2) = c2
         count(3) = c3
         count(4) = c4

         return
         end

! --------------------------------------------------------------------------------------------
                                     
!        SUBROUTINE loader ( u, work_f, nxs, nxe, nys, nye,  nzs, nze )

!             integer nxs, nxe, nzs, nze, nys, nye, i, j, k, l
!             real * 4 work_f ( ( nxe - nxs + 1 ) * ( nye - nys + 1 ) * ( nze - nzs + 1 ) )
!             real  u ( nxs: nxe, nys: nye, nzs: nze ), missing_value
!             data  missing_value / -1.0e+34 /

!                 l = 0

!             do k = nzs, nze
!             do j = nys, nye
!             do i = nxs, nxe

!                 l = l + 1
!		 work_f( l ) =  u( i, j, k )  

!             enddo
!             enddo
!             enddo

!          return

!          end

! ----------------------------------------------------------------------------------
        subroutine define_grid ( grid_id, xdim_id, ydim_id, zdim_id, tdim_id )

                integer grid_id(4), xdim_id, ydim_id, zdim_id, tdim_id

                grid_id(1) = xdim_id
                grid_id(2) = ydim_id
                grid_id(3) = zdim_id
                grid_id(4) = tdim_id

                return
        end

! ----------------------------------------------------------------

         subroutine load_axis (cdf_id, axis_id, nx, xs, xinc)

         integer cdf_id, axis_id, nx, i, rcode
         real  xs, xinc,  axis_data (1:nx)

             do 10 i = 1,nx
  10              axis_data(i) = xinc * ( i - 1 )  + xs

             call ncvpt( cdf_id, axis_id, 1, nx, real(axis_data,kind=4), rcode )

             return
             end