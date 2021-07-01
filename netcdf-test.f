      module check_mod
         real, public:: a(0:10,0:10,0:10),pi
      end module check_mod

      program check

!pgf95 -r8 -Mextend -O3 netcdf-test.f -L/usr/local/netcdf/4.1.3/pgi123/lib -lnetcdf -L/usr/bin -lcurl

         use check_mod
         integer  cid1, cid2
         pi = 4.0*atan2(1.0,1.0)

         call write_snapshot(1,10,1,10,1,10, 'bob1.nc', 'open', cid1)
         call write_snapshot(1,10,1,10,1,10, 'bob2.nc', 'open', cid2)

         a = -12.0
         call write_snapshot(1,10,1,10,1,10, 'bob1.nc', 'write', cid1)

         do k= 0,10;do j=0,10;do i=0,10
              a(i,j,k) = sin(real(j)/3.0*pi)*cos(real(i)/5.0*pi)*real(k)
         enddo;enddo;enddo
         
         !Note only writing the 1:10 range of A not the full 0:10 range)
         call write_snapshot(1,10,1,10,1,10, 'bob2.nc', 'write', cid2)

         call write_snapshot(1,10,1,10,1,10, 'bob1.nc', 'close', cid1)
         call write_snapshot(1,10,1,10,1,10, 'bob2.nc', 'close', cid2)

         write(6,*) 'I created two netcdf files: bob1.nc and bob1.nc in this directory'

        end
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         subroutine write_snapshot (nxs, nxe, nys, nye, nzs, nze, file_name, action, cdf_id)

         use check_mod
         include '/usr/local/netcdf/4.1.3/pgi123/include/netcdf.inc' 

         character(len=*) action, file_name
         integer nxs, nxe, nys, nye, nzs, nze
         integer cdf_id, grid_id_a(4)
         real work(1:10), time

         integer dim_id_x, dim_id_y,dim_id_z,dim_id_t
         integer axis_id_x, axis_id_y,axis_id_z,axis_id_t, var_id_a, rcode

         if  (action .eq. 'WRITE'. or. action .eq. 'write') then
              write(6,*) 'write netcdf data'
              call write
          elseif (action .eq. 'OPEN' . or. action .eq. 'open' ) then
              write(6,*) 'open netcdf file'
              call open
          elseif (action .eq. 'CLOSE' .or. action .eq. 'close') then
              write(6,*) 'close netcdf file'
              call close

          else
                stop
          endif

        contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subroutine OPEN

           integer i

           do i=1, 10
              work(i)=i*2.5
           enddo
 
           write(6,*) file_name 

           cdf_id = nccre (file_name,ncclob,rcode)


           dim_id_t = ncddef(cdf_id,'time', ncunlim, rcode)
           dim_id_x = ncddef(cdf_id,'x_centers', 10, rcode)
           dim_id_y = ncddef(cdf_id,'y_centers', 10, rcode)
           dim_id_z = ncddef(cdf_id,'z_centers', 10, rcode)

!          write(6,*) 'rcode,cdf_id ', rcode,cdf_id, dim_id_x, dim_id_y, dim_id_z, dim_id_t

           grid_id_a = (/dim_id_x, dim_id_y, dim_id_z, dim_id_t/)

           axis_id_t = ncvdef(cdf_id,'time', ncdouble, 1, dim_id_t, rcode)

           axis_id_x  = ncvdef(cdf_id,'xlocations',  ncfloat, 1, dim_id_x, rcode)
           axis_id_y  = ncvdef(cdf_id,'ylocations',  ncfloat, 1, dim_id_y, rcode)
           axis_id_z  = ncvdef(cdf_id,'zlocations',  ncfloat, 1, dim_id_z, rcode)
!          axis_id_t  = ncvdef(cdf_id,'tlocations',  ncfloat, 1, dim_id_t, rcode)

           var_id_a = ncvdef(cdf_id,'a',ncfloat,4,grid_id_a,rcode)

!           write(6,*) ' grid_id_a',  grid_id_a

           call ncendf(cdf_id,rcode)

        end subroutine open

        subroutine CLOSE
                                                                                                                                    
          call  ncclos(cdf_id,rcode)
                                                                                                                                    
          return
                                                                                                                                    
        end subroutine close

        subroutine write

           time = 1.0d0  

           call ncvpt1(cdf_id,axis_id_t,1,time,rcode)         
           call ncvpt( cdf_id, axis_id_x, 1,10, real(work, kind=4 ), rcode )          
           call ncvpt( cdf_id, axis_id_y, 1,10, real(work, kind=4 ), rcode )
           call ncvpt( cdf_id, axis_id_z, 1,10, real(work, kind=4 ), rcode )

         call ncvpt (cdf_id, var_id_a, (/1,1,1,1/), (/10, 10, 10, 1/),
     &       real (a( 1: 10, 1: 10, 1:10), kind =4), rcode)

        end subroutine write

       end subroutine write_snapshot
