                                                                                                                   
        subroutine stat_3d(s,nx_s,nx_e,ny_s,ny_e,nz_s,nz_e)!,field_name)
                                                                                                                                    
!  subroutine to stat a  3-d array with ferret
                                                                                                                                    
!       character field_name(20)
                                                                                                                            
        integer nx_s,nx_e,ny_s,ny_e,nz_s,nz_e
        integer i,j,k

        real s(nx_s:nx_e,ny_s:ny_e,nz_s:nz_e), s_ave, sum , std, missing_value
        integer number

        missing_value = -1.0E+34

          sum = 0.0d0
          number= 0

       do k = nz_s,nz_e
       do j = ny_s,ny_e
       do i = nx_s,nx_e

          if ( s(i,j,k) .gt. 0.98 * missing_value) then 

             sum = s(i,j,k) + sum
             number= number+1

          endif

       enddo;enddo;enddo

!         write(6,*) 'in stat', sum
         
          s_ave = 0.0d0                                                                                                           
          if( number .gt. 0) s_ave = sum / real(number)
          sum = 0.0d0
                                                                                                                       
       do k = nz_s, nz_e
       do j = ny_s, ny_e
       do i = nx_s, nx_e

          if ( s(i,j,k) .gt. 0.99* missing_value)  then
             sum = (s(i,j,k)-s_ave)**2 + sum
          endif                                                                                            
                                                                                                                        
       enddo;enddo;enddo
          
          std = 0.0d0                                                                                         
          if( number .gt. 0.0) std = sqrt( sum /  real(number) )
 
!         write (6,*) 'number of bad values;', (nx_e-nx_s+1)*(ny_e-ny_s+1)*(nz_e-nz_s+1)-number

          write(6,'(2d18.10,i10)')   s_ave, std , (nx_e-nx_s+1)*(ny_e-ny_s+1)*(nz_e-nz_s+1)-number!, 's_ave, std, # bad:',       

                 
        return
        end
