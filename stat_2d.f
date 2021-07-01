      subroutine stat_2d (s, nx_s, nx_e, ny_s, ny_e)
                                                                                                                               
!  subroutine to stat a  3-d array with ferret
                                                                                                                                 
!       character field_name(20)
                                                                                                                                  
        integer nx_s,nx_e,ny_s,ny_e,i,j 
        real s(nx_s:nx_e,ny_s:ny_e)

        real sum,save,std    
                                                                                                                        
          sum = 0.0

       do j = ny_s+1,ny_e-1
       do i = nx_s+1,nx_e-1
          sum = s(i,j) + sum
       enddo;enddo
                                                                                                            
          save = sum / real(nx_e- nx_s-1) / real(ny_e- ny_s-1) 
                                                                                                                                    
          sum = 0.0
                                                                                                                                    
       do j = ny_s+1,ny_e-1
       do i = nx_s+1,nx_e-1
                                                                                                                                    
          sum = (s(i,j)-save)**2 + sum
                                                                                                                                    
       enddo;enddo
                                                                                                                                    
          std = sqrt( sum / ( real(nx_e- nx_s-1) * real(ny_e- ny_s-1) -1 ))                                                         
          write(6,*) 'save,std:',real(save),real(std)                                                     
                           
        return
        end
