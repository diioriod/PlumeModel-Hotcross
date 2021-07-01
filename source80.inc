! pointed to by source.inc

            q_s =  0.0d0
            q_t =  0.0d0
            q_c =  0.0d0
            
            nx_offset =  0
            ramp_hrs = 2.0
                                                                                                                    
            nnx =  nint(source_length/dx)
            nny =  nint(source_width /dy)

            if (mod(nnx,2) .eq. 0) then
                nxs = nx/3 - (nnx/2 - 1)
                nxe = nx/3 + (nnx/2 - 1) + 1
            else
                nxs = nx/3 - (nnx/2 - 1) - 1
                nxe = nx/3 + (nnx/2 - 1) + 1
            endif

            if (mod(nny,2) .eq. 0) then
                nys = ny/2 - (nny/2 - 1)
                nye = ny/2 + (nny/2 - 1) + 1
            else
                nys = ny/2 - (nny/2 - 1) - 1
                nye = ny/2 + (nny/2 - 1) + 1
            endif

              isource_map = 0
              

              wee = .025                  !the randomness is diluted by background
!              write(6,*) 'nxs,nxe,nys,nye: ',  nxs,nxe,nys,nye

            do  j = nys, nye
            do  i = nxs, nxe

                isource_map(i,j) = 1

                call random_number(randnum1) 
                call random_number(randnum2)
                call random_number(randnum3)

               !write(6,*) 'random numbers:',randnum1,randnum2,randnum3

               q_t(i,j,nz) = dQtdt *  min (1.0, mtime * dt / (ramp_hrs * 3600.0)  ) * ( 1.0 + wee * (2.0 * (randnum1 -.5)) )
               q_s(i,j,nz) = dQSdt *  min (1.0, mtime * dt / (ramp_hrs * 3600.0)  ) * ( 1.0 + wee * (2.0 * (randnum2 -.5)) )
               q_c(i,j,nz) = dqcdt *  min (1.0, mtime * dt / (ramp_hrs * 3600.0)  ) * ( 1.0 + wee * (2.0 * (randnum3 -.5)) )
               q_c(i,j,nz) = 1.0

               !write(1776,'(4f15.7)') ( 1.0 + wee * (2.0 * (randnum1 -.5)) ), randnum1,randnum2,randnum3 
                
            enddo;enddo

!                 write(6,*) 'stat_3d', mtime, dQtdt,  min (1.0, mtime * dt / (ramp_hrs * 3600.0)  )
                 !write(6,*)  dQtdt, dQsdt, dQcdt
!                 call view_2d(isource_map*1.0,0,nx+1,0,ny+1,0)
!                 stop 89












