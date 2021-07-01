
         c0 = 0.0d0 

         t_bkg ( k ) =  1.6341 + 3.08e-4 * ( zend - zz)               !between 1000 amd 2300 at Cleft
         s_bkg ( k ) =  34.6155 - 5.5e-5 * ( zend - zz)               !bewteen 1000 and cc2300 at Cleft

!        t_bkg ( k ) =  1.66d0 + 9.7d-4 *  ( zend - zz)                  !from Scott  
!        s_bkg ( k ) =  34.614d0 - 1.0d-4 * ( zend - zz)                !from scott

          t_bkg ( k ) =  1.673 + 6.666e-4 * ( zend - zz)               !above 2175 at Dante
          s_bkg ( k ) =  34.612 -1.135e-4 * ( zend - zz)               !above 2175 at Dante 

          !These values are a fit to at1536ba.cnv data in the depth range 2175:1800

!following are bogus profiles

!         t_bkg ( k ) =  1.6339 + 3.0921d-4 * ( zend - zz)               !above 2175 at Dante
!         s_bkg ( k ) =  34.616 - 5.4523d-5 * ( zend - zz)

!         t_bkg ( k ) =  1.66d0 + 2*9.7d-4 *  ( zend - zz)              !test only


         c_bkg (k) = c0

         azzmin(k) = ( azmin + ( azmax - azmin ) * exp ( - ( ( zend -zstart ) - k * dz ) / mixing_height_scale ) )  
						                    !background vertical mixing -- non-constant profile!
							            !allows for more realistic near-bottom T and S profiles
                                                                    !and for Ekman bottom boundary layer
                                                                    !if azzmin is not constant, the t_bkg and s_bkg
                                                                    !profiles stated above are changed to make them consistent
                                                                    !with vertical diffusive mixing

         spongez = 0.0d0

!        azzmin(k) = azmin                                           
!        t0 = 2.0
!        s0 = 35.0
!        t_bkg ( k ) = t0
!        s_bkg ( k ) = s0

!!       s_bkg ( k) = s0 + N2 / g / alpha * zz
!!       t_bkg (k) =   0.372d-3 * zz + 7.74 * exp ( - zz / 1020.0 )
!!       s_bkg (k) = ( 0.372E-3 * zz + 7.74 * exp ( - zz / 1020.0 ) ) * (-0.16967) + 34.91   ! ridge.f profiles

!!       t_bkg ( k) = 23.0 - (23.0 -20.0) * (zstart - dz/2 + k * dz)/ (zstart + zend + dz)
!!       s_bkg ( k) = 0.0 - ( 0.0  - 2.0 ) * ( zstart - dz/2 + k * dz )/ (zstart + zend + dz )
!!       N2 = g ( drho/dz)/rho0 = g*alpha*dS/dz so dS/dz = N2/g/alpha --->  s(z) = s(0) + n2/g/alpha*dz
!!       s_bkg ( k) = s0 + N2 / g / alpha * ( zstart - dz/2 + k * dz )
!!       s_bkg ( k) = s0 + N2 / g / alpha * zz
