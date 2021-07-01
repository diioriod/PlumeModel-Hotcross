!       e.g. call  view_3d(temp,0,21,0,33,0,55)

        subroutine view_3d (a,nx_s,nx_e,ny_s,ny_e,nz_s,nz_e)!,field_name)

!  subroutine to view 3-d array with ferret
!  ferret> go view.jnl
!  ferret> exit

!       character field_name(20)
        integer nx_s,nx_e,ny_s,ny_e,i,j,k, nz_s,nz_e
        real a(nx_s:nx_e,ny_s:ny_e,nz_s:nz_e)

        rewind(71)
        write(71,1801) nx_s,nx_e,ny_s,ny_e,nz_s,nz_e
1801    format(6i10)

!        rewind(72)
!        write(72,1802) field_name
!1802    format (a20)         
!  	write(6,*) field_name

        rewind (70)
        write(70,1800) ( ( ( ( a(i,j,k) ), i=nx_s,nx_e), j=ny_s,ny_e), k = nz_s,nz_e)
1800    format(1x,e15.7)
 
        rewind(70) ; rewind (71)

        write(6,*) 'go view_3d'
!       call system ("/usr/local/ferret/ferret/bin/ferret")
        call stat_3d(a,nx_s,nx_e,ny_s,ny_e,nz_s,nz_e)
        call paws
        return
        end
