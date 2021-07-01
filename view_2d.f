       subroutine view_2d (a,nx_s,nx_e,ny_s,ny_e)!,field_name)

!  subroutine to view 2-d array with ferret
!  ferret> go view.jnl
!  ferret> exit
!       use constants_module
!       character field_name(20)

        integer nx_s,nx_e,ny_s,ny_e,i,j
        real a(nx_s:nx_e,ny_s:ny_e)

        rewind(71)
        write(71,1801) nx_s,nx_e,ny_s,ny_e
1801    format(4i10)

!        rewind(72)
!        write(72,1802) field_name
!1802    format (a20)         
!	write(6,*) field_name

        rewind (70)
        write(70,1800) ( ( ( a(i,j) ), i=nx_s,nx_e), j=ny_s,ny_e)
1800    format(1x,e15.7) 

        rewind (70); rewind(71)

         write(6,*) 'go view_2d'
!        call system ("rsh rumba; cd /home/rumba/lavelle/Tmp; cat message.txt;logout)
!        call system (/usr/local/ferret/ferret/bin/ferret ")

        call paws
        return
        end
