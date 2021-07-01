!       HotCross Distribution - Version 1.0
!       tag: subroutine ../vents2/lavelle/NewCross/source.f - version 6.4.97

!-------------------------------------------------------------------------------------------------
! --------------------- SUBROUTINE SOURCE ---------------------------------------------------------
! -------------------------------------------------------------------------------------------------

        subroutine source_terms
!
!       import heat, salt, and tracer fluxes, and body forces
!       q_t, q_s, q_c, q_cc1 ...., fbx, fby, fbz 
!

        use constants_module
        use grid_module
        use tracer_module       !q-t,q_s linked here
        use signals_module
        use force_module
        use velocity_module
        use source_module
        include 'array_size.inc'

        real heat_rate, heat_area, heat_flux, rand, source_radius, r, rx, ry, ramp_hrs, wee
        real randnum1, randnum2, randnum3
        integer nx_offset, i, j, k, l, m, ny_offset
        integer ksource, nsource, nsource_halflength, nsource_halfwidth
        integer clock(1:8),nnx,nny,nxs,nxe,nys,nye

        include 'source80.inc'

        return
        end
