!       HotCross Distribution - Version 1.0
!       tag: subroutine ../vents2/lavelle/NewCross/put_full_block.f - version 6.4.97

    	subroutine put_snapshot_header(nxs,nxe,nys,nye,nzs,nze,output_file_shot)

          use grid_module
          use velocity_module
          use viscosity_module
          use tracer_module
          use signals_module

          include 'array_size.inc'
          include '/apps/eb/netCDF/4.1.3-PGI-17.9/include/netcdf.inc'
          integer, save:: dim_id_sx,dim_id_sy,dim_id_sz,dim_id_ux,dim_id_vy,dim_id_wz,dim_id_time
          integer, save:: ind_s,ind_e,start(4),count(4),rcode,cdf_id
          integer, save:: axis_id_sx,axis_id_sy,axis_id_sz
          integer, save:: axis_id_ux,axis_id_vy,axis_id_wz,axis_id_time
          integer, save:: grid_id_s(4),grid_id_u(4),grid_id_v(4),grid_id_w(4),grid_id_p(4),
     &                    grid_id_kxt(3),grid_id_kyt(3),grid_id_kzt(3)

          integer, save:: var_id_u,var_id_v,var_id_w,var_id_p,var_id_ss,var_id_temp,var_id_s,var_id_c,
     &           var_id_cc1,var_id_cc2,var_id_cc3,var_id_cc4,var_id_az,var_id_alphaee, var_id_alphau, var_id_alphav,
     &           var_id_ss_v
          integer, save:: var_id_axu,var_id_axv,var_id_axw,var_id_ayu,var_id_ayv,var_id_ayw,
     &                    var_id_azu,var_id_azv,var_id_azw,var_id_kxt,var_id_kyt,var_id_kzt,var_id_Ri 
          integer l,mct,i,j,k,m
          integer, save::  grid_id_alphaee(3),grid_id_alphau(2),grid_id_alphav(2)
          integer  nxs,nxe,nys,nye,nzs,nze
 
          real dil(0:nx+1,0:ny+1,0:nz+1),sumd,div_p, div_m
          real  ww(0:nx+1,0:ny+1,-1: nz+1 )
          real * 4  missing_value
          !real * 4 work_f ( (nx+3)*(ny+3)*(nz+3) )
          !real * 4 ux_axis(nx+3), wz_axis(nz+3), vy_axis(ny+3)
          !real * 4 sx_axis(nx+2),sy_axis(ny+2),sz_axis(nz+2)                               !axes data
          real  sx_start, sy_start, sz_start, ux_start, vy_start, wz_start, dx_s, dy_s, dz_s

          common/start_info1/output_file, title, history0
          character * 90 input_file, output_file, output_file_shot, title, message, history0, history1, history2

          data message  /'HotCross Version 1.0 - January, 1998, created by J.W.Lavelle, NOAA/PMEL, Seattle  '/               
          data history1 /'Non-hydrostatic model of convection occuring in a cross flow                       '/
          data history2 /'Ref: Buoyancy Driven Plumes in Rotating, Strat. Cross Flows..., J. G. R., 1996     '/

          data missing_value, mct /-1.0e+34, 0/

        entry PUT_snapshot_OPEN(nxs, nxe, nys, nye, nzs, nze, output_file_shot)

          cdf_id = nccre (output_file_shot, ncclob, rcode)     

          call ncaptc(cdf_id, ncglobal, 'title'  ,ncchar,80, title,rcode)
          call ncaptc(cdf_id, ncglobal, 'message',ncchar,80, message,rcode)
          call ncaptc(cdf_id, ncglobal, 'history0', ncchar, 80, history0, rcode)
          call ncaptc(cdf_id, ncglobal, 'history1', ncchar, 80, history1, rcode)
          call ncaptc(cdf_id, ncglobal, 'history2', ncchar, 80, history2, rcode)


          dim_id_time   = ncddef(cdf_id,'time', ncunlim, rcode)
          dim_id_sx     = ncddef(cdf_id,'sx', nx+2, rcode)
          dim_id_sy     = ncddef(cdf_id,'sy', ny+2, rcode)
          dim_id_sz     = ncddef(cdf_id,'sz', nz+2, rcode)
          dim_id_ux     = ncddef(cdf_id,'ux', nx+3, rcode)
          dim_id_vy     = ncddef(cdf_id,'vy', ny+3, rcode)
          dim_id_wz     = ncddef(cdf_id,'wz', nz+3, rcode)

          grid_id_s = (/ dim_id_sx,  dim_id_sy, dim_id_sz, dim_id_time/)
          grid_id_u = (/ dim_id_ux,  dim_id_sy, dim_id_sz, dim_id_time/)
          grid_id_v = (/ dim_id_sx,  dim_id_vy, dim_id_sz, dim_id_time/)
          grid_id_w = (/ dim_id_sx,  dim_id_sy, dim_id_wz, dim_id_time/)
          grid_id_p = (/ dim_id_sx,  dim_id_sy, dim_id_sz, dim_id_time/)

          grid_id_kxt = (/dim_id_ux,  dim_id_sy, dim_id_sz/)
          grid_id_kyt = (/dim_id_sx,  dim_id_vy, dim_id_sz/)
          grid_id_kzt = (/dim_id_sx,  dim_id_sy, dim_id_wz/)          

          grid_id_alphau = (/ dim_id_ux,  dim_id_sy/)
          grid_id_alphav = (/ dim_id_sx,  dim_id_vy/)
          grid_id_alphaee = (/ dim_id_sx,  dim_id_sy,dim_id_sz/)

!         call define_grid (grid_id_az,  dim_id_point,  dim_id_point, dim_id_sz, dim_id_point)
                                                                                                    
          axis_id_time = ncvdef(cdf_id,'time', ncdouble, 1, dim_id_time, rcode)
                 call ncaptc(cdf_id, axis_id_time,'units',ncchar,3,'sec',rcode)
!                call ncaptc(cdf_id,axis_id_time,'time_origin', ncchar,11,'15-JAN-1901',rcode)

          axis_id_sx   = ncvdef(cdf_id,'sx',  ncfloat, 1, dim_id_sx, rcode)
                 call ncaptc(cdf_id, axis_id_sx,'units',ncchar,1,'m',rcode)
!                call ncaptc(cdf_id, axis_id_sx,'point_spacing',ncchar,4,'even',rcode)

          axis_id_sy   = ncvdef(cdf_id,'sy',  ncfloat, 1, dim_id_sy, rcode)
                 call ncaptc(cdf_id, axis_id_sy,'units',ncchar,1,'m',rcode)
!                call ncaptc(cdf_id, axis_id_sy,'point_spacing',ncchar,4,'even',rcode)

          axis_id_sz   = ncvdef(cdf_id,'sz',  ncfloat, 1, dim_id_sz, rcode)
                call ncaptc(cdf_id, axis_id_sz, 'units', ncchar, 1, 'm', rcode)
                call ncaptc(cdf_id, axis_id_sz, 'positive', ncchar, 4, 'down',rcode)
!               call ncaptc(cdf_id, axis_id_sz,'point_spacing',ncchar,4,'even',rcode)

          axis_id_ux   = ncvdef(cdf_id,'ux',  ncfloat,1,dim_id_ux,rcode)
                 call ncaptc(cdf_id, axis_id_ux,'units',ncchar,1,'m',rcode)
!                call ncaptc(cdf_id, axis_id_ux,'point_spacing',ncchar,4,'even',rcode)

          axis_id_vy   = ncvdef(cdf_id,'vy',  ncfloat,1,dim_id_vy,rcode)
                 call ncaptc(cdf_id, axis_id_vy,'units',ncchar,1,'m',rcode)
!                call ncaptc(cdf_id, axis_id_vy,'point_spacing',ncchar,4,'even',rcode)

          axis_id_wz   = ncvdef(cdf_id,'wz',  ncfloat,1,dim_id_wz,rcode)
                 call ncaptc(cdf_id, axis_id_wz,'units',ncchar,1,'m',rcode)
                 call ncaptc(cdf_id, axis_id_wz, 'positive',ncchar,4,'down',rcode)
!                 call ncaptc(cdf_id, axis_id_wz,'point_spacing',ncchar,4,'even',rcode)


          if (toggle(1) .eq. 1) then
          var_id_u = ncvdef(cdf_id,'u',ncfloat,4,grid_id_u,rcode)
               call ncaptc(cdf_id, var_id_u,'units',ncchar,4,'ms-1',rcode)
               call ncaptc(cdf_id, var_id_u,'long_name',ncchar,10,'u-velocity',rcode)
               call ncapt(cdf_id,  var_id_u,'missing_value',ncfloat,1,missing_value,rcode)
          endif

          if (toggle(2) .eq. 1) then
          var_id_v = ncvdef(cdf_id,'v',ncfloat,4,grid_id_v,rcode)
               call ncaptc(cdf_id, var_id_v,'units',ncchar,4,'ms-1',rcode)
               call ncaptc(cdf_id, var_id_v,'long_name',ncchar,10,'v-velocity',rcode)
               call ncapt(cdf_id,  var_id_v,'missing_value',ncfloat,1,missing_value,rcode)
          endif

          if (toggle(3) .eq. 1) then
          var_id_w = ncvdef(cdf_id,'w',ncfloat,4,grid_id_w,rcode)
               call ncaptc(cdf_id, var_id_w,'units',ncchar,4,'ms-1',rcode)
               call ncaptc(cdf_id, var_id_w,'long_name',ncchar,17,'vertical velocity',rcode)
               call ncapt(cdf_id,  var_id_w,'missing_value',ncfloat,1,missing_value,rcode)
          endif

          if (toggle(4) .eq. 1) then
          var_id_p = ncvdef(cdf_id,'p',ncfloat,4,grid_id_p,rcode)
               call ncaptc(cdf_id, var_id_p,'units',ncchar,10,'newtons/m2',rcode)
               call ncaptc(cdf_id, var_id_p,'long_name',ncchar,16,'Pressure Anomaly',rcode)
               call ncapt(cdf_id,  var_id_p,'missing_value',ncfloat,1,missing_value,rcode)
          endif

          if (toggle(5) .eq. 1) then
          var_id_ss = ncvdef(cdf_id,'ss',ncfloat,4,grid_id_p,rcode)
               call ncaptc(cdf_id, var_id_ss,'units',ncchar,5,'m2s-1',rcode)
               call ncaptc(cdf_id, var_id_ss,'long_name',ncchar,17,'Horz Mixing Coef.',rcode)
               call ncapt(cdf_id,  var_id_ss,'missing_value',ncfloat,1,missing_value,rcode)
          var_id_ss_v = ncvdef(cdf_id,'ss_v',ncfloat,4,grid_id_p,rcode)
               call ncaptc(cdf_id, var_id_ss_v,'units',ncchar,5,'m2s-1',rcode)
               call ncaptc(cdf_id, var_id_ss_v,'long_name',ncchar,18,'Vert. Shear Mixing',rcode)
               call ncapt(cdf_id,  var_id_ss_v,'missing_value',ncfloat,1,missing_value,rcode)
          endif

          var_id_ri = ncvdef(cdf_id,'ri',ncfloat,4,grid_id_p,rcode)
               call ncaptc(cdf_id, var_id_ri,'units',ncchar,1,' ',rcode)
               call ncaptc(cdf_id, var_id_ri,'long_name',ncchar,17,'Richardson Number',rcode)
               call ncapt(cdf_id,  var_id_ri,'missing_value',ncfloat,1,missing_value,rcode)

          if (toggle(6) .eq. 1) then
          var_id_temp = ncvdef(cdf_id,'temp',ncfloat,4,grid_id_s,rcode)
               call ncaptc(cdf_id, var_id_temp,'units',ncchar,2,'oC',rcode)
               call ncaptc(cdf_id, var_id_temp,'long_name',ncchar,21,'Potential Temperature',rcode)
               call ncapt(cdf_id,  var_id_temp,'missing_value',ncfloat,1,missing_value,rcode)
          endif

          if (toggle(7) .eq. 1) then
          var_id_s = ncvdef(cdf_id,'s',ncfloat,4,grid_id_s,rcode)
               call ncaptc(cdf_id, var_id_s,'units',ncchar,4,'o/oo',rcode)
               call ncaptc(cdf_id, var_id_s,'long_name',ncchar,8,'Salinity',rcode)
               call ncapt(cdf_id,  var_id_s,'missing_value',ncfloat,1,missing_value,rcode)
          endif

          if (toggle(8) .eq. 1) then
          var_id_c = ncvdef(cdf_id,'c',ncfloat,3,grid_id_s,rcode)
               call ncaptc(cdf_id, var_id_c,'units',ncchar,13,'mass_units/m3',rcode)
               call ncaptc(cdf_id, var_id_c,'long_name',ncchar,8,'Tracer-0',rcode)
               call ncapt(cdf_id,  var_id_c,'missing_value',ncfloat,1,missing_value,rcode)
          endif

          var_id_kxt = ncvdef(cdf_id,'kxt',ncfloat,3,grid_id_kxt,rcode)
               call ncaptc(cdf_id, var_id_kxt,'units',ncchar,4,'m2/s',rcode)
               call ncaptc(cdf_id, var_id_kxt,'long_name',ncchar,3,'KXT',rcode)
               call ncapt(cdf_id,  var_id_kxt,'missing_value',ncfloat,1,missing_value,rcode)

          var_id_kyt = ncvdef(cdf_id,'kyt',ncfloat,3,grid_id_kyt,rcode)
               call ncaptc(cdf_id, var_id_kyt,'units',ncchar,4,'m2/s',rcode)
               call ncaptc(cdf_id, var_id_kyt,'long_name',ncchar,3,'KYT',rcode)
               call ncapt(cdf_id,  var_id_kyt,'missing_value',ncfloat,1,missing_value,rcode)

          var_id_kzt = ncvdef(cdf_id,'kzt',ncfloat,3,grid_id_kzt,rcode)
               call ncaptc(cdf_id, var_id_kzt,'units',ncchar,4,'m2/s',rcode)
               call ncaptc(cdf_id, var_id_kzt,'long_name',ncchar,3,'KZT',rcode)
               call ncapt(cdf_id,  var_id_kzt,'missing_value',ncfloat,1,missing_value,rcode)



          if (toggle(9) .eq. 1) then                                                                                              
           var_id_cc1 = ncvdef(cdf_id,'cc1',ncfloat,4,grid_id_s,rcode)
               call ncaptc(cdf_id, var_id_cc1,'units',ncchar,13,'mass_units/m3',rcode)
               call ncaptc(cdf_id, var_id_cc1,'long_name',ncchar,8,'Tracer-1',rcode)
               call ncapt(cdf_id,  var_id_cc1,'missing_value',ncfloat,1,missing_value,rcode)
          endif

          if (toggle(10) .eq. 1) then
           var_id_cc2 = ncvdef(cdf_id,'cc2',ncfloat,4,grid_id_s,rcode)
               call ncaptc(cdf_id, var_id_cc2,'units',ncchar,13,'mass_units/m3',rcode)
               call ncaptc(cdf_id, var_id_cc2,'long_name',ncchar,8,'Tracer-2',rcode)
               call ncapt(cdf_id,  var_id_cc2,'missing_value',ncfloat,1,missing_value,rcode)
          endif

          if (toggle(11) .eq. 1) then                                                          
           var_id_cc3 = ncvdef(cdf_id,'cc3',ncfloat,4,grid_id_s,rcode)
               call ncaptc(cdf_id, var_id_cc3,'units',ncchar,13,'mass_units/m3',rcode)
               call ncaptc(cdf_id, var_id_cc3,'long_name',ncchar,8,'Tracer-3',rcode)
               call ncapt(cdf_id,  var_id_cc3,'missing_value',ncfloat,1,missing_value,rcode)
          endif

          if (toggle(12) .eq. 1) then
          var_id_cc4 = ncvdef(cdf_id,'cc4',ncfloat,4,grid_id_s,rcode)
               call ncaptc(cdf_id, var_id_cc4,'units',ncchar,13,'mass_units/m3',rcode)
               call ncaptc(cdf_id, var_id_cc4,'long_name',ncchar,8,'Tracer-4',rcode)
               call ncapt(cdf_id,  var_id_cc4,'missing_value',ncfloat,1,missing_value,rcode)
          endif                                                                                 

          if (toggle(13) .eq. 1) then
          var_id_az = ncvdef(cdf_id,'Az_bkg',ncfloat,1,dim_id_sz,rcode)
               call ncaptc(cdf_id, var_id_az,'units',ncchar,4,'m2/s',rcode)
               call ncaptc(cdf_id, var_id_az,'long_name',ncchar,13,'Background Az',rcode)
               call ncapt(cdf_id,  var_id_az,'missing_value',ncfloat,1,missing_value,rcode)
          endif

         var_id_alphaee = ncvdef(cdf_id,'alphaee',ncfloat,3,grid_id_alphaee,rcode)
               call ncaptc(cdf_id, var_id_alphaee,'units',ncchar,6,'s^(-1)',rcode)
               call ncaptc(cdf_id, var_id_alphaee,'long_name',ncchar,10,'eta_sponge',rcode)
               call ncapt(cdf_id,  var_id_alphaee,'missing_value',ncfloat,1,missing_value,rcode)

          var_id_alphau = ncvdef(cdf_id,'alphau',ncfloat,2,grid_id_alphau,rcode)
               call ncaptc(cdf_id, var_id_alphau,'units',ncchar,6,'s^(-1)',rcode)
               call ncaptc(cdf_id, var_id_alphau,'long_name',ncchar,8,'u_sponge',rcode)
               call ncapt(cdf_id,  var_id_alphau,'missing_value',ncfloat,1,missing_value,rcode)

          var_id_alphav = ncvdef(cdf_id,'alphav',ncfloat,2,grid_id_alphav,rcode)
               call ncaptc(cdf_id, var_id_alphav,'units',ncchar,6,'s^(-1)',rcode)
               call ncaptc(cdf_id, var_id_alphav,'long_name',ncchar,8,'v_sponge',rcode)
               call ncapt(cdf_id,  var_id_alphav,'missing_value',ncfloat,1,missing_value,rcode)
                                                                                              
          call ncendf(cdf_id,rcode)

!         ------------------------------------------------------------------------------------------
           sx_start = xstart -.5d0 * dx
           sy_start = ystart -.5d0 * dy
           sz_start = zstart -.5d0 * dz
           ux_start = xstart - dx
           vy_start = ystart - dy
           wz_start = zstart - dz
           dx_s = dx
           dy_s = dy
           dz_s = dz


           call load_axis ( cdf_id, axis_id_sx,  nx +2,  sx_start,   dx_s )
           call load_axis ( cdf_id, axis_id_sy,  ny +2,  sy_start,   dy_s )
           call load_axis ( cdf_id, axis_id_sz,  nz +2,  sz_start,   dz_s )
           call load_axis ( cdf_id, axis_id_ux,  nx +3,  ux_start,   dx_s )
           call load_axis ( cdf_id, axis_id_vy,  ny +3,  vy_start,   dy_s )
           call load_axis ( cdf_id, axis_id_wz,  nz +3,  wz_start,   dz_s )

          return
          !end subroutine PUT_snapshot_open
    
      entry PUT_snapshot(nxs, nxe, nys, nye, nzs, nze, output_file_shot)

          mct = mct + 1
          mct = 1

          call ncvpt1(cdf_id,axis_id_time,mct,time,rcode)

          if (toggle(1) .eq. 1) then
!         call start_n_count ( start, 1, 1, 1, mct, count, nx+3, ny+2, nz+2, 1)
          call ncvpt (cdf_id, var_id_u, (/1,1,1,mct/) , (/nx+3,ny+2,nz+2,1/), 
     &                       real ( u(-1: nx+1,0: ny+1,0: nz+1 ), kind = 4 ), rcode)

          endif

          if (toggle(2) .eq. 1) then
!          call start_n_count ( start, 1, 1, 1, mct, count, nx+2, ny+3, nz+2, 1)
          call ncvpt (cdf_id, var_id_v,(/1,1,1,mct/) , (/nx+2,ny+3,nz+2,1/), real ( v(0: nx+1,-1: ny+1,0: nz+1 ), kind = 4 ), rcode)
          endif

          if (toggle(3) .eq. 1) then
!          call start_n_count ( start, 1, 1, 1, mct, count, nx+2, ny+2, nz+3, 1)
          ww = - w
          call fill_value(ww,0,nx+1,0,ny+1,-1,nz+1)
          call ncvpt (cdf_id, var_id_w, (/1,1,1,mct/),(/nx+2,ny+2,nz+3,1/), real ( ww(0:nx+1,0:ny+1,-1: nz+1 ), kind = 4 ), rcode)
          endif

          if (toggle(4) .eq. 1) then
!          call start_n_count ( start, 1, 1, 1, mct, count, nx+2, ny+2, nz+2, 1)
          call ncvpt (cdf_id, var_id_p, (/1,1,1,mct/) , (/nx+2,ny+2,nz+2,1/) , real( p(0:nx+1,0:ny+1,0:nz+1 ), kind = 4 ), rcode)

          endif

          if (toggle(5) .eq. 1) then
!          call start_n_count ( start, 1, 1, 1, mct, count, nx+2, ny+2, nz+2, 1)
          call ncvpt (cdf_id, var_id_ss,  (/1,1,1,mct/) , (/nx+2,ny+2,nz+2,1/) , 
     &                 real ( ss(0:nx+1,0:ny+1,0: nz+1 ), kind = 4 ), rcode)
          call ncvpt (cdf_id, var_id_ss_v, (/1,1,1,mct/) , (/nx+2,ny+2,nz+2,1/) , 
     &                 real ( ss_v(0:nx+1,0:ny+1,0: nz+1 ), kind = 4 ), rcode)
          call ncvpt (cdf_id, var_id_ri,(/1,1,1,mct/) , (/nx+2,ny+2,nz+2,1/) , 
     &                 real ( ri(0:nx+1,0:ny+1,0: nz+1 ), kind = 4 ), rcode)
          endif



          if (toggle(6) .eq. 1) then
!          call start_n_count ( start, 1, 1, 1, mct, count, nx+2, ny+2, nz+2, 1)
          call ncvpt (cdf_id, var_id_temp, (/1,1,1,mct/) , (/nx+2,ny+2,nz+2,1/) , 
     &                       real ( t(0:nx+1,0:ny+1,0: nz+1 ), kind = 4), rcode)
          endif

          if (toggle(7) .eq. 1) then
!          call start_n_count ( start, 1, 1, 1, mct, count, nx+2, ny+2, nz+2, 1)
          call ncvpt (cdf_id, var_id_s,  (/1,1,1,mct/) , (/nx+2,ny+2,nz+2,1/), 
     &                  real ( s(0:nx+1,0:ny+1,0: nz+1 ), kind = 4), rcode)
          endif

          call ncvpt (cdf_id, var_id_kxt, (/1,1,1/),(/nx+3,ny+2,nz+2/),real ( kxt(-1:nx+1,0:ny+1,0: nz+1 ), kind = 4), rcode)
          call ncvpt (cdf_id, var_id_kyt, (/1,1,1/),(/nx+2,ny+3,nz+2/),real ( kyt(0:nx+1,-1:ny+1,0: nz+1 ), kind = 4), rcode) 
          call ncvpt (cdf_id, var_id_kzt, (/1,1,1/),(/nx+2,ny+2,nz+3/),real ( kzt(0:nx+1,0:ny+1,-1: nz+1 ), kind = 4), rcode)

          if (toggle(8) .eq. 1) then
!    	  call start_n_count ( start, 1, 1, 1, mct, count, nx+2, ny+2, nz+2, 1)
          call ncvpt (cdf_id, var_id_c,  (/1,1,1,mct/) , (/nx+2,ny+2,nz+2,1/) , 
     &                     real ( c( 0: nx+1,0: ny+1,0: nz+1 ), kind = 4), rcode)
          endif

          if (toggle(9) .eq. 1) then
!          call start_n_count ( start, 1, 1, 1, mct, count, nx+2, ny+2, nz+2, 1)
          call ncvpt (cdf_id, var_id_cc1,  (/1,1,1,mct/) , (/nx+2,ny+2,nz+2,1/) , 
     &                    real ( cc1(0:nx+1,0:ny+1,0: nz+1 ), kind =4 ), rcode)
          endif

          if (toggle(10) .eq. 1) then                                                                                      
!          call start_n_count ( start, 1, 1, 1, mct, count, nx+2, ny+2, nz+2, 1)
          call ncvpt (cdf_id, var_id_cc1,  (/1,1,1,mct/) , (/nx+2,ny+2,nz+2,1/) , 
     &                           real ( cc2(0:nx+1,0:ny+1,0: nz+1 ), kind =4 ), rcode)
          endif

          if (toggle(11) .eq. 1) then
!          call start_n_count ( start, 1, 1, 1, mct, count, nx+2, ny+2, nz+2, 1)
          call ncvpt (cdf_id, var_id_cc1, (/1,1,1,mct/) , (/nx+2,ny+2,nz+2,1/) , 
     &                              real ( cc3(0:nx+1,0:ny+1,0: nz+1 ), kind =4 ), rcode)
          endif

          if (toggle(12) .eq. 1) then
!          call start_n_count ( start, 1, 1, 1, mct, count, nx+2, ny+2, nz+2, 1)
          call ncvpt (cdf_id, var_id_cc1,  (/1,1,1,mct/) , (/nx+2,ny+2,nz+2,1/) , 
     &                        real ( cc4(0:nx+1,0:ny+1,0: nz+1 ), kind =4 ), rcode)
          endif

          if (toggle(13) .eq. 1 .and. mct .eq. 1 ) then
!          call start_n_count ( start, 1, 1, 1, 1, count, 1, 1, nz+2, 1)
          call ncvpt (cdf_id, var_id_az, 1, nz+2,  real ( azzmin(0: nz+1 ), kind =4 ), rcode)
          endif

          call ncvpt (cdf_id, var_id_alphaee, (/1,1,1/), (/nxe - nxs + 1, nye - nys + 1, nze - nzs + 1/),
     &       real ( alphaee( nxs: nxe, nys: nye, nzs: nze), kind =4), rcode)

          call ncvpt (cdf_id, var_id_alphau, (/1,1/), (/nxe - nxs + 2, nye - nys + 1/),
     &       real ( alphau( nxs-1: nxe, nys: nye), kind =4), rcode)

          call ncvpt (cdf_id, var_id_alphav, (/1,1/), (/nxe - nxs + 1, nye - nys + 2/),
     &       real ( alphav( nxs: nxe, nys-1: nye), kind =4), rcode)

          return
          !end subroutine put_snapshot
 
      entry PUT_snapshot_CLOSE(nxs, nxe, nys, nye, nzs, nze,output_file_shot)
       
          call  ncclos(cdf_id,rcode)
          return

      contains 

         subroutine load_axis (cdf_id, axis_id, nx, xs, xinc)

         integer cdf_id, axis_id, nx, i, rcode
         real  xs, xinc,  axis_data (1:nx)

             do 10 i = 1,nx
  10              axis_data(i) = xinc * ( i - 1 )  + xs

             call ncvpt( cdf_id, axis_id, 1, nx, real(axis_data,kind=4), rcode )

             return
             end subroutine load_axis

      end subroutine put_snapshot_header

