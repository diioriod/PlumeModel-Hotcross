        ! a copy of  put_block2.f except block is changed to block2
!       HotCross Distribution - Version 1.0
!       tag: subroutine ../vents2/lavelle/NewCross/put_full_block2.f - version 6.4.97

    	subroutine push_full_block2_header(nxs,nxe,nys,nye,nzs,nze,output_file_shot,cdf_id)


          use grid_module
          use velocity_module
          use viscosity_module
          use tracer_module
          use signals_module

          include 'array_size.inc'
          include '/apps/eb/netCDF/4.1.3-PGI-17.9/include/netcdf.inc'

          integer cdf_id
          integer, save::  dim_id_sx,dim_id_sy,dim_id_sz,dim_id_ux,dim_id_vy,dim_id_wz,dim_id_time
          integer, save:: axis_id_sx,axis_id_sy,axis_id_sz
          integer, save:: axis_id_ux,axis_id_vy,axis_id_wz,axis_id_time
          integer, save::  grid_id_s(4),grid_id_u(4),grid_id_v(4),grid_id_w(4),grid_id_p(4)
          integer  ind_s,ind_e,start(4),count(4),rcode

          integer, save:: var_id_u,var_id_v,var_id_w,var_id_p,var_id_ss,var_id_temp,var_id_s,var_id_c,
     &           var_id_cc1,var_id_cc2,var_id_cc3,var_id_cc4,var_id_az, var_id_alphae, var_id_alphau, var_id_alphav
     &           , var_id_epsh,var_id_epsv,var_id_chi_v, var_id_chi_h, var_id_Ri, var_id_rho
          integer l,mct,i,j,k,m
          integer, save::  grid_id_alphae(2),grid_id_alphau(2),grid_id_alphav(2)
          integer nxs,nxe,nys,nye,nzs,nze
 
          real   ww(0:nx+1,0:ny+1,-1: nz+1 ), axis_data(nx*ny*nz)
          real * 4  missing_value
          !real * 4 work_f ((nxe-nxs+2)*(nye-nys+2)*(nze-nzs+2) )
          !real * 4 ux_axis(nxe-nxs+2), wz_axis(nze-nzs+2), vy_axis(nye-nys+2)
          !real * 4 sx_axis(nxe-nxs+1),sy_axis(nye-nys+1),sz_axis(nze-nzs+1)                               !axes data
          real sx_start, sy_start, sz_start, ux_start, vy_start, wz_start
          

          common/start_info1/output_file, title, history0
          character * 90 input_file, output_file, output_file_shot, title, message, history0, history1, history2
          data message  /'Code HotCross created by J.W.Lavelle, NOAA/PMEL, Seattle  '/               
          data history1 /'Non-hydrostatic 3-d model of convection occuring in a time_variable cross flow        '/
          data history2 /'Ref: Buoyancy Driven Plumes in Rotating, Strat. Cross Flows...,Lavelle, J. G. R., 1996'/

          data missing_value, mct /-1.0e+34, 0/

      ENTRY PUT_full_block2_OPEN(nxs, nxe, nys, nye, nzs, nze, output_file_shot,cdf_id)

          cdf_id = nccre (output_file_shot,ncclob,rcode)     

          call ncaptc(cdf_id, ncglobal, 'title'  ,ncchar,80, title,rcode)
          call ncaptc(cdf_id, ncglobal, 'message',ncchar,80, message,rcode)
          call ncaptc(cdf_id, ncglobal, 'history0', ncchar, 80, history0, rcode)
          call ncaptc(cdf_id, ncglobal, 'history1', ncchar, 80, history1, rcode)
          call ncaptc(cdf_id, ncglobal, 'history2', ncchar, 80, history2, rcode)

          dim_id_time   = ncddef(cdf_id,'time', ncunlim, rcode)

          dim_id_sx     = ncddef(cdf_id,'sx', nxe-nxs+1, rcode)
          dim_id_sy     = ncddef(cdf_id,'sy', nye-nys+1, rcode)
          dim_id_sz     = ncddef(cdf_id,'sz', nze-nzs+1, rcode)

          dim_id_ux     = ncddef(cdf_id,'ux', nxe-nxs+2, rcode)
          dim_id_vy     = ncddef(cdf_id,'vy', nye-nys+2, rcode)
          dim_id_wz     = ncddef(cdf_id,'wz', nze-nzs+2, rcode)


           call define_grid (grid_id_s,  dim_id_sx,  dim_id_sy, dim_id_sz, dim_id_time)
          call define_grid (grid_id_u,  dim_id_ux,  dim_id_sy, dim_id_sz, dim_id_time)
          call define_grid (grid_id_v,  dim_id_sx,  dim_id_vy, dim_id_sz, dim_id_time)

!           grid_id_u = (/dim_id_time/)
!           grid_id_v = (/dim_id_time/)

           call define_grid (grid_id_w,  dim_id_sx,  dim_id_sy, dim_id_wz, dim_id_time)
!          call define_grid (grid_id_p,  dim_id_sx,  dim_id_sy, dim_id_sz, dim_id_time)

!          grid_id_alphau = (/ dim_id_ux,  dim_id_sy/)
!          grid_id_alphav = (/ dim_id_sx,  dim_id_vy/)
!          grid_id_alphae = (/ dim_id_sx,  dim_id_sy/)
                                                                                                    
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
                call ncaptc(cdf_id, axis_id_sz,'point_spacing',ncchar,4,'even',rcode)

          axis_id_ux   = ncvdef(cdf_id,'ux',  ncfloat,1,dim_id_ux,rcode)
                 call ncaptc(cdf_id, axis_id_ux,'units',ncchar,1,'m',rcode)
                call ncaptc(cdf_id, axis_id_ux,'point_spacing',ncchar,4,'even',rcode)

          axis_id_vy   = ncvdef(cdf_id,'vy',  ncfloat,1,dim_id_vy,rcode)
                 call ncaptc(cdf_id, axis_id_vy,'units',ncchar,1,'m',rcode)
                call ncaptc(cdf_id, axis_id_vy,'point_spacing',ncchar,4,'even',rcode)

          axis_id_wz   = ncvdef(cdf_id,'wz',  ncfloat,1,dim_id_wz,rcode)
                 call ncaptc(cdf_id, axis_id_wz,'units',ncchar,1,'m',rcode)
                 call ncaptc(cdf_id, axis_id_wz, 'positive',ncchar,4,'down',rcode)
!                call ncaptc(cdf_id, axis_id_wz,'point_spacing',ncchar,4,'even',rcode)


!         if (toggle(1) .eq. -1) then
          var_id_u = ncvdef(cdf_id,'u',ncfloat,4,grid_id_u,rcode)
               call ncaptc(cdf_id, var_id_u,'units',ncchar,4,'ms-1',rcode)
               call ncaptc(cdf_id, var_id_u,'long_name',ncchar,19,'bkground u-velocity',rcode)
               call ncapt(cdf_id,  var_id_u,'missing_value',ncfloat,1,missing_value,rcode)
!          endif

!          if (toggle(2) .eq. -1) then
          var_id_v = ncvdef(cdf_id,'v',ncfloat,4,grid_id_v,rcode)
               call ncaptc(cdf_id, var_id_v,'units',ncchar,4,'ms-1',rcode)
               call ncaptc(cdf_id, var_id_v,'long_name',ncchar,19,'bkground v-velocity',rcode)
               call ncapt(cdf_id,  var_id_v,'missing_value',ncfloat,1,missing_value,rcode)
!          endif

          if (toggle(3) .eq. 1) then
          var_id_w = ncvdef(cdf_id,'w',ncfloat,4,grid_id_w,rcode)
               call ncaptc(cdf_id, var_id_w,'units',ncchar,4,'ms-1',rcode)
               call ncaptc(cdf_id, var_id_w,'long_name',ncchar,17,'vertical velocity',rcode)
               call ncapt(cdf_id,  var_id_w,'missing_value',ncfloat,1,missing_value,rcode)
          endif

          if (toggle(4) .eq. -1) then
          var_id_p = ncvdef(cdf_id,'p',ncfloat,4,grid_id_p,rcode)
               call ncaptc(cdf_id, var_id_p,'units',ncchar,10,'newtons/m2',rcode)
               call ncaptc(cdf_id, var_id_p,'long_name',ncchar,16,'Pressure Anomaly',rcode)
               call ncapt(cdf_id,  var_id_p,'missing_value',ncfloat,1,missing_value,rcode)
          endif

          if (toggle(5) .eq. -1) then
          var_id_ss = ncvdef(cdf_id,'ss',ncfloat,4,grid_id_p,rcode)
               call ncaptc(cdf_id, var_id_ss,'units',ncchar,5,'m2s-1',rcode)
               call ncaptc(cdf_id, var_id_ss,'long_name',ncchar,18,'Mixing coefficient',rcode)
               call ncapt(cdf_id,  var_id_ss,'missing_value',ncfloat,1,missing_value,rcode)
          endif

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

!          if (toggle(8) .eq. -1) then
          var_id_c = ncvdef(cdf_id,'c',ncfloat,4,grid_id_s,rcode)
               call ncaptc(cdf_id, var_id_c,'units',ncchar,13,'mass_units/m3',rcode)
               call ncaptc(cdf_id, var_id_c,'long_name',ncchar,8,'Tracer-0',rcode)
               call ncapt(cdf_id,  var_id_c,'missing_value',ncfloat,1,missing_value,rcode)
!          endif

          if (toggle(9) .eq. -1) then                                                                                              
           var_id_cc1 = ncvdef(cdf_id,'cc1',ncfloat,4,grid_id_s,rcode)
               call ncaptc(cdf_id, var_id_cc1,'units',ncchar,13,'mass_units/m3',rcode)
               call ncaptc(cdf_id, var_id_cc1,'long_name',ncchar,8,'Tracer-1',rcode)
               call ncapt(cdf_id,  var_id_cc1,'missing_value',ncfloat,1,missing_value,rcode)
          endif

          if (toggle(10) .eq. -1) then
           var_id_cc2 = ncvdef(cdf_id,'cc2',ncfloat,4,grid_id_s,rcode)
               call ncaptc(cdf_id, var_id_cc2,'units',ncchar,13,'mass_units/m3',rcode)
               call ncaptc(cdf_id, var_id_cc2,'long_name',ncchar,8,'Tracer-2',rcode)
               call ncapt(cdf_id,  var_id_cc2,'missing_value',ncfloat,1,missing_value,rcode)
          endif

          if (toggle(11) .eq. -1) then                                                          
           var_id_cc3 = ncvdef(cdf_id,'cc3',ncfloat,4,grid_id_s,rcode)
               call ncaptc(cdf_id, var_id_cc3,'units',ncchar,13,'mass_units/m3',rcode)
               call ncaptc(cdf_id, var_id_cc3,'long_name',ncchar,8,'Tracer-3',rcode)
               call ncapt(cdf_id,  var_id_cc3,'missing_value',ncfloat,1,missing_value,rcode)
          endif

          if (toggle(12) .eq. -1) then
          var_id_cc4 = ncvdef(cdf_id,'cc4',ncfloat,4,grid_id_s,rcode)
               call ncaptc(cdf_id, var_id_cc4,'units',ncchar,13,'mass_units/m3',rcode)
               call ncaptc(cdf_id, var_id_cc4,'long_name',ncchar,8,'Tracer-4',rcode)
               call ncapt(cdf_id,  var_id_cc4,'missing_value',ncfloat,1,missing_value,rcode)
          endif                                                                                 

          if (toggle(13) .eq. -1) then
          var_id_az = ncvdef(cdf_id,'Az_bkg',ncfloat,1,dim_id_sz,rcode)
               call ncaptc(cdf_id, var_id_az,'units',ncchar,4,'m2/s',rcode)
               call ncaptc(cdf_id, var_id_az,'long_name',ncchar,13,'Background Az',rcode)
               call ncapt(cdf_id,  var_id_az,'missing_value',ncfloat,1,missing_value,rcode)
          endif

!         var_id_alphae = ncvdef(cdf_id,'alphae',ncfloat,2,grid_id_alphae,rcode)
!               call ncaptc(cdf_id, var_id_alphae,'units',ncchar,6,'s^(-1)',rcode)
!               call ncaptc(cdf_id, var_id_alphae,'long_name',ncchar,10,'eta_sponge',rcode)
!               call ncapt(cdf_id,  var_id_alphae,'missing_value',ncfloat,1,missing_value,rcode)

!          var_id_alphau = ncvdef(cdf_id,'alphau',ncfloat,2,grid_id_alphau,rcode)
!               call ncaptc(cdf_id, var_id_alphau,'units',ncchar,6,'s^(-1)',rcode)
!               call ncaptc(cdf_id, var_id_alphau,'long_name',ncchar,8,'u_sponge',rcode)
!               call ncapt(cdf_id,  var_id_alphau,'missing_value',ncfloat,1,missing_value,rcode)

!          var_id_alphav = ncvdef(cdf_id,'alphav',ncfloat,2,grid_id_alphav,rcode)
!               call ncaptc(cdf_id, var_id_alphav,'units',ncchar,6,'s^(-1)',rcode)
!               call ncaptc(cdf_id, var_id_alphav,'long_name',ncchar,8,'v_sponge',rcode)
!               call ncapt(cdf_id,  var_id_alphav,'missing_value',ncfloat,1,missing_value,rcode)

         var_id_epsh = ncvdef(cdf_id,'epsh',ncfloat,4,grid_id_s,rcode)
               call ncaptc(cdf_id, var_id_epsh,'units',ncchar,5,'m2/s3',rcode)
               call ncaptc(cdf_id, var_id_epsh,'long_name',ncchar,12,'Horiz. Diss.',rcode)
               call ncapt (cdf_id, var_id_epsh,'missing_value',ncfloat,1,missing_value,rcode)

         var_id_epsv = ncvdef(cdf_id,'epsv',ncfloat,4,grid_id_s,rcode)
               call ncaptc(cdf_id, var_id_epsv,'units',ncchar,5,'m2/s3',rcode)
               call ncaptc(cdf_id, var_id_epsv,'long_name',ncchar,11,'Vert. Diss.',rcode)
               call ncapt (cdf_id, var_id_epsv,'missing_value',ncfloat,1,missing_value,rcode)

         var_id_chi_h = ncvdef(cdf_id,'chih',ncfloat,4,grid_id_s,rcode)
               call ncaptc(cdf_id, var_id_chi_h,'units',ncchar,5,'m2/s3',rcode)
               call ncaptc(cdf_id, var_id_chi_h,'long_name',ncchar,18,'horz temp var diss',rcode)
               call ncapt (cdf_id, var_id_chi_h,'missing_value',ncfloat,1,missing_value,rcode)

         var_id_chi_v = ncvdef(cdf_id,'chiv',ncfloat,4,grid_id_s,rcode)
               call ncaptc(cdf_id, var_id_chi_v,'units',ncchar,5,'m2/s3',rcode)
               call ncaptc(cdf_id, var_id_chi_v,'long_name',ncchar,18,'vert temp var diss',rcode)
               call ncapt (cdf_id, var_id_chi_v,'missing_value',ncfloat,1,missing_value,rcode)

         var_id_Ri = ncvdef(cdf_id,'Ri',ncfloat,4,grid_id_s,rcode)
               call ncaptc(cdf_id, var_id_Ri,'units',ncchar,1,' ',rcode)
               call ncaptc(cdf_id, var_id_Ri,'long_name',ncchar,17,'Richardson number',rcode)
               call ncapt (cdf_id, var_id_Ri,'missing_value',ncfloat,1,missing_value,rcode)

         var_id_rho = ncvdef(cdf_id,'rrho',ncfloat,4,grid_id_s,rcode)
               call ncaptc(cdf_id, var_id_rho,'units',ncchar,5,'kg/m3',rcode)
               call ncaptc(cdf_id, var_id_rho,'long_name',ncchar,7,'Density',rcode)
               call ncapt (cdf_id, var_id_rho,'missing_value',ncfloat,1,missing_value,rcode)


          call ncendf(cdf_id,rcode)

!         ------------------------------------------------------------------------------------------
           sx_start = xstart -.5d0 * dx + real(nxs) * dx
           sy_start = ystart -.5d0 * dy + real(nys) * dy 
           sz_start = zstart -.5d0 * dz + real(nzs) * dz
           
           ux_start = xstart - dx + real(nxs) * dx
           vy_start = ystart - dy + real(nys) * dy
           wz_start = zstart - dz + real(nzs) * dz

           call load_axis ( cdf_id, axis_id_sx,  nxe-nxs +1,  sx_start,   dx )
           call load_axis ( cdf_id, axis_id_sy,  nye-nys +1,  sy_start,   dy )
           call load_axis ( cdf_id, axis_id_sz,  nze-nzs +1,  sz_start,   dz )
           call load_axis ( cdf_id, axis_id_ux,  nxe-nxs +2,  ux_start,   dx )
           call load_axis ( cdf_id, axis_id_vy,  nye-nys +2,  vy_start,   dy )
           call load_axis ( cdf_id, axis_id_wz,  nze-nzs +2,  wz_start,   dz )

!          call ncvpt (cdf_id, var_id_alphae, (/1,1/), (/nxe - nxs + 1, nye - nys + 1/),
!     &       real ( alphae( nxs: nxe, nys: nye), kind =4), rcode)

!          call ncvpt (cdf_id, var_id_alphau, (/1,1/), (/nxe - nxs + 2, nye - nys + 1/),
!     &       real ( alphau( nxs-1: nxe, nys: nye), kind =4), rcode)

!          call ncvpt (cdf_id, var_id_alphav, (/1,1/), (/nxe - nxs + 1, nye - nys + 2/),
!     &       real ( alphav( nxs: nxe, nys-1: nye), kind =4), rcode)                                                                                              

          return
          
      ENTRY PUT_full_block2(nxs, nxe, nys, nye, nzs, nze, output_file_shot,cdf_id)

          mct = mct + 1
          call ncvpt1(cdf_id,axis_id_time,mct,time,rcode)



          call ncvpt (cdf_id, var_id_u,(/1,1,1,mct/),(/nxe-nxs+2, nye-nys+1, nze-nzs+1, 1/),
     &              real ( u(nxs-1:nxe,nys:nye,nzs:nze ), kind = 4 ), rcode)

          call ncvpt (cdf_id, var_id_v, (/1,1,1,mct/),(/nxe-nxs+1, nye-nys+2, nze-nzs+1, 1/), 
     &         real ( v(nxs:nxe,nys-1:nye,nzs:nze ), kind = 4 ), rcode)

          if (toggle(3) .eq. 1) then
          ww = - w
          call fill_value(ww,0,nx+1,0,ny+1,-1,nz+1)
          call ncvpt (cdf_id, var_id_w, (/1, 1, 1, mct/),(/ nxe-nxs+1, nye-nys+1, nze-nzs+2, 1/), 
     &         real ( ww(nxs:nxe,nys:nye, nzs-1: nze ), kind = 4 ), rcode)
          endif

          if (toggle(4) .eq. -1) then
          call ncvpt (cdf_id, var_id_p, (/1, 1, 1, mct/),(/nxe-nxs+1, nye-nys+1, nze-nzs+1, 1/), 
     &              real( p(0:nx+1,nys:nye,nzs:nze ), kind = 4 ), rcode)

          endif

          if (toggle(5) .eq. 1) then
!          call ncvpt (cdf_id, var_id_ss, (/1, 1, 1, mct/),(/nxe-nxs+1, nye-nys+1, nze-nzs+1, 1/), 
!     &              real ( ss(0:nx+1,nys:nye,nzs:nze ), kind = 4 ), rcode)

          call ncvpt (cdf_id, var_id_epsh, (/1, 1, 1, mct/),(/nxe-nxs+1, nye-nys+1, nze-nzs+1, 1/),
     &              real ( epsh(nxs:nxe,nys:nye,nzs:nze ), kind = 4 ), rcode)

          call ncvpt (cdf_id, var_id_epsv, (/1, 1, 1, mct/),(/nxe-nxs+1, nye-nys+1, nze-nzs+1, 1/),
     &              real ( epsv(nxs:nxe, nys:nye,nzs:nze ), kind = 4 ), rcode)
          call ncvpt (cdf_id, var_id_chi_v, (/1, 1, 1, mct/),(/nxe-nxs+1, nye-nys+1, nze-nzs+1, 1/),
     &              real ( chi_v(nxs:nxe, nys:nye,nzs:nze ), kind = 4 ), rcode)
          call ncvpt (cdf_id, var_id_chi_h, (/1, 1, 1, mct/),(/nxe-nxs+1, nye-nys+1, nze-nzs+1, 1/),
     &              real ( chi_h(nxs:nxe,nys:nye,nzs:nze ), kind = 4 ), rcode)
          call ncvpt (cdf_id, var_id_Ri, (/1, 1, 1, mct/),(/nxe-nxs+1, nye-nys+1, nze-nzs+1, 1/),
     &              real ( Ri(nxs:nxe,nys:nye,nzs:nze ), kind = 4 ), rcode)
          call ncvpt (cdf_id, var_id_rho, (/1, 1, 1, mct/),(/nxe-nxs+1, nye-nys+1, nze-nzs+1, 1/),
     &              real ( rrho(nxs:nxe,nys:nye,nzs:nze ), kind = 4 ), rcode)
          endif

!          if (toggle(6) .eq. -1) then
          call ncvpt (cdf_id, var_id_temp, (/1, 1, 1, mct/),(/nxe-nxs+1, nye-nys+1, nze-nzs+1, 1/),
     &               real ( t(nxs:nxe,nys:nye,nzs:nze ), kind = 4), rcode)
!          endif


!          if (toggle(7) .eq. -1) then
          call ncvpt (cdf_id, var_id_s,(/1, 1, 1, mct/),(/nxe-nxs+1, nye-nys+1, nze-nzs+1, 1/), 
     &               real ( s(nxs:nxe,nys:nye,nzs:nze  ), kind = 4), rcode)
!          endif

          if (toggle(8) .eq. -1) then
          call ncvpt (cdf_id, var_id_c, (/1, 1, 1, mct/),(/nxe-nxs+1, nye-nys+1, nze-nzs+1, 1/), 
     &               real ( c( nxs:nxe,nys: nye,nzs:nze  ), kind = 4), rcode)
          endif

          if (toggle(9) .eq. -1) then
          call ncvpt (cdf_id, var_id_cc1, (/1, 1, 1, mct/),(/nxe-nxs+1, nye-nys+1, nze-nzs+1, 1/),
     &           real ( cc1(0:nx+1,nys:nye,nzs:nze  ), kind =4 ), rcode)
          endif

          if (toggle(10) .eq. -1) then                                                                                      
          call ncvpt (cdf_id, var_id_cc2, (/1, 1, 1, mct/),(/nxe-nxs+1, nye-nys+1, nze-nzs+1, 1/), 
     &                 real ( cc2(0:nx+1,nys:nye,nzs:nze  ), kind =4 ), rcode)
          endif

          if (toggle(11) .eq. -1) then
          call ncvpt (cdf_id, var_id_cc3, (/1, 1, 1, mct/),(/nxe-nxs+1, nye-nys+1, nze-nzs+1, 1/), 
     &               real ( cc3(0:nx+1,nys:nye,nzs:nze  ), kind =4 ), rcode)
          endif

          if (toggle(12) .eq. -1) then
          call ncvpt (cdf_id, var_id_cc4,(/1, 1, 1, mct/),(/nxe-nxs+1, nye-nys+1, nze-nzs+1, 1/), 
     &            real ( cc4(0:nx+1,0:ny+1,nzs:nze  ), kind =4 ), rcode)
          endif

          return

      ENTRY PUT_full_block2_CLOSE(nxs, nxe, nys, nye, nzs, nze,output_file_shot,cdf_id)
       
          call  ncclos(cdf_id,rcode)
          return
    
      end

