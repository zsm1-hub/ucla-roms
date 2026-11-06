      subroutine def_bry_file(fname, ncid, xi_rho,eta_rho,N,
     &                OBC_WEST, OBC_EAST, OBC_SOUTH, OBC_NORTH,
     &                        theta_s, theta_b, hc, Cs_w,Cs_r)
      implicit none
      character(len=*) fname
      integer(kind=4)  ncid, xi_rho,eta_rho,N
      logical OBC_WEST, OBC_EAST, OBC_SOUTH, OBC_NORTH
      real(kind=8) theta_s,theta_b,hc, Cs_w(0:N),Cs_r(N)
      integer(kind=4) r2d_EW(2), u2d_EW(2), r3d_EW(3), u3d_EW(3),
     &        r2d_NS(2), v2d_NS(2), r3d_NS(3), v3d_NS(3),
     &        vert_dim, time_dim, varid, ierr,
     &        old_fill_mode, lenstr
      real(kind=4), parameter :: FillValue=1.D+33
      real(kind=8), parameter :: spv=1.D+9, spv_set=1.D+33
      integer(kind=4), parameter :: nf_byte=1, nf_char=2, nf_short=3, 
     &                             nf_int=4,
     &   nf_float=5,  nf_double=6,  nf_int2=nf_short,  nf_int1=nf_byte,
     &   nf_real=nf_float,   nf_clobber=0,   nf_64bit_offset=512,
     &   nf_netcdf4=4096,   nf_nowrite=0,   nf_write=1,   nf_chunked=0,
     &   nf_nofill=256, nf_unlimited=0, nf_global=0,
     &   nf_format_classic=1,  nf_format_64bit=2, nf_format_netcdf4=3,
     &   nf_format_netcdf4_classic=4,  nf_enameinuse=-42,  nf_noerr=0
      integer(kind=4), external :: nf_open, nf_inq, nf_inq_format,
     &   nf_inq_ndims,    nf_inq_nvars, nf_inq_unlimdim, nf_inq_dim,
     &   nf_inq_dimid,    nf_inq_dimname, nf_inq_dimlen, nf_inq_att,
     &   nf_inq_attname,  nf_inq_var, nf_inq_varid,  nf_inq_varname,
     &   nf_inq_varndims, nf_inq_varnatts, nf_inq_vardimid,
     &   nf_inq_var_chunking,
     &   nf_create, nf_set_fill, nf_def_dim, nf_def_var, nf_copy_att,
     &   nf_def_var_chunking, nf_def_var_deflate, nf_enddef,
     &   nf_redef, nf_rename_dim, nf_rename_att, nf_rename_var,
     &   nf_get_att_text,  nf_get_att_int1,  nf_get_att_int2,
     &   nf_get_att_int,   nf_get_att_real,  nf_get_att_double,
     &   nf_put_att_text,  nf_put_att_int1,  nf_put_att_int2,
     &   nf_put_att_int,   nf_put_att_real,  nf_put_att_double,
     &   nf_get_var_text,  nf_get_var_int1,  nf_get_var_int2,
     &   nf_get_var_int,   nf_get_var_real,  nf_get_var_double,
     &   nf_get_var1_text, nf_get_var1_int1, nf_get_var1_int2,
     &   nf_get_var1_int,  nf_get_var1_real, nf_get_var1_double,
     &   nf_get_vara_text, nf_get_vara_int1, nf_get_vara_int2,
     &   nf_get_vara_int,  nf_get_vara_real, nf_get_vara_double,
     &   nf_put_var_text,  nf_put_var_int1,  nf_put_var_int2,
     &   nf_put_var_int,   nf_put_var_real,  nf_put_var_double,
     &   nf_put_var1_text, nf_put_var1_int1, nf_put_var1_int2,
     &   nf_put_var1_int,  nf_put_var1_real, nf_put_var1_double,
     &   nf_put_vara_text, nf_put_vara_int1, nf_put_vara_int2,
     &   nf_put_vara_int,  nf_put_vara_real, nf_put_vara_double,
     &   nf_sync,          nf_close
      character(len=80), external :: nf_strerror
      write(*,'(1x,3A)',advance='no') 'def_bry_file :: creating ''',
     &                            fname(1:lenstr(fname)),  ''' ...'
      ierr=nf_set_fill(ncid, nf_nofill, old_fill_mode)
      ierr=nf_def_dim(ncid, 'xi_rho',  xi_rho,      r2d_EW(1))
      ierr=nf_def_dim(ncid, 'xi_u',    xi_rho-1,    u2d_EW(1))
      ierr=nf_def_dim(ncid, 'eta_rho', eta_rho,     r2d_NS(1))
      ierr=nf_def_dim(ncid, 'eta_v',   eta_rho-1,   v2d_NS(1))
      ierr=nf_def_dim(ncid, 's_rho',   N,            vert_dim)
      ierr=nf_def_dim(ncid, 'rec_time',nf_unlimited, time_dim)
      r3d_EW(1)=r2d_EW(1) ; r2d_EW(2)=time_dim
      u3d_EW(1)=u2d_EW(1) ; u2d_EW(2)=time_dim
      r3d_NS(1)=r2d_NS(1) ; r2d_NS(2)=time_dim
      v3d_NS(1)=v2d_NS(1) ; v2d_NS(2)=time_dim
      r3d_EW(2)=vert_dim  ; r3d_EW(3)=time_dim
      u3d_EW(2)=vert_dim  ; u3d_EW(3)=time_dim
      r3d_NS(2)=vert_dim  ; r3d_NS(3)=time_dim
      v3d_NS(2)=vert_dim  ; v3d_NS(3)=time_dim
      ierr=nf_put_att_text(ncid, nf_global, 'VertCoordType',4,'SM09')
      ierr=nf_put_att_double(ncid, nf_global,'theta_s', nf_double, 1,
     &                                                       theta_s)
      ierr=nf_put_att_double(ncid, nf_global,'theta_b', nf_double, 1,
     &                                                       theta_b)
      ierr=nf_put_att_double(ncid, nf_global, 'hc', nf_double, 1, hc)
      ierr=nf_put_att_double(ncid, nf_global, 'Cs_w', nf_double, N+1,
     &                                                          Cs_w)
      ierr=nf_put_att_double(ncid, nf_global, 'Cs_r', nf_double, N,
     &                                                         Cs_r)
      ierr=nf_def_var (ncid, 'bry_time', nf_double, 1, time_dim, varid)
      ierr=nf_put_att_text (ncid, varid, 'units', 4, 'days')
      if (OBC_WEST) then
        ierr=nf_def_var(ncid, 'zeta_west',nf_real, 2,r2d_NS,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 5, 'meter')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)
        ierr=nf_def_var(ncid, 'ubar_west',nf_real, 2,r2d_NS,varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real,  1,
     &                                                    FillValue)
        ierr=nf_def_var(ncid, 'vbar_west',nf_real, 2,v2d_NS,varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                    FillValue)
        ierr=nf_def_var(ncid, 'u_west',nf_real, 3, r3d_NS, varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                    FillValue)
        ierr=nf_def_var(ncid, 'v_west',nf_real, 3, v3d_NS, varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                    FillValue)
        ierr=nf_def_var(ncid, 'temp_west',nf_real, 3,r3d_NS,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 7, 'Celsius')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                   FillValue)
        ierr=nf_def_var(ncid, 'salt_west',nf_real, 3,r3d_NS,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 3, 'PSU')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)
      endif
      if (OBC_EAST) then
        ierr=nf_def_var(ncid, 'zeta_east',nf_real, 2,r2d_NS,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 5, 'meter')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)
        ierr=nf_def_var(ncid, 'ubar_east',nf_real, 2,r2d_NS,varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                    FillValue)
        ierr=nf_def_var(ncid, 'vbar_east',nf_real, 2,v2d_NS,varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                    FillValue)
        ierr=nf_def_var(ncid, 'u_east',nf_real, 3, r3d_NS, varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                    FillValue)
        ierr=nf_def_var(ncid, 'v_east',nf_real, 3, v3d_NS, varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                   FillValue)
        ierr=nf_def_var(ncid, 'temp_east',nf_real, 3,r3d_NS,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 7, 'Celsius')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                   FillValue)
        ierr=nf_def_var(ncid, 'salt_east',nf_real, 3,r3d_NS,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 3, 'PSU')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)
      endif
      if (OBC_SOUTH) then
        ierr=nf_def_var(ncid, 'zeta_south',nf_real,2,r2d_EW,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 5, 'meter')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)
        ierr=nf_def_var(ncid, 'ubar_south',nf_real,2,u2d_EW,varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                   FillValue)
        ierr=nf_def_var(ncid, 'vbar_south',nf_real,2,r2d_EW,varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                    FillValue)
        ierr=nf_def_var(ncid, 'u_south',nf_real, 3, u3d_EW, varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                    FillValue)
        ierr=nf_def_var(ncid, 'v_south',nf_real, 3, r3d_EW, varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                    FillValue)
        ierr=nf_def_var(ncid, 'temp_south',nf_real,3,r3d_EW,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 7, 'Celsius')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)
        ierr=nf_def_var(ncid, 'salt_south',nf_real,3,r3d_EW,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 3, 'PSU')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)
      endif
      if (OBC_NORTH) then
        ierr=nf_def_var(ncid, 'zeta_north',nf_real,2,r2d_EW,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 5, 'meter')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)
        ierr=nf_def_var(ncid, 'ubar_north',nf_real,2,u2d_EW,varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)
        ierr=nf_def_var(ncid, 'vbar_north',nf_real,2,r2d_EW,varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)
        ierr=nf_def_var(ncid, 'u_north',nf_real, 3, u3d_EW, varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                   FillValue)
        ierr=nf_def_var(ncid, 'v_north',nf_real, 3, r3d_EW, varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                    FillValue)
        ierr=nf_def_var(ncid, 'temp_north',nf_real,3,r3d_EW,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 7, 'Celsius')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)
        ierr=nf_def_var(ncid, 'salt_north',nf_real,3,r3d_EW,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 3, 'PSU')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)
      endif
      write(*,*) '...done.'
      end
