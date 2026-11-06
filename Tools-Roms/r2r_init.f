      program r2r_init
      use mod_io_size_acct
      implicit none
      character(len=16) VertCoordType, str
      character(len=80) chldgrd, prntgrd, prnt_data
      character(len=12), parameter :: roms_init='roms_init.nc'
      character(len=10), parameter :: time_var_name='ocean_time'
      integer(kind=4), parameter :: ntrc=2
      character(len=4) tracer(ntrc)
      real(kind=8) theta_s,theta_b, hc,hcp, xcmin,xcmax, cff,cff1,cff2
      real(kind=8), allocatable, dimension(:) :: Cs_w,Cs_r,Csp_w,Csp_r,
     &                                           srX,srY,sXY,sYX, kprnt
      real(kind=8), allocatable, dimension(:,:) :: h,xc,yc,  hp,xp,yp,
     &      xi,eta, xpu,ypu,xiu,etau, xpv,ypv,xiv,etav, hprnt, csA,snA
      integer(kind=4), allocatable, dimension(:,:) :: ip,jp, ipu,jpu,
     &                                                       ipv,jpv
      integer(kind=2), allocatable, dimension(:,:) :: mskp, umsp, vmsp,
     &                                                mask, umask,vmask
      real(kind=4), allocatable, dimension(:) :: wrk1,wrk2,wrk3,wrmx
      integer(kind=4) net_alloc_mem, nargs, rec, ncpgrd,nccgrd, 
     &                           ncsrc,nctarg,
     &        nx,ny,Np, ncx,ncy,N, i,j,k, itrc,isrc,itrg, size, varid,
     &        tvar_in,tvar_out, natts, ierr, lstr, lpgrd,lprnt,lcgrd
      integer(kind=4) read_str_int, read_str_dble, lenstr
      integer(kind=4), parameter :: btm_slp=-1
     &                    , btm_trc=+1
      integer(kind=4) iwest,jsouth, imin,imax,jmin,jmax
      integer(kind=4), parameter :: margn=8
      real(kind=4) tstart,tend
      integer(kind=4) nclk, iclk(2), clk_rate, clk_max
      integer(kind=8) :: inc_clk, net_clk=0
      real(kind=8), parameter  ::  pi=3.14159265358979323846D0,
     &     Eradius=6371315.D0, SiderealDay=23.9344699D0*3600.D0,
     &                   day2sec=86400.D0, sec2day=1.D0/day2sec,
     &                     deg2rad=pi/180.D0, rad2deg=180.D0/pi
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
      call cpu_time(tstart)
      nclk=1
      call system_clock(iclk(nclk), clk_rate, clk_max)
      tracer(1)='temp' ; tracer(2)='salt'
      nargs=iargc()
      if (nargs == 8) then
        call getarg(1,chldgrd) ; lcgrd=lenstr(chldgrd)
        ierr=nf_open(chldgrd(1:lcgrd), nf_nowrite, nccgrd)
        if (ierr == nf_noerr) then
          call getarg(2,str) ; ierr=read_str_dble(str, theta_s)
          if (ierr == 0) then
            call getarg(3,str) ; ierr=read_str_dble(str,theta_b)
            if (ierr == 0) then
              call getarg(4,str) ; ierr=read_str_dble(str, hc)
              if (ierr == 0) then
                call getarg(5,str) ; ierr=read_str_int(str, N)
                if (ierr == 0) then
                  call getarg(6,prntgrd) ; lpgrd=lenstr(prntgrd)
                  ierr=nf_open(prntgrd(1:lpgrd), nf_nowrite, ncpgrd)
                  if (ierr == nf_noerr) then
                    call getarg(7,prnt_data) ; lprnt=lenstr(prnt_data)
                    ierr=nf_open(prnt_data(1:lprnt), nf_nowrite, ncsrc)
                    if (ierr == nf_noerr) then
                      call getarg(8,str) ; ierr=read_str_int(str, rec)
                      if (ierr == 0) then
                        write(*,'(/1x,3A,T52,A,I4,T62,A,F9.5)')
     &                               'r2r_init args: child grid = ''',
     &                                         chldgrd(1:lcgrd), '''',
     &                                 'N =',N,  'theta_s =', theta_s
                        write(*,'(15x,3A,T62,A,F9.5)')
     &                           'parent grid = ''', prntgrd(1:lpgrd),
     &                                    '''',  'theta_b =', theta_b
                        write(*,'(15x,3A,T52,A,I4,T67,A,F9.3/)')
     &                        'parent data = ''',  prnt_data(1:lprnt),
     &                                '''', 'rec =', rec,  'hc =', hc
                        goto 1
                      else
                        write(*,'(/1x,3A/)')     '### ERROR: Wrong ',
     &                     'argument #8, rec = ', str(1:lenstr(str))
                      endif
                    else
                      write(*,'(/1x,4A/22x,A/)') '### ERROR: arg #7 ',
     &                        ':: Cannot open ''', prnt_data(1:lprnt),
     &                                       '''.', nf_strerror(ierr)
                    endif
                  else
                    write(*,'(/1x,4A/22x,A/)') '### ERROR: arg #6 ::',
     &                     ' Cannot open ''', prntgrd(1:lpgrd), '''.',
     &                                              nf_strerror(ierr)
                  endif
                else
                  write(*,'(/1x,3A/)') '### ERROR: Wrong argument #5',
     &                                   ', N = ', str(1:lenstr(str))
                endif
              else
                write(*,'(/1x,3A/)') '### ERROR: Wrong argument #4, ',
     &                                   'hc = ',  str(1:lenstr(str))
              endif
            else
              write(*,'(/1x,3A/)')   '### ERROR: Wrong argument #3, ',
     &                               'theta_b = ', str(1:lenstr(str))
            endif
          else
            write(*,'(/1x,3A/)')     '### ERROR: Wrong argument #2, ',
     &                               'theta_s = ', str(1:lenstr(str))
          endif
        else
          write(*,'(/1x,4A/22x,A/)') '### ERROR: arg #1 :: Cannot ',
     &        'open ''', chldgrd(1:lcgrd), '''.', nf_strerror(ierr)
        endif
      else
        write(*,'(2(/1x,2A)/1x,A/)')    'This command interpolates ',
     &           'user-specified record of the parent-grid solution',
     &  'into the child grid and creates an initial condition file ',
     &  'to start the child', 'grid model.  It should be invoked as'
        write(*,'(/2x,2A/)')       'r2r_init child_grid.nc theta_s ',
     &              'theta_b hc N parent_grid.nc parent_data.nc rec'
        write(*,'(3(/1x,2A)/)')    'where theta_s theta_b hc N are ',
     &              'for vertical coordinate of the intended child-',
     &     'grid configuration (parent-grid vertical coordinate is ',
     &     'determined automatically',  'from the parent-grid data ',
     &         'file);  rec is record number within parent_data.nc.'
      endif
      stop
  1   call roms_find_dims(nccgrd, chldgrd, i,j,k)
      ncx=i+2 ;  ncy=j+2
      allocate(Cs_r(N), Cs_w(0:N))
      call set_scoord(theta_s,theta_b, N, Cs_r,Cs_w)
      call roms_find_dims(ncsrc,  prnt_data, i,j,k)
      call roms_check_dims(ncpgrd, prntgrd, i,j,0)
      nx=i+2 ; ny=j+2 ; Np=k
      allocate(Csp_w(0:Np),Csp_r(Np))
      call read_scoord(ncsrc, Np, Csp_r,Csp_w,hcp, VertCoordType)
      allocate( csA(ncx,ncy),snA(ncx,ncy), h(ncx,ncy), hprnt(ncx,ncy),
     &               mask(ncx,ncy), umask(ncx-1,ncy), vmask(ncx,ncy-1),
     &                                                     stat=ierr )
      if (ierr == 0) then
        net_alloc_mem=2*(2*N+1) +2*(2*Np+1) +4*2*ncx*ncy +3*ncx*ncy/2
        allocate( xi(ncx,ncy),eta(ncx,ncy), xiu(ncx,ncy),etau(ncx,ncy),
     &                                      xiv(ncx,ncy),etav(ncx,ncy),
     &                                                     stat=ierr )
        if (ierr == 0) then
          net_alloc_mem=net_alloc_mem +6*2*ncx*ncy
          allocate( ip(ncx,ncy),jp(ncx,ncy), ipu(ncx,ncy),jpu(ncx,ncy),
     &                                       ipv(ncx,ncy),jpv(ncx,ncy),
     &                                                     stat=ierr )
          if (ierr == 0) then
            net_alloc_mem = net_alloc_mem + 6*ncx*ncy
            allocate( xc(ncx,ncy),yc(ncx,ncy), stat=ierr)
            if (ierr == 0) then
              net_alloc_mem = net_alloc_mem + 2*2*ncx*ncy
            endif
          endif
        endif
      endif
      if (ierr == 0) then
        write(*,'(10x,A,F9.2,1x,A)') 'allocated', dble(net_alloc_mem)
     &                                          /262144.D0, 'MB total'
      else
        write(*,'(/1x,A/)')  '### ERROR: Memory allocation failure 1.'
        stop
      endif
      ierr=nf_inq_varid(nccgrd, 'mask_rho', varid)
      if (ierr == nf_noerr) then
        ierr=nf_get_var_double(nccgrd, varid, xc)
        if (ierr == nf_noerr) then
          call set_contour_mask(ncx,ncy, xc, mask,umask,vmask)
        else
          write(*,'(/1x,5A/)') '### ERROR: Cannot read ''mask_rho'' ',
     &         'from ''', chldgrd(1:lcgrd), ''', ', nf_strerror(ierr)
          stop
        endif
      else
        mask=1 ; umask=1 ; vmask=1
        write(*,'(9x,4A)')   'No land mask ''mask_rho'' is present ',
     &  'in ''', chldgrd(1:lcgrd), ''', assuming mask=1 everywhere.'
      endif
      call get_var_by_name_double(nccgrd, 'h',  h)
      call get_var_by_name_double(nccgrd, 'lon_rho',  xc)
      call get_var_by_name_double(nccgrd, 'lat_rho',  yc)
      write(*,'(1x,A)',advance='no') 'child grid longitude '
      call compute_min_max(ncx,ncy, xc, xcmin,xcmax)
      allocate(xp(nx,ny), yp(nx,ny), stat=ierr)
      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem + 2*2*nx*ny
      endif
      call get_var_by_name_double(ncpgrd, 'lon_rho',  xp)
      call get_var_by_name_double(ncpgrd, 'lat_rho',  yp)
      write(*,'(1x,A)',advance='no') 'parent grid longitude'
      call adjust_lon_into_range(nx,ny, xp, xcmin,xcmax)
C$OMP PARALLEL SHARED(nx,ny, xp,yp, ncx,ncy, xc,yc, ip,jp)
      call r2r_interp_search_thread(  nx,ny, xp,yp, ncx,ncy, xc,yc,
     &                                                        ip,jp)
C$OMP END PARALLEL
      call check_search_indices(nx,ny, xp,yp, ncx,ncy, xc,yc, ip,jp)
      call compute_index_bounds( ncx,ncy, ip,jp, mask, imin,imax,
     &                                                 jmin,jmax)
      write(*,'(1x,2A/2(4x,A,2I5,1x,A,I5))') 'minimal parent-grid ',
     &                    'index bounds to accommodate child grid:',
     &                     'imin,imax =', imin,imax, 'of nx =', nx,
     &                     'jmin,jmax =', jmin,jmax, 'of ny =', ny
      imin=max(imin-margn,1)  ; jmin=max(jmin-margn,1)
      imax=min(imax+margn,nx) ; jmax=min(jmax+margn,ny)
      write(*,'(1x,A/4x,A,2I5,17x,A,2I5)')  'adjusted to',
     &   'imin,imax =',imin,imax, 'jmin,jmax =',jmin,jmax
      iwest=imin ; nx=imax-imin+1 ; jsouth=jmin ; ny=jmax-jmin+1
      write(*,'(/2(2x,A,I5)/)') 'setting subdomain sizes  nx =', nx,
     &                                                   'ny =', ny
      deallocate(xp,yp)
      net_alloc_mem=net_alloc_mem -2*2*nx*ny
      write(*,'(8x,A,F9.2,1x,A)') 'deallocated',
     &           dble(2*2*nx*ny)/262144.D0,'MB'
      allocate( srX(nx*ny), srY(nx*ny), sXY(nx*ny), sYX(nx*ny),
     &                                              stat=ierr )
      if (ierr == 0) then
        allocate(hp(nx,ny), mskp(nx,ny), umsp(nx-1,ny), vmsp(nx,ny-1),
     &           xp(nx,ny), yp(nx,ny),    xpu(nx-1,ny),  xpv(nx,ny-1),
     &                                    ypu(nx-1,ny),  ypv(nx,ny-1),
     &                                                    stat=ierr )
      endif
      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem +(4+7)*2*nx*ny +3*(nx*ny)/2
        write(*,'(10x,A,F9.2,1x,A)') 'allocated',dble(net_alloc_mem)
     &                                         /262144.D0, 'MB total'
      else
        write(*,'(/1x,A/)') '### ERROR: Memory allocation failure 2.'
        stop
      endif
      call get_patch_by_name_double(ncpgrd, prntgrd, 'lon_rho',
     &                            iwest,jsouth, nx,ny,0,0, xp)
      call get_patch_by_name_double(ncpgrd, prntgrd, 'lat_rho',
     &                            iwest,jsouth, nx,ny,0,0, yp)
      write(*,'(1x,A)',advance='no') 'parent grid longitude'
      call adjust_lon_into_range(nx,ny, xp, xcmin,xcmax)
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      write(*,'(/F8.2,1x,A)') dble(net_clk)/dble(clk_rate),
     &           'initializing horizontal interpolation...'
C$OMP PARALLEL SHARED(nx,ny, xp,yp, ncx,ncy, xc,yc, ip,jp, xi,eta)
      call r2r_interp_init_thread(nx,ny, xp,yp, ncx,ncy, xc,yc,
     &                                               ip,jp, xi,eta)
C$OMP END PARALLEL
      call check_search_indices(nx,ny,xp,yp, ncx,ncy, xc,yc, ip,jp)
      call check_offsets(nx,ny,xp,yp, ncx,ncy,xc,yc, ip,jp, xi,eta)
      do j=1,ny
        do i=1,nx-1
          xpu(i,j)=0.5D0*(xp(i,j)+xp(i+1,j))
          ypu(i,j)=0.5D0*(yp(i,j)+yp(i+1,j))
        enddo
      enddo
C$OMP PARALLEL SHARED(nx,ny, xpu,ypu, ncx,ncy,xc,yc, ipu,jpu,xiu,etau)
      call r2r_interp_init_thread(nx-1,ny, xpu,ypu, ncx,ncy, xc,yc,
     &                                               ipu,jpu, xiu,etau)
C$OMP END PARALLEL
      call check_search_indices(nx-1,ny,xpu,ypu, ncx,ncy,xc,yc,ipu,jpu)
      call check_offsets(nx-1,ny,xpu,ypu, ncx,ncy,xc,yc,
     &                                               ipu,jpu,xiu,etau)
      do j=1,ny-1
        do i=1,nx
          xpv(i,j)=0.5D0*(xp(i,j)+xp(i,j+1))
          ypv(i,j)=0.5D0*(yp(i,j)+yp(i,j+1))
        enddo
      enddo
C$OMP PARALLEL SHARED(nx,ny, xpv,ypv, ncx,ncy,xc,yc, ipv,jpv,xiv,etav)
      call r2r_interp_init_thread(nx,ny-1, xpv,ypv, ncx,ncy, xc,yc,
     &                                               ipv,jpv, xiv,etav)
C$OMP END PARALLEL
      call check_search_indices(nx,ny-1,xpv,ypv, ncx,ncy,xc,yc,ipv,jpv)
      call check_offsets(nx,ny-1,xpv,ypv, ncx,ncy,xc,yc,
     &                                               ipv,jpv,xiv,etav)
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      write(*,'(F8.2,1x,A,F8.2,1x,A/)') dble(net_clk)/dble(clk_rate),
     &    'initialization complete in', dble(inc_clk)/dble(clk_rate),
     &    'sec'
      ierr=nf_inq_varid(ncpgrd, 'mask_rho', varid)
      if (ierr == nf_noerr) then
        call get_patch_by_name_double(ncpgrd, prntgrd, 'mask_rho',
     &                               iwest,jsouth, nx,ny,0,0, hp)
        call set_mask(nx,ny, hp, mskp,umsp,vmsp)
      else
        mskp=1 ; umsp=1 ; vmsp=1
        write(*,'(9x,4A)')   'No land mask ''mask_rho'' is present ',
     &  'in ''', prntgrd(1:lpgrd), ''', assuming mask=1 everywhere.'
      endif
      call get_patch_by_name_double(ncpgrd, prntgrd, 'h',
     &                      iwest,jsouth, nx,ny,0,0, hp)
      call spln2d_double(nx,ny, hp, srX,srY,sXY,sYX)
C$OMP PARALLEL SHARED(nx,ny, hp,  srX,srY,sXY, ncx,ncy,
C$OMP&                            ip,jp, xi,eta, hprnt)
      call spln2d_interp_double(nx,ny, hp,  srX,srY,sXY,
     &                   ncx,ncy, ip,jp, xi,eta, hprnt)
C$OMP END PARALLEL
      do j=1,ncy
        do i=1,ncx
          if (hprnt(i,j) < 0.0001D0) hprnt(i,j)=h(i,j)
        enddo
      enddo
      call read_angle(nccgrd, chldgrd, 1,1,    ncx,ncy,   csA,snA)
      call read_angle(ncpgrd, prntgrd, iwest,jsouth, nx,ny, xp,yp)
      call spln2d_double(nx,ny, xp, srX,srY,sXY,sYX)
C$OMP PARALLEL SHARED(nx,ny, xp, srX,srY,sXY, ncx,ncy,
C$OMP&                              ip,jp, xi,eta, xc)
      call spln2d_interp_double(nx,ny, xp, srX,srY,sXY,
     &                     ncx,ncy, ip,jp, xi,eta, xc)
C$OMP END PARALLEL
      call spln2d_double(nx,ny, yp, srX,srY,sXY,sYX)
C$OMP PARALLEL SHARED(nx,ny, yp, srX,srY,sXY, ncx,ncy,
C$OMP&                              ip,jp, xi,eta, yc)
      call spln2d_interp_double(nx,ny, yp, srX,srY,sXY,
     &                     ncx,ncy, ip,jp, xi,eta, yc)
C$OMP END PARALLEL
      do j=1,ncy
        do i=1,ncx
          cff=1.D0/sqrt(xc(i,j)*xc(i,j)+yc(i,j)*yc(i,j))
          cff1=csA(i,j)*xc(i,j) +snA(i,j)*yc(i,j)
          cff2=snA(i,j)*xc(i,j) -csA(i,j)*yc(i,j)
          csA(i,j)=cff*cff1  ;  snA(i,j)=cff*cff2
        enddo
      enddo
      deallocate(hp, xp,yp, xpu,ypu, xpv,ypv, xc,yc)
      net_alloc_mem=net_alloc_mem-(7*2*nx*ny+2*2*ncy*ncy)
      write(*,'(8x,A,F10.2,1x,A)') 'deallocated',
     &      dble(7*2*nx*ny+2*2*ncy*ncy)/262144.D0,'MB'
      ierr=nf_close(nccgrd) ; ierr=nf_close(ncpgrd)
      size=ncx*ncy*N ; allocate(kprnt(size), stat=ierr)
      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem + 2*size
        write(*,'(/8x,A,F9.2,1x,A/)')    'allocated kprnt, reaching',
     &                      dble(net_alloc_mem)/262144.D0, 'MB total'
      else
        write(*,'(/1x,A/)') '### ERROR: Memory allocation failure 3.'
        stop
      endif
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      write(*,'(/F8.2,1x,A)') dble(net_clk)/dble(clk_rate),
     &          'initializing vertical interpolation...'
C$OMP PARALLEL SHARED(ncx,ncy, hprnt,Np,hcp,Csp_r, h,N,hc,Cs_r, kprnt)
      call r2r_init_vertint_thread(ncx,ncy,  hprnt, Np,hcp,Csp_r,
     &                                          h, N, hc, Cs_r, kprnt)
      call r2r_check_vertint_thread(ncx,ncy, hprnt, Np,hcp,Csp_r,
     &                                          h, N, hc, Cs_r, kprnt)
C$OMP END PARALLEL
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      write(*,'(F8.2,1x,A,F8.2,1x,A/)') dble(net_clk)/dble(clk_rate),
     &    'initialization complete in', dble(inc_clk)/dble(clk_rate),
     &    'sec'
      call r2r_init_diag_file(ncx,ncy, N, ip,jp,xi,eta,
     &               ipu,jpu,xiu,etau, ipv,jpv,xiv,etav,
     &                        csA,snA, h, hprnt, kprnt)
      write(*,'(1x,A/)') 'initialization complete'
      size=max(nx*ny*(Np+1), ncx*ncy*(N+1))
      allocate(wrk1(size), stat=ierr)
      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem+size
        write(*,'(8x,A,F10.2,1x,A)')    'allocated wrk1, reaching',
     &                    dble(net_alloc_mem)/262144.D0, 'MB total'
        allocate(wrk2(size), stat=ierr)
        if (ierr == 0) then
          net_alloc_mem=net_alloc_mem+size
          write(*,'(8x,A,F10.2,1x,A)')   'allocated wrk2, reaching',
     &                     dble(net_alloc_mem)/262144.D0, 'MB total'
          allocate(wrk3(size), stat=ierr)
          if (ierr == 0) then
            net_alloc_mem=net_alloc_mem+size
            write(*,'(8x,A,F10.2,1x,A)')  'allocated wrk3, reaching',
     &                      dble(net_alloc_mem)/262144.D0, 'MB total'
            size=max(size, ncx*ncy*(Np+1))
            allocate(wrmx(size), stat=ierr)
            if (ierr == 0) then
              net_alloc_mem=net_alloc_mem+size
              write(*,'(8x,A,F10.2,1x,A/)')'allocated wrmx, reaching',
     &                       dble(net_alloc_mem)/262144.D0, 'MB total'
            endif
          endif
        endif
      endif
      if (ierr /= 0) then
        write(*,'(/1x,A/)') '### ERROR: Memory allocation failure 4.'
        stop
      endif
      ierr=nf_create(roms_init, nf_64bit_offset, nctarg)
      if (ierr == nf_noerr) then
        call def_roms_file(nctarg, roms_init, ncx,ncy,N,
     &                     theta_s,theta_b, hc, Cs_w,Cs_r)
        ierr=nf_put_att_text(nctarg, nf_global, 'memo', 19,
     &                              'created by r2r_init')
        if (ierr == nf_noerr) then
          ierr=nf_put_att_text(nctarg, nf_global, 'parent_grid',
     &                                 lpgrd, prntgrd(1:lpgrd))
          if (ierr == nf_noerr) then
            ierr=nf_put_att_text(nctarg, nf_global,
     &            'parent_data_file', lprnt, prnt_data(1:lprnt))
            if (ierr == nf_noerr) then
              ierr=nf_put_att_int(nctarg, nf_global,
     &             'parent_data_file_record', nf_int, 1, rec)
              if (ierr == nf_noerr) then
                write(*,*) 'added r2r_init attributes'
              else
                write(*,*) '### ERROR 5: put_att_int'
              endif
            else
              write(*,*) '### ERROR 4: put_att_text'
            endif
          else
            write(*,*) '### ERROR 3: put_att_text'
          endif
        else
          write(*,*) '### ERROR 2: put_att_text'
        endif
      else
        write(*,*) '### ERROR 1: nf_create'
      endif
      if (ierr /= nf_noerr) stop
      ierr=nf_inq_varid(nctarg, time_var_name, tvar_out)
      if (ierr == nf_noerr) then
        ierr=nf_inq_varid(ncsrc, time_var_name, tvar_in)
        if (ierr == nf_noerr) then
          ierr=nf_inq_varnatts(ncsrc, tvar_in, natts)
          if (ierr == nf_noerr) then
            do i=1,natts
              ierr=nf_inq_attname(ncsrc, tvar_in, i, str)
              if (ierr == nf_noerr) then
                lstr=lenstr(str)
                ierr=nf_copy_att(ncsrc, tvar_in, str(1:lstr),
     &                                     nctarg, tvar_out)
                if (ierr == nf_noerr) then
                  write(*,*) 'copied attribute ''',str(1:lstr),''''
                else
                  write(*,*) '### ERROR 10: copy_att' ; stop
                endif
              else
                write(*,*) '### ERROR 9: inq_attname' ; stop
              endif
            enddo
          else
            write(*,*) '### ERROR 8: inq_varnatts'
          endif
        else
          write(*,*) '### ERROR 7: inq_varid'
        endif
      else
        write(*,*) '### ERROR 6: inq_varid'
      endif
      if (ierr /= nf_noerr) stop
      ierr=nf_enddef(nctarg)
      if (ierr == nf_noerr) then
        ierr=nf_get_vara_double(ncsrc, tvar_in, rec,1, cff1)
        if (ierr == nf_noerr) then
          ierr=nf_put_var_double(nctarg, tvar_out, cff1)
          if (ierr == nf_noerr) then
            write(*,*) 'wrote ''', time_var_name, ''''
          else
            write(*,*) '### ERROR 13: put_var'
          endif
        else
          write(*,*) '### ERROR 12: get_vara'
        endif
      else
        write(*,*) '### ERROR 11: nf_enddef'
      endif
      if (ierr /= nf_noerr) stop
      call get_patch_by_name_real(ncsrc, prnt_data, 'zeta',
     &                  iwest,jsouth, nx,ny,0,rec,  wrk1)
C$OMP PARALLEL SHARED(nx,ny, mskp, wrk1)
      call etch_into_land_thread(nx,ny, mskp, wrk1)
C$OMP END PARALLEL
      call spln2d_real(nx,ny, wrk1, srX,srY,sXY,sYX)
C$OMP PARALLEL SHARED(nx,ny, wrk1, srX,srY,sXY, ncx,ncy,
C$OMP&                              ip,jp, xi,eta, wrk2)
      call spln2d_interp_real(nx,ny, wrk1,  srX,srY,sXY,
     &                     ncx,ncy, ip,jp, xi,eta, wrk2)
C$OMP END PARALLEL
      call put_rec_by_name_real(nctarg,roms_init, 'zeta',
     &                                ncx,ncy,0,1, wrk2)
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      write(*,'(F8.2,1x,A,F8.2,1x,A/)') dble(net_clk)/dble(clk_rate),
     &            'complete in', dble(inc_clk)/dble(clk_rate), 'sec'
      call get_patch_by_name_real(ncsrc, prnt_data, 'u', iwest,jsouth,
     &                                         nx-1,ny,Np, rec, wrk1)
      do k=1,Np
        isrc=1+(k-1)*(nx-1)*ny ; itrg=1+(k-1)*ncx*ncy
C$OMP PARALLEL SHARED(nx,ny, umsp, wrk1, isrc)
        call etch_into_land_thread(nx-1,ny, umsp, wrk1(isrc))
C$OMP END PARALLEL
        call spln2d_real(nx-1,ny, wrk1(isrc), srX,srY,sXY,sYX)
C$OMP PARALLEL SHARED(nx,ny, isrc, wrk1, srX,srY,sXY,
C$OMP&                   ncx,ncy, itrg, ipu,jpu, xiu,etau, wrmx)
        call spln2d_interp_real(nx-1,ny, wrk1(isrc), srX,srY,sXY,
     &                   ncx,ncy, ipu,jpu, xiu,etau, wrmx(itrg))
C$OMP END PARALLEL
        write(6,'(A)',advance='no') '.' ; flush(6)
      enddo
      write(6,*)
C$OMP PARALLEL SHARED(ncx,ncy, mask,  Np,wrmx,  N,kprnt, wrk1)
      call r2r_vrtint_thread(ncx,ncy, 0,mask, btm_slp, Np,wrmx,
     &                                          N,kprnt, wrk1)
C$OMP END PARALLEL
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      write(*,'(F8.2,1x,A,F8.2,1x,A/)') dble(net_clk)/dble(clk_rate),
     &            'complete in', dble(inc_clk)/dble(clk_rate), 'sec'
      call get_patch_by_name_real(ncsrc, prnt_data, 'v', iwest,jsouth,
     &                                         nx,ny-1,Np, rec, wrk2)
      do k=1,Np
        isrc=1+(k-1)*nx*(ny-1) ; itrg=1+(k-1)*ncx*ncy
C$OMP PARALLEL SHARED(nx,ny, vmsp, wrk2, isrc)
        call etch_into_land_thread(nx,ny-1, vmsp, wrk2(isrc))
C$OMP END PARALLEL
        call spln2d_real(nx,ny-1, wrk2(isrc), srX,srY,sXY,sYX)
C$OMP PARALLEL SHARED(nx,ny, isrc,  wrk2, srX,srY,sXY,
C$OMP&                   ncx,ncy, itrg, ipv,jpv, xiv,etav, wrmx)
        call spln2d_interp_real(nx,ny-1, wrk2(isrc), srX,srY,sXY,
     &                   ncx,ncy, ipv,jpv, xiv,etav, wrmx(itrg))
C$OMP END PARALLEL
        write(6,'(A)',advance='no') '.' ; flush(6)
      enddo
      write(6,*)
C$OMP PARALLEL SHARED(ncx,ncy, mask,  Np,wrmx,  N,kprnt, wrk2)
      call r2r_vrtint_thread(ncx,ncy, 0,mask, btm_slp, Np,wrmx,
     &                                          N,kprnt, wrk2)
C$OMP END PARALLEL
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      write(*,'(F8.2,1x,A,F8.2,1x,A/)') dble(net_clk)/dble(clk_rate),
     &            'complete in', dble(inc_clk)/dble(clk_rate), 'sec'
C$OMP PARALLEL SHARED(ncx,ncy, N, csA,snA, wrk1,wrk2, wrk3,wrmx)
      call r2r_rotate_shift_thread(ncx,ncy,N, csA,snA, wrk1,wrk2,
     &                                                 wrk3,wrmx)
!$OMP END PARALLEL
      call put_rec_by_name_real(nctarg, roms_init, 'u',
     &                              ncx-1,ncy,N,1,  wrk3)
      call put_rec_by_name_real(nctarg, roms_init, 'v',
     &                              ncx,ncy-1,N,1,  wrmx)
C$OMP PARALLEL SHARED(ncx,ncy, N, hc,Cs_w, h, wrk3,wrmx, wrk1,wrk2)
      call compute_uvbar_thread(ncx,ncy, N, hc,Cs_w, h,  wrk3,wrmx,
     &                                                   wrk1,wrk2)
C$OMP END PARALLEL
      call put_rec_by_name_real(nctarg, roms_init, 'ubar',
     &                              ncx-1,ncy,0,1,  wrk1)
      call put_rec_by_name_real(nctarg, roms_init, 'vbar',
     &                              ncx,ncy-1,0,1,  wrk2)
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      write(*,'(F8.2,1x,A,F8.2,1x,A/)') dble(net_clk)/dble(clk_rate),
     &            'complete in', dble(inc_clk)/dble(clk_rate), 'sec'
      do itrc=1,ntrc
        call get_patch_by_name_real(ncsrc, prnt_data, tracer(itrc),
     &                          iwest,jsouth, nx,ny,Np,rec, wrk1)
        do k=1,Np
          isrc=1+(k-1)*nx*ny ; itrg=1+(k-1)*ncx*ncy
C$OMP PARALLEL SHARED(nx,ny, wrk1, isrc)
          call etch_into_land_thread(nx,ny, mskp, wrk1(isrc))
C$OMP END PARALLEL
          call spln2d_real(nx,ny, wrk1(isrc), srX,srY,sXY,sYX)
C$OMP PARALLEL SHARED(nx,ny, isrc,  wrk1, srX,srY,sXY,
C$OMP&                    ncx,ncy, itrg, ip,jp, xi,eta, wrmx)
          call spln2d_interp_real(nx,ny, wrk1(isrc), srX,srY,sXY,
     &                      ncx,ncy, ip,jp, xi,eta, wrmx(itrg))
C$OMP END PARALLEL
          write(6,'(A)',advance='no') '.' ; flush(6)
        enddo
        write(6,*)
C$OMP PARALLEL SHARED(ncx,ncy, mask, Np,wrmx, N,kprnt,wrk1)
      call r2r_vrtint_thread(ncx,ncy, 1,mask, btm_trc, Np,wrmx,
     &                                        N,kprnt,wrk1)
C$OMP END PARALLEL
        call put_rec_by_name_real(nctarg,roms_init, tracer(itrc),
     &                                        ncx,ncy,N,1, wrk1)
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      write(*,'(F8.2,1x,A,F8.2,1x,A/)') dble(net_clk)/dble(clk_rate),
     &            'complete in', dble(inc_clk)/dble(clk_rate), 'sec'
      enddo
      ierr=nf_close(nctarg)
      write(*,'(/4x,A,F12.3,1x,A,F10.1,1x,A)')    'total data read ',
     &                  dble(sz_read_acc)/dble(1024**2), 'MBytes in',
     &                          dble(read_clk)/dble(clk_rate), 'sec'
      write(*,'(1x,A,F12.3,1x,A,F10.1,1x,A/)') 'total data written ',
     &                 dble(sz_write_acc)/dble(1024**2), 'MBytes in',
     &                         dble(write_clk)/dble(clk_rate), 'sec'
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate, clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      call cpu_time(tend)
      write(*,'(1x,A,F8.2,1x,A,4x,A,F8.2,1x,A,F8.1,1x,A)')
     &   'Wall Clock time:', dble(net_clk)/dble(clk_rate), 'sec',
     &                           'CPU time:', tend-tstart, 'sec',
     & (tend-tstart)*dble(clk_rate)/dble(net_clk)*100.D0,'% CPU'
      lstr=lenstr(roms_init)
      write(*,'(/1x,5A/)')  'Files ''', roms_init(1:lstr),
     &     ''' and ''', 'r2r_init_diag.nc', ''' are ready.'
      end
      subroutine def_roms_file(ncid, fname, xi_rho,eta_rho,s_rho,
     &                            theta_s,theta_b, hc, Cs_w,Cs_r)
      implicit none
      character(len=*) fname
      integer(kind=4) ncid, xi_rho,xi_u, eta_rho,eta_v, s_rho, ierr,
     &                        old_fill_mode, varid,  lfnm, lenstr
      real(kind=8) theta_s,theta_b, hc, Cs_w(0:s_rho),Cs_r(s_rho)
      integer(kind=4), parameter :: n2d=2, n3d=3
      integer(kind=4) r2dgrd(n2d), u2dgrd(n2d), v2dgrd(n2d),
     &        r3dgrd(n3d), u3dgrd(n3d), v3dgrd(n3d)
      real(kind=8), parameter :: cycle_length=360.D0
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
      lfnm=lenstr(fname) ;  write(*,'(1x,3A)', advance='no')
     &    'creating netCDF file ''', fname(1:lfnm), ''' ...'
      ierr=nf_set_fill(ncid, nf_nofill, old_fill_mode)
      if (ierr /= nf_noerr) then
        write(*,'(/1x,A/)') '### WARNING: Cannot set nofill mode.'
      endif
      xi_u=xi_rho-1 ; eta_v=eta_rho-1
      ierr=nf_def_dim(ncid, 'xi_rho',  xi_rho,    r2dgrd(1))
      ierr=nf_def_dim(ncid, 'xi_u',    xi_u,      u2dgrd(1))
      ierr=nf_def_dim(ncid, 'eta_rho', eta_rho,   r2dgrd(2))
      ierr=nf_def_dim(ncid, 'eta_v',   eta_v,     v2dgrd(2))
      ierr=nf_def_dim(ncid, 's_rho',   s_rho,     r3dgrd(3))
      v2dgrd(1)=r2dgrd(1) ; u2dgrd(2)=r2dgrd(2)
      r3dgrd(1)=r2dgrd(1) ; r3dgrd(2)=r2dgrd(2)
      u3dgrd(1)=u2dgrd(1) ; u3dgrd(2)=u2dgrd(2)
      v3dgrd(1)=v2dgrd(1) ; v3dgrd(2)=v2dgrd(2)
      u3dgrd(3)=r3dgrd(3) ; v3dgrd(3)=r3dgrd(3)
      ierr=nf_put_att_text(ncid,nf_global,'VertCoordType',4,'SM09')
      ierr=nf_put_att_double(ncid, nf_global,'theta_s', nf_double,
     &                                                 1, theta_s)
      ierr=nf_put_att_double(ncid, nf_global,'theta_b', nf_double,
     &                                                 1, theta_b)
      ierr=nf_put_att_double(ncid, nf_global, 'hc',nf_double,1,hc)
      ierr=nf_put_att_double(ncid, nf_global, 'Cs_w', nf_double,
     &                                              s_rho+1, Cs_w)
      ierr=nf_put_att_double(ncid, nf_global, 'Cs_r', nf_double,
     &                                                s_rho, Cs_r)
      ierr=nf_def_var (ncid, 'ocean_time',  nf_double, n2d-2,
     &                                      r3dgrd(n3d), varid)
      ierr=nf_put_att_text (ncid, varid, 'long_name', 22,
     &                             'time since initialization')
      ierr=nf_put_att_text (ncid, varid, 'units', 4, 'days')
      ierr=nf_def_var(ncid, 'zeta', nf_real, n2d, r2dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name', 22,
     &                                  'free-surface elevation')
      ierr=nf_put_att_text(ncid, varid, 'units', 5, 'meter')
      ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                 FillValue)
      ierr=nf_def_var(ncid, 'ubar', nf_real, n2d, u2dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name', 22,
     &                                   'barotropic XI-velocity')
      ierr=nf_put_att_text(ncid, varid, 'units',12,'meter/second')
      ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)
      ierr=nf_def_var(ncid, 'vbar', nf_real, n2d, v2dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name', 23,
     &                                  'barotropic ETA-velocity')
      ierr=nf_put_att_text(ncid, varid, 'units',12,'meter/second')
      ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)
      ierr=nf_def_var(ncid, 'u', nf_real, n3d, u3dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name', 21,
     &                                    'XI-velocity component')
      ierr=nf_put_att_text(ncid, varid, 'units',12,'meter/second')
      ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)
      ierr=nf_def_var(ncid, 'v', nf_real, n3d, v3dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name', 22,
     &                                   'ETA-velocity component')
      ierr=nf_put_att_text(ncid, varid, 'units',12,'meter/second')
      ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     &                                                  FillValue)
      ierr=nf_def_var(ncid, 'temp', nf_real, n3d, r3dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name', 21,
     &                                   'potential temperature')
      ierr=nf_put_att_text(ncid, varid, 'units', 7,    'Celsius')
      ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real,
     &                                              1, FillValue)
      ierr=nf_def_var(ncid, 'salt', nf_real, n3d, r3dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name', 8,'salinity')
      ierr=nf_put_att_text(ncid, varid, 'units',     3,     'PSU')
      ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real,
     &                                              1, FillValue)
      end
