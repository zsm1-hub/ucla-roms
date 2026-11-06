      program zslice
      use mod_io_size_acct
      implicit none
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
      real(kind=8)  hc, cff
      real(kind=8), dimension(:), allocatable :: Cs_r, Cs_w,
     &                                    h, rmask, z_r, z_w
      real(kind=4), dimension(:), allocatable :: z_lev, bffr,
     &                                              bffr_out
      real(kind=4), parameter :: FillValue=1.D+33
      character(len=64) string,  trg_fname, VertCoordType
      character(len=64), dimension(:), allocatable :: fname
      character(len=32), dimension(:), allocatable :: dname,vname,
     &                                         sel_vars, tg_vname
      logical, dimension(:), allocatable :: Ld_in_use, lmask
      integer(kind=4), dimension(:), allocatable :: dsize,  vtype, 
     &                               vdims,
     &                     tg_dimid,  dimv
      character(len=4), dimension(:), allocatable :: var_grd_type,
     &                                            tg_var_grd_type
      integer(kind=4), dimension(:,:), allocatable :: dimids
      integer(kind=4) nargs,  iarg,  first_file, last_file, indx,  ierr,
     &        Lm, Mm, N, nz, max_horiz_pts, imin,jmin,kmin,i,j,k,
     &        ncsrc,  ndims, max_dims,   tg_nvars,  tsize,  rec,
     &        nctarg, nvars, nsel_vars,  max_vars,  ngatts, vatts,
     &        varid,  ivar,  is,ie, lsrc, ltrg, lstr, lvar, ldim,
     &        id_xi_rho,   id_eta_rho,   id_xi_u,     id_eta_v,
     &        id_xi_r_tgt, id_et_r_tgt,  id_xi_u_tgt, id_et_v_tgt,
     &        id_s_rho,    id_s_w,       id_time,     ndim_tgt,
     &        id_nz_tgt,   id_time_tgt,  unlimdimid,  net_alloc_size
      integer(kind=4) lenstr, read_str_dble
      real(kind=4) tstart, run_time
      integer(kind=4) iclk_start, iclk_end, iclk_rate, iclk_max
      call cpu_time(tstart)
      call system_clock(iclk_start, iclk_rate, iclk_max)
      net_alloc_size=0
      max_horiz_pts=0
      max_dims=0
      max_vars=0
      N=0
      nargs=iargc()
      if (nargs < 2) then
        write(*,'(3(/1x,2A)/2(/6x,2A/1x,A)/6x,2A/5(/1x,2A)/)')
     &     'Command-line operator "zslice" takes one or more ROMS-',
     &    'conforming netCDF files',     'and performs Sigma-to-Z ',
     &    'interpolation for either all(default) or user-selected',
     &   '3D fields available in the data.  It should be invoked ',
     &   'as', 'zslice z1 [z2 z3...] file1.nc [file2.nc file3.nc ',
     &   '...]', 'or', 'zslice z1 [z2 z3 ...] [--vars=v1,v2,v3] ',
     &            'file1.nc [file2.nc file3.nc...]',   'or',
     &           'zslice z1 [z2 z3...] [--vars=v1,v2,v3] ',
     &           'grid.nc file1.nc [file2.nc file3.nc ...]',
     &  'where z1...zN is the set of desired z-levels; optionally ',
     &  '--vars=v1,v2,v3 is', 'a comma-separated list of selected ',
     &              'variables (must match the actual names of',
     &     'netCDF variables; omitting this argument makes zslice ',
     &      'interpolate all 3D', 'fields available in the file). ',
     &         'The first file, either file1.nc or grid.nc,',
     &     'must contain bottom topography "h", and [optionally] ',
     &     'land mask data "rmask".'
        stop
      endif
      allocate(z_lev(nargs))
      iarg=0 ; nz=0 ; nsel_vars=0 ; ierr=0
      do while(ierr == 0 .and. iarg < nargs)
        iarg=iarg+1
        call getarg(iarg,string)
        lstr=lenstr(string)
        ierr=read_str_dble(string, cff)
        if (ierr == 0) then
          nz=nz+1
          z_lev(nz)=-abs(cff)
          write(*,'(4x,A,I2,A,F11.4)') 'z_lev(',nz,') =', z_lev(nz)
        elseif (iarg == 1) then
          write(*,'(/1x,4A/)')    '### ERROR: Cannot read z-level ',
     &          'from the first argument ''', string(1:lstr), '''.'
        endif
      enddo
      iarg=nz
      if (iarg < nargs+1) then
        if (string(1:7) == '--vars=') then
          nsel_vars=1
          do i=8,lstr
            if (string(i:i) == ',') nsel_vars=nsel_vars+1
          enddo
          allocate(sel_vars(nsel_vars))
          j=0 ; is=8
          do while(is < lstr+1)
            ie=is
            do while(string(ie:ie) /= ',' .and. ie < lstr)
              ie=ie+1
            enddo
            if (string(ie:ie) == ',') ie=ie-1
            j=j+1 ; sel_vars(j)=string(is:ie) ; is=ie+2
          enddo
          iarg=iarg+1
        endif
        last_file=nargs-iarg-1
        allocate(fname(0:last_file))
        do indx=0,last_file
          iarg=iarg+1
          call getarg(iarg,fname(indx))
          lsrc=lenstr(fname(indx))
          write(*,'(4x,I3,1x,3A)') indx, 'file =''',
     &                    fname(indx)(1:lsrc), ''''
        enddo
      else
        write(*,'(/1x,A/)') '### ERROR: No file names are given.'
        stop
      endif
      indx=0 ; lsrc=lenstr(fname(indx))
      ierr=nf_open(fname(indx), nf_nowrite, ncsrc)
      if (ierr == nf_noerr) then
        call roms_find_dims(ncsrc, fname(indx), Lm,Mm,N)
        max_horiz_pts=(Lm+2)*(Mm+2)
        allocate (h(max_horiz_pts))
        allocate (rmask(max_horiz_pts))
        net_alloc_size=net_alloc_size + 2*8*max_horiz_pts
        call get_var_by_name_double(ncsrc, 'h', h)
        if (ierr == nf_noerr) then
          ierr=nf_inq_varid(ncsrc, 'mask_rho', varid)
          if (ierr == nf_noerr) then
            ierr=nf_get_var_double(ncsrc, varid, rmask)
            if (ierr /= nf_noerr) then
              write(*,*)  '### ERROR: Cannot read ''mask_rho'' ',
     &                 'from ROMS file ''',  fname(indx)(1:lsrc),
     &                                 ''': ', nf_strerror(ierr)
            endif
          else
            write(*,'(/1x,4A/14x,2A/)')   '### WARNING: Cannot ',
     &                        'find land mask from grid file ''',
     &             fname(indx)(1:lsrc), ''': Assuming rmask = 1',
     &                        'everywhere:  ', nf_strerror(ierr)
            do i=1,max_horiz_pts
              rmask(i)=1.D0
            enddo
            ierr=nf_noerr
          endif
          first_file=min(1,last_file)
          if (N > 0) then
            allocate(bffr(2*(N+1)))
            ierr=nf_get_att_double(ncsrc, nf_global,'Cs_r', bffr)
            if (ierr == nf_noerr) then
              ierr=nf_get_att_double(ncsrc, nf_global,'Cs_w',bffr)
              if (ierr == nf_noerr) then
                ierr=nf_get_att_double(ncsrc, nf_global,'hc',bffr)
                if (ierr == nf_noerr) then
                  first_file=0
                endif
              endif
            endif
            deallocate(bffr)
            ierr=nf_noerr
          endif
          ierr=nf_close(ncsrc)
          if (ierr/=nf_noerr) then
            write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot close ROMS ',
     &                      'grid file ''', fname(indx)(1:lsrc), '''.',
     &                                               nf_strerror(ierr)
          endif
        else
          write(*,'(/1x,4A/12x,A/)')  '### ERROR: Cannot read bottom ',
     &                'topography from ''', fname(indx)(1:lsrc), '''.',
     &                                               nf_strerror(ierr)
        endif
      else
        write(*,'(1x,4A/12x,A)') '### ERROR: Cannot open grid file ''',
     &                  fname(indx)(1:lsrc), ''': ', nf_strerror(ierr)
      endif
      if (ierr /= nf_noerr) stop
      do indx=first_file, last_file
        lsrc=lenstr(fname(indx))
        ierr=nf_open(fname(indx), nf_nowrite, ncsrc)
        if (ierr == nf_noerr) then
          write(*,*) 'Open ''', fname(indx)(1:lsrc), '''...'
          call roms_check_dims(ncsrc, fname(indx), Lm,Mm,N)
          call  roms_find_dims(ncsrc, fname(indx), Lm,Mm,N)
          ierr=nf_inq(ncsrc, ndims, nvars, ngatts, unlimdimid)
          if (ierr == nf_noerr) then
            if (ndims > max_dims) then
              if (allocated(dname)) then
                deallocate(dname)     ; deallocate(tg_dimid)
                deallocate(dimv)      ; deallocate(Ld_in_use)
                deallocate(dsize)
              endif
              max_dims=ndims          ; allocate(dname(ndims))
              allocate(dsize(ndims))  ; allocate(Ld_in_use(ndims))
              allocate(dimv(ndims))   ; allocate(tg_dimid(ndims))
            endif
            tsize=0    ; id_xi_rho=-1  ; id_xi_u=-1  ; id_s_rho=-1
            id_time=-1 ; id_eta_rho=-1 ; id_eta_v=-1 ; id_s_w=-1
            do j=1,ndims
              Ld_in_use(j)=.false.
              ierr=nf_inq_dim(ncsrc, j, dname(j), dsize(j))
              if (ierr == nf_noerr) then
                ldim=lenstr(dname(j))
                if (ldim==6 .and. dname(j)(1:ldim) == 'xi_rho') then
                  id_xi_rho=j
                elseif (ldim==4 .and. dname(j)(1:ldim) == 'xi_u') then
                  id_xi_u=j
                elseif (ldim==7 .and. dname(j)(1:ldim)=='eta_rho') then
                  id_eta_rho=j
                elseif (ldim==5 .and. dname(j)(1:ldim) == 'eta_v') then
                  id_eta_v=j
                elseif (ldim==5 .and. dname(j)(1:ldim) == 's_rho') then
                  id_s_rho=j
                elseif (ldim==3 .and. dname(j)(1:ldim) == 's_w') then
                  id_s_w=j
                elseif (ldim==4  .and. dname(j)(1:ldim) == 'time') then
                  id_time=j
                  tsize=dsize(j)
                  dsize(j)=nf_unlimited
                endif
              else
                write(*,'(/1x,2A,I4/12x,3A/12x,A/)')  '### ERROR: ',
     &             'Cannot determine name and size for dimension #',
     &              j,   'in netCDF file ''',   fname(indx)(1:lsrc),
     &                                   '''.',   nf_strerror(ierr)
                stop
              endif
            enddo
            if (nvars > max_vars) then
              if (allocated(vname)) then
                deallocate(lmask)
                deallocate(dimids)
                deallocate(vdims)
                deallocate(vtype)
                deallocate(tg_var_grd_type)
                deallocate(var_grd_type)
                deallocate(tg_vname)
                deallocate(vname)
              endif
              max_vars=nvars
              allocate(vname(nvars))
              allocate(lmask(nvars))
              allocate(tg_vname(nvars+1))
              allocate(var_grd_type(nvars))
              allocate(tg_var_grd_type(nvars+1))
              allocate(vtype(nvars))
              allocate(vdims(nvars))
              allocate(dimids(ndims,nvars))
            endif
            do i=1,nvars
              vname(i)='                '
              var_grd_type(i)='    '
              ierr=nf_inq_var(ncsrc,  i, vname(i),  vtype(i),
     &                        vdims(i), dimids(1,i),  vatts)
              if (ierr == nf_noerr) then
                lstr=lenstr(vname(i))
                if (nsel_vars > 0) then
                  lmask(i)=.false.
                  do j=1,nsel_vars
                    if (vname(i) == sel_vars(j)) lmask(i)=.true.
                  enddo
                else
                  lmask(i)=.true.
                endif
                if (lmask(i)) then
                  do j=1,vdims(i)
                    if (dimids(j,i) == id_xi_rho) then
                      var_grd_type(i)(1:1)='R'
                    elseif (dimids(j,i) == id_xi_u) then
                      var_grd_type(i)(1:1)='U'
                    elseif (dimids(j,i) == id_eta_rho) then
                      var_grd_type(i)(2:2)='R'
                    elseif (dimids(j,i) == id_eta_v) then
                      var_grd_type(i)(2:2)='V'
                    elseif (dimids(j,i) == id_s_rho) then
                      var_grd_type(i)(3:3)='R'
                    elseif (dimids(j,i) == id_s_w) then
                      var_grd_type(i)(3:3)='W'
                    elseif (dimids(j,i) == id_time) then
                      var_grd_type(i)(4:4)='T'
                    endif
                  enddo
                  if (var_grd_type(i)(3:3) /= ' ' .or.
     &                              nsel_vars > 0) then
                    do j=1,vdims(i)
                      Ld_in_use(dimids(j,i))=.true.
                    enddo
                  endif
                endif
              else
                write(*,'(/1x,2A,I3/12x,3A/12x,A/)')  '### ERROR: ',
     &            'Cannot make general inquiry about variable ID =',
     &             i,   'in netCDF file ''',    fname(indx)(1:lsrc),
     &                                     ''',', nf_strerror(ierr)
                stop
              endif
            enddo
            j=1
            do i=1,lsrc
              if (fname(indx)(i:i) == '/') j=i+1
            enddo
            trg_fname='z_'/  /fname(indx)(j:lsrc)
            ltrg=lenstr(trg_fname) ; tg_nvars=0
            ierr=nf_create(trg_fname(1:ltrg), nf_64bit_offset, nctarg)
            if (ierr == nf_noerr) then
              write(*,'(/1x,3A/)')    'Created netCDF file ''',
     &                                 trg_fname(1:ltrg), '''.'
              j=0
              if (id_xi_rho > 0) then
                if (Ld_in_use(id_xi_rho)) then
                  j=j+1
                  ierr=nf_def_dim(nctarg,'xi_rho', dsize(id_xi_rho),
     &                                                 tg_dimid(j))
                  id_xi_r_tgt=tg_dimid(j)
                endif
              endif
              if (id_xi_u > 0 .and. ierr == nf_noerr) then
                if (Ld_in_use(id_xi_u)) then
                  j=j+1
                  ierr=nf_def_dim(nctarg, 'xi_u',   dsize(id_xi_u),
     &                                                tg_dimid(j))
                  id_xi_u_tgt=tg_dimid(j)
                endif
              endif
              if (id_eta_rho > 0 .and. ierr == nf_noerr) then
                if (Ld_in_use(id_eta_rho)) then
                  j=j+1
                  ierr=nf_def_dim(nctarg,'eta_rho',dsize(id_eta_rho),
     &                                                  tg_dimid(j))
                  id_et_r_tgt=tg_dimid(j)
                endif
              endif
              if (id_eta_v > 0 .and. ierr == nf_noerr) then
                if (Ld_in_use(id_eta_v)) then
                  j=j+1
                  ierr=nf_def_dim(nctarg, 'eta_v', dsize(id_eta_v),
     &                                                tg_dimid(j))
                  id_et_v_tgt=tg_dimid(j)
                endif
              endif
              if (nz > 1 .and. ierr == nf_noerr) then
                j=j+1
                ierr=nf_def_dim(nctarg, 'depth', nz, tg_dimid(j))
                id_nz_tgt=tg_dimid(j)
              endif
              if (tsize > 1 .and. ierr == nf_noerr) then
                j=j+1
                ierr=nf_def_dim(nctarg, 'time', nf_unlimited,
     &                                          tg_dimid(j))
                id_time_tgt=tg_dimid(j)
              endif
              if (ierr == nf_noerr) then
                ierr=nf_inq_att(ncsrc, nf_global, 'partition',
     &                                               i, lstr)
                if (ierr == nf_noerr) then
                  ierr=nf_copy_att(ncsrc, nf_global, 'partition',
     &                                        nctarg, nf_global)
                  if (ierr == nf_noerr) then
                    write(*,*) 'copied global attribute ',
     &                                    '''partition'''
                  endif
                endif
                ndim_tgt=j
                if (nz > 1) ierr=nf_def_var(nctarg, 'depth', nf_real,
     &                                         1, id_nz_tgt, varid)
                do i=1,nvars+1
                  tg_vname(i)='                '
                  tg_var_grd_type(i)='    '
                enddo
                do i=1,nvars
                  lvar=lenstr(vname(i))
                  if ( (var_grd_type(i)(3:3)==' ' .and. nsel_vars==0)
     &                     .and. ( var_grd_type(i)(1:1) /= ' ' .or.
     &                             var_grd_type(i)(2:2) /= ' ' )) then
                    var_grd_type(i)='    '
                  elseif (lmask(i)) then
                    j=0
                    if (var_grd_type(i)(1:1) == 'R') then
                      j=j+1 ; dimv(j)=id_xi_r_tgt
                    elseif (var_grd_type(i)(1:1) == 'U') then
                      j=j+1 ; dimv(j)=id_xi_u_tgt
                    elseif (var_grd_type(i)(1:1) /= ' ') then
                      write(*,'(/1x,A/)') '### ERROR: XI-grid type?'
                      stop
                    endif
                    if (var_grd_type(i)(2:2) == 'R') then
                      j=j+1 ;; dimv(j)=id_et_r_tgt
                    elseif (var_grd_type(i)(2:2) == 'V') then
                      j=j+1
                      dimv(j)=id_et_v_tgt
                    elseif (var_grd_type(i)(2:2) /= ' ') then
                      write(*,'(/1x,A/)') '### ERROR: ETA-grid type?'
                      stop
                    endif
                    if (var_grd_type(i)(3:3) /= ' ' .and. nz > 1) then
                      j=j+1 ; dimv(j)=id_nz_tgt
                    endif
                    if (var_grd_type(i)(4:4) /=' ' .and. tsize>1) then
                      j=j+1 ; dimv(j)=id_time_tgt
                    endif
                    ierr=-1
                    if (j < 2 .and. lvar > 3) then
                      if (vname(i)(lvar-3:lvar) == 'time') ierr=1
                    endif
                    if (ierr < 0 .and.  vtype(i) == nf_double)
     &                                        vtype(i)=nf_real
                    ierr=nf_def_var(nctarg,  vname(i)(1:lvar),
     &                              vtype(i), j, dimv, varid)
                    tg_nvars=varid
                    if (ierr == nf_noerr) then
                      tg_vname(varid)=vname(i)(1:lvar)
                      tg_var_grd_type(varid)=var_grd_type(i)
                      k=0
                      do while(ierr == nf_noerr)
                        k=k+1
                        ierr=nf_inq_attname(ncsrc, i, k, string)
                        if (ierr == nf_noerr) then
                          lstr=lenstr(string)
                          ierr=nf_copy_att(ncsrc, i, string(1:lstr),
     &                                               nctarg, varid)
                        endif
                      enddo
                      if (j > 1) then
                        ierr=nf_put_att_real(nctarg,varid,
     &                     '_FillValue', vtype(i), 1, FillValue)
                      endif
                    else
                      write(*,'(/1x,4A/12x,A/)')       '### ERROR: ',
     &                 'Cannot define variable ''', vname(i)(1:lvar),
     &                ''' in output netCDF file.', nf_strerror(ierr)
                    endif
                  endif
                enddo
                ierr=nf_enddef(nctarg)
                if (nz > 1) then
                  call put_var_by_name_real(nctarg, 'depth', z_lev)
                endif
              else
                write(*,'(/1x,4A/12x,A/)')      '### ERROR: Cannot ',
     &               'define dimension(s) in output netCDF file. ''',
     &                   trg_fname(1:ltrg), '''.', nf_strerror(ierr)
              endif
            else
              write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot create ',
     &                    'netCDF file ''', trg_fname(1:ltrg), '''.',
     &                                             nf_strerror(ierr)
            endif
          else
            write(*,'(/1x,4A/12x,A/)')     '### ERROR: Cannot make ',
     &                         'general inquiry into netCDF file ''',
     &                 fname(indx)(1:lsrc), '''.', nf_strerror(ierr)
          endif
        else
          write(*,'(1x,4A/12x,A)')  '### ERROR: Cannot open netCDF ',
     &             'file ''', fname(indx)(1:lsrc), ''' for reading.',
     &                                             nf_strerror(ierr)
        endif
        if (ierr /= nf_noerr) stop
        if (.not.allocated(z_r) ) then
          allocate(Cs_r(N),   z_r(max_horiz_pts*N))
          allocate(Cs_w(0:N), z_w(max_horiz_pts*(N+1)))
          net_alloc_size=net_alloc_size +8*(2*N+1)*(1+max_horiz_pts)
          allocate(bffr(max_horiz_pts*(N+1)))
          allocate(bffr_out(max_horiz_pts*(nz+1)))
          net_alloc_size=net_alloc_size + 4*(N+1 + nz+1)*max_horiz_pts
          call read_scoord(ncsrc, N, Cs_r,Cs_w,hc, VertCoordType)
          call set_depth  (Lm,Mm, N, Cs_r,Cs_w,hc, VertCoordType,
     &                                                h, z_r,z_w)
          write(*,'(/1x,2A,F8.2,1x,A/)')    'Aggregate allocated ',
     &              'memory:', 1.D-6*dble(net_alloc_size), 'MBytes'
        endif
        do rec=1,max(tsize,1)
          do ivar=1,tg_nvars
            lvar=lenstr(tg_vname(ivar))
            if ( lvar > 0 .and. ( tg_var_grd_type(ivar)(4:4) == 'T'
     &                                         .or. rec == 1 )) then
              if (tg_var_grd_type(ivar)(1:2) /= '  ') then
                if (tg_var_grd_type(ivar)(1:1) == 'R') then
                  imin=0
                elseif (tg_var_grd_type(ivar)(1:1) == 'U') then
                  imin=1
                else
                  write(*,'(/1x,4A/)')     '### ERROR: illegal ',
     &           'var_grd_type =''', tg_var_grd_type(ivar), '''.'
                  stop
                endif
                if (tg_var_grd_type(ivar)(2:2) == 'R') then
                  jmin=0
                elseif (tg_var_grd_type(ivar)(2:2) == 'V') then
                  jmin=1
                else
                  write(*,'(/1x,4A/)')     '### ERROR: illegal ',
     &           'var_grd_type =''', tg_var_grd_type(ivar), '''.'
                  stop
                endif
                if (tg_var_grd_type(ivar)(3:3) == 'W') then
                  kmin=0
                elseif (tg_var_grd_type(ivar)(3:3) == 'R') then
                  kmin=1
                elseif (nsel_vars == 0) then
                  write(*,'(/1x,4A/)')     '### ERROR: illegal ',
     &           'var_grd_type =''', tg_var_grd_type(ivar), '''.'
                  stop
                endif
                if ( tg_var_grd_type(ivar)(3:3) == 'W'  .or.
     &               tg_var_grd_type(ivar)(3:3) == 'R' ) then
                  call get_rec_by_name_real(ncsrc, fname(indx)(1:lsrc),
     &                                          tg_vname(ivar)(1:lvar),
     &                        Lm+2-imin,Mm+2-jmin,N+1-kmin, rec, bffr)
C$OMP PARALLEL SHARED(Lm,Mm,N, z_r,z_w, rmask,  nz, z_lev,
C$OMP&                        imin,jmin,kmin, bffr, bffr_out)
                  call sigma_to_z_thread(Lm,Mm,N, z_r, z_w,  rmask,
     &                                 nz, z_lev, imin,jmin,kmin,
     &                                 bffr, bffr_out, FillValue)
C$OMP END PARALLEL
                  k=nz ; if (k == 1) k=0
                  call put_rec_by_name_real(nctarg, trg_fname(1:ltrg),
     &                                         tg_vname(ivar)(1:lvar),
     &                        Lm+2-imin,Mm+2-jmin, k,  rec, bffr_out)
                elseif ( tg_var_grd_type(ivar)(3:3) == ' ' .and.
     &                                       nsel_vars > 0 ) then
                  call get_rec_by_name_real(ncsrc, fname(indx)(1:lsrc),
     &                                          tg_vname(ivar)(1:lvar),
     &                              Lm+2-imin,Mm+2-jmin, 1, rec, bffr)
                  call put_rec_by_name_real(nctarg, trg_fname(1:ltrg),
     &                                          tg_vname(ivar)(1:lvar),
     &                             Lm+2-imin,Mm+2-jmin, 1,  rec, bffr)
                endif
              else
                if (tg_vname(ivar)(1:lvar) /= 'spherical') then
                  call read_by_name(ncsrc,   tg_vname(ivar)(1:lvar),
     &                                                   rec, bffr)
                  call write_by_name(nctarg, tg_vname(ivar)(1:lvar),
     &                                                   rec, bffr)
                endif
              endif
            endif
          enddo
          if (mod(rec,16) == 0) then
            write(*,*) 'synchronizing rec =', rec
            ierr=nf_sync(nctarg)
          endif
        enddo
        ierr=nf_close(ncsrc)
        ierr=nf_close(nctarg)
      enddo
      call cpu_time(run_time) ; run_time=run_time-tstart
      write(*,'(/3x,A,F11.2,1x,A)') 'cpu time =', run_time, 'sec'
      call system_clock (iclk_end, iclk_rate, iclk_max)
      if (iclk_rate > 0) then
        write(*,'(/1x,A,F20.3,1x,A/)')  'Elapsed wall-clock time =',
     &           float(iclk_end-iclk_start)/float(iclk_rate), 'sec'
      endif
      write(*,'(1x,A,F16.2,1x,2A,F15.2,1x,A)')   'Total data read:',
     &                                 sz_read_acc*1.D-6, 'MBytes,',
     &                   ' written:', sz_write_acc*1.D-6, 'MBytes.'
      end
