      program partit
      use netcdf
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
      integer(kind=4) :: deflate_level=1
      logical :: shuffle=.true.
      character(len=128) src_name, string
      integer(kind=4) nargs, nnodes, NP_XI,  xi_rho,  id_xi_rho,  
     &                             id_xi_psi,
     &        arg,   node,   NP_ETA, eta_rho, id_eta_rho, id_eta_psi,
     &        ncsrc, ndims,  ngatts, xi_u,    id_xi_u,    id_xi_v,
     &        nvars, varatts,        eta_v,   id_eta_v,   id_eta_u,
     &        i,j,k, ierr,   tsize,  rec,     unlimdimid, ierr_all,
     &        lstr,  lsrc,   lfnm,   lvar,    lenstr,     chk_alloc
      character(len=128), allocatable, dimension(:) :: fname
      character(len=32),  allocatable, dimension(:) :: dimname, vname
      integer(kind=4), allocatable, dimension(:,:) :: vid, dimids
      integer(kind=4), allocatable, dimension(:) :: ncid,  dimid,  
     &                              dimsize,
     &                     xi_start, xi_size,   eta_start, eta_size,
     &              vtype, vdims,    part_type, start, count, start1
     &                                                      , count1
      logical, allocatable, dimension(:) ::  series, western_edge,
     &                 eastern_edge, southern_edge, northern_edge
      integer(kind=4) :: maxdims=0, maxvars=0
      integer(kind=4) size, size_cbfr,  size_ibfr1, size_ibfr2,
     &              size_ibffr, size_bffr4, size_bffr8
      integer(kind=4) :: alloc_cbfr=0,  alloc_ibfr1=0, alloc_ibfr2=0,
     &           alloc_ibffr=0, alloc_bffr4=0, alloc_bffr8=0
      character(len=1), allocatable, dimension(:) :: cbfr
      integer(kind=1),  allocatable, dimension(:) :: ibfr1
      integer(kind=2),  allocatable, dimension(:) :: ibfr2
      integer(kind=4),  allocatable, dimension(:) :: ibffr
      real(kind=4), allocatable, dimension(:) :: bffr4
      real(kind=8), allocatable, dimension(:) :: bffr8
      integer(kind=4) size1, size_cbf_in,  size_ibf1_in, size_ibf2_in,
     &               size_ibff_in, size_bff4_in, size_bff8_in
      integer(kind=4) :: allc_cbf_in=0,  allc_ibf1_in=0, allc_ibf2_in=0,
     &           allc_ibff_in=0, allc_bff4_in=0, allc_bff8_in=0
      character(len=1), allocatable, dimension(:) :: cbf_in
      integer(kind=1),  allocatable, dimension(:) :: ibf1_in
      integer(kind=2),  allocatable, dimension(:) :: ibf2_in
      integer(kind=4),  allocatable, dimension(:) :: ibff_in
      real(kind=4), allocatable, dimension(:) :: bff4_in
      real(kind=8), allocatable, dimension(:) :: bff8_in
      real(kind=4) tstart, run_time
      real(kind=8) InitTime, ReadSize, ReadTime, WrtSize,  WrtTime,
     &                                           SyncTime, GrayTime
      integer(kind=4) iclk(2), nclk, clk_rate, clk_max, vtypsize
      integer(kind=8) inc_clk
      integer(kind=8) :: net_clk=0, net_read_size=0, net_wrt_size=0,
     &              net_init_clk=0,  net_read_clk=0,  net_wrt_clk=0,
     &                               net_sync_clk=0, net_gray_clk=0
      real(kind=8)  AssmTime
      integer(kind=8) :: net_assm_clk=0
      call cpu_time(tstart)
      nclk=1
      call system_clock(iclk(nclk), clk_rate, clk_max)
      nargs=iargc()
      if (nargs < 3) then
        write(*,'(/1x,A/30x,A/)')      'Usage of partit should be:',
     &                     'partit NP_X NP_E file1.nc file2.nc ...'
        stop
      endif
      call getarg(1,string) ; lstr=lenstr(string) ; NP_XI=0
      do i=1,lstr
        j=ichar(string(i:i))-48
        if (j < 0 .or. j > 9) then
          write(*,'(/8x,3A/)')     '### ERROR: First argument ',
     &           string(1:lstr), ', must be an integer number.'
          stop
        else
          NP_XI=10*NP_XI+j
        endif
      enddo
      call getarg(2,string) ; lstr=lenstr(string) ; NP_ETA=0
      do i=1,lstr
        j=ichar(string(i:i))-48
        if (j < 0 .or. j > 9) then
          write(*,'(/8x,3A/)')    '### ERROR: Second argument ',
     &           string(1:lstr), ', must be an integer number.'
          stop
        else
          NP_ETA=10*NP_ETA+j
        endif
      enddo
      nnodes=NP_XI*NP_ETA
      write(*,'(/1x,2(4x,A,I3)/)') 'NP_XI =',NP_XI, 'NP_ETA =',NP_ETA
      allocate( ncid(0:nnodes-1),          fname(0:nnodes-1),
     &          xi_size(0:nnodes-1),       eta_size(0:nnodes-1),
     &          xi_start(0:nnodes-1),      eta_start(0:nnodes-1),
     &          western_edge(0:nnodes-1),  eastern_edge(0:nnodes-1),
     &          southern_edge(0:nnodes-1), northern_edge(0:nnodes-1) )
      do arg=3,nargs
        ierr_all=0
        call getarg(arg,src_name)
        lsrc=lenstr(src_name)
        ierr=nf_open(src_name(1:lsrc), nf_nowrite, ncsrc)
        if (ierr == nf_noerr) then
          ierr=nf_inq_att (ncsrc, nf_global, 'partition', i,j)
          if (ierr == nf_noerr) then
            write(*,'(/1x,4A/14x,2A/)')'### WARNING: netCDF file ''',
     &        src_name(1:lsrc),  ''' is already ', 'a partial file ',
     &       'and cannot be partitioned any further ==> ignoring it.'
            goto 97
          endif
        else
          write(*,'(/1x,4A/14x,A)')      '### WARNING: Cannot open ',
     &    'netCDF file ''',src_name(1:lsrc),'''.', nf_strerror(ierr)
          goto 97
        endif
        write(*,'(1x,3A)') 'Processing netCDF file ''',
     &                        src_name(1:lsrc), '''...'
        ierr=nf_inq(ncsrc, ndims, nvars, ngatts, unlimdimid)
        if (ierr == nf_noerr) then
          if (ndims > maxdims .or. nvars > maxvars) then
            if (allocated(dimname)) then
              deallocate( dimname, dimid, dimsize, start, count,
     &                                                       start1,
     &                                                       count1,
     &        vname, vtype, vdims, part_type, series, dimids, vid )
            endif
            maxdims=max(ndims,4) ; maxvars=nvars
            allocate( dimname(maxdims), dimid(maxdims),
     &                dimsize(maxdims), start(maxdims), count(maxdims),
     &                                 start1(maxdims)
     &                                               , count1(maxdims)
     &                                                               )
            allocate(vname(maxvars),  vtype(maxvars),  vdims(maxvars),
     &                            part_type(maxvars), series(maxvars))
            allocate(dimids(maxdims,maxvars), vid(maxvars,0:nnodes-1))
            write(*,'(1x,2A,I3,1x,A,I3)')   'allocated book-keeping ',
     &           'arrays, maxdims =', maxdims, 'maxvars =', maxvars
          endif
        else
          write(*,'(/1x,4A/12x,A/)')       '### ERROR: Cannot make ',
     &       'general inquiry into netCDF file ''', src_name(1:lsrc),
     &                                     '''.', nf_strerror(ierr)
        endif
        if (ierr /= nf_noerr) stop
        tsize=-1
        do j=1,ndims
          ierr=nf_inq_dim(ncsrc, j, dimname(j), dimsize(j))
          if (ierr == nf_noerr) then
            if (j == unlimdimid) then
              tsize=dimsize(j)
              dimsize(j)=nf_unlimited
            endif
          else
            write(*,'(/1x,2A,I4/12x,3A/12x,A/)')       '### ERROR: ',
     &           'Cannot determine name and size for dimension #', j,
     &                  'in netCDF file ''', src_name(1:lsrc), '''.',
     &                                             nf_strerror(ierr)
             goto 97
          endif
        enddo
        xi_rho=0   ; id_xi_rho=0
        xi_u=0     ; id_xi_u=0
        eta_rho=0  ; id_eta_rho=0
        eta_v=0    ; id_eta_v=0
        id_xi_v=0  ; id_xi_psi=0
        id_eta_u=0 ; id_eta_psi=0
        do j=1,ndims
          lvar=lenstr(dimname(j))
          if (lvar==6 .and. dimname(j)(1:lvar)=='xi_rho') then
            id_xi_rho=j  ; xi_rho=dimsize(j)
          elseif (lvar==4 .and. dimname(j)(1:lvar)=='xi_u') then
            id_xi_u=j    ; xi_u=dimsize(j)
          elseif (lvar==7 .and. dimname(j)(1:lvar)=='eta_rho') then
            id_eta_rho=j ; eta_rho=dimsize(j)
          elseif (lvar==5 .and. dimname(j)(1:lvar)=='eta_v') then
            id_eta_v=j   ; eta_v=dimsize(j)
          elseif (lvar==6 .and. dimname(j)(1:lvar)=='xi_psi') then
            id_xi_psi=j
          elseif (lvar==4 .and. dimname(j)(1:lvar)=='xi_v') then
            id_xi_v=j
          elseif (lvar==7 .and. dimname(j)(1:lvar)=='eta_psi') then
            id_eta_psi=j
          elseif (lvar==5 .and. dimname(j)(1:lvar)=='eta_u') then
            id_eta_u=j
          endif
        enddo
        if (id_xi_rho /= 0 .and.  id_xi_u == 0) then
          xi_u=xi_rho-1
        elseif (id_xi_rho == 0 .and.  id_xi_u /= 0) then
          xi_rho=xi_u+1
        endif
        if (id_eta_rho /= 0 .and. id_eta_v == 0) then
          eta_v=eta_rho-1
        elseif (id_eta_rho == 0 .and. id_eta_v == 0) then
          eta_rho=eta_v+1
        endif
        if (xi_rho == 0  .or. xi_u == 0  .or.
     &      eta_rho == 0 .or. eta_v == 0) then
          write(*,'(/8x,2A/15x,3A/)')  '### ERROR: not all ',
     &       'partitionable dimensions are found in netCDF ',
     &                    'file ''', src_name(1:lsrc), '''.'
          goto 97
        endif
        call mpi_setup (NP_XI,NP_ETA, xi_rho,eta_rho,
     &           xi_start,xi_size, eta_start,eta_size,
     &                    western_edge, eastern_edge,
     &                   southern_edge, northern_edge)
        j=1
        do i=1,lsrc
          if (src_name(i:i)=='/') j=i+1
        enddo
        fname(nnodes-1)=src_name(j:lsrc)
        do node=0,nnodes-1
          fname(node)=fname(nnodes-1) ; lfnm=lenstr(fname(node))
          call insert_node(fname(node), lfnm, node, nnodes, ierr)
          if (ierr == 0) then
            ierr=nf90_create(fname(node)(1:lfnm), NF90_NETCDF4, 
     &                            ncid(node))
            if (ierr == nf_noerr) then
              if ( node < 16 .or. ( nnodes > 16 .and.
     &                              node == nnodes-1 )) then
                write(*,'(4x,3A)') 'Created partitioned file ''',
     &                                fname(node)(1:lfnm), '''.'
              elseif (nnodes > 16 .and. node < 18) then
                write(*,'(16x,A)') '.................................'
              endif
              ierr=nf_set_fill(ncid(node), nf_nofill, i)
              if (ierr /= nf_noerr) then
                 write(*,'(/8x,A,1x,3A/8x,A)')   '### ERROR: cannot ',
     &                      'set "nf_nofill" mode for netCDF file ''',
     &                  fname(node)(1:lfnm), '''.', nf_strerror(ierr)
              endif
            else
              write(*,'(/1x,4A/8x,A)' )   '### ERROR: cannot create ',
     &                          'netCDF file ''', fname(node)(1:lfnm),
     &                                      '''.',  nf_strerror(ierr)
            endif
          endif
          if (ierr /=  nf_noerr) stop
          do j=1,ndims
            if (j == id_xi_rho) then
              size=xi_size(node)
              if (western_edge(node)) size=size+1
              if (eastern_edge(node)) size=size+1
            elseif (j == id_xi_u) then
              size=xi_size(node)
              if (eastern_edge(node)) size=size+1
            elseif (j == id_eta_rho) then
              size=eta_size(node)
              if (southern_edge(node)) size=size+1
              if (northern_edge(node)) size=size+1
            elseif (j == id_eta_v) then
              size=eta_size(node)
              if (northern_edge(node)) size=size+1
            else
              size=dimsize(j)
            endif
            dimid(j)=0 ; lvar=lenstr(dimname(j))
            if (j /= unlimdimid .and. size == 0) then
              if (node == 0) write(*,'(4x,4A)')   'Suppressing ',
     &       'zero-size dimension ''', dimname(j)(1:lvar), '''.'
            elseif  (j == unlimdimid .and. tsize == 0) then
              if (node == 0) write(*,'(4x,4A)')   'Suppressing ',
     &                        'zero-size unlimited dimension ''',
     &                                 dimname(j)(1:lvar), '''.'
            elseif (j == id_xi_psi) then
              if (node == 0) write(*,'(4x,2A)')   'Suppressing ',
     &                          'obsolete dimension ''xi_psi''.'
            elseif (j == id_xi_v) then
              if (node == 0) write(*,'(4x,2A)')   'Suppressing ',
     &                            'obsolete dimension ''xi_v''.'
            elseif (j == id_eta_psi) then
              if (node == 0) write(*,'(4x,2A)')   'Suppressing ',
     &                         'obsolete dimension ''eta_psi''.'
            elseif (j == id_eta_u)  then
              if (node == 0) write(*,'(4x,2A)')   'Suppressing ',
     &                            'obsolete dimension ''xi_u''.'
            else
              ierr=nf_def_dim(ncid(node), dimname(j)(1:lvar),
     &                                        size, dimid(j))
              if (ierr /= nf_noerr) then
                ierr_all=ierr_all+1
                write(*,'(/1x,4A,I8,A,I4/12x,A)')    '### ERROR: ',
     &            'Cannot define dimension ''', dimname(j)(1:lvar),
     &                       ''' of size =', size,',  node=', node,
     &                                           nf_strerror(ierr)
              endif
            endif
          enddo
          start(1)=node    ; start(3)=xi_start(node)
          start(2)=nnodes  ; start(4)=eta_start(node)
          ierr=nf_put_att_int(ncid(node), nf_global, 'partition',
     &                                       nf_int,   4, start)
        enddo
        if (ierr_all /= 0) stop
        do i=1,ngatts
          ierr=nf_inq_attname (ncsrc, nf_global, i, string)
          if (ierr == nf_noerr) then
            lvar=lenstr(string)
            do node=0,nnodes-1
              ierr=nf_copy_att (ncsrc, nf_global, string(1:lvar),
     &                                     ncid(node), nf_global)
              if (ierr. ne. nf_noerr) then
                ierr_all=ierr_all+1
                write(*,'(/1x,7A/12x,A)')  '### ERROR: Cannot copy ',
     &             'global attribute ''', string(1:lvar), ''' into ',
     &             'netCDF file ''',   fname(node)(1:lfnm),   '''.',
     &                                             nf_strerror(ierr)
                goto 97
              endif
            enddo
          else
            write(*,'(/1x,2A,I3,1x,4A/12x,A/)') '### ERROR: Cannot ',
     &       'determine name of global attribute #', i, 'in netCDF ',
     &       'file ''', src_name(1:lsrc), '''.',  nf_strerror(ierr)
            goto 97
          endif
        enddo
        do i=1,nvars
          ierr=nf_inq_var(ncsrc,  i,    vname(i),  vtype(i),
     &                    vdims(i), dimids(1,i),   varatts)
          if (ierr /= nf_noerr) then
             write(*,'(/1x,2A,I3/12x,3A/12x,A/)')     '### ERROR: ',
     &            'Cannot make general inquiry about variable ID =',
     &             i, 'in netCDF file ''',  src_name(1:lsrc), ''',',
     &                                            nf_strerror(ierr)
            goto 97
          endif
          do j=1,vdims(i)
            if (dimids(j,i) == id_xi_psi) then
              dimids(j,i)=id_xi_u
            elseif (dimids(j,i) == id_xi_v) then
              dimids(j,i)=id_xi_rho
            elseif (dimids(j,i) == id_eta_psi) then
              dimids(j,i)=id_eta_v
            elseif (dimids(j,i) == id_eta_u) then
              dimids(j,i)=id_eta_rho
            endif
          enddo
          series(i)=.false.
          part_type(i)=0
          do j=1,vdims(i)
            if ( dimids(j,i) == id_xi_rho .or.
     &           dimids(j,i) == id_xi_u  ) then
              part_type(i)=part_type(i)+1
            elseif ( dimids(j,i) == id_eta_rho .or.
     &               dimids(j,i) == id_eta_v ) then
              part_type(i)=part_type(i)+2
            elseif (dimids(j,i)==unlimdimid) then
              series(i)=.true.
            endif
          enddo
          if (tsize > 0 .or. .not.series(i)) then
            do j=1,vdims(i)
              do k=1,ndims
                if (dimids(j,i) == k) start(j)=dimid(k)
              enddo
            enddo
            lvar=lenstr(vname(i))
            do node=0,nnodes-1
              vid(i,node)=0
              if (lvar > 5) then
                if (vname(i)(lvar-4:lvar) == '_west' .and.
     &              .not.western_edge(node)) vid(i,node)=-1000
                if (vname(i)(lvar-4:lvar) == '_east'  .and.
     &              .not.eastern_edge(node)) vid(i,node)=-1000
              endif
              if (lvar > 6) then
                if (vname(i)(lvar-5:lvar) == '_south' .and.
     &              .not.southern_edge(node)) vid(i,node)=-1000
                if (vname(i)(lvar-5:lvar) == '_north'  .and.
     &              .not.northern_edge(node)) vid(i,node)=-1000
              endif
              if (vid(i,node) == 0) then
                ierr=nf90_def_var (ncid(node), 
     &                     vname(i)(1:lvar),vtype(i),
     &                             start(1:vdims(i)), vid(i,node),
     &                             deflate_level=deflate_level, 
     &                          shuffle=shuffle)
                if (ierr /= nf_noerr) then
                  write(*,'(/1x,4A,I4/12x,A/)') '### ERROR: Cannot ',
     &                    'create variable ''',   vname(i)(1:lvar),
     &                  ''', node =',   node,    nf_strerror(ierr)
                endif
              endif
            enddo
            do j=1,varatts
              ierr=nf_inq_attname(ncsrc, i, j, string)
              lvar=lenstr(string)
              do node=0,nnodes-1
                if (vid(i,node)>0) then
                  ierr=nf_copy_att(ncsrc, i, string(1:lvar),
     &                              ncid(node), vid(i,node))
                  if (ierr /= nf_noerr) then
                    write(*,'(/1x,4A,I4/12x,A/)') '### ERROR: ',
     &               'Cannot copy attribute ''', string(1:lvar),
     &               ''' for variable ''',   vname(i)(1:lvar),
     &               ''', node =',   node,    nf_strerror(ierr)
                  endif
                endif
              enddo
            enddo
          endif
        enddo
        do node=0,nnodes-1
          ierr=nf_enddef(ncid(node))
          if (ierr /= nf_noerr) then
            write(*,'(/1x,5A/)')  '### ERROR: Cannot switch partial ',
     &                          'netCDF file ''', fname(node)(1:lfnm),
     &                      ''' into input mode, ', nf_strerror(ierr)
            ierr_all=ierr_all+1
          endif
        enddo
        if (ierr_all /= 0) stop
  1     size_cbfr=0     ; size_ibfr1=0    ; size_ibfr2=0
        size_ibffr=0    ; size_bffr4=0    ; size_bffr8=0
        size_cbf_in=0   ; size_ibf1_in=0  ; size_ibf2_in=0
        size_ibff_in=0  ; size_bff4_in=0  ; size_bff8_in=0
        do i=1,nvars
          do node=0,nnodes-1
            if (vid(i,node) > 0) then
              size=1
              do j=1,vdims(i)
                if ( dimids(j,i) == id_xi_rho .or.
     &               dimids(j,i) == id_xi_u ) then
                  count(j)=xi_size(node)
                  if ( western_edge(node) .and.
     &                 dimids(j,i) == id_xi_rho ) then
                    count(j)=count(j)+1
                  endif
                  if (eastern_edge(node)) then
                    count(j)=count(j)+1
                  endif
                elseif ( dimids(j,i) == id_eta_rho .or.
     &                   dimids(j,i) == id_eta_v ) then
                  count(j)=eta_size(node)
                  if ( southern_edge(node) .and.
     &                 dimids(j,i) == id_eta_rho) then
                    count(j)=count(j)+1
                  endif
                  if (northern_edge(node)) then
                    count(j)=count(j)+1
                  endif
                elseif (dimids(j,i) == unlimdimid) then
                  count(j)=1
                else
                  count(j)=dimsize(dimids(j,i))
                endif
                size=size*count(j)
              enddo
              if (vtype(i) == nf_char) then
                size_cbfr=max(size_cbfr,size)
              elseif (vtype(i) == nf_byte) then
                size_ibfr1=max(size_ibfr1,size)
              elseif (vtype(i) == nf_short) then
                size_ibfr2=max(size_ibfr2,size)
              elseif (vtype(i) == nf_int) then
                size_ibffr=max(size_ibffr,size)
              elseif (vtype(i) == nf_float) then
                size_bffr4=max(size_bffr4,size)
              elseif (vtype(i) == nf_double) then
                size_bffr8=max(size_bffr8,size)
              else
                write(*,'(/8x,3A/)')  '### ERROR: Variable ''',
     &              vname(i)(1:lvar), ''' is of unknown type.'
                stop
             endif
            endif
          enddo
          if (part_type(i) > 0) then
            size1=1
            do j=1,vdims(i)
              if (dimids(j,i) /= unlimdimid) then
                size1=size1*dimsize(dimids(j,i))
              endif
            enddo
            if (vtype(i) == nf_char) then
              size_cbf_in=max(size_cbf_in,size1)
            elseif (vtype(i) == nf_byte) then
              size_ibf1_in=max(size_ibf1_in,size1)
            elseif (vtype(i) == nf_short) then
              size_ibf2_in=max(size_ibf2_in,size1)
            elseif (vtype(i) == nf_int) then
              size_ibff_in=max(size_ibff_in,size1)
            elseif (vtype(i) == nf_float) then
              size_bff4_in=max(size_bff4_in,size1)
            elseif (vtype(i) == nf_double) then
              size_bff8_in=max(size_bff8_in,size1)
            endif
          endif
        enddo
        if (size_cbfr > alloc_cbfr) then
          if (allocated(cbfr)) deallocate(cbfr)
          allocate(cbfr(size_cbfr)) ;  alloc_cbfr=size_cbfr
          write(*,*) 'allocated nf_char workspace', size_cbfr
        endif
        if (size_ibfr1 > alloc_ibfr1) then
          if (allocated(ibfr1)) deallocate(ibfr1)
          allocate(ibfr1(size_ibfr1)) ; alloc_ibfr1=size_ibfr1
          write(*,*) 'allocated nf_byte workspace', size_ibfr1
        endif
        if (size_ibfr2 > alloc_ibfr2) then
          if (allocated(ibfr2)) deallocate(ibfr2)
          allocate(ibfr2(size_ibfr2)) ; alloc_ibfr2=size_ibfr2
          write(*,*) 'allocated nf_int2 workspace', size_ibfr2
        endif
        if (size_ibffr > alloc_ibffr) then
          if (allocated(ibffr)) deallocate(ibffr)
          allocate(ibffr(size_ibffr)) ; alloc_ibffr=size_ibffr
          write(*,*)  'allocated nf_int workspace', size_ibffr
        endif
        if (size_bffr4 > alloc_bffr4) then
          if (allocated(bffr4)) deallocate(bffr4)
          allocate(bffr4(size_bffr4)) ; alloc_bffr4=size_bffr4
          write(*,*) 'allocated nf_real workspace', size_bffr4
        endif
        if (size_bffr8 > alloc_bffr8) then
          if (allocated(bffr8)) deallocate(bffr8)
          allocate(bffr8(size_bffr8)) ; alloc_bffr8=size_bffr8
          write(*,*)'allocated nf_double workspace',size_bffr8
        endif
        if (size_cbf_in > allc_cbf_in) then
          if (allocated(cbf_in)) deallocate(cbf_in)
          allocate(cbf_in(size_cbf_in)) ; allc_cbf_in=size_cbf_in
          write(*,*) 'allocated nf_char workspace', size_cbf_in
        endif
        if (size_ibf1_in > allc_ibf1_in) then
          if (allocated(ibf1_in)) deallocate(ibf1_in)
          allocate(ibf1_in(size_ibf1_in)) ; allc_ibf1_in=size_ibf1_in
          write(*,*) 'allocated nf_byte workspace',size_ibf1_in
        endif
        if (size_ibf2_in > allc_ibf2_in) then
          if (allocated(ibf2_in)) deallocate(ibf2_in)
          allocate(ibf2_in(size_ibf2_in)) ; allc_ibf2_in=size_ibf2_in
          write(*,*) 'allocated nf_int2 workspace', size_ibf2_in
        endif
        if (size_ibff_in > allc_ibff_in) then
          if (allocated(ibff_in)) deallocate(ibff_in)
          allocate(ibff_in(size_ibff_in)) ; allc_ibff_in=size_ibff_in
          write(*,*) 'allocated nf_int workspace',size_ibff_in
        endif
        if (size_bff4_in > allc_bff4_in) then
          if (allocated(bff4_in)) deallocate(bff4_in)
          allocate(bff4_in(size_bff4_in)) ; allc_bff4_in=size_bff4_in
          write(*,*) 'allocated nf_real workspace',size_bff4_in
        endif
        if (size_bff8_in > allc_bff8_in) then
          if (allocated(bff8_in)) deallocate(bff8_in)
          allocate(bff8_in(size_bff8_in)) ; allc_bff8_in=size_bff8_in
            write(*,*) 'allocated nf_double workspace', size_bff8_in
        endif
        nclk=3-nclk
        call system_clock(iclk(nclk), clk_rate,clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        if (inc_clk < 0) inc_clk=inc_clk+clk_max
        net_clk=net_clk+inc_clk
        net_init_clk=net_init_clk+inc_clk
        do rec=1,max(tsize,1)
          if (tsize > 1) then
            nclk=3-nclk
            call system_clock(iclk(nclk), clk_rate,clk_max)
            inc_clk=iclk(nclk)-iclk(3-nclk)
            if (inc_clk < 0) inc_clk=inc_clk+clk_max
            net_clk=net_clk+inc_clk
            net_gray_clk=net_gray_clk+inc_clk
            write(*,'(F8.1,1x,2(A,I8,1x),A)')
     &                         dble(net_clk)/dble(clk_rate),
     &       'rec', rec, 'out of', tsize,  '...'
          endif
          do i=1,nvars
            lvar=lenstr(vname(i))
            if (rec==1 .or. series(i)) then
              if (vtype(i) == nf_char .or.
     &            vtype(i) == nf_byte) then
                vtypsize=1
              elseif (vtype(i) == nf_short) then
                vtypsize=2
              elseif (vtype(i) == nf_int .or.
     &                vtype(i) == nf_float) then
                vtypsize=4
              elseif (vtype(i) == nf_double) then
                vtypsize=8
              endif
              if (part_type(i) == 0 .and. .not.series(i)) then
                write(*,'(21x,A,I4,1x,A,I4,1x,3A)') 'var', i, 'out of',
     &           nvars,'non-part.  non-rec. ''',vname(i)(1:lvar),'''.'
                size=1
                do j=1,vdims(i)
                  size=size * dimsize(dimids(j,i))
                enddo
                if (vtype(i) == nf_char) then
                  ierr=chk_alloc(size, alloc_cbfr, 'nf_char')
                elseif (vtype(i) == nf_byte) then
                  ierr=chk_alloc(size, alloc_ibfr1, 'nf_byte')
                elseif (vtype(i) == nf_short) then
                  ierr=chk_alloc(size, alloc_ibfr2, 'nf_short')
                elseif (vtype(i) == nf_int) then
                  ierr=chk_alloc(size, alloc_ibffr,   'nf_int')
                elseif (vtype(i) == nf_float) then
                  ierr=chk_alloc(size, alloc_bffr4, 'nf_float')
                elseif (vtype(i) == nf_double) then
                  ierr=chk_alloc(size,alloc_bffr8, 'nf_double')
                else
                  write(*,'(/1x,4A/)')  '### ERROR: Unknown type ',
     &                 'of variable ''', vname(i)(1:lvar), '''.'
                  stop
                endif
                if (ierr /= 0) stop
                if (vtype(i) == nf_char) then
                  ierr=nf_get_var_text (ncsrc, i, cbfr)
                elseif (vtype(i) == nf_byte) then
                  ierr=nf_get_var_int1   (ncsrc, i, ibfr1)
                elseif (vtype(i) == nf_short) then
                  ierr=nf_get_var_int2   (ncsrc, i, ibfr2)
                elseif (vtype(i) == nf_int) then
                  ierr=nf_get_var_int    (ncsrc, i, ibffr)
                elseif (vtype(i) == nf_float) then
                  ierr=nf_get_var_real   (ncsrc, i, bffr4)
                elseif (vtype(i) == nf_double) then
                  ierr=nf_get_var_double (ncsrc, i, bffr8)
                endif
                if (ierr /= nf_noerr) then
                  write(*,'(/1x,6A/12x,A/)')      '### ERROR: ',
     &              'Cannot read variable ''', vname(i)(1:lvar),
     &              ''' from netCDF file ''',  src_name(1:lsrc),
     &                                 '''.', nf_strerror(ierr)
                  goto 97
                else
                  nclk=3-nclk
                  call system_clock(iclk(nclk),clk_rate,clk_max)
                  inc_clk=iclk(nclk)-iclk(3-nclk)
                  if (inc_clk < 0) inc_clk=inc_clk+clk_max
                  net_clk=net_clk+inc_clk
                  net_read_clk=net_read_clk+inc_clk
                  net_read_size=net_read_size+size*vtypsize
                endif
                do node=0,nnodes-1
                  if (vid(i,node) > 0) then
                    if (vtype(i) == nf_char) then
                      ierr=nf_put_var_text (ncid(node), vid(i,node),
     &                                                        cbfr)
                    elseif (vtype(i) == nf_byte) then
                      ierr=nf_put_var_int1 (ncid(node), vid(i,node),
     &                                                       ibfr1)
                    elseif (vtype(i) == nf_short) then
                      ierr=nf_put_var_int2 (ncid(node), vid(i,node),
     &                                                       ibfr2)
                    elseif (vtype(i) == nf_int) then
                      ierr=nf_put_var_int  (ncid(node), vid(i,node),
     &                                                       ibffr)
                    elseif (vtype(i) == nf_float) then
                      ierr=nf_put_var_real (ncid(node), vid(i,node),
     &                                                       bffr4)
                    elseif (vtype(i) == nf_double) then
                      ierr=nf_put_var_double(ncid(node),vid(i,node),
     &                                                       bffr8)
                    endif
                    if (ierr /= nf_noerr) then
                      write(*,'(/1x,6A/12x,A/)')       '### ERROR:',
     &                'Cannot write variable ''',  vname(i)(1:lvar),
     &                ''' into netCDF file ''', fname(node)(1:lfnm),
     &                                     '''.', nf_strerror(ierr)
                      goto 97
                    else
                      nclk=3-nclk
                      call system_clock(iclk(nclk), clk_rate,clk_max)
                      inc_clk=iclk(nclk)-iclk(3-nclk)
                      if (inc_clk < 0) inc_clk=inc_clk+clk_max
                      net_clk=net_clk+inc_clk
                      net_wrt_clk=net_wrt_clk+inc_clk
                      net_wrt_size=net_wrt_size+size*vtypsize
                    endif
                  endif
                enddo
              elseif (part_type(i) == 0) then
                if (rec < 5) then
                  write(*,'(21x,A,I4,1x,A,I4,1x,3A)') 'var',i,'out of',
     &                 nvars, 'non-part.  ''', vname(i)(1:lvar), '''.'
                elseif (rec < 6) then
                   write(*,'(24x,A)') '.........................'
                endif
                size=1
                do j=1,vdims(i)
                  if (dimids(j,i)==unlimdimid) then
                    start(j)=rec ; count(j)=1
                  else
                    start(j)=1 ; count(j)=dimsize(dimids(j,i))
                  endif
                  size=size*count(j)
                enddo
                if (vtype(i) == nf_char) then
                  ierr=chk_alloc(size, alloc_cbfr,  'nf_char')
                elseif (vtype(i) == nf_byte) then
                  ierr=chk_alloc(size, alloc_ibfr1,  'nf_byte')
                elseif (vtype(i) == nf_short) then
                  ierr=chk_alloc(size, alloc_ibfr2,  'nf_short')
                elseif (vtype(i) == nf_int) then
                  ierr=chk_alloc(size, alloc_ibffr,    'nf_int')
                elseif (vtype(i) == nf_float) then
                  ierr=chk_alloc(size, alloc_bffr4,  'nf_float')
                elseif (vtype(i) == nf_double) then
                  ierr=chk_alloc(size, alloc_bffr8, 'nf_double')
                else
                  write(*,'(/1x,4A/)')   '### ERROR: Unknown ',
     &            'type of variable ''',vname(i)(1:lvar),'''.'
                  stop
                endif
                if (ierr /= 0) stop
                if (vtype(i) == nf_char) then
                  ierr=nf_get_vara_text (ncsrc, i, start,count, cbfr)
                elseif (vtype(i) == nf_byte) then
                  ierr=nf_get_vara_int1  (ncsrc, i, start,count, ibfr1)
                elseif (vtype(i) == nf_short) then
                  ierr=nf_get_vara_int2  (ncsrc, i, start,count, ibfr2)
                elseif (vtype(i) == nf_int) then
                  ierr=nf_get_vara_int   (ncsrc, i, start,count, ibffr)
                elseif (vtype(i) == nf_float) then
                  ierr=nf_get_vara_real  (ncsrc, i, start,count, bffr4)
                elseif (vtype(i) == nf_double) then
                  ierr=nf_get_vara_double(ncsrc, i, start,count, bffr8)
                endif
                if (ierr /= nf_noerr) then
                  write(*,'(/1x,A,I4,1x,3A/12x,3A/12x,A/)')
     &              '### ERROR: Cannot read time record =',    rec,
     &              'for nonpartitionable array ''', vname(i)(1:lvar),
     &              ''' from netCDF',   'file ''',   src_name(1:lsrc),
     &                                       '''.', nf_strerror(ierr)
                  goto 97
                else
                  nclk=3-nclk
                  call system_clock(iclk(nclk),clk_rate,clk_max)
                  inc_clk=iclk(nclk)-iclk(3-nclk)
                  if (inc_clk < 0) inc_clk=inc_clk+clk_max
                  net_clk=net_clk+inc_clk
                  net_read_clk=net_read_clk+inc_clk
                  net_read_size=net_read_size+size*vtypsize
                endif
                do node=0,nnodes-1
                  if (vid(i,node) > 0) then
                    if (vtype(i) == nf_char) then
                      ierr=nf_put_vara_text (ncid(node), vid(i,node),
     &                                             start,count, cbfr)
                    elseif (vtype(i) == nf_byte) then
                      ierr=nf_put_vara_int1  (ncid(node), vid(i,node),
     &                                            start,count, ibfr1)
                    elseif (vtype(i) == nf_short) then
                      ierr=nf_put_vara_int2  (ncid(node), vid(i,node),
     &                                            start,count, ibfr2)
                    elseif (vtype(i) == nf_int) then
                      ierr=nf_put_vara_int   (ncid(node), vid(i,node),
     &                                            start,count, ibffr)
                    elseif (vtype(i) == nf_float) then
                      ierr=nf_put_vara_real  (ncid(node), vid(i,node),
     &                                            start,count, bffr4)
                    elseif (vtype(i) == nf_double) then
                      ierr=nf_put_vara_double(ncid(node), vid(i,node),
     &                                            start,count, bffr8)
                  endif
                    if (ierr /= nf_noerr) then
                      write(*,'(/1x,A,I4,1x,3A/12x,3A/12x,A/)')
     &               '### ERROR: Cannot write time record =',     rec,
     &               'for nonpartitionable array ''',vname(i)(1:lvar),
     &               '''into netCDF',  'file ''', fname(node)(1:lfnm),
     &                                       '''.', nf_strerror(ierr)
                      goto 97
                    else
                      nclk=3-nclk
                      call system_clock(iclk(nclk), clk_rate,clk_max)
                      inc_clk=iclk(nclk)-iclk(3-nclk)
                      if (inc_clk < 0) inc_clk=inc_clk+clk_max
                      net_clk=net_clk+inc_clk
                      net_wrt_clk=net_wrt_clk+inc_clk
                      net_wrt_size=net_wrt_size+size*vtypsize
                    endif
                  endif
                enddo
              elseif (part_type(i) > 0) then
                if (rec < 5) then
                  write(*,'(21x,A,I4,1x,A,I4,1x,3A)') 'var',i,'out of',
     &              nvars,'partitioned array ''',vname(i)(1:lvar),''''
                elseif (rec < 6) then
                   write(*,'(24x,A)') '.........................'
                endif
                do j=1,vdims(i)
                  if (dimids(j,i)==unlimdimid) then
                    start1(j)=rec
                  else
                    start1(j)=1
                  endif
                enddo
                size1=1
                do j=1,vdims(i)
                  if (dimids(j,i)==unlimdimid) then
                    count1(j)=1
                  else
                    count1(j)=dimsize(dimids(j,i))
                  endif
                  size1=size1*count1(j)
                enddo
                if (vtype(i) == nf_char) then
                  ierr=chk_alloc(size1, allc_cbf_in, 'nf_char')
                elseif (vtype(i) == nf_byte) then
                  ierr=chk_alloc(size1, allc_ibf1_in,  'nf_byte')
                elseif (vtype(i) == nf_short) then
                  ierr=chk_alloc(size1, allc_ibf2_in,  'nf_short')
                elseif (vtype(i) == nf_int) then
                  ierr=chk_alloc(size1, allc_ibff_in,    'nf_int')
                elseif (vtype(i) == nf_float) then
                  ierr=chk_alloc(size1, allc_bff4_in,  'nf_float')
                elseif (vtype(i) == nf_double) then
                  ierr=chk_alloc(size1, allc_bff8_in, 'nf_double')
                else
                  write(*,'(/1x,4A/)')   '### ERROR: Unknown ',
     &            'type of variable ''',vname(i)(1:lvar),'''.'
                  stop
                endif
                if (ierr /= 0) stop
                if (vtype(i) == nf_char) then
                  ierr=nf_get_vara_text (ncsrc, i, start1,count1,
     &                                                     cbf_in)
                elseif (vtype(i) == nf_byte) then
                  ierr=nf_get_vara_int1   (ncsrc, i, start1,count1,
     &                                                    ibf1_in)
                elseif (vtype(i) == nf_short) then
                  ierr=nf_get_vara_int2   (ncsrc, i, start1,count1,
     &                                                    ibf2_in)
                elseif (vtype(i) == nf_int) then
                  ierr=nf_get_vara_int    (ncsrc, i, start1,count1,
     &                                                    ibff_in)
                elseif (vtype(i) == nf_float) then
                  ierr=nf_get_vara_real   (ncsrc, i, start1,count1,
     &                                                    bff4_in)
                elseif (vtype(i) == nf_double) then
                  ierr=nf_get_vara_double (ncsrc, i, start1,count1,
     &                                                    bff8_in)
                endif
                if (ierr /= nf_noerr) then
                  write(*,'(/1x,A,I4,1x,3A/12x,3A/12x,A)')
     &                '### ERROR: Cannot read time record =',  rec,
     &               'of partitioned array ''',   vname(i)(1:lvar),
     &               '''', 'from netCDF file ''', src_name(1:lsrc),
     &                                   '''.',  nf_strerror(ierr)
                  write(*,*) 'start1 =', (start1(j), j=1,vdims(i))
                  write(*,*) 'count1 =', (count1(j), j=1,vdims(i))
                  goto 97
                else
                  nclk=3-nclk
                  call system_clock(iclk(nclk),clk_rate,clk_max)
                  inc_clk=iclk(nclk)-iclk(3-nclk)
                  if (inc_clk < 0) inc_clk=inc_clk+clk_max
                  net_clk=net_clk+inc_clk
                  net_read_clk=net_read_clk+inc_clk
                  net_read_size=net_read_size+size1*vtypsize
                endif
                do node=0,nnodes-1
                  if (vid(i,node) > 0) then
                    size=1
                    do j=1,vdims(i)
                      if (dimids(j,i) == id_xi_rho) then
                        start(j)=xi_start(node)
                        count(j)=xi_size(node)
                        if (western_edge(node)) then
                          count(j)=count(j)+1
                        endif
                        if (eastern_edge(node)) then
                          count(j)=count(j)+1
                        endif
                      elseif (dimids(j,i) == id_xi_u) then
                        start(j)=xi_start(node)
                        count(j)=xi_size(node)
                        if (.not.western_edge(node)) then
                          start(j)=start(j)-1
                        endif
                        if (eastern_edge(node)) then
                          count(j)=count(j)+1
                        endif
                      elseif (dimids(j,i) == id_eta_rho) then
                        start(j)=eta_start(node)
                        count(j)=eta_size(node)
                        if (southern_edge(node)) then
                          count(j)=count(j)+1
                        endif
                        if (northern_edge(node)) then
                          count(j)=count(j)+1
                        endif
                      elseif (dimids(j,i) == id_eta_v) then
                        start(j)=eta_start(node)
                        count(j)=eta_size(node)
                        if (.not.southern_edge(node)) then
                          start(j)=start(j)-1
                        endif
                        if (northern_edge(node)) then
                          count(j)=count(j)+1
                        endif
                      elseif (dimids(j,i) == unlimdimid) then
                        start(j)=rec  ; count(j)=1
                      else
                        start(j)=1 ; count(j)=dimsize(dimids(j,i))
                      endif
                      size=size*count(j)
                    enddo
                    if (vtype(i) == nf_char) then
                      ierr=chk_alloc(size, alloc_cbfr,  'nf_char')
                    elseif (vtype(i) == nf_byte) then
                      ierr=chk_alloc(size, alloc_ibfr1,  'nf_byte')
                    elseif (vtype(i) == nf_short) then
                      ierr=chk_alloc(size, alloc_ibfr2,  'nf_short')
                    elseif (vtype(i) == nf_int) then
                      ierr=chk_alloc(size, alloc_ibffr,    'nf_int')
                    elseif (vtype(i) == nf_float) then
                      ierr=chk_alloc(size, alloc_bffr4,  'nf_float')
                    elseif (vtype(i) == nf_double) then
                      ierr=chk_alloc(size, alloc_bffr8, 'nf_double')
                    else
                      write(*,'(/1x,4A/)')   '### ERROR: Unknown ',
     &                'type of variable ''',vname(i)(1:lvar),'''.'
                      stop
                    endif
                    if (ierr /= 0) stop
                    if (vtype(i) == nf_char) then
                      call extract_text (cbf_in, cbfr, start, count,
     &                                              count1, vdims(i))
                    elseif (vtype(i) == nf_byte) then
                      call extract_int1  (ibf1_in, ibfr1, start,count,
     &                                              count1, vdims(i))
                    elseif (vtype(i) == nf_short) then
                      call extract_int2  (ibf2_in, ibfr2, start,count,
     &                                              count1, vdims(i))
                    elseif (vtype(i) == nf_int) then
                      call extract_int   (ibff_in, ibffr, start,count,
     &                                              count1, vdims(i))
                    elseif (vtype(i) == nf_float) then
                      call extract_real  (bff4_in, bffr4, start,count,
     &                                              count1, vdims(i))
                    elseif (vtype(i) == nf_double) then
                      call extract_double(bff8_in, bffr8, start,count,
     &                                              count1, vdims(i))
                    endif
                    nclk=3-nclk
                    call system_clock(iclk(nclk), clk_rate, clk_max)
                    inc_clk=iclk(nclk)-iclk(3-nclk)
                    if (inc_clk < 0) inc_clk=inc_clk+clk_max
                    net_clk=net_clk+inc_clk
                    net_assm_clk=net_assm_clk+inc_clk
                    if (vtype(i) == nf_char) then
                      ierr=nf_put_vara_text (ncid(node), vid(i,node),
     &                                           start1, count, cbfr)
                    elseif (vtype(i) == nf_byte) then
                      ierr=nf_put_vara_int1  (ncid(node), vid(i,node),
     &                                         start1,  count, ibfr1)
                    elseif (vtype(i) == nf_short) then
                      ierr=nf_put_vara_int2  (ncid(node), vid(i,node),
     &                                          start1, count, ibfr2)
                    elseif (vtype(i) == nf_int) then
                      ierr=nf_put_vara_int   (ncid(node), vid(i,node),
     &                                        start1, count, ibffr)
                    elseif (vtype(i) == nf_float) then
                      ierr=nf_put_vara_real  (ncid(node), vid(i,node),
     &                                          start1, count, bffr4)
                    elseif (vtype(i) == nf_double) then
                      ierr=nf_put_vara_double(ncid(node), vid(i,node),
     &                                          start1, count, bffr8)
                    endif
                    if (ierr /= nf_noerr) then
                      write(*,'(/1x,A,I4,1x,3A/12x,3A/12x,A)')
     &                   '### ERROR: Cannot write time record =',  rec,
     &                     'of partitioned array ''', vname(i)(1:lvar),
     &                '''', 'into netCDF file ''', fname(node)(1:lfnm),
     &                                        '''.', nf_strerror(ierr)
                      write(*,*) 'start1 =', (start1(j), j=1,vdims(i))
                      write(*,*) ' count =',  (count(j), j=1,vdims(i))
                      goto 97
                    else
                      nclk=3-nclk
                      call system_clock(iclk(nclk), clk_rate,clk_max)
                      inc_clk=iclk(nclk)-iclk(3-nclk)
                      if (inc_clk < 0) inc_clk=inc_clk+clk_max
                      net_clk=net_clk+inc_clk
                      net_wrt_clk=net_wrt_clk+inc_clk
                      net_wrt_size=net_wrt_size+size*vtypsize
                    endif
                  endif
                enddo
              endif
            endif
          enddo
        enddo
  97    write(*,*) 'closing files...'
        ierr=nf_close(ncsrc)
        nclk=3-nclk
        call system_clock(iclk(nclk), clk_rate, clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        if (inc_clk < 0) inc_clk=inc_clk+clk_max
        net_clk=net_clk+inc_clk
        net_gray_clk=net_gray_clk+inc_clk
        write(*,*) '...........input'
        do node=0,nnodes-1
         ierr=nf_close(ncid(node))
        enddo
        nclk=3-nclk
        call system_clock(iclk(nclk), clk_rate, clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        if (inc_clk < 0) inc_clk=inc_clk+clk_max
        net_clk=net_clk+inc_clk
        net_sync_clk=net_sync_clk+inc_clk
        write(*,*) '..........output'
      enddo
      call cpu_time(run_time) ; run_time=run_time-tstart
      write(*,'(/1x,A,F11.2,1x,A/)') 'cpu time =', run_time, 'sec'
      if (clk_rate > 0) then
        ReadSize=1.0D-6*dble(net_read_size)
        WrtSize=1.0D-6*dble(net_wrt_size)
        ReadTime=net_read_clk/dble(clk_rate)
        WrtTime = net_wrt_clk/dble(clk_rate)
        InitTime=net_init_clk/dble(clk_rate)
        SyncTime=net_sync_clk/dble(clk_rate)
        write(*,'(/1x,A,22x,F12.2,1x,A)') 'Analysis/file creation :',
     &                                               InitTime, 'sec'
        write(*,'(8x,A,F12.2,1x,A,F12.2,1x,A,F8.2,1x,A)')
     &         'Total data read :', ReadSize, 'MBytes in',  ReadTime,
     &                          'sec (', ReadSize/ReadTime, 'MB/sec)'
        write(*,'(5x,A,F12.2,1x,A,F12.2,1x,A,F8.2,1x,A)')
     &      'Total data written :', WrtSize,  'MBytes in',   WrtTime,
     &                          'sec (',  WrtSize/WrtTime,  'MB/sec)'
        AssmTime=net_assm_clk/dble(clk_rate)
        write(*,'(3x,A,22x,F12.2,1x,A)')    'Array splitting time :',
     &                                               AssmTime, 'sec'
        write(*,'(2x,A,22x,F12.2,1x,A)')   'Output file sync time :',
     &                                               SyncTime, 'sec'
        nclk=3-nclk
        call system_clock(iclk(nclk), clk_rate, clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        if (inc_clk < 0) inc_clk=inc_clk+clk_max
        net_clk=net_clk+inc_clk
        net_gray_clk=net_gray_clk+inc_clk
        GrayTime=dble(net_gray_clk)/dble(clk_rate)
        write(*,'(14x,A,22x,F12.2,1x,A)') 'Gray time :',GrayTime,'sec'
        write(*,'(47x,A/12x,A,11x,F12.2,1x,A/)') '------------------',
     &   'Elapsed wall-clock time:',dble(net_clk)/dble(clk_rate),'sec'
      endif
      stop
      end
      subroutine mpi_setup (NP_XI,NP_ETA, xi_rho,eta_rho,
     &               xi_start,xi_size, eta_start,eta_size,
     &                        western_edge, eastern_edge,
     &                       southern_edge, northern_edge)
      implicit none
      integer(kind=4), intent(in) :: NP_XI,NP_ETA, xi_rho,eta_rho
      integer(kind=4), intent(out), dimension(0:NP_XI*NP_ETA-1) ::
     &              xi_start, xi_size, eta_start, eta_size
      logical, intent(out), dimension(0:NP_XI*NP_ETA-1) ::
     & western_edge, eastern_edge, southern_edge, northern_edge
      integer(kind=4) LLm,Lm, MMm,Mm, nnodes,node, inode,jnode,
     &        iwest,ieast,jsouth,jnorth, iSW_corn,jSW_corn,
     &                                   off_XI,off_ETA
      nnodes=NP_XI*NP_ETA
      LLm=xi_rho-2           ; MMm=eta_rho-2
      Lm=(LLm+NP_XI-1)/NP_XI ; Mm=(MMm+NP_ETA-1)/NP_ETA
      do node=0,nnodes-1
        jnode=node/NP_XI
        inode=node-jnode*NP_XI
        off_XI=NP_XI*Lm-LLm
        iSW_corn=inode*Lm-off_XI/2
        if (inode == 0) then
          iwest=1+off_XI/2
        else
          iwest=1
        endif
        if (inode < NP_XI-1) then
          ieast=Lm
        else
          ieast=Lm -(off_XI+1)/2
        endif
        off_ETA=NP_ETA*Mm-MMm
        jSW_corn=jnode*Mm-off_ETA/2
        if (jnode == 0) then
          jsouth=1+off_ETA/2
        else
          jsouth=1
        endif
        if (jnode < NP_ETA-1) then
          jnorth=Mm
        else
          jnorth=Mm -(off_ETA+1)/2
        endif
        xi_size(node)=ieast-iwest+1
        eta_size(node)=jnorth-jsouth+1
        if (inode == 0) then
          xi_start(node)=iSW_corn+iwest
        else
          xi_start(node)=iSW_corn+iwest+1
        endif
        if (jnode == 0) then
          eta_start(node)=jSW_corn+jsouth
        else
          eta_start(node)=jSW_corn+jsouth+1
        endif
        if (inode == 0) then
          western_edge(node)=.true.
        else
          western_edge(node)=.false.
        endif
        if (inode < NP_XI-1) then
          eastern_edge(node)=.false.
        else
          eastern_edge(node)=.true.
        endif
        if (jnode == 0) then
          southern_edge(node)=.true.
        else
          southern_edge(node)=.false.
        endif
        if (jnode < NP_ETA-1) then
          northern_edge(node)=.false.
        else
          northern_edge(node)=.true.
        endif
      enddo
      end
      function chk_alloc(size, alloc_size, type_str)
      implicit none
      integer(kind=4) chk_alloc, size, alloc_size
      character(len=*) type_str
      if (size > alloc_size) then
        write(*,'(/1x,4A/12x,A,I10,2x,A,I10/)')          '### ERROR: ',
     & 'Insufficient size of allocated workspace of type ''', type_str,
     &     ''':', 'needed', size, 'but actually allocated', alloc_size
        chk_alloc=1
      else
       chk_alloc=0
      endif
      end
      subroutine extract_text (bfr_in, bffr, start,count,
     &                                      count1,vdims)
      implicit none
      character(len=1) bfr_in(*), bffr(*)
      integer(kind=4) vdims, start(vdims), count(vdims), count1(vdims)
      integer(kind=4) ndims, i, istr,imax,imax1, j,js,js1, 
     &                          jstr,jmax,jmax1,
     &        k,ks,ks1, kstr,kmax,kmax1, l,ls,ls1, lstr,lmax,lmax1
      if (count1(vdims)==1) then
        ndims=vdims-1
      else
        ndims=vdims
      endif
      istr=start(1)
      imax=count(1)
      imax1=count1(1)
      if (ndims > 1) then
        jstr=start(2)
        jmax=count(2)
        jmax1=count1(2)
      else
        jstr=1 ; jmax=1 ; jmax1=1
      endif
      if (ndims > 2) then
        kstr=start(3) ; kmax=count(3) ; kmax1=count1(3)
      else
        kstr=1 ; kmax=1 ; kmax1=1
      endif
      if (ndims > 3) then
        lstr=start(4) ; lmax=count(4) ; lmax1=count1(4)
      else
        lstr=1 ; lmax=1 ; lmax1=1
      endif
      if (ndims > 4) then
        write(*,'(/1x,2A/12x,A/)')   '### ERROR: Exceeding limit of ',
     &                           '4 dimensions for partitioned array',
     &                        '[unlimited dimension does not count].'
        stop
      endif
      do l=1,lmax
        ls=l-1           ; ls1=l+lstr-2
        do k=1,kmax
          ks=k-1 +ls*kmax  ; ks1=k+kstr-2 +ls1*kmax1
          do j=1,jmax
            js=j-1 +ks*jmax  ; js1=j+jstr-2 +ks1*jmax1
            do i=1,imax
              bffr(i +js*imax)  = bfr_in(i+istr-1 +js1*imax1)
            enddo
          enddo
        enddo
      enddo
      end
      subroutine extract_int1 (bfr_in, bffr, start,count,
     &                                     count1, vdims)
      implicit none
      integer(kind=1) ::  bfr_in(*), bffr(*)
      integer(kind=4) vdims, start(vdims), count(vdims), count1(vdims)
      integer(kind=4) ndims, i, istr,imax,imax1, j,js,js1, 
     &                          jstr,jmax,jmax1,
     &        k,ks,ks1, kstr,kmax,kmax1, l,ls,ls1, lstr,lmax,lmax1
      if (count1(vdims)==1) then
        ndims=vdims-1
      else
        ndims=vdims
      endif
      istr=start(1)
      imax=count(1)
      imax1=count1(1)
      if (ndims > 1) then
        jstr=start(2)
        jmax=count(2)
        jmax1=count1(2)
      else
        jstr=1 ; jmax=1 ; jmax1=1
      endif
      if (ndims > 2) then
        kstr=start(3) ; kmax=count(3) ; kmax1=count1(3)
      else
        kstr=1 ; kmax=1 ; kmax1=1
      endif
      if (ndims > 3) then
        lstr=start(4) ; lmax=count(4) ; lmax1=count1(4)
      else
        lstr=1 ; lmax=1 ; lmax1=1
      endif
      if (ndims > 4) then
        write(*,'(/1x,2A/12x,A/)')   '### ERROR: Exceeding limit of ',
     &                           '4 dimensions for partitioned array',
     &                        '[unlimited dimension does not count].'
        stop
      endif
      do l=1,lmax
        ls=l-1           ; ls1=l+lstr-2
        do k=1,kmax
          ks=k-1 +ls*kmax  ; ks1=k+kstr-2 +ls1*kmax1
          do j=1,jmax
            js=j-1 +ks*jmax  ; js1=j+jstr-2 +ks1*jmax1
            do i=1,imax
              bffr(i +js*imax)  = bfr_in(i+istr-1 +js1*imax1)
            enddo
          enddo
        enddo
      enddo
      end
      subroutine extract_int2 (bfr_in, bffr, start,count,
     &                                     count1, vdims)
      implicit none
      integer(kind=2) :: bfr_in(*), bffr(*)
      integer(kind=4) vdims, start(vdims), count(vdims), count1(vdims)
      integer(kind=4) ndims, i, istr,imax,imax1, j,js,js1, 
     &                          jstr,jmax,jmax1,
     &        k,ks,ks1, kstr,kmax,kmax1, l,ls,ls1, lstr,lmax,lmax1
      if (count1(vdims)==1) then
        ndims=vdims-1
      else
        ndims=vdims
      endif
      istr=start(1)
      imax=count(1)
      imax1=count1(1)
      if (ndims > 1) then
        jstr=start(2)
        jmax=count(2)
        jmax1=count1(2)
      else
        jstr=1 ; jmax=1 ; jmax1=1
      endif
      if (ndims > 2) then
        kstr=start(3) ; kmax=count(3) ; kmax1=count1(3)
      else
        kstr=1 ; kmax=1 ; kmax1=1
      endif
      if (ndims > 3) then
        lstr=start(4) ; lmax=count(4) ; lmax1=count1(4)
      else
        lstr=1 ; lmax=1 ; lmax1=1
      endif
      if (ndims > 4) then
        write(*,'(/1x,2A/12x,A/)')   '### ERROR: Exceeding limit of ',
     &                           '4 dimensions for partitioned array',
     &                        '[unlimited dimension does not count].'
        stop
      endif
      do l=1,lmax
        ls=l-1           ; ls1=l+lstr-2
        do k=1,kmax
          ks=k-1 +ls*kmax  ; ks1=k+kstr-2 +ls1*kmax1
          do j=1,jmax
            js=j-1 +ks*jmax  ; js1=j+jstr-2 +ks1*jmax1
            do i=1,imax
              bffr(i +js*imax)  = bfr_in(i+istr-1 +js1*imax1)
            enddo
          enddo
        enddo
      enddo
      end
      subroutine extract_int   (bfr_in, bffr, start,count,
     &                                      count1, vdims)
      implicit none
      integer(kind=4) :: bfr_in(*), bffr(*)
      integer(kind=4) vdims, start(vdims), count(vdims), count1(vdims)
      integer(kind=4) ndims, i, istr,imax,imax1, j,js,js1, 
     &                          jstr,jmax,jmax1,
     &        k,ks,ks1, kstr,kmax,kmax1, l,ls,ls1, lstr,lmax,lmax1
      if (count1(vdims)==1) then
        ndims=vdims-1
      else
        ndims=vdims
      endif
      istr=start(1)
      imax=count(1)
      imax1=count1(1)
      if (ndims > 1) then
        jstr=start(2)
        jmax=count(2)
        jmax1=count1(2)
      else
        jstr=1 ; jmax=1 ; jmax1=1
      endif
      if (ndims > 2) then
        kstr=start(3) ; kmax=count(3) ; kmax1=count1(3)
      else
        kstr=1 ; kmax=1 ; kmax1=1
      endif
      if (ndims > 3) then
        lstr=start(4) ; lmax=count(4) ; lmax1=count1(4)
      else
        lstr=1 ; lmax=1 ; lmax1=1
      endif
      if (ndims > 4) then
        write(*,'(/1x,2A/12x,A/)')   '### ERROR: Exceeding limit of ',
     &                           '4 dimensions for partitioned array',
     &                        '[unlimited dimension does not count].'
        stop
      endif
      do l=1,lmax
        ls=l-1           ; ls1=l+lstr-2
        do k=1,kmax
          ks=k-1 +ls*kmax  ; ks1=k+kstr-2 +ls1*kmax1
          do j=1,jmax
            js=j-1 +ks*jmax  ; js1=j+jstr-2 +ks1*jmax1
            do i=1,imax
              bffr(i +js*imax)  = bfr_in(i+istr-1 +js1*imax1)
            enddo
          enddo
        enddo
      enddo
      end
      subroutine extract_real   (bfr_in, bffr, start,count,
     &                                      count1, vdims)
      implicit none
      real(kind=4) :: bfr_in(*), bffr(*)
      integer(kind=4) vdims, start(vdims), count(vdims), count1(vdims)
      integer(kind=4) ndims, i, istr,imax,imax1, j,js,js1, 
     &                          jstr,jmax,jmax1,
     &        k,ks,ks1, kstr,kmax,kmax1, l,ls,ls1, lstr,lmax,lmax1
      if (count1(vdims)==1) then
        ndims=vdims-1
      else
        ndims=vdims
      endif
      istr=start(1)
      imax=count(1)
      imax1=count1(1)
      if (ndims > 1) then
        jstr=start(2)
        jmax=count(2)
        jmax1=count1(2)
      else
        jstr=1 ; jmax=1 ; jmax1=1
      endif
      if (ndims > 2) then
        kstr=start(3) ; kmax=count(3) ; kmax1=count1(3)
      else
        kstr=1 ; kmax=1 ; kmax1=1
      endif
      if (ndims > 3) then
        lstr=start(4) ; lmax=count(4) ; lmax1=count1(4)
      else
        lstr=1 ; lmax=1 ; lmax1=1
      endif
      if (ndims > 4) then
        write(*,'(/1x,2A/12x,A/)')   '### ERROR: Exceeding limit of ',
     &                           '4 dimensions for partitioned array',
     &                        '[unlimited dimension does not count].'
        stop
      endif
      do l=1,lmax
        ls=l-1           ; ls1=l+lstr-2
        do k=1,kmax
          ks=k-1 +ls*kmax  ; ks1=k+kstr-2 +ls1*kmax1
          do j=1,jmax
            js=j-1 +ks*jmax  ; js1=j+jstr-2 +ks1*jmax1
            do i=1,imax
              bffr(i +js*imax)  = bfr_in(i+istr-1 +js1*imax1)
            enddo
          enddo
        enddo
      enddo
      end
      subroutine extract_double (bfr_in, bffr, start,count,
     &                                      count1, vdims)
      implicit none
      real(kind=8) :: bfr_in(*), bffr(*)
      integer(kind=4) vdims, start(vdims), count(vdims), count1(vdims)
      integer(kind=4) ndims, i, istr,imax,imax1, j,js,js1, 
     &                          jstr,jmax,jmax1,
     &        k,ks,ks1, kstr,kmax,kmax1, l,ls,ls1, lstr,lmax,lmax1
      if (count1(vdims)==1) then
        ndims=vdims-1
      else
        ndims=vdims
      endif
      istr=start(1)
      imax=count(1)
      imax1=count1(1)
      if (ndims > 1) then
        jstr=start(2)
        jmax=count(2)
        jmax1=count1(2)
      else
        jstr=1 ; jmax=1 ; jmax1=1
      endif
      if (ndims > 2) then
        kstr=start(3) ; kmax=count(3) ; kmax1=count1(3)
      else
        kstr=1 ; kmax=1 ; kmax1=1
      endif
      if (ndims > 3) then
        lstr=start(4) ; lmax=count(4) ; lmax1=count1(4)
      else
        lstr=1 ; lmax=1 ; lmax1=1
      endif
      if (ndims > 4) then
        write(*,'(/1x,2A/12x,A/)')   '### ERROR: Exceeding limit of ',
     &                           '4 dimensions for partitioned array',
     &                        '[unlimited dimension does not count].'
        stop
      endif
      do l=1,lmax
        ls=l-1           ; ls1=l+lstr-2
        do k=1,kmax
          ks=k-1 +ls*kmax  ; ks1=k+kstr-2 +ls1*kmax1
          do j=1,jmax
            js=j-1 +ks*jmax  ; js1=j+jstr-2 +ks1*jmax1
            do i=1,imax
              bffr(i +js*imax)  = bfr_in(i+istr-1 +js1*imax1)
            enddo
          enddo
        enddo
      enddo
      end
