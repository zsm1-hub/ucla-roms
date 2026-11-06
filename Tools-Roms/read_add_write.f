      subroutine write_by_name (ncid, vname, rec, var)
      use mod_io_size_acct
      implicit none
      character(len=*) vname
      character(len=32) name
      integer(kind=4) ncid, rec, indx, varid, vtype, ndims, dimid(8), 
     &                               natts,
     &        ierr, size,  start(8),  count(8), i, lvar, lnm, lenstr
      real(kind=4), dimension(*) :: var
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
      lvar=lenstr(vname)
      do i=1,8
        start(i)=0
        count(i)=0
        dimid(i)=0
      enddo
      indx=0
      ierr=nf_inq_varid (ncid, vname, varid)
      if (ierr==nf_noerr) then
        ierr=nf_inq_var (ncid, varid, name,  vtype,
     &                         ndims, dimid, natts)
        if (ierr == nf_noerr) then
          do i=1,ndims
            ierr=nf_inq_dim (ncid, dimid(i), name, size)
            if (ierr == nf_noerr) then
              lnm=lenstr(name)
              start(i)=1
              count(i)=size
              if (i == ndims  .and. lnm > 3) then
                if (name(lnm-3:lnm) == 'time') then
                  start(i)=rec
                  count(i)=1
                  indx=rec
                endif
              endif
            else
              write(*,'(/1x,3A,I3,A/12x,A/)') '### ERROR: Cannot ',
     &                'determine name and size of dimension with ',
     &                'ID =',  dimid(i),  ':',  nf_strerror(ierr)
              stop
            endif
          enddo
          size=1
          do i=1,ndims
            size=size*count(i)
          enddo
          if (vtype == nf_real) then
            ierr=nf_put_vara_real (ncid, varid, start,count, var)
            sz_write_acc=sz_write_acc + 4*size
          elseif (vtype == nf_double) then
            ierr=nf_put_vara_double(ncid, varid, start,count, var)
            sz_write_acc=sz_write_acc + 8*size
          elseif (vtype == nf_int) then
            ierr=nf_put_vara_int   (ncid, varid, start,count, var)
            sz_write_acc=sz_write_acc + 4*size
          elseif (vtype == nf_char) then
            ierr=nf_put_vara_text  (ncid, varid, start,count, var)
            sz_write_acc=sz_write_acc + size
          else
            write(*,'(/1x,3A/)')   '### ERROR: variable ''',
     &              vname(1:lvar), ''' is of unknown type.'
            stop
          endif
          if (ierr == nf_noerr) then
            write(*,'(5x,A,I5,1x,3A)')  'wrote rec', indx, 'of ''',
     &                                         vname(1:lvar), ''''
          else
            write(*,'(/1x,4A/12x,A/)')  '### ERROR: Cannot write ',
     &                         'netCDF variable ''', vname(1:lvar),
     &                                   ''':',  nf_strerror(ierr)
          endif
        else
          write(*,'(/1x,2A,I3,2x,3A/12x,A/)')   '### ERROR: Cannot ',
     &      'make general inquiry for variable #', varid, 'named ''',
     &                       vname(1:lvar), ''':', nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,4A/12x,A/)')   '### ERROR: Cannot get netCDF ',
     &                         'ID for variable ''', vname(1:lvar),
     &                                      ''':', nf_strerror(ierr)
      endif
      if (ierr/=nf_noerr) stop
      return
      end
      subroutine read_by_name  (ncid, vname, rec, var)
      use mod_io_size_acct
      implicit none
      character(len=*) vname
      character(len=32) name
      integer(kind=4) ncid, rec, indx, varid, vtype, ndims, dimid(8), 
     &                               natts,
     &        ierr, size,  start(8),  count(8), i, lvar, lnm, lenstr
      real(kind=4), dimension(*) :: var
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
      lvar=lenstr(vname)
      do i=1,8
        start(i)=0
        count(i)=0
        dimid(i)=0
      enddo
      indx=0
      ierr=nf_inq_varid (ncid, vname, varid)
      if (ierr==nf_noerr) then
        ierr=nf_inq_var (ncid, varid, name,  vtype,
     &                         ndims, dimid, natts)
        if (ierr == nf_noerr) then
          do i=1,ndims
            ierr=nf_inq_dim (ncid, dimid(i), name, size)
            if (ierr == nf_noerr) then
              lnm=lenstr(name)
              start(i)=1
              count(i)=size
              if (i == ndims  .and. lnm > 3) then
                if (name(lnm-3:lnm) == 'time') then
                  start(i)=rec
                  count(i)=1
                  indx=rec
                endif
              endif
            else
              write(*,'(/1x,3A,I3,A/12x,A/)') '### ERROR: Cannot ',
     &                'determine name and size of dimension with ',
     &                'ID =',  dimid(i),  ':',  nf_strerror(ierr)
              stop
            endif
          enddo
          size=1
          do i=1,ndims
            size=size*count(i)
          enddo
          if (vtype == nf_real) then
            ierr=nf_get_vara_real(ncid, varid, start,count, var)
            sz_read_acc=sz_read_acc + 4*size
          elseif (vtype == nf_double) then
            ierr=nf_get_vara_double(ncid, varid, start,count, var)
            sz_read_acc=sz_read_acc + 8*size
          elseif (vtype == nf_int) then
            ierr=nf_get_vara_int  (ncid, varid, start,count, var)
            sz_read_acc=sz_read_acc + 4*size
          elseif (vtype == nf_char) then
            ierr=nf_get_vara_text (ncid, varid, start,count, var)
            sz_read_acc=sz_read_acc + size
          else
            write(*,'(/1x,3A/)')   '### ERROR: netCDF variable ''',
     &          vname(1:lvar), ''' is of unknown or illegal type.'
            stop
          endif
          if (ierr == nf_noerr) then
            write(*,'(6x,A,I5,1x,3A)')  'read rec', indx, 'of ''',
     &                                       vname(1:lvar), '''.'
          else
            write(*,'(/1x,4A/12x,A)')    '### ERROR: Cannot read ',
     &              'netCDF variable ''',    vname(1:lvar),  ''':',
     &                                           nf_strerror(ierr)
            write(*,'(12x,A,4I5)') 'start =', (start(i),i=1,ndims)
            write(*,'(12x,A,4I5)') 'count =', (count(i),i=1,ndims)
            stop
          endif
        else
          write(*,'(/1x,2A,I3,2x,3A/12x,A/)')   '### ERROR: Cannot ',
     &      'make general inquiry for variable #', varid, 'named ''',
     &                       vname(1:lvar), ''':', nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,4A/12x,A/)')   '### ERROR: Cannot get netCDF ',
     &                         'ID for variable ''', vname(1:lvar),
     &                                      ''':', nf_strerror(ierr)
      endif
      if (ierr/=nf_noerr) stop
      return
      end
      subroutine read_and_add  (ncid, vname, rec, var,  buff)
      use mod_io_size_acct
      implicit none
      character(len=*) vname
      character(len=32) name
      integer(kind=4) ncid, rec, indx, varid, vtype, ndims, dimid(8), 
     &                               natts,
     &        ierr, size,  start(8),  count(8), i, lvar, lnm, lenstr
      real(kind=4), dimension(*) :: var
     &                             , buff
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
      lvar=lenstr(vname)
      do i=1,8
        start(i)=0
        count(i)=0
        dimid(i)=0
      enddo
      indx=0
      ierr=nf_inq_varid (ncid, vname, varid)
      if (ierr==nf_noerr) then
        ierr=nf_inq_var (ncid, varid, name,  vtype,
     &                         ndims, dimid, natts)
        if (ierr == nf_noerr) then
          do i=1,ndims
            ierr=nf_inq_dim (ncid, dimid(i), name, size)
            if (ierr == nf_noerr) then
              lnm=lenstr(name)
              start(i)=1
              count(i)=size
              if (i == ndims  .and. lnm > 3) then
                if (name(lnm-3:lnm) == 'time') then
                  start(i)=rec
                  count(i)=1
                  indx=rec
                endif
              endif
            else
              write(*,'(/1x,3A,I3,A/12x,A/)') '### ERROR: Cannot ',
     &                'determine name and size of dimension with ',
     &                'ID =',  dimid(i),  ':',  nf_strerror(ierr)
              stop
            endif
          enddo
          size=1
          do i=1,ndims
            size=size*count(i)
          enddo
          if (vtype == nf_real) then
            ierr=nf_get_vara_real(ncid, varid, start,count, buff)
            sz_read_acc=sz_read_acc + 4*size
          elseif (vtype == nf_double) then
            ierr=nf_get_vara_double(ncid, varid, start,count, buff)
            sz_read_acc=sz_read_acc + 8*size
          else
            write(*,'(/1x,3A/)')   '### ERROR: netCDF variable ''',
     &          vname(1:lvar), ''' is of unknown or illegal type.'
            stop
          endif
          if (ierr == nf_noerr) then
            write(*,'(6x,A,I5,1x,3A)')  'read rec', indx, 'of ''',
     &                                       vname(1:lvar), '''.'
          else
            write(*,'(/1x,4A/12x,A)')    '### ERROR: Cannot read ',
     &              'netCDF variable ''',    vname(1:lvar),  ''':',
     &                                           nf_strerror(ierr)
            write(*,'(12x,A,4I5)') 'start =', (start(i),i=1,ndims)
            write(*,'(12x,A,4I5)') 'count =', (count(i),i=1,ndims)
            stop
          endif
            call add_buff_to_var_real (var, buff, size)
        else
          write(*,'(/1x,2A,I3,2x,3A/12x,A/)')   '### ERROR: Cannot ',
     &      'make general inquiry for variable #', varid, 'named ''',
     &                       vname(1:lvar), ''':', nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,4A/12x,A/)')   '### ERROR: Cannot get netCDF ',
     &                         'ID for variable ''', vname(1:lvar),
     &                                      ''':', nf_strerror(ierr)
      endif
      if (ierr/=nf_noerr) stop
      return
      end
       subroutine add_buff_to_var_real (var, buff, size)
       implicit none
       integer(kind=4) size, i
       real(kind=4), dimension(size) :: var, buff
       do i=1,size
         var(i)=var(i) + buff(i)
       enddo
       end
       subroutine add_buff_to_var_double (var, buff, size)
       implicit none
       integer(kind=4) size, i
       real(kind=8), dimension(size) :: var, buff
       do i=1,size
         var(i)=var(i) + buff(i)
       enddo
       end
      subroutine quad_var (ncid, vname, rec, var, buff, average)
      use mod_io_size_acct
      implicit none
      character(len=*) vname
      character(len=32) name
      integer(kind=4) ncid, rec, indx, varid, vtype, ndims, dimid(8), 
     &                               natts,
     &        ierr, size,  start(8),  count(8), i, lvar, lnm, lenstr
      real(kind=4), dimension(*) :: var
     &                             , buff
     &                             , average
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
      lvar=lenstr(vname)
      do i=1,8
        start(i)=0
        count(i)=0
        dimid(i)=0
      enddo
      indx=0
      ierr=nf_inq_varid (ncid, vname, varid)
      if (ierr==nf_noerr) then
        ierr=nf_inq_var (ncid, varid, name,  vtype,
     &                         ndims, dimid, natts)
        if (ierr == nf_noerr) then
          do i=1,ndims
            ierr=nf_inq_dim (ncid, dimid(i), name, size)
            if (ierr == nf_noerr) then
              lnm=lenstr(name)
              start(i)=1
              count(i)=size
              if (i == ndims  .and. lnm > 3) then
                if (name(lnm-3:lnm) == 'time') then
                  start(i)=rec
                  count(i)=1
                  indx=rec
                endif
              endif
            else
              write(*,'(/1x,3A,I3,A/12x,A/)') '### ERROR: Cannot ',
     &                'determine name and size of dimension with ',
     &                'ID =',  dimid(i),  ':',  nf_strerror(ierr)
              stop
            endif
          enddo
          size=1
          do i=1,ndims
            size=size*count(i)
          enddo
          if (vtype == nf_real) then
            ierr=nf_get_vara_real(ncid, varid, start,count, buff)
            sz_read_acc=sz_read_acc + 4*size
          elseif (vtype == nf_double) then
            ierr=nf_get_vara_double(ncid, varid, start,count, buff)
            sz_read_acc=sz_read_acc + 8*size
          else
            write(*,'(/1x,3A/)')   '### ERROR: netCDF variable ''',
     &          vname(1:lvar), ''' is of unknown or illegal type.'
            stop
          endif
          if (ierr == nf_noerr) then
            write(*,'(6x,A,I5,1x,3A)')  'read rec', indx, 'of ''',
     &                                       vname(1:lvar), '''.'
          else
            write(*,'(/1x,4A/12x,A)')    '### ERROR: Cannot read ',
     &              'netCDF variable ''',    vname(1:lvar),  ''':',
     &                                           nf_strerror(ierr)
            write(*,'(12x,A,4I5)') 'start =', (start(i),i=1,ndims)
            write(*,'(12x,A,4I5)') 'count =', (count(i),i=1,ndims)
            stop
          endif
            call add_quad_var_real (var, buff, average, size)
        else
          write(*,'(/1x,2A,I3,2x,3A/12x,A/)')   '### ERROR: Cannot ',
     &      'make general inquiry for variable #', varid, 'named ''',
     &                       vname(1:lvar), ''':', nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,4A/12x,A/)')   '### ERROR: Cannot get netCDF ',
     &                         'ID for variable ''', vname(1:lvar),
     &                                      ''':', nf_strerror(ierr)
      endif
      if (ierr/=nf_noerr) stop
      return
      end
       subroutine add_quad_var_real (var, buff, average, size)
       implicit none
       integer(kind=4) size, i
       real(kind=4), dimension(size) :: var, buff, average
       do i=1,size
         var(i)=var(i) + (buff(i)-average(i))**2
       enddo
       return
       end
       subroutine add_quad_var_double (var, buff, average, size)
       implicit none
       integer(kind=4) size, i
       real(kind=8), dimension(size) :: var, buff, average
       do i=1,size
         var(i)=var(i) + (buff(i)-average(i))**2
       enddo
       return
       end
