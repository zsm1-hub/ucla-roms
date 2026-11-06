      subroutine set_mask(nx,ny, rmsk, mask,umask,vmask)
      implicit none
      integer(kind=4) nx,ny, i,j
      real(kind=8) rmsk(nx,ny)
      integer(kind=2) mask(nx,ny), umask(2:nx,ny), vmask(nx,2:ny)
      do j=1,ny
        do i=1,nx
          if (rmsk(i,j) < 0.5D0) then
            mask(i,j)=0
          else
            mask(i,j)=1
          endif
        enddo
        do i=2,nx
          if (rmsk(i,j) < 0.5D0 .and. rmsk(i-1,j) < 0.5D0) then
            umask(i,j)=0
          else
            umask(i,j)=1
          endif
        enddo
        if (j > 1) then
          do i=1,nx
            if (rmsk(i,j) < 0.5D0 .and. rmsk(i,j-1) < 0.5D0) then
              vmask(i,j)=0
            else
              vmask(i,j)=1
            endif
          enddo
        endif
      enddo
      end
      subroutine set_contour_mask(nx,ny, rmsk, mask,umask,vmask)
      implicit none
      integer(kind=4) nx,ny, i,j, iw2,iw,ie,js2,js,jn
      real(kind=8) rmsk(nx,ny)
      integer(kind=2) mask(nx,ny), umask(2:nx,ny), vmask(nx,2:ny)
      do j=1,ny
        js2=max(1,j-2) ; js=max(1,j-1) ; jn=min(j+1,ny)
        do i=1,nx
          iw=max(1,i-1) ; ie=min(i+1,nx)
          if (rmsk(i,j) < 0.5D0 .and. ( rmsk(ie,j) > 0.5D0
     &       .or. rmsk(i,jn) > 0.5D0 .or. rmsk(iw,j) > 0.5D0
     &                       .or. rmsk(i,js) > 0.5D0 )) then
            mask(i,j)=0
          else
            mask(i,j)=1
          endif
        enddo
        do i=2,nx
          iw2=max(1,i-2) ; iw=max(1,i-1) ; ie=min(i+1,nx)
          if ( rmsk(i-1,j) < 0.5D0 .and. rmsk(i,j) < 0.5D0 .and.
     &        (  rmsk(i-1,jn) > 0.5D0 .or. rmsk(i,jn) > 0.5D0  .or.
     &       rmsk(iw2,j) > 0.5D0      .or.     rmsk(ie,j) > 0.5D0 .or.
     &           rmsk(i-1,js) > 0.5D0 .or. rmsk(i,js) > 0.5D0 )) then
            umask(i,j)=0
          else
            umask(i,j)=1
          endif
        enddo
        if (j > 1) then
          do i=1,nx
            iw=max(1,i-1) ; ie=min(i+1,nx)
            if ( rmsk(i,j-1) < 0.5D0 .and. rmsk(i,j) < 0.5D0 .and.
     &          (             rmsk(i,jn) > 0.5D0 .or.
     &           rmsk(iw,j  ) > 0.5D0 .or. rmsk(ie,j  ) > 0.5D0 .or.
     &           rmsk(iw,j-1) > 0.5D0 .or. rmsk(ie,j-1) > 0.5D0 .or.
     &                        rmsk(i,js2) > 0.5D0 )) then
              vmask(i,j)=0
            else
              vmask(i,j)=1
            endif
          enddo
        endif
      enddo
      end
      subroutine apply_mask(nx,ny,N, mask,q)
      integer(kind=4) nx,ny,N, i,j,k
      integer(kind=2) mask(nx,ny)
      real (kind=4) q(nx,ny,N)
      real(kind=4), parameter :: FillValue=1.D+33
      real(kind=8), parameter :: spv=1.D+9, spv_set=1.D+33
      if (N > 1) then
        do j=1,ny
          do i=1,nx
            if (mask(i,j) == 0) then
              do k=1,N
                q(i,j,k)=FillValue
              enddo
            endif
          enddo
        enddo
      else
        k=1
        do j=1,ny
          do i=1,nx
            if (mask(i,j) == 0) q(i,j,k)=FillValue
          enddo
        enddo
      endif
      end
      subroutine bry_apply_mask(nx,N, mask,q)
      integer(kind=4) nx,N, i,k
      integer(kind=2) mask(nx)
      real (kind=4) q(nx,N)
      real(kind=4), parameter :: FillValue=1.D+33
      real(kind=8), parameter :: spv=1.D+9, spv_set=1.D+33
      if (N > 1) then
        do i=1,nx
          if (mask(i) == 0) then
            do k=1,N
              q(i,k)=FillValue
            enddo
          endif
        enddo
      else
        k=1
        do i=1,nx
          if (mask(i) == 0) q(i,k)=FillValue
        enddo
      endif
      end
      subroutine compute_min_max(nx,ny, x, xmin,xmax)
      implicit none
      integer(kind=4) nx,ny, i,j
      real(kind=8) x(nx,ny), xmin,xmax
      xmin=x(1,1) ; xmax=x(1,1)
      do j=1,ny
        do i=1,nx
          if (x(i,j) < xmin) then
            xmin=x(i,j)
          elseif (x(i,j) > xmax) then
            xmax=x(i,j)
          endif
        enddo
      enddo
      write(*,*) 'min,max =', xmin,xmax
      end
      subroutine adjust_lon_into_range(nx,ny, x, west,east)
      implicit none
      integer(kind=4) nx,ny, i,j
      real(kind=8) x(nx,ny), west,east, xmin,xmax
      xmin=x(1,1) ; xmax=x(1,1)
      do j=1,ny
        do i=1,nx
          if (x(i,j) < xmin) then
            xmin=x(i,j)
          elseif (x(i,j) > xmax) then
            xmax=x(i,j)
          endif
        enddo
      enddo
      write(*,*) 'min,max =', xmin,xmax
      if (xmin > east) then
        do j=1,ny
          do i=1,nx
            x(i,j)=x(i,j)-360.D0
          enddo
        enddo
        write(*,*) 'adjusted to', xmin-360.D0,xmax-360.D0
      elseif (xmax < west) then
        do j=1,ny
          do i=1,nx
            x(i,j)=x(i,j)+360.D0
          enddo
        enddo
        write(*,*) 'adjusted to', xmin+360.D0, xmax+360.D0
      else
        write(*,*) 'no need to adjust parent-grid longitude'
      endif
      end
      subroutine compute_index_bounds(ncx,ncy, ip,jp, mask, imin,imax,
     &                                                      jmin,jmax)
      implicit none
      integer(kind=4) ncx,ncy,ip(ncx,ncy),jp(ncx,ncy), 
     &                      imin,imax,jmin,jmax,i,j
      integer(kind=2) mask(ncx,ncy)
      imin=10000000 ; imax=-1
      jmin=10000000 ; jmax=-1
      do j=1,ncy
        do i=1,ncx
          if (mask(i,j) > 0 .and. ip(i,j) > 0 .and. jp(i,j) > 0) then
            imin=min(imin,ip(i,j)) ; jmin=min(jmin,jp(i,j))
            imax=max(imax,ip(i,j)) ; jmax=max(jmax,jp(i,j))
          endif
        enddo
      enddo
      imax=imax+1 ; jmax=jmax+1
      end
      subroutine r2r_bry_index_bounds(ncx, ip,jp, mask, imin,imax,
     &                                                  jmin,jmax)
      implicit none
      integer(kind=4) ncx, ip(ncx),jp(ncx), imin,imax,jmin,jmax, ic
      integer(kind=2) mask(ncx)
      do ic=1,ncx
        if (mask(ic) > 0 .and. ip(ic) > 0 .and. jp(ic) > 0) then
          imin=min(imin,ip(ic)) ; jmin=min(jmin,jp(ic))
          imax=max(imax,ip(ic)) ; jmax=max(jmax,jp(ic))
        endif
      enddo
      end
      subroutine r2r_init_diag_file(ncx,ncy, N, ip,jp,xi,eta,
     &                     ipu,jpu,xiu,etau, ipv,jpv,xiv,etav,
     &                              csA,snA, h, hprnt, kprnt)
      implicit none
      integer(kind=4) ncx,ncy, N
      integer(kind=4), dimension(ncx,ncy) :: ip,jp, ipu,jpu, ipv,jpv
      real(kind=8), dimension(ncx,ncy)  :: csA,snA, h,hprnt, xi,eta,
     &                                             xiu,etau, xiv,etav
      real(kind=8), dimension(ncx,ncy,N) :: kprnt
      integer(kind=4) ncid, ierr,i,j,dimids(4), ipvar, jpvar,  xivar, 
     &                              etavar,
     &            old_fill_mode, csvar, ipuvar,jpuvar, xiuvar,etauvar,
     &                           snvar, ipvvar,jpvvar, xivvar,etavvar,
     &                                  hvar,  hpvar,  dhvar, kpvar
      real(kind=8), allocatable, dimension(:,:) :: dh
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
      write(*,'(1x,2A)',advance='no') 'creating diagnostic file ',
     &                                  '''r2r_init_diag.nc''...'
      ierr=nf_create('r2r_init_diag.nc', nf_64bit_offset, ncid)
      ierr=nf_set_fill(ncid, nf_nofill, old_fill_mode)
      ierr=nf_def_dim(ncid, 'xi_rho',  ncx, dimids(1))
      ierr=nf_def_dim(ncid, 'eta_rho', ncy, dimids(2))
      ierr=nf_def_dim(ncid, 's_rho',   N,   dimids(3))
      ierr=nf_def_var(ncid, 'ip', nf_int, 2, dimids, ipvar)
      ierr=nf_put_att_text(ncid, ipvar, 'long_name', 37,
     &                    'parent-to-child interpolation i-index')
      ierr=nf_put_att_int(ncid, ipvar, '_FillValue', nf_int, 1,-1)
      ierr=nf_def_var(ncid, 'jp', nf_int, 2, dimids, jpvar)
      ierr=nf_put_att_text(ncid, jpvar, 'long_name', 37,
     &                    'parent-to-child interpolation j-index')
      ierr=nf_put_att_int(ncid, jpvar, '_FillValue', nf_int, 1,-1)
      ierr=nf_def_var(ncid, 'xi',  nf_double, 2, dimids,  xivar)
      ierr=nf_put_att_text(ncid, xivar, 'long_name', 38,
     &                  'parent-to-child interpolation offset X')
      ierr=nf_put_att_double(ncid, xivar, '_FillValue', nf_double,
     &                                                  1, -1.D0)
      ierr=nf_def_var(ncid, 'eta', nf_double, 2, dimids, etavar)
      ierr=nf_put_att_text(ncid, etavar, 'long_name', 38,
     &                  'parent-to-child interpolation offset Y')
      ierr=nf_put_att_double(ncid, etavar,'_FillValue', nf_double,
     &                                                  1, -1.D0)
      ierr=nf_def_var(ncid, 'ipu', nf_int, 2, dimids, ipuvar)
      ierr=nf_put_att_text(ncid, ipuvar, 'long_name', 39,
     &                   'parent-to-child U-interpolation i-index')
      ierr=nf_put_att_int(ncid, ipuvar, '_FillValue', nf_int, 1,-1)
      ierr=nf_def_var(ncid, 'jpu', nf_int, 2, dimids, jpuvar)
      ierr=nf_put_att_text(ncid, jpuvar, 'long_name', 39,
     &                   'parent-to-child U-interpolation j-index')
      ierr=nf_put_att_int(ncid, jpuvar, '_FillValue', nf_int, 1,-1)
      ierr=nf_def_var(ncid, 'xiu',  nf_double, 2, dimids,  xiuvar)
      ierr=nf_put_att_text(ncid, xiuvar, 'long_name', 40,
     &                 'parent-to-child U-interpolation offset X')
      ierr=nf_put_att_double(ncid, xiuvar, '_FillValue', nf_double,
     &                                                   1, -1.D0)
      ierr=nf_def_var(ncid, 'etau', nf_double, 2, dimids, etauvar)
      ierr=nf_put_att_text(ncid, etauvar, 'long_name', 40,
     &                  'parent-to-child U-interpolation offset Y')
      ierr=nf_put_att_double(ncid, etauvar,'_FillValue', nf_double,
     &                                                    1, -1.D0)
      ierr=nf_def_var(ncid, 'ipv', nf_int, 2, dimids, ipvvar)
      ierr=nf_put_att_text(ncid, ipvvar, 'long_name', 39,
     &                   'parent-to-child V-interpolation i-index')
      ierr=nf_put_att_int(ncid, ipvvar, '_FillValue', nf_int, 1,-1)
      ierr=nf_def_var(ncid, 'jpv', nf_int, 2, dimids, jpvvar)
      ierr=nf_put_att_text(ncid, jpvvar, 'long_name', 39,
     &                   'parent-to-child V-interpolation j-index')
      ierr=nf_put_att_int(ncid, jpvvar, '_FillValue', nf_int, 1,-1)
      ierr=nf_def_var(ncid, 'xiv',  nf_double, 2, dimids,  xivvar)
      ierr=nf_put_att_text(ncid, xivvar, 'long_name', 40,
     &                 'parent-to-child V-interpolation offset X')
      ierr=nf_put_att_double(ncid, xivvar, '_FillValue', nf_double,
     &                                                   1, -1.D0)
      ierr=nf_def_var(ncid, 'etav', nf_double, 2, dimids, etavvar)
      ierr=nf_put_att_text(ncid, etauvar, 'long_name', 40,
     &                  'parent-to-child V-interpolation offset Y')
      ierr=nf_put_att_double(ncid, etavvar,'_FillValue', nf_double,
     &                                                    1, -1.D0)
      ierr=nf_def_var(ncid, 'csA', nf_double, 2, dimids, csvar)
      ierr=nf_put_att_text(ncid, csvar, 'long_name', 44,
     &          'cosine of child-parent grid angle difference')
      ierr=nf_def_var(ncid, 'snA', nf_double, 2, dimids, snvar)
      ierr=nf_put_att_text(ncid, snvar, 'long_name', 42,
     &            'sine of child-parent grid angle difference')
      ierr=nf_def_var(ncid, 'h', nf_double, 2, dimids, hvar)
      ierr=nf_put_att_text(ncid, hvar, 'long_name', 35,
     &                     'child grid native bottom topography')
      ierr=nf_def_var(ncid, 'hprnt', nf_double, 2, dimids, hpvar)
      ierr=nf_put_att_text(ncid, hpvar, 'long_name', 35,
     &                     'interpolated parent grid topography')
      ierr=nf_def_var(ncid, 'dh', nf_double, 2, dimids, dhvar)
      ierr=nf_put_att_text(ncid, dhvar, 'long_name', 34,
     &                     'child-parent topography difference')
      ierr=nf_def_var(ncid, 'kprnt', nf_real, 3, dimids, kpvar)
      ierr=nf_put_att_text(ncid, kpvar, 'long_name', 56,
     & 'index coordinate for parent-child vertical interpolation')
      ierr=nf_enddef(ncid)
      ierr=nf_put_var_int(ncid, ipvar, ip)
      ierr=nf_put_var_int(ncid, jpvar, jp)
      ierr=nf_put_var_double(ncid, xivar,  xi)
      ierr=nf_put_var_double(ncid, etavar, eta)
      ierr=nf_put_var_int(ncid, ipuvar, ipu)
      ierr=nf_put_var_int(ncid, jpuvar, jpu)
      ierr=nf_put_var_double(ncid, xiuvar,  xiu)
      ierr=nf_put_var_double(ncid, etauvar, etau)
      ierr=nf_put_var_int(ncid, ipvvar, ipv)
      ierr=nf_put_var_int(ncid, jpvvar, jpv)
      ierr=nf_put_var_double(ncid, xivvar,  xiv)
      ierr=nf_put_var_double(ncid, etavvar, etav)
      ierr=nf_put_var_double(ncid, csvar, csA)
      ierr=nf_put_var_double(ncid, snvar, snA)
      ierr=nf_put_var_double(ncid, hvar,  h)
      ierr=nf_put_var_double(ncid, hpvar, hprnt)
      ierr=nf_put_var_double(ncid, kpvar, kprnt)
      allocate(dh(ncx,ncy))
      do j=1,ncy
        do i=1,ncx
          dh(i,j)=h(i,j)-hprnt(i,j)
        enddo
      enddo
      ierr=nf_put_var_double(ncid, dhvar, dh)
      deallocate(dh)
      ierr=nf_close(ncid)
      write(*,'(2x,A)') '...done'
      end
