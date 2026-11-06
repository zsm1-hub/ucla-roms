      module r2r_vertint_vars
       integer(kind=4), save :: trd_count=0
       real(kind=8), save :: vert_int_error=0._8
       integer(kind=4), save :: alloc_zc_size=0, allc_zpr_size=0
C$OMP THREADPRIVATE(allc_zpr_size, alloc_zc_size)
       real(kind=8), allocatable, dimension(:) :: zc,zp_r,drv
C$OMP THREADPRIVATE(zc,zp_r,drv)
      end module r2r_vertint_vars
      subroutine r2r_init_vertint_thread(ncx,ncy, hprnt, Np,hcp,Csp_r,
     &                                            h, N,hc,Cs_r, kprnt)
      use r2r_vertint_vars
      implicit none
      integer(kind=4) ncx,ncy, Np,N
      real(kind=8) hprnt(ncx,ncy),h(ncx,ncy),hcp,hc,Csp_r(Np),Cs_r(N),
     &                                                kprnt(ncx,ncy,N)
      integer(kind=4) ntrds,trd, nsub_x,nsub_y, ntls, isize,
     &        my_first,my_last,tile, istr,iend,jstr,jend
C$    integer omp_get_thread_num, omp_get_num_threads
      ntrds=1 ; trd=0
C$    ntrds=omp_get_num_threads() ; trd=omp_get_thread_num()
      call set_tiles(ncx,ncy, nsub_x,nsub_y)
      ntls=(nsub_x*nsub_y+ntrds-1)/ntrds
      my_first=trd*ntls -(ntls*ntrds-nsub_x*nsub_y)/2
      my_last=min(my_first+ntls-1, nsub_x*nsub_y-1)
      my_first=max(my_first, 0)
      isize=(ncx+nsub_x-1)/nsub_x
      if (alloc_zc_size < isize*(N+1) .or.
     &    allc_zpr_size < isize*(Np+2)) then
        alloc_zc_size=isize*(N+1); allc_zpr_size=isize*(Np+2)
        if (allocated(zp_r)) deallocate(zp_r,drv,zc)
        allocate( zc(alloc_zc_size), zp_r(allc_zpr_size),
     &                                drv(allc_zpr_size) )
C$OMP CRITICAL(r2r_vert_crgn)
        write(*,'(1x,2A,F8.4,1x,A,I3)')  'r2r_init_vertint_thread: ',
     & 'allocated', dble(2*allc_zpr_size+alloc_zc_size)/dble(262144),
     &                           'MB private workspace, trd =', trd
C$OMP END CRITICAL(r2r_vert_crgn)
      elseif (trd == 0) then
        write(*,'(2x,A)',advance='no')
     & 'entering r2r_init_vertint_thread ...'
      endif
      do tile=my_first,my_last,+1
        call comp_tile_bounds( tile, ncx,ncy, nsub_x,nsub_y,
     &                                    istr,iend,jstr,jend )
        call r2r_init_vrtint_tile( istr,iend,jstr,jend, ncx,ncy,
     &                              hprnt,Np,hcp,Csp_r,zp_r,drv,
     &                                 h, N,hc,Cs_r,zc, kprnt)
      enddo
      if (trd == 0) write(*,'(2x,A)')
     & 'leaving r2r_init_vertint_thread'
C$OMP BARRIER
      end
      subroutine r2r_check_vertint_thread(ncx,ncy, hprnt, Np,hcp,Csp_r,
     &                                             h, N,hc,Cs_r, kprnt)
      use r2r_vertint_vars
      implicit none
      integer(kind=4) ncx,ncy, Np,N
      real(kind=8) hprnt(ncx,ncy),h(ncx,ncy), hcp,hc,Csp_r(Np),Cs_r(N),
     &                                        kprnt(ncx,ncy,N),my_error
      integer(kind=4) ntrds,trd, nsub_x,nsub_y, ntls, isize,
     &        my_first,my_last, tile, istr,iend,jstr,jend
C$    integer omp_get_thread_num, omp_get_num_threads
      ntrds=1 ; trd=0
C$    ntrds=omp_get_num_threads() ; trd=omp_get_thread_num()
      call set_tiles(ncx,ncy, nsub_x,nsub_y)
      ntls=(nsub_x*nsub_y+ntrds-1)/ntrds
      my_first=trd*ntls -(ntls*ntrds-nsub_x*nsub_y)/2
      my_last=min(my_first+ntls-1, nsub_x*nsub_y-1)
      my_first=max(my_first, 0)
      isize=(ncx+nsub_x-1)/nsub_x
      if (alloc_zc_size < isize*(N+1) .or.
     &    allc_zpr_size < isize*(Np+2)) then
        alloc_zc_size=isize*(N+1); allc_zpr_size=isize*(Np+2)
        if (allocated(zp_r)) deallocate(zp_r,drv,zc)
        allocate(zc(alloc_zc_size), zp_r(allc_zpr_size),
     &                               drv(allc_zpr_size))
      endif
      my_error=0.D0
      do tile=my_last,my_first,-1
        call comp_tile_bounds( tile, ncx,ncy, nsub_x,nsub_y,
     &                                    istr,iend,jstr,jend )
        call r2r_check_vrtint_tile(istr,iend,jstr,jend, ncx,ncy,
     &                             hprnt, Np,hcp,Csp_r,zp_r,drv,
     &                        h, N,hc,Cs_r,zc, kprnt, my_error)
      enddo
C$OMP CRITICAL(r2r_vert_crgn)
      if (trd_count == 0) vert_int_error=0.D0
      trd_count=trd_count+1
      vert_int_error=max(vert_int_error, my_error)
      if (trd_count == ntrds) then
        trd_count=0
        write(*,*)
        write(*,*) 'maximum vert_int_error =', vert_int_error
        write(*,*)
      endif
C$OMP END CRITICAL(r2r_vert_crgn)
      end
      subroutine r2r_init_vrtint_tile( istr,iend,jstr,jend, ncx,ncy,
     &                                 hprnt,Np,hcp,Csp_r, zp_r,drv,
     &                                     h,N, hc, Cs_r,z_r, kprnt)
      implicit none
      integer(kind=4) istr,iend,jstr,jend, ncx,ncy, Np,N
      real(kind=8) hprnt(ncx,ncy), h(ncx,ncy), Csp_r(Np), Cs_r(N),
     &             zp_r(istr:iend,0:Np+1), z_r(istr:iend,N), hcp,hc,
     &              drv(istr:iend,0:Np+1), kprnt(ncx,ncy,N), p
      integer(kind=4) i,j,k, kp
      integer(kind=4) iter
      real(kind=8) q,pq, zerr, dZds
     &                 , dZs2, d2Zds2
      do j=jstr,jend
        call r2r_set_depth_tile(istr,iend, j,j, ncx,ncy,
     &                               h, N, hc,Cs_r, z_r)
        call r2r_set_depth_tile(istr,iend, j,j, ncx,ncy,
     &               hprnt, Np,hcp,Csp_r, zp_r(istr,1) )
        do i=istr,iend
          zp_r(i,Np+1)=0.D0
          zp_r(i,0)=-hprnt(i,j)
        enddo
        call comp_zspline_tile(istr,iend, Np, zp_r,drv)
        do i=istr,iend
          kp=Np
          do k=N,1,-1
            do while(z_r(i,k) < zp_r(i,kp) .and. kp > 0)
              kp=kp-1
            enddo
            kprnt(i,j,k)=dble(kp)+0.5D0
          enddo
        enddo
        do k=1,N
          do i=istr,iend
            kp=int(kprnt(i,j,k))
            p=(z_r(i,k)-zp_r(i,kp))/(zp_r(i,kp+1)-zp_r(i,kp))
            kprnt(i,j,k)=dble(kp)+p
          enddo
          do i=istr,iend
            if (kprnt(i,j,k) > 0.499999999999D0) then
              kp=int(kprnt(i,j,k)) ; p=kprnt(i,j,k)-dble(kp)
              do iter=1,8
                q=1.D0-p ; pq=p*q
                zerr=p*((p+2.D0*pq)*zp_r(i,kp+1) -pq*drv(i,kp+1))
     &              +q*((q+2.D0*pq)*zp_r(i,kp  ) +pq*drv(i,kp  ))
     &                                                  -z_r(i,k)
                dZds=6.D0*pq*(zp_r(i,kp+1)-zp_r(i,kp))
     &                           +(p*p-2.D0*pq)*drv(i,kp+1)
     &                           +(q*q-2.D0*pq)*drv(i,kp  )
                d2Zds2=(2.D0*p-q)*drv(i,kp+1) -(2.D0*q-p)*drv(i,kp)
     &                        -3.D0*(p-q)*(zp_r(i,kp+1)-zp_r(i,kp))
                dZs2=dZds*dZds
                p=p-zerr*(dZs2+d2Zds2*zerr)/(dZs2*dZds)
              enddo
              kprnt(i,j,k)=dble(kp)+p
            else
              kprnt(i,j,k)=0.D0
            endif
          enddo
        enddo
      enddo
      end
      subroutine r2r_check_vrtint_tile(istr,iend,jstr,jend, ncx,ncy,
     &                                 hprnt, Np,hcp,Csp_r, zp_r,drv,
     &                            h, N,hc,Cs_r,z_r, kprnt, my_error)
      implicit none
      integer(kind=4) istr,iend,jstr,jend, ncx,ncy, Np,N, i,j,k, kp
      real(kind=8) hprnt(ncx,ncy), h(ncx,ncy), Csp_r(Np), Cs_r(N),
     &             zp_r(istr:iend,0:Np+1), z_r(istr:iend,N), hcp,hc,
     &             drv(istr:iend,0:Np+1), kprnt(ncx,ncy,N), my_error
     &           , p
     &           , q, pq
      do j=jstr,jend
        call r2r_set_depth_tile(istr,iend, j,j, ncx,ncy,
     &                               h, N, hc,Cs_r, z_r)
        call r2r_set_depth_tile(istr,iend, j,j, ncx,ncy,
     &               hprnt, Np,hcp,Csp_r, zp_r(istr,1) )
        do i=istr,iend
          zp_r(i,Np+1)=0.D0
          zp_r(i,0)=-hprnt(i,j)
        enddo
        call comp_zspline_tile(istr,iend, Np, zp_r,drv)
        do i=istr,iend
          do k=1,N
            kp=max(int(kprnt(i,j,k)), 0) ; p=kprnt(i,j,k)-float(kp)
            if (kprnt(i,j,k) > 0.499999999999D0) then
              q=1.D0-p ; pq=p*q
              my_error=max( my_error, abs(
     &             p*((p+2.D0*pq)*zp_r(i,kp+1) -pq*drv(i,kp+1))
     &            +q*((q+2.D0*pq)*zp_r(i,kp  ) +pq*drv(i,kp  ))
     &                                            -z_r(i,k) ) )
            endif
          enddo
        enddo
      enddo
      end
      subroutine comp_zspline_tile(istr,iend, Np, zp_r,drv)
      implicit none
      integer(kind=4) istr,iend, Np
      real(kind=8) zp_r(istr:iend,0:Np+1), drv(istr:iend,0:Np+1)
      integer(kind=4) i,k
      real(kind=8) CF(1:Np-1), cff,cff1
      CF(1)=1.D0/3.D0 ; cff=8.D0/9.D0
      do i=istr,iend
        drv(i,1)=cff*(zp_r(i,2) -zp_r(i,0))
      enddo
      do k=2,Np-1,+1
        cff=1.D0/(4.D0-CF(k-1)) ; CF(k)=cff
        do i=istr,iend
          drv(i,k)=cff*( 3.D0*(zp_r(i,k+1)-zp_r(i,k-1)) -drv(i,k-1))
        enddo
      enddo
      cff=1.D0/(3.D0-CF(Np-1)) ; cff1=8.D0/3.D0
      do i=istr,iend
        drv(i,Np)=cff*( cff1*(zp_r(i,Np+1)-zp_r(i,Np-1)) -drv(i,Np-1) )
      enddo
      do k=Np-1,1,-1
        cff=CF(k)
        do i=istr,iend
          drv(i,k)=drv(i,k)-cff*drv(i,k+1)
        enddo
      enddo
      do i=istr,iend
        drv(i,0)=2.D0*drv(i,2)-drv(i,1) -8.D0*zp_r(i,0)+12.D0*zp_r(i,1)
     &                                                  -4.D0*zp_r(i,2)
        zp_r(i,0)=4.D0*zp_r(i,0)-4.D0*zp_r(i,1)+zp_r(i,2)
     &                         -0.5D0*(drv(i,2)-drv(i,1))
        drv(i,Np+1)=2.D0*drv(i,Np-1)-drv(i,Np)  +8.D0*zp_r(i,Np+1)
     &                          -12.D0*zp_r(i,Np)+4.D0*zp_r(i,Np-1)
        zp_r(i,Np+1)=4.D0*zp_r(i,Np+1)-4.D0*zp_r(i,Np)+zp_r(i,Np-1)
     &                               -0.5D0*(drv(i,Np)-drv(i,Np-1))
      enddo
      end
      subroutine r2r_vrtint_thread(ncx,ncy, lmsk,mask, btm_bc,
     &                                  Np,qsrc, N,kprnt,qtr)
      use r2r_vertint_vars
      implicit none
      integer(kind=4) ncx,ncy, lmsk, btm_bc, Np,N
      integer(kind=2) mask(ncx,ncy)
      real(kind=8) kprnt(ncx,ncy,N)
      real(kind=4) qsrc(ncx,ncy,Np), qtr(ncx,ncy,N)
      integer(kind=4) ntrds,trd, nsub_x,nsub_y, ntls, my_first,my_last,
     &                                tile, istr,iend,jstr,jend
C$    integer omp_get_thread_num, omp_get_num_threads
      ntrds=1 ; trd=0
C$    ntrds=omp_get_num_threads() ; trd=omp_get_thread_num()
      call set_tiles(ncx,ncy, nsub_x,nsub_y)
      ntls=(nsub_x*nsub_y+ntrds-1)/ntrds
      my_first=trd*ntls -(ntls*ntrds-nsub_x*nsub_y)/2
      my_last=min(my_first+ntls-1, nsub_x*nsub_y-1)
      my_first=max(my_first, 0)
      do tile=my_first,my_last,+1
        call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     &                                  istr,iend,jstr,jend)
        call r2r_vsplnint_tile( istr,iend,jstr,jend, ncx,ncy,
     &           lmsk,mask, btm_bc, Np,qsrc,drv, N,kprnt,qtr)
      enddo
      end
      subroutine r2r_vretint_tile(istr,iend,jstr,jend, ncx,ncy,
     &                         lmsk,mask, Np,qsrc, N,kprnt,qtr)
      implicit none
      integer(kind=4) istr,iend,jstr,jend, ncx,ncy, lmsk,Np,N, i,j,k,kp
      integer(kind=2) mask(ncx,ncy)
      real(kind=8) kprnt(ncx,ncy,N), p
      real(kind=4) qsrc(ncx,ncy,Np), qtr(ncx,ncy,N)
      real(kind=4), parameter :: FillValue=1.D+33
      real(kind=8), parameter :: spv=1.D+9, spv_set=1.D+33
      do j=jstr,jend
        do k=1,N
          do i=istr,iend
            kp=int(kprnt(i,j,k)) ; p=kprnt(i,j,k)-float(kp)
            if (kp > Np-1) then
              qtr(i,j,k)=qsrc(i,j,Np)
            elseif (kp > 0) then
              qtr(i,j,k)=p*qsrc(i,j,kp+1)+(1.D0-p)*qsrc(i,j,kp)
            else
              qtr(i,j,k)=qsrc(i,j,1)
            endif
          enddo
        enddo
        if (lmsk == 1) then
          do i=istr,iend
            if (mask(i,j) == 0) then
              do k=1,N
                qtr(i,j,k)=FillValue
              enddo
            endif
          enddo
        endif
      enddo
      end
      subroutine r2r_vsplnint_tile(istr,iend,jstr,jend, ncx,ncy, lmsk,
     &                          mask, btm_bc,Np,qsrc,drv, N,kprnt,qtr)
      implicit none
      integer(kind=4) istr,iend,jstr,jend, ncx,ncy, lmsk,btm_bc,Np,N, 
     &                              i,j,k,kp
      integer(kind=2) mask(ncx,ncy)
      real(kind=4) qsrc(ncx,ncy,Np), qtr(ncx,ncy,N)
      real(kind=8) kprnt(ncx,ncy,N), CF(Np),  drv(istr:iend,Np),
     &                                              cff, p,q,pq
      real(kind=4), parameter :: FillValue=1.D+33
      real(kind=8), parameter :: spv=1.D+9, spv_set=1.D+33
      do j=jstr,jend
        if (btm_bc == -1) then
          CF(1)=1.D0/3.D0 ; cff=8.D0/9.D0
          do i=istr,iend
            drv(i,1)=cff*qsrc(i,j,2)
          enddo
        elseif (btm_bc == +1) then
          CF(1)=7.D0/15.D0 ; cff=1.2D0
          do i=istr,iend
            drv(i,1)=cff*(qsrc(i,j,2)-qsrc(i,j,1))
          enddo
        else
          CF(1)=0.5D0
          do i=istr,iend
            drv(i,1)=1.5D0*(qsrc(i,j,2)-qsrc(i,j,1))
          enddo
        endif
        do k=2,Np-1,+1
          CF(k)=1.D0/(4.D0-CF(k-1)) ; cff=CF(k)
          do i=istr,iend
            drv(i,k)=cff*( 3.D0*(qsrc(i,j,k+1)-qsrc(i,j,k-1))
     &                                           -drv(i,k-1) )
          enddo
        enddo
        cff=1.D0/(2.D0-CF(Np-1))
        do i=istr,iend
          drv(i,Np)=cff*( 3.D0*(qsrc(i,j,Np)-qsrc(i,j,Np-1))
     &                                         -drv(i,Np-1) )
        enddo
        do k=Np-1,1,-1
          cff=CF(k)
          do i=istr,iend
            drv(i,k)=drv(i,k)-cff*drv(i,k+1)
          enddo
        enddo
        do k=1,N
          do i=istr,iend
            kp=int(kprnt(i,j,k))
            kp=max( 1, int(kprnt(i,j,k)) )
            p=max(-0.5D0, kprnt(i,j,k)-float(kp) )
            if (kp > Np-1) then
              qtr(i,j,k)=qsrc(i,j,Np)+p*drv(i,Np)
            elseif (kp > 0) then
              q=1.D0-p ; pq=p*q
              qtr(i,j,k)=p*((p+2.D0*pq)*qsrc(i,j,kp+1) -pq*drv(i,kp+1))
     &                  +q*((q+2.D0*pq)*qsrc(i,j,kp  ) +pq*drv(i,kp  ))
            endif
          enddo
        enddo
        if (lmsk == 1) then
          do i=istr,iend
            if (mask(i,j) == 0) then
              do k=1,N
                qtr(i,j,k)=FillValue
              enddo
            endif
          enddo
        endif
      enddo
      end
      subroutine r2r_set_depth_tile(istr,iend,jstr,jend, ncx,ncy,
     &                                             h, N,hc,Cs, z)
      implicit none
      integer(kind=4) istr,iend,jstr,jend, ncx,ncy,N, i,j,k
      real(kind=8) h(ncx,ncy), hc, Cs(N), z(istr:iend,jstr:jend,N),
     &                                           hcds,bias,cf1,cf2
      if (Cs(1) > -0.999999D0) then
        hcds=hc/dble(N) ; bias=0.5D0
      else
        hcds=hc/dble(N-1) ; bias=0.D0
      endif
      do j=jstr,jend
        do i=istr,iend
          z(i,j,N)=h(i,j)/(hc+h(i,j))
        enddo
        do k=1,N
          cf1=hcds*(dble(k-N)-bias) ; cf2=Cs(k)
          do i=istr,iend
            z(i,j,k)= z(i,j,N)*(cf1 + cf2*h(i,j))
          enddo
        enddo
      enddo
      end
      module comp_uvbar_vars
      integer(kind=4), save :: alloc_zw_size=0
C$OMP THREADPRIVATE(alloc_zw_size)
      real(kind=8), allocatable, dimension(:) :: D,z_w,Hz
C$OMP THREADPRIVATE(D,z_w,Hz)
      end module comp_uvbar_vars
      subroutine compute_uvbar_thread(ncx,ncy, N, hc,Cs_w, h, u,v,
     &                                                  ubar,vbar)
      use comp_uvbar_vars
      implicit none
      integer(kind=4) ncx,ncy, N
      real(kind=8) hc, Cs_w(0:N),  h(ncx,ncy)
      real(kind=4) u(2:ncx,ncy,N), ubar(2:ncx,ncy),
     &             v(ncx,2:ncy,N), vbar(ncx,2:ncy)
      integer(kind=4) ntrds,trd, nsub_x,nsub_y, tile,my_first,my_last,
     &        ntls,  isize,  istr,iend,jstr,jend, istrU,jstrV
C$    integer omp_get_thread_num, omp_get_num_threads
      ntrds=1 ; trd=0
C$    ntrds=omp_get_num_threads() ; trd=omp_get_thread_num()
      call set_tiles(ncx,ncy, nsub_x,nsub_y)
      ntls=(nsub_x*nsub_y+ntrds-1)/ntrds
      my_first=trd*ntls -(ntls*ntrds-nsub_x*nsub_y)/2
      my_last=min(my_first+ntls-1, nsub_x*nsub_y-1)
      my_first=max(my_first, 0)
      isize=(ncx+nsub_x-1)/nsub_x
      if (alloc_zw_size < (isize+2)*(N+1)) then
        alloc_zw_size=(isize+2)*(N+1)
        if (allocated(z_w)) deallocate(Hz,z_w,D)
        allocate(D(isize+2), z_w(alloc_zw_size), Hz(2*alloc_zw_size))
C$OMP CRITICAL(comp_uvbar_crgn)
        write(*,'(1x,2A,F8.4,1x,A,I3)')  'compute_uvbar_thread: ',
     &   'allocated', (2*(isize+1)+6*alloc_zw_size)/float(262144),
     &                        'MB private workspace, trd =', trd
C$OMP END CRITICAL(comp_uvbar_crgn)
      endif
      do tile=my_last,my_first,-1
        call comp_tile_bounds( tile, ncx,ncy, nsub_x,nsub_y,
     &                                    istr,iend,jstr,jend)
        istrU=max(istr,2) ; jstrV=max(jstr,2)
        call compute_uvbar_tile(istr,iend,jstr,jend, istrU,jstrV,
     &         ncx,ncy, N, hc,Cs_w, h, u,v, ubar,vbar, D,z_w,Hz)
      enddo
      end
      subroutine compute_uvbar_tile(istr,iend,jstr,jend, istrU,jstrV,
     &              ncx,ncy, N, hc,Cs_w, h, u,v, ubar,vbar, D,z_w,Hz)
      implicit none
      integer(kind=4) istr,iend,jstr,jend, istrU,jstrV, ncx,ncy, N, 
     &                              i,j,k,j1
      real(kind=8) hc, Cs_w(0:N),  h(ncx,ncy),  cff,  D(istrU-1:iend),
     &                   z_w(istrU-1:iend,0:N),  Hz(istrU-1:iend,2,N)
      real(kind=4) u(2:ncx,ncy,N), ubar(2:ncx,ncy),
     &             v(ncx,2:ncy,N), vbar(ncx,2:ncy)
      j1=2
      do j=jstrV-1,jend,+1
        call r2r_set_depth_tile(istrU-1,iend, j,j, ncx,ncy,
     &                                h, N+1,hc,Cs_w, z_w)
        j1=3-j1
        do k=1,N
          do i=istrU-1,iend
            Hz(i,j1,k)=z_w(i,k)-z_w(i,k-1)
          enddo
        enddo
        if (j > jstr-1) then
          do i=istrU,iend
            ubar(i,j)=0._8 ; D(i)=0.D0
          enddo
          do k=N,1,-1
            do i=istrU,iend
              cff=Hz(i,j1,k)+Hz(i-1,j1,k)
              D(i)=D(i)+cff
              ubar(i,j)=ubar(i,j)+cff*u(i,j,k)
            enddo
          enddo
          do i=istrU,iend
            ubar(i,j)=ubar(i,j)/D(i)
          enddo
        endif
        if (j > jstrV-1) then
          do i=istr,iend
            vbar(i,j)=0._8 ; D(i)=0.D0
          enddo
          do k=N,1,-1
            do i=istr,iend
              cff=Hz(i,j1,k)+Hz(i,3-j1,k)
              D(i)=D(i)+cff
              vbar(i,j)=vbar(i,j)+cff*v(i,j,k)
            enddo
          enddo
          do i=istr,iend
            vbar(i,j)=vbar(i,j)/D(i)
          enddo
        endif
      enddo
      end
