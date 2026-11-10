      subroutine fix_outside_map(ncx,ncy,mask,xc,yc,ip,jp,xi,eta,
     &                           label)
      implicit none
      integer ncx,ncy
      integer(kind=2) mask(ncx,ncy)
      integer(kind=4) ip(ncx,ncy),jp(ncx,ncy)
      real(kind=8) xc(ncx,ncy),yc(ncx,ncy),xi(ncx,ncy),eta(ncx,ncy)
      character*(*) label
      integer i,j,ii,jj,di,dj,iter,nleft,nfixed,nwet,nland
      logical found

      nwet=0
      nland=0
      do j=1,ncy
        do i=1,ncx
          if (ip(i,j) <= 0 .or. jp(i,j) <= 0) then
            if (mask(i,j) == 1) then
              nwet=nwet+1
            else
              nland=nland+1
            endif
          endif
        enddo
      enddo
      write(*,'(1x,3A,2(A,I8))') 'fix_outside_map(',trim(label),
     &          '):',' wet outside=',nwet,' land outside=',nland

      do iter=1,ncx+ncy
        nfixed=0
        do j=1,ncy
          do i=1,ncx
            if (ip(i,j) <= 0 .or. jp(i,j) <= 0) then
              found=.false.
              do dj=-1,1
                do di=-1,1
                  if (.not.found .and. (di /= 0 .or. dj /= 0)) then
                    ii=i+di
                    jj=j+dj
                    if (ii >= 1 .and. ii <= ncx .and.
     &                  jj >= 1 .and. jj <= ncy) then
                      if (ip(ii,jj) > 0 .and. jp(ii,jj) > 0) then
                        ip(i,j)=ip(ii,jj)
                        jp(i,j)=jp(ii,jj)
                        xi(i,j)=xi(ii,jj)
                        eta(i,j)=eta(ii,jj)
                        xc(i,j)=xc(ii,jj)
                        yc(i,j)=yc(ii,jj)
                        found=.true.
                        nfixed=nfixed+1
                      endif
                    endif
                  endif
                enddo
              enddo
            endif
          enddo
        enddo
        if (nfixed == 0) exit
      enddo

      nleft=0
      do j=1,ncy
        do i=1,ncx
          if (ip(i,j) <= 0 .or. jp(i,j) <= 0) nleft=nleft+1
        enddo
      enddo
      if (nleft > 0) then
        write(*,*) '### ERROR: fix_outside_map failed for ',
     &             trim(label),', points left=',nleft
        stop
      endif
      end


      subroutine bilin2d_interp_real(nx,ny,src,ncx,ncy,
     &                               ip,jp,xi,eta,targ)
      implicit none
      integer nx,ny,ncx,ncy,ic,jc,i,j,jcmin,jcmax
      integer numthreads,trd,chunk_size
      integer(kind=4) ip(ncx,ncy),jp(ncx,ncy)
      real(kind=8) xi(ncx,ncy),eta(ncx,ncy),p,q
      real(kind=4) src(nx,ny),targ(ncx,ncy)
C$    integer omp_get_num_threads,omp_get_thread_num

      numthreads=1
      trd=0
C$    numthreads=omp_get_num_threads()
C$    trd=omp_get_thread_num()
      chunk_size=(ncy+numthreads-1)/numthreads
      jcmin=1+trd*chunk_size
      jcmax=min(ncy,jcmin+chunk_size-1)

      do jc=jcmin,jcmax
        do ic=1,ncx
          i=ip(ic,jc)
          j=jp(ic,jc)
          if (i > 0 .and. i < nx .and. j > 0 .and. j < ny) then
            p=max(0.D0,min(1.D0,xi(ic,jc)))
            q=max(0.D0,min(1.D0,eta(ic,jc)))
            targ(ic,jc)=real(
     &        (1.D0-p)*(1.D0-q)*dble(src(i  ,j  ))+
     &               p *(1.D0-q)*dble(src(i+1,j  ))+
     &        (1.D0-p)*       q *dble(src(i  ,j+1))+
     &               p *       q *dble(src(i+1,j+1)),kind=4)
          else
            write(*,*) '### ERROR: invalid interpolation index ',
     &                  ic,jc,i,j
            stop
          endif
        enddo
      enddo
      end


      subroutine r2r_vlinint_thread(ncx,ncy,lmsk,mask,Np,qsrc,
     &                               N,kprnt,qtr)
      implicit none
      integer ncx,ncy,lmsk,Np,N
      integer(kind=2) mask(ncx,ncy)
      real(kind=8) kprnt(ncx,ncy,N)
      real(kind=4) qsrc(ncx,ncy,Np),qtr(ncx,ncy,N)
      integer i,j,k,kp,jmin,jmax,numthreads,trd,chunk_size
      real(kind=8) p
#include "spval.h"
C$    integer omp_get_num_threads,omp_get_thread_num

      numthreads=1
      trd=0
C$    numthreads=omp_get_num_threads()
C$    trd=omp_get_thread_num()
      chunk_size=(ncy+numthreads-1)/numthreads
      jmin=1+trd*chunk_size
      jmax=min(ncy,jmin+chunk_size-1)

      do j=jmin,jmax
        do k=1,N
          do i=1,ncx
            if (lmsk == 1 .and. mask(i,j) == 0) then
              qtr(i,j,k)=FillValue
            else
              kp=int(kprnt(i,j,k))
              if (kp < 1) then
                qtr(i,j,k)=qsrc(i,j,1)
              elseif (kp >= Np) then
                qtr(i,j,k)=qsrc(i,j,Np)
              else
                p=max(0.D0,min(1.D0,kprnt(i,j,k)-dble(kp)))
                qtr(i,j,k)=real((1.D0-p)*dble(qsrc(i,j,kp))+
     &                                p *dble(qsrc(i,j,kp+1)),kind=4)
              endif
            endif
          enddo
        enddo
      enddo
      end


      subroutine enforce_nonnegative(nx,ny,N,mask,q,vname)
      implicit none
      integer nx,ny,N,i,j,k,nneg
      integer(kind=2) mask(nx,ny)
      real(kind=4) q(nx,ny,N),qmin
      character*(*) vname

      nneg=0
      qmin=0.
      do k=1,N
        do j=1,ny
          do i=1,nx
            if (mask(i,j) == 1 .and. q(i,j,k) < 0.) then
              qmin=min(qmin,q(i,j,k))
              q(i,j,k)=0.
              nneg=nneg+1
            endif
          enddo
        enddo
      enddo
      if (nneg > 0) write(*,'(1x,3A,I12,A,1PE12.4)')
     &    'positivity fix for ',trim(vname),': count=',nneg,
     &    ' original minimum=',qmin
      end
