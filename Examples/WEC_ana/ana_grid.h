      real    :: SizeX,SizeY
      real    :: f0,beta
      real    :: x0,y0,dx,dy
      integer :: i,j

      SizeX = 1400.0
      SizeY = 2.0
      f0 = 0.0
      beta = 0.0

      xl = SizeX
      el = SizeY
      dx = SizeX/gnx
      dy = SizeY/gny

# ifdef MPI
      x0=dx*dble(iSW_corn)
      y0=dy*dble(jSW_corn)
# else
      x0=0.0
      y0=0.0
# endif

      do j=-1,ny+2
        do i=-1,nx+2
          xp(i,j)=x0+dx*dble(i-1)
          xr(i,j)=x0+dx*(dble(i)-0.5D0)
          yp(i,j)=y0+dy*dble(j-1)
          yr(i,j)=y0+dy*(dble(j)-0.5D0)
          pm(i,j)=1.0/dx
          pn(i,j)=1.0/dy
          f(i,j)=f0+beta*(yr(i,j)-0.5*SizeY)
          h(i,j)=10.0
# ifdef MASKING
          rmask(i,j)=1.0
# endif
        enddo
      enddo
