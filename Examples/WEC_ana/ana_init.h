      integer :: i,j,k
      real    :: Clw,Cg,cff

      Clw=sqrt(g*h(1,1))
      Cg=0.5*sqrt(g/(2.0*pi))
      cff=-(Clw**2)/(Clw**2-Cg**2)

      do j=-1,ny+2
        do i=-1,nx+2
          zeta(i,j,1)=0.0
          zeta(i,j,2)=0.0
          ubar(i,j,1)=cff*ust2d(i,j)
          ubar(i,j,2)=ubar(i,j,1)
          vbar(i,j,1)=cff*vst2d(i,j)
          vbar(i,j,2)=vbar(i,j,1)
          do k=1,nz
            t(i,j,k,1,itemp)=18.0
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            u(i,j,k,1)=ubar(i,j,1)
            u(i,j,k,2)=u(i,j,k,1)
            v(i,j,k,1)=vbar(i,j,1)
            v(i,j,k,2)=v(i,j,k,1)
          enddo
        enddo
      enddo

      do j=-1,ny+2
        do i=-1,nx+1
          zeta(i,j,1)=cff*0.5*(ust2d(i,j)+ust2d(i+1,j))*Cg/g
          zeta(i,j,2)=zeta(i,j,1)
        enddo
        zeta(nx+2,j,1)=zeta(nx+1,j,1)
        zeta(nx+2,j,2)=zeta(nx+2,j,1)
      enddo
