      subroutine def_bry_bec2_file(fname,ncid,xi_rho,eta_rho,N,
     &                OBC_WEST,OBC_EAST,OBC_SOUTH,OBC_NORTH,
     &                        theta_s,theta_b,hc,Cs_w,Cs_r)
      implicit none
      character(len=*) fname
      integer ncid,xi_rho,eta_rho,N
      logical OBC_WEST,OBC_EAST,OBC_SOUTH,OBC_NORTH
      real(kind=8) theta_s,theta_b,hc,Cs_w(0:N),Cs_r(N)
      integer dim_xi,dim_eta,dim_s,dim_time
      integer varid,ierr,old_fill_mode,lenstr
#include "netcdf.inc"

      write(*,'(1x,3A)',advance='no')
     & 'def_bry_bec2_file :: creating ''',
     & fname(1:lenstr(fname)),''' ...'

      ierr=nf_set_fill(ncid,nf_nofill,old_fill_mode)
      ierr=nf_def_dim(ncid,'xi_rho',xi_rho,dim_xi)
      ierr=nf_def_dim(ncid,'eta_rho',eta_rho,dim_eta)
      ierr=nf_def_dim(ncid,'s_rho',N,dim_s)
      ierr=nf_def_dim(ncid,'rec_time',nf_unlimited,dim_time)

      ierr=nf_put_att_text(ncid,nf_global,'VertCoordType',4,'SM09')
      ierr=nf_put_att_double(ncid,nf_global,'theta_s',nf_double,1,
     &                                                       theta_s)
      ierr=nf_put_att_double(ncid,nf_global,'theta_b',nf_double,1,
     &                                                       theta_b)
      ierr=nf_put_att_double(ncid,nf_global,'hc',nf_double,1,hc)
      ierr=nf_put_att_double(ncid,nf_global,'Cs_w',nf_double,N+1,
     &                                                          Cs_w)
      ierr=nf_put_att_double(ncid,nf_global,'Cs_r',nf_double,N,
     &                                                          Cs_r)
      ierr=nf_put_att_text(ncid,nf_global,'boundary_content',9,
     &                                                    'BEC2 only')

      ierr=nf_def_var(ncid,'bry_time',nf_double,1,dim_time,varid)
      ierr=nf_put_att_text(ncid,varid,'units',4,'days')
      if (ierr /= nf_noerr) then
        write(*,*) '### ERROR defining BEC2 boundary dimensions'
        stop
      endif
      write(*,*) '...done.'
      end
