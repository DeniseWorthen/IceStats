module regmask_regrid_north

  use param
  use netcdf
  use cdf

  implicit none

  integer, parameter :: n_a = ietop*jetop
  integer, parameter :: n_b = igrid*jgrid
#ifdef use_cfsv2
  integer, parameter :: n_s = 48334
#endif
#ifdef use_m6c5  
  integer, parameter :: n_s = 1555200 
#endif
#ifdef use_sis2  
  integer, parameter :: n_s = 1555200 
#endif
#ifdef use_cpc  
  integer, parameter :: n_s = 295200 
#endif

  integer(kind=8), dimension(n_s), private :: col, row
     real(kind=8), dimension(n_s), private :: S

  contains
 
!---------------------------------------------------------------------
 
  subroutine get_weights_regmask(cdffile)

  implicit none

  character(len=*), intent(in) :: cdffile
 
  integer :: rc, ncid, datid, ilen
 
  !---------------------------------------------------------------------

  rc = nf90_open(trim(cdffile), nf90_nowrite, ncid)
  !print *,'get weights ',trim(cdffile)
  if(rc .ne. 0)stop

  rc = nf90_inq_dimid(ncid, 'n_s', datid)
  rc = nf90_inquire_dimension(ncid, datid, len=ilen)
  if(ilen .ne. n_s)stop
  rc = nf90_inq_dimid(ncid, 'n_a', datid)
  rc = nf90_inquire_dimension(ncid, datid, len=ilen)
  if(ilen .ne. n_a)stop
  rc = nf90_inq_dimid(ncid, 'n_b', datid)
  rc = nf90_inquire_dimension(ncid, datid, len=ilen)
  if(ilen .ne. n_b)stop

  rc = nf90_inq_varid(ncid, 'col', datid)
  rc = nf90_get_var(ncid, datid, col)
  rc = nf90_inq_varid(ncid, 'row', datid)
  rc = nf90_get_var(ncid, datid, row)
  rc = nf90_inq_varid(ncid,   'S', datid)
  rc = nf90_get_var(ncid, datid,  S)

  rc = nf90_close(ncid)

  end subroutine get_weights_regmask

!---------------------------------------------------------------------

  subroutine regmask_2model(ain,aout)

  implicit none

  real(kind=4), dimension(ietop,jetop), intent( in) :: ain
  real(kind=4), dimension(igrid,jgrid), intent(out) :: aout

  integer :: i,j,ii,jj
  real(kind=4), dimension(ietop*jetop) :: src_field
  real(kind=4), dimension(igrid*jgrid) :: dst_field

  !---------------------------------------------------------------------

    aout = 0.0
      ii = 0
    do j = 1,jetop
     do i = 1,ietop
                  ii = ii + 1
       src_field(ii) = ain(i,j)
     enddo
    enddo

    do i = 1,n_b
     dst_field(i) = 0.0
    enddo

    do i = 1,n_s
      ii = row(i); jj = col(i)
      dst_field(ii) = dst_field(ii) + S(i)*src_field(jj)
    enddo

      ii = 0
    do j = 1,jgrid
     do i = 1,igrid
             ii = ii + 1
      aout(i,j) = dst_field(ii)
     enddo
    enddo

  end subroutine regmask_2model
end module regmask_regrid_north

