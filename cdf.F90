module cdf

  use param
  use netcdf

  implicit none

  integer, parameter :: ndim1 = 1, ndim2 = 2, ndim3 = 3, ndim4 = 4

  integer, dimension(ndim4) :: dim4, corner4, edge4
  integer, dimension(ndim3) :: dim3, corner3, edge3
  integer, dimension(ndim2) :: dim2, corner2, edge2
  integer, dimension(ndim1) :: dim1, corner1, edge1

  integer :: xtdim, ytdim, ztdim, zrdim, tdim, tmdim
  integer :: xtid,  ytid, ztid, zrid

  integer :: rc, ncid, timid, tbid, tbdim, datid
  integer :: clen

  real(kind=4) :: mval = -9999.

  real(kind=8)                      :: tstamp
  real(kind=8), dimension(2)        :: tsbound
  character(len=31)    :: torg = 'hours since 1991-01-01 00:00:00'
  character(len=31)    :: tcal = 'gregorian'

  real(kind=8), allocatable, dimension(:) :: taxis

  character(len=64)    :: vname,vunit,vlong,fname
  contains

!---------------------------------------------------------------------

  subroutine set_taxis(ytmp,mtmp,dtmp)

  implicit none
  real*8 tm_secs_from_bc

  integer, intent(in) :: ytmp, mtmp, dtmp

  real(kind = 8) :: tday

  tday = tm_secs_from_bc(ytmp, mtmp, dtmp,  &
                           12,    0,    0) - &
         tm_secs_from_bc(1991,    1,    1,   &
                            0,    0,    0)
    tday = tday/(60.d0*60.d0)
  tstamp = tday

  tsbound(1) = tstamp - 12.0
  tsbound(2) = tstamp + 12.0
  end subroutine set_taxis

!---------------------------------------------------------------------

  subroutine get_cdate(yr,mon,day,cdate1,cdate2)

            integer, intent( in) :: yr,mon,day
  character(len= *), intent(out) :: cdate1,cdate2

  character(len=4) :: cyear
  character(len=2) :: cmon
  character(len=2) :: cday

  character(len=8) :: i2fmt = '(i2.2)'
  character(len=8) :: i4fmt = '(i4.4)'

   write(cyear, i4fmt)yr
   write( cmon, i2fmt)mon
   write( cday, i2fmt)day

   write(cdate1,'(a4,a1,a2,a1,a2)')cyear,'-',cmon,'-',cday
   cdate2 = trim(cyear)//trim(cmon)//trim(cday)

  end subroutine get_cdate

!---------------------------------------------------------------------

  subroutine get_pfld(cdffile,varname,aout, lstep)

  implicit none

                 character(len=*), intent( in) :: cdffile
                 character(len=*), intent( in) :: varname
                          integer, intent( in) :: lstep
     real, dimension(igrid,jgrid), intent(out) :: aout

  real :: lmval

  aout = mval
    rc = nf90_open(trim(cdffile), nf90_nowrite, ncid)
  !if(rc .ne. 0)print *,'looking for file ',trim(cdffile)
  if(rc .ne. 0)then
   rc = nf90_close(ncid)
   return
  endif

  rc = nf90_inq_varid(ncid, trim(varname), datid)
  if(rc .ne. 0)then
   rc = nf90_close(ncid)
   return
  endif

   corner3(3) = lstep
   corner3(2) = 1
   corner3(1) = 1
     edge3(3) = 1
     edge3(2) = jgrid
     edge3(1) = igrid
  rc = nf90_get_var(ncid, datid,  aout, corner3, edge3)
  rc = nf90_get_att(ncid, datid, '_FillValue', lmval)
  rc = nf90_close(ncid)

  ! change the missing value
  where(aout .eq. lmval)aout = mval

  end subroutine get_pfld
end module cdf
