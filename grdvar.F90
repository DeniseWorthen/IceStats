module grdvar

  use param
  use charstrings
  use cdf

  implicit none

  integer(kind=4), dimension(nreg,ijmax)    :: indx,jndx
  integer(kind=4), dimension(nreg)          :: ijsize

  real(kind=4), dimension(igrid,jgrid)      :: pmask
  real(kind=4), dimension(igrid,jgrid)      :: parea
  real(kind=4), dimension(igrid,jgrid)      :: plat, plon

  real(kind=4), dimension(igrid,jgrid,nreg) :: csm

  real(kind=4), dimension(ietop,jetop) :: emask
  real(kind=4), dimension(igrid,jgrid) :: rmask

  contains

  !---------------------------------------------------------------------

  subroutine get_grid(cdffile)

  character(len=*), intent(in) :: cdffile

       integer :: i,j
#ifdef use_cfsv2
   integer :: ij
   real(kind=8), parameter :: RADIUS = 6.3712e+6           !< Radius of the Earth [m]

   ! from intermediate t126_SCRIP.nc file created
   ! using NCO
      real(kind=8), dimension(igrid*jgrid) :: dlon
      real(kind=8), dimension(igrid*jgrid) :: dlat
      real(kind=8), dimension(igrid*jgrid) :: darea

#else
   integer(kind=4), dimension(igrid,jgrid) :: imask
#endif

  rc = nf90_open(trim(cdffile), nf90_nowrite, ncid)
  if(rc .ne. 0)print *,trim(nf90_strerror(rc))

#ifdef use_cfsv2
  rc = nf90_open(trim(cdffile), nf90_nowrite, ncid)
  if(rc .ne. 0)print *,trim(nf90_strerror(rc))

  rc = nf90_inq_varid(ncid, 'grid_center_lat', datid)
  rc = nf90_get_var(ncid,               datid,  dlat)
  rc = nf90_inq_varid(ncid, 'grid_center_lon', datid)
  rc = nf90_get_var(ncid,               datid,  dlon)
  rc = nf90_inq_varid(ncid,       'grid_area', datid)
  rc = nf90_get_var(ncid,               datid, darea)
  rc = nf90_inq_varid(ncid,      'grid_imask', datid)
  rc = nf90_get_var(ncid,               datid, imask)
  rc = nf90_close(ncid)

    ij = 0
  do j = 1,jgrid
   do i = 1,igrid
            ij = ij+1
     plon(i,j) =  dlon(ij)
     plat(i,j) =  dlat(ij)
    parea(i,j) = darea(ij)*radius*radius
    pmask(i,j) = imask(ij)
   enddo
  enddo
  !print '(10f12.5)',plon(1:10,10)
  !print '(10f12.5)',plat(10,1:10)
  print '(10f12.5)',pmask(10,1:10)
#else
  rc = nf90_inq_varid(ncid, 'area', datid)
  rc = nf90_get_var(ncid,    datid, parea)

  rc = nf90_inq_varid(ncid,  'wet', datid)
  print *,trim(nf90_strerror(rc))
  rc = nf90_get_var(ncid,    datid, imask)
  
  rc = nf90_inq_varid(ncid, 'latCt', datid)
  rc = nf90_get_var(ncid,     datid,  plat)

  rc = nf90_inq_varid(ncid, 'lonCt', datid)
  rc = nf90_get_var(ncid,     datid,  plon)
  rc = nf90_close(ncid)

  pmask = real(imask,4)
#endif
  where(plon .lt. 0.0)plon = plon + 360.0

  end subroutine get_grid

  !---------------------------------------------------------------------

  subroutine get_regmask(cdffile)

  character(len=*), intent(in) :: cdffile

       integer :: i,j,imval
  integer(kind=4), dimension(ietop,jetop) :: i2d

  rc = nf90_open(trim(cdffile), nf90_nowrite, ncid)
  if(rc .ne. 0)print *,trim(nf90_strerror(rc))

  rc = nf90_inq_varid(ncid, 'allmask_fix',    datid)
  rc = nf90_get_var(ncid,           datid,      i2d)
  rc = nf90_get_att(ncid, datid, '_FillValue',imval)
  rc = nf90_close(ncid)
  where(i2d .eq. imval)i2d = 0
 ! print *,minval(i2d),maxval(i2d)

  emask = real(i2d,4)
  print *,minval(emask),maxval(emask)

  end subroutine get_regmask

  !---------------------------------------------------------------------

  subroutine setup_masks

  use regmask_regrid_north

  implicit none

!---------------------------------------------------------------------
! local variables
!---------------------------------------------------------------------

  integer :: i,j,nr,nr1,maxreg,rnum
     real :: rlon

  character(len=240) :: cdffile

  !---------------------------------------------------------------------

   cdffile = trim(wgtsrc)//'etopo12_'//trim(wsrc)//'.nc'
   print *,'using weights file ',trim(cdffile)
   call get_weights_regmask(trim(cdffile))

  ! get the shapefile mask
   cdffile = trim(wgtsrc)//trim('etopo12_regions.nc')
   call get_regmask(trim(cdffile))

   call regmask_2model(emask,rmask)
   where(pmask .eq. 0.0)rmask = 0.0
   where(rmask .lt. 0.0)rmask = mval
   cdffile=trim(dirout)//'rmask.nc'
   call check_regmask(trim(cdffile),rmask)

   ! maximum number of regions from shapefile 
   maxreg = maxval(rmask)
   print *,' found ',maxreg,' regions in etopo mask file'

   csm = 0.0
   do j = 1,jgrid
    do i = 1,igrid
     do nr = 1,maxreg
      if(rmask(i,j) .eq. real(nr,4))csm(i,j,nr) = real(nr,4)
     enddo
    enddo
   enddo

   ! make special regions
    nr1 = maxreg
   do j = 1,jgrid
    do i = 1,igrid
     if(pmask(i,j) .eq. 1.0)then
      ! N 50N
      if(plat(i,j) .ge.  50.0)csm(i,j,nr1+1) = float(nr1+1)
      ! Antarctic Wedges
      if(plat(i,j) .le. -50.0)then
                      rlon = plon(i,j)
       if(rlon .lt. 0)rlon = rlon + 360

       if((rlon .ge. 300.0) .or. &
           (rlon .lt.  20.0))csm(i,j,nr1+2) = float(nr1+2)  ! Weddell
       if((rlon .ge.  20.0) .and. &
           (rlon .lt.  90.0))csm(i,j,nr1+3) = float(nr1+3)  ! Indian
       if((rlon .ge.  90.0) .and. &
           (rlon .lt. 160.0))csm(i,j,nr1+4) = float(nr1+4)  ! Pacific
       if((rlon .ge. 160.0) .and. &
           (rlon .lt. 230.0))csm(i,j,nr1+5) = float(nr1+5)  ! Ross
       if((rlon .ge. 230.0) .and. &
           (rlon .lt. 300.0))csm(i,j,nr1+6) = float(nr1+6)  ! Belling-Amund
      endif
      if(plat(i,j) .le. -50.0)csm(i,j,nr1+7) = float(nr1+7)  ! S 50S
     endif
    enddo
   enddo
   do nr = 1,nreg
    where(pmask(:,:) .eq. 0.0)csm(:,:,nr) = mval
    print *,nr,maxval(csm(:,:,nr))
   enddo

  ! check the final region mask
  cdffile=trim(dirout)//'csm.nc'
  call check_csm(trim(cdffile),csm)

  do nr = 1,nreg
   print *,'region ',nr,trim(regname(nr))
  enddo

end subroutine setup_masks

!---------------------------------------------------------------------

  subroutine check_regmask(cdffile,ain)

  implicit none

                       character(len=*), intent(in) :: cdffile
           real, dimension(igrid,jgrid), intent(in) :: ain

            integer :: n
  character(len=64) :: varname, varunit, varlong

  !---------------------------------------------------------------------

   rc = nf90_create(trim(cdffile), nf90_clobber, ncid)
   print *,'setting up ',trim(cdffile)

   rc = nf90_def_dim(ncid,    'Xt',         igrid,     xtdim)
   rc = nf90_def_dim(ncid,    'Yt',         jgrid,     ytdim)

   dim2(2) = ytdim
   dim2(1) = xtdim
   rc = nf90_def_var(ncid, 'plat',   nf90_float,          dim2, ytid)
   rc = nf90_put_att(ncid, ytid,        'units',     'degrees_north')
   rc = nf90_put_att(ncid, ytid,    'long_name',          'Latitude')

   dim2(2) = ytdim
   dim2(1) = xtdim
   rc = nf90_def_var(ncid, 'plon',   nf90_float,          dim2,  xtid)
   rc = nf90_put_att(ncid, xtid,        'units',       'degrees_east')
   rc = nf90_put_att(ncid, xtid,    'long_name',          'Longitude')

   dim2(2) = ytdim
   dim2(1) = xtdim
    varname = 'regmask'
    varunit = '  '
    varlong = 'IHO regions regridded to model grid'
    rc = nf90_def_var(ncid, trim(varname), nf90_float, dim2, datid)
    rc = nf90_put_att(ncid, datid,     'units', trim(varunit))
    rc = nf90_put_att(ncid, datid, 'long_name', trim(varlong))
    rc = nf90_put_att(ncid, datid, 'missing_value', mval)
    rc = nf90_put_att(ncid, datid,    '_FillValue', mval)
    rc = nf90_enddef(ncid)

    rc = nf90_inq_varid(ncid,     'regmask',  datid)
    rc = nf90_put_var(ncid,           datid,    ain)

    rc = nf90_put_var(ncid,  xtid,  plon)
    rc = nf90_put_var(ncid,  ytid,  plat)

    rc = nf90_close(ncid)

  end subroutine check_regmask

!---------------------------------------------------------------------

  subroutine check_csm(cdffile,ain)

  implicit none

                      character(len=*), intent(in) :: cdffile
     real, dimension(igrid,jgrid,nreg), intent(in) :: ain

            integer :: n
  character(len=64) :: varname, varunit, varlong

  !---------------------------------------------------------------------

   rc = nf90_create(trim(cdffile), nf90_clobber, ncid)
   print *,'setting up ',trim(cdffile)

   rc = nf90_def_dim(ncid,    'Xt',         igrid,     xtdim)
   rc = nf90_def_dim(ncid,    'Yt',         jgrid,     ytdim)
   rc = nf90_def_dim(ncid,    'Zt',          nreg,     ztdim)
#ifdef test
   dim2(2) = ytdim
   dim2(1) = xtdim
   rc = nf90_def_var(ncid, 'plat',   nf90_float,          dim2, ytid)
   rc = nf90_put_att(ncid, ytid,        'units',     'degrees_north')
   rc = nf90_put_att(ncid, ytid,    'long_name',          'Latitude')

   dim2(2) = ytdim
   dim2(1) = xtdim
   rc = nf90_def_var(ncid, 'plon',   nf90_float,          dim2,  xtid)
   rc = nf90_put_att(ncid, xtid,        'units',       'degrees_east')
   rc = nf90_put_att(ncid, xtid,    'long_name',          'Longitude')
#endif
   dim2(2) = ytdim
   dim2(1) = xtdim
    varname = 'parea'
    varunit = 'm2'
    rc = nf90_def_var(ncid, trim(varname), nf90_float, dim2, datid)
    rc = nf90_put_att(ncid, datid,     'units', trim(varunit))

   dim3(3) = ztdim
   dim3(2) = ytdim
   dim3(1) = xtdim
    varname = 'csm'
    varunit = '  '
    varlong = ' '
    rc = nf90_def_var(ncid, trim(varname), nf90_float, dim3, datid)
    rc = nf90_put_att(ncid, datid,     'units', trim(varunit))
    rc = nf90_put_att(ncid, datid, 'long_name', trim(varlong))
    rc = nf90_put_att(ncid, datid, 'missing_value', mval)
    rc = nf90_put_att(ncid, datid,    '_FillValue', mval)
    rc = nf90_enddef(ncid)

    !rc = nf90_put_var(ncid,  xtid,  plon)
    !rc = nf90_put_var(ncid,  ytid,  plat)

    rc = nf90_inq_varid(ncid,   'parea',  datid)
    rc = nf90_put_var(ncid,       datid,  parea)
    rc = nf90_inq_varid(ncid,     'csm',  datid)
    rc = nf90_put_var(ncid,       datid,    ain)
    rc = nf90_close(ncid)

  end subroutine check_csm

end module grdvar
