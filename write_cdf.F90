subroutine write_cdf(cdffile,lstep,kstep)

  use param
  use runparams
  use stats
  use variablelist
  use netcdf
  use charstrings
  use cdf

  implicit none

                       character(len=*), intent(in) :: cdffile
                                integer, intent(in) :: kstep, lstep

  !-----------------------------------------------------------------------------

            integer :: n,nlast
  character(len=64) :: varname, varunit, varlong

  !-----------------------------------------------------------------------------

  if(lstep .eq. 0)then
   rc = nf90_create(trim(cdffile), nf90_clobber, ncid)
   !print *,'setting up ',trim(cdffile)

   rc = nf90_def_dim(ncid,    'Yt',           nreg,     ytdim)
   rc = nf90_def_dim(ncid,    'Zt',         nelast,     ztdim)
   rc = nf90_def_dim(ncid,  'time', nf90_unlimited,      tdim)

   dim1(1) =  tdim
   rc = nf90_def_var(ncid, 'time', nf90_double,  dim1,      timid)
   rc = nf90_put_att( ncid, timid,            'units', trim(torg))
   rc = nf90_put_att( ncid, timid,         'calendar', trim(tcal))
   rc = nf90_put_att( ncid, timid,             'axis',        'T')

   dim3(3) =  tdim
   dim3(2) = ztdim
   dim3(1) = ytdim
   do n = 1,nvars
    varname = icefields(n)%field_name
    varunit = icefields(n)%unit_name
    varlong = icefields(n)%long_name
    rc = nf90_def_var(ncid, trim(varname), nf90_float, dim3, datid)
    rc = nf90_put_att(ncid, datid,     'units', trim(varunit))
    rc = nf90_put_att(ncid, datid, 'long_name', trim(varlong))
    rc = nf90_put_att(ncid, datid, 'missing_value', mval)
    rc = nf90_put_att(ncid, datid,    '_FillValue', mval)
   enddo
    rc = nf90_enddef(ncid)

    rc = nf90_close(ncid)
  endif

  !-----------------------------------------------------------------------------

  rc = nf90_open(trim(cdffile), nf90_write, ncid)

   corner1(1) = lstep
     edge1(1) = 1
   rc = nf90_inq_varid(ncid, 'time',  timid)
   rc = nf90_put_var(ncid, timid, tstamp, corner1)

   corner3(3) = lstep
   corner3(2) = kstep
   corner3(1) = 1
     edge3(3) = 1
     edge3(2) = 1
     edge3(1) = nreg
   rc = nf90_inq_varid(ncid,    'exmod',  datid)
   rc = nf90_put_var(ncid, datid, exmod, corner3, edge3)
   rc = nf90_inq_varid(ncid,    'armod',  datid)
   rc = nf90_put_var(ncid, datid, armod, corner3, edge3)
   rc = nf90_inq_varid(ncid,     'ar15',  datid)
   rc = nf90_put_var(ncid, datid,  ar15, corner3, edge3)

   rc = nf90_inq_varid(ncid,    'himod',  datid)
   rc = nf90_put_var(ncid, datid, himod, corner3, edge3)
   rc = nf90_inq_varid(ncid,    'hsmod',  datid)
   rc = nf90_put_var(ncid, datid, hsmod, corner3, edge3)
   rc = nf90_close(ncid)

end subroutine write_cdf
