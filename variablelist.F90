module variablelist

  use param
  use cdf

  implicit none

  type IceFieldDefs
    character(len=64)                           :: field_name
    character(len=64)                           ::  long_name
    character(len=64)                           ::  unit_name
  end type IceFieldDefs

  type(IceFieldDefs) :: icefields(nvars)

  contains
  !---------------------------------------------------------------------

  subroutine def_vartypes

            integer :: ii
  character(len=40) :: vname, uname, lname

  ii = 0

  ii = ii + 1
  icefields(ii)%field_name    = 'exmod'
  icefields(ii)%long_name     = 'Model ice extent'
  icefields(ii)%unit_name     = 'm2'

  ii = ii + 1
  icefields(ii)%field_name    = 'armod'
  icefields(ii)%long_name     = 'Model ice area, ice>0.15'
  icefields(ii)%unit_name     = '10^6 km2'

  ii = ii + 1
  icefields(ii)%field_name    = 'ar15'
  icefields(ii)%long_name     = 'Model ice area, ice<0.15'
  icefields(ii)%unit_name     = '10^6 km2'

  ii = ii + 1
  icefields(ii)%field_name    = 'himod'
  icefields(ii)%long_name     = 'Model ice thickness'
  icefields(ii)%unit_name     = 'm'

  ii = ii + 1
  icefields(ii)%field_name    = 'hsmod'
  icefields(ii)%long_name     = 'Model snow thickness'
  icefields(ii)%unit_name     = 'm'

  end subroutine def_vartypes
end module variablelist
