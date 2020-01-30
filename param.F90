module param

  implicit none
 
#ifdef use_cfsv2
  integer, parameter :: igrid =  384, jgrid = 190
#endif
#ifdef use_m6c5
  integer, parameter :: igrid = 1440, jgrid = 1080
#endif
#ifdef use_sis2
  integer, parameter :: igrid = 1440, jgrid = 1080
#endif
#ifdef use_cpc
  integer, parameter :: igrid =  720, jgrid = 410
#endif
  ! region mask on 12min etopo grid
  integer, parameter :: ietop = 1800, jetop = 901
  integer, parameter :: ijmax = igrid*jgrid

  integer, parameter ::    yrbeg = 2011, yrend = 2019
  !integer, parameter ::    yrbeg = 2013, yrend = 2016
  integer, parameter ::   nyears = (yrend-yrbeg)+1
  integer, parameter ::     nmon = 12
  integer, parameter ::    ndays = 35
  !integer, parameter ::    ndays = 5
  integer, parameter :: maxsteps = nyears*366
  integer, parameter ::  maxexps = 2*nmon*nyears

  integer, parameter ::     nreg = 14 &    !IHO regions
                                 + 1  &    ! NOcn (lat>50N)
                                 + 5  &    ! GSFC regions
                                 + 1       ! SOcn (lat<50S)

  character(len=18), dimension(nreg) :: regname = (/ "Gulf of Alaska   ", & 
                                                     "Bering Sea       ", & 
                                                     "Chukchi Sea      ", &
                                                     "Beaufort Sea     ", & 
                                                     "Baffin Bay       ", &
                                                     "Lincoln Sea      ", &
                                                     "White Sea        ", &
                                                     "EastSib Sea      ", &
                                                     "NW Passages      ", &
                                                     "Central Arctic   ", &
                                                     "Barents Sea      ", &
                                                     "Greenland Sea    ", &
                                                     "Kara Sea         ", &
                                                     "Laptev Sea       ", &
                                                     "Nocn >50N        ", &
                                                     "Weddell Sea      ", & !60W:20E
                                                     "Indian Ocean     ", & !20E:90E
                                                     "Pacific Ocean    ", & !90:160E
                                                     "Ross Sea         ", & !160E:130W
                                                     "Belling-Amund Sea", & !130W:60W
                                                     "Socn <50S        "  &
                                                   /)

  integer, parameter :: nvars = 1  & !extent
                              + 1  & !area
                              + 1  & !area 15
                              + 1  & !hi
                              + 1  & !hs
                              + 2    !melt,prod

end module param
