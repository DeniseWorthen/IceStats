module charstrings

  use param

  implicit none

   character(len=120) ::  wgtsrc="/scratch2/NCEPDEV/climate/Denise.Worthen/NEMS_INPUT0.1/regrids/"
   character(len=120) ::  obssrc="/scratch2/NCEPDEV/climate/Denise.Worthen/IceData/"
   character(len=120) ::  grdsrc="/scratch2/NCEPDEV/climate/Denise.Worthen/TTout/"
#ifdef use_cfsv2
   character(len= 80) ::   rtsrc="/scratch2/NCEPDEV/climate/Denise.Worthen/CFSv2/"
   character(len= 80) ::  dirout="/scratch2/NCEPDEV/climate/Denise.Worthen/CFSv2/"
   character(len= 10) ::    fsrc="cfs."
   character(len= 20) ::    wsrc="t126"
#endif
#ifdef use_m6c5
   character(len= 80) ::   rtsrc="/scratch2/NCEPDEV/stmp3/Denise.Worthen/BM3_ice/"
   character(len= 80) ::  dirout="/scratch2/NCEPDEV/stmp3/Denise.Worthen/BM3_ice/"

   character(len= 10) ::    fsrc=""
   character(len= 20) ::    wsrc="tripole.mx025"
#endif
#ifdef use_sis2
   character(len= 80) ::   rtsrc="/scratch2/NCEPDEV/stmp3/Denise.Worthen/SIS2_ice/"
   character(len= 80) ::  dirout="/scratch2/NCEPDEV/stmp3/Denise.Worthen/SIS2_ice/"

   character(len= 10) ::    fsrc="sis2."
   character(len= 20) ::    wsrc="tripole.mx025"
#endif
#ifdef use_cpc
   character(len= 80) ::   rtsrc="/scratch2/NCEPDEV/stmp3/Denise.Worthen/CPC_ice/"
   character(len= 80) ::  dirout="/scratch2/NCEPDEV/stmp3/Denise.Worthen/CPC_ice/"

   character(len= 10) ::    fsrc="cpc."
   character(len= 20) ::    wsrc="tripole.mx05"
#endif
   character(len= 80) :: rtname

   character(len=8) :: i2fmt = '(i2.2)'
   character(len=8) :: i4fmt = '(i4.4)'

end module charstrings
