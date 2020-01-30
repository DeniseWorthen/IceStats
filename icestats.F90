program icestats

  use param
  use runparams
  use cdf
  use variablelist
  use stats
  use grdvar
  use charstrings
  use netcdf

  implicit none

!---------------------------------------------------------------------
! local variables
!---------------------------------------------------------------------

  character(len=240) :: dirsrc, cdffile, outcdf

  integer :: i,j,ij,ne,nf,nh,lll,nex,nr
  integer :: nt, ny, nm, nd, year, mon, day, lstep
  integer :: varid

  character(len=10) :: cdate1, initdate00
  character(len= 8) :: cdate2, initdate

     real(kind=8) :: tstamp0
     real(kind=4) :: rnum, area, extent, area15, mvar

     real(kind=4) :: modmin,modmax

     real, dimension(igrid,jgrid)      :: ai, aihi, aihs
     ! generic field
     real, dimension(igrid,jgrid)      :: icevar, icesum
     real, dimension(igrid,jgrid)      :: icemelt, iceprod

  !---------------------------------------------------------------------
  !
  !---------------------------------------------------------------------

   call setup_runs
   allocate(taxis(1:nsteps))

   call def_vartypes
   !stop

  !---------------------------------------------------------------------
  ! get the grids
  !---------------------------------------------------------------------
#ifdef use_cfsv2
   cdffile = trim(rtsrc)//'t126_SCRIP.nc'
#else
   cdffile = trim(grdsrc)//'tripole.mx025.nc'
#endif

   print *,'getting grid from ',trim(cdffile)
   ! get the model grid variables
   call get_grid(trim(cdffile))
   print *,'plat ',minval(plat),maxval(plat)
   print *,'plon ',minval(plon),maxval(plon)
   print *,'pmask ',minval(pmask),maxval(pmask)

   call setup_masks

   print *,ijmax,' total grid points'
   do nr = 1,nreg
      ij = 0
    rnum = float(nr)
    do j = 1,jgrid
     do i = 1,igrid
      if(csm(i,j,nr) .eq. rnum)then
                ij = ij + 1
       indx(nr,ij) = i
       jndx(nr,ij) = j
      endif
     enddo
    enddo
    print *,'region ',nr,' ij points ',ij
    ijsize(nr) = ij
   enddo

  !---------------------------------------------------------------------
  !
  !---------------------------------------------------------------------

  outcdf = trim(rtsrc)//'stats.nc'
  print *,'working on ',trim(outcdf)
  call write_cdf(trim(outcdf),0,0)

  do nex = 1,nelast
  !do nex = 1,10
   year = ymd(1,lbeg(nex)); mon = ymd(2,lbeg(nex)); day = ymd(3,lbeg(nex)) 
   call get_cdate(year,mon,day,cdate1,cdate2)
   initdate = trim(cdate2)
     rtname = trim(fsrc)//trim(initdate)//'/'

   armod = mval; exmod = mval; himod = mval; ar15 = mval
    lstep = 0
    do lll = lbeg(nex),lend(nex)
        lstep = lstep + 1
         year = ymd(1,lll); mon = ymd(2,lll); day = ymd(3,lll) 
     if(nex .ge. nefrst)then
       call get_cdate(year,mon,day,cdate1,cdate2)
      cdffile = trim(rtsrc)//trim(rtname)//'ice'//trim(cdate2)//'.01.'//trim(initdate)//'00.nc'
      call get_pfld(trim(cdffile),     trim('aice_h'),       ai,     1)
      call get_pfld(trim(cdffile),     trim('hi_h'  ),     aihi,     1)
      call get_pfld(trim(cdffile),     trim('hs_h'  ),     aihs,     1)
      do nr = 1,nreg
       rnum = float(nr)
       call areaextent(rnum,    ai,area,extent,area15)
       armod(nr) =   area
       exmod(nr) = extent
        ar15(nr) = area15
       call    meanvar(rnum, aihi,ai,mvar,0.70,.false.)
       himod(nr) = mvar
       call    meanvar(rnum, aihs,ai,mvar,0.00,.false.)
       hsmod(nr) = mvar

       !call    meanvar(rnum, icemelt,ai,mvar,0.00,.true.)
       !tmelt(nr) = mvar
       !call    meanvar(rnum, iceprod,ai,mvar,0.00,.true.)
       !tprod(nr) = mvar

      enddo !nr
      modmin = minval(ai,mask=ai .ne. mval)
      modmax = maxval(ai,mask=ai .ne. mval)
      if(lstep .eq. 1 .or. lstep .eq. 35)print '(5i5,a2,a,e12.5)',nex,lstep,year,mon,day,'  ',trim(cdffile),exmod(15)
      if(modmin .eq. 0.0 .and. modmax .eq. 0.0)print '(4i5,2e12.5)',lstep,year,mon,day,modmin,modmax
     end if
 
     call set_taxis(year,mon,day)
     call write_cdf(trim(outcdf),lll,nex)

   enddo !lll
   print * 
  enddo !nex

end program icestats
