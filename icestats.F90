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
  integer :: lij

  character(len=10) :: cdate1, initdate00
  character(len= 8) :: cdate2, initdate

     real(kind=8) :: tstamp0
     real(kind=4) :: rnum, area, extent, area15, mvar

     real(kind=4) :: modmin,modmax

     real, dimension(igrid,jgrid)         :: ai, aihi, aihs
     real, allocatable, dimension(:,:,:)  :: aisum, hisum, hssum
     real, dimension(nmon,31)             :: cnt

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

    lij = 0
   do j = 1,jgrid
    do i = 1,igrid
     lij = lij + pmask(i,j)
     if(i.eq.400.and.j.eq.1000)print *,'i=400,j=1000',lij
    enddo
   enddo

   allocate(aisum(1:lij,1:nmon,1:31))
   allocate(hisum(1:lij,1:nmon,1:31))
   allocate(hssum(1:lij,1:nmon,1:31))

  !---------------------------------------------------------------------
  !
  !---------------------------------------------------------------------

  outcdf = trim(dirout)//trim(ssrc)//'.ice.stats.nc'
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
      !cdffile = trim(rtsrc)//trim(rtname)//'ice'//trim(cdate2)//'.01.'//trim(initdate)//'00.nc'
      cdffile = trim(rtsrc)//trim(rtname)//'ice'//trim(cdate2)//'.01.'//trim(initdate)//'00.subset.nc'
      call get_pfld(trim(cdffile),     trim('aice_h'),       ai,     1)
      call get_pfld(trim(cdffile),     trim('hi_h'  ),     aihi,     1)
      call get_pfld(trim(cdffile),     trim('hs_h'  ),     aihs,     1)

        ij = 0
      do j = 1,jgrid
       do i = 1,igrid
        if(pmask(i,j) .eq. 1.0)then
         ij = ij + 1
         if(  ai(i,j) .ne. mval)aisum(ij,mon,day) = aisum(ij,mon,day) +   ai(i,j)
         if(aihi(i,j) .ne. mval)hisum(ij,mon,day) = hisum(ij,mon,day) + aihi(i,j)
         if(aihs(i,j) .ne. mval)hssum(ij,mon,day) = hssum(ij,mon,day) + aihs(i,j)
        endif
       enddo
      enddo
      cnt(mon,day) = cnt(mon,day) + 1.0
      write(20,'(4i5,3f12.5)')nex,year,mon,day,cnt(mon,day),ai(400,1000),aisum(887140,mon,day)

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
  write(20,*)

  !---------------------------------------------------------------------
  !
  !---------------------------------------------------------------------

  do ij = 1,lij
   where(cnt(:,:) .ne. 0.0)aisum(ij,:,:) = aisum(ij,:,:)/cnt(:,:)
   where(cnt(:,:) .ne. 0.0)hisum(ij,:,:) = hisum(ij,:,:)/cnt(:,:)
   where(cnt(:,:) .ne. 0.0)hssum(ij,:,:) = hssum(ij,:,:)/cnt(:,:)
  end do
  mon = 9; day = 25
  write(20,'(3i5,2f12.5)')nex,mon,day,cnt(mon,day),aisum(887140,mon,day)
  write(20,*)

     ai = mval; aihi = mval; aihs = mval
   year = 2012
    lll = 0
  do nm = 1,nmon
   do nd = 1,mnend(nm)
     lll = lll + 1

      ij = 0
    do j = 1,jgrid
     do i = 1,igrid
      if(pmask(i,j) .eq. 1.0)then
               ij = ij + 1
          ai(i,j) = aisum(ij,nm,nd)
        aihi(i,j) = hisum(ij,nm,nd)
        aihs(i,j) = hssum(ij,nm,nd)
      endif
     enddo
    enddo

    call set_taxis(year,nm,nd)
    call write_fieldcdf(trim(dirout)//trim(ssrc)//'.ai.dm.nc','aice_h',  ai,lll)
    call write_fieldcdf(trim(dirout)//trim(ssrc)//'.hi.dm.nc',  'hi_h',aihi,lll)
    call write_fieldcdf(trim(dirout)//trim(ssrc)//'.hs.dm.nc',  'hs_h',aihs,lll)
    write(20,'(3i5,2f12.5)')nex,nm,nd,cnt(nm,nd),ai(400,1000)

   enddo
  enddo

  
end program icestats
