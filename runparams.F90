module runparams

  use param
  use caldata

  implicit none

  integer, dimension(3,maxsteps) :: ymd

  integer :: nsteps, nexps, nefrst, nelast
  integer :: mnend(nmon)
  integer :: lbeg(maxexps)
  integer :: lend(maxexps)
  integer :: lybeg(nyears), lyend(nyears)

  contains

  subroutine setup_runs

  implicit none

  integer :: ny, nm, nd, ne
  integer :: lll, nex, year, mon, day

  !---------------------------------------------------------------------

     lbeg = -1;  lend = -1
    lybeg = -1; lyend = -1
      lll = 0
      nex = 0
      ymd = 0
    do ny = 1,nyears
     year = yrbeg + (ny-1)
                           mnend = mnendn
     if(mod(year,4) .eq. 0)mnend = mnendl

     do nm = 1,nmon
      do nd = 1, mnend(nm)
              lll = lll +1
       ymd(1,lll) = year
       ymd(2,lll) = nm
       ymd(3,lll) = nd
       if(nd .eq. 1 .or. nd .eq. 15)then
              nex = nex+1
        lbeg(nex) = lll
        write(*,*)year,nm,nd,lll,nex,lbeg(nex)
       endif
      enddo
     enddo
   enddo!ny
    nexps = nex
   where(lbeg .ne. -1)lend = lbeg + (ndays-1)

   do ne = 1,nexps
    year = ymd(1,lbeg(ne)); mon = ymd(2,lbeg(ne)); day = ymd(3,lbeg(ne))
    !if((year .eq. 2018) .and.  (mon .eq. 3))nelast = ne
    if((year .eq. yrend) .and.  (mon .eq. 3))nelast = ne
   enddo
   do ne = nexps,1,-1
    year = ymd(1,lbeg(ne)); mon = ymd(2,lbeg(ne)); day = ymd(3,lbeg(ne))
    if((year .eq. 2011) .and.  (mon .eq. 4))nefrst = ne
    !if((year .eq. 2013) .and.  (mon .eq. 1))nefrst = ne
    !if((year .eq. 2011) .and.  (mon .eq. 8))nefrst = ne
   end do

   do ne = nefrst,nelast
    year = ymd(1,lbeg(ne)); mon = ymd(2,lbeg(ne)); day = ymd(3,lbeg(ne))
    print '(3i8,6i6)',ne,lbeg(ne),lend(ne),ymd(:,lbeg(ne)),ymd(:,lend(ne))
   enddo
   nsteps = lend(nelast)
   print *,'nelast = ',nelast,' nsteps = ',nsteps
 
   do ne = 1,nelast
    year = ymd(1,lbeg(ne)); mon = ymd(2,lbeg(ne)); day = ymd(3,lbeg(ne))
   if((mon .eq. 1) .and. &
       (day .eq. 1))then
          year = ymd(1,lbeg(ne))
            ny = (year-yrbeg) + 1
     lybeg(ny) = lbeg(ne)
    endif
   enddo
   do ne = 1,nelast
     year = ymd(1,lbeg(ne))
       ny = (year-yrbeg) + 1
    if(mod(ne,24) .eq. 0)lyend(ny) = lend(ne)
   enddo
   do ny = 1,nyears
    print *,lybeg(ny),lyend(ny),24*(ny-1)+1,ny*24
   enddo

  end subroutine setup_runs
end module runparams
