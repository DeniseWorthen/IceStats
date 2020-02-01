module stats

  use param
  use grdvar
 
  implicit none

  real, dimension(nreg)         :: exmod, armod, himod, hsmod
  real, dimension(nreg)         :: rmean
  real, dimension(nreg)         :: ar15
  real, dimension(nreg)         :: tmelt, tprod
  contains

  !---------------------------------------------------------------------

   subroutine areaextent(reg,ain,asum,esum,a15)

       real, intent( in) :: reg
       real, intent( in), dimension(igrid,jgrid) :: ain
       real, intent(out) :: asum, esum, a15

    integer :: nr,i,j,ij

      nr = int(reg)
    asum = 0.0; esum = 0.0; a15 = 0.0
    do ij = 1,ijsize(nr)
        i = indx(nr,ij)
        j = jndx(nr,ij)
     if(ain(i,j) .ne. mval)then
      if(ain(i,j)    .ge.  0.15)then
        asum = asum + parea(i,j)*ain(i,j)
        esum = esum + parea(i,j)
      endif !
      if(ain(i,j) .lt. 0.15)then
        a15 = a15 + parea(i,j)*ain(i,j)
      endif 
     endif
    enddo
    !m2
    if(asum .eq. 0.0)asum = mval
    if(esum .eq. 0.0)esum = mval
    if( a15 .eq. 0.0) a15 = mval
    ! areas in km2
    if(asum .ne. mval)asum = asum/1.0e12
    if(esum .ne. mval)esum = esum/1.0e12
    if( a15 .ne. mval) a15 =  a15/1.0e12

   end subroutine areaextent

  !---------------------------------------------------------------------

   subroutine meanvar(reg,avar,ai,amean,aimin,weight)

       real, intent( in) :: reg,aimin
       real, intent( in), dimension(igrid,jgrid) :: avar, ai
    logical, intent( in) :: weight
       real, intent(out) :: amean

    integer :: nr,i,j,ij
       real :: asum, acnt

    ! sum(hi*ai*area)/sum(ai*area) = total ice vol/total ice area 
    ! => mean ice thick
       nr = int(reg)
    amean = mval
     asum = 0.0; acnt = 0.0
    do ij = 1,ijsize(nr)
        i = indx(nr,ij)
        j = jndx(nr,ij)
      if((avar(i,j) .ne. mval) .and. &
          (  ai(i,j) .ge. aimin))then
        if(weight)then
         asum = asum +  ai(i,j)*avar(i,j)*parea(i,j)
        else
         asum = asum +          avar(i,j)*parea(i,j)
        end if
        acnt = acnt +    ai(i,j)*parea(i,j)
       !print *,j,avar(i,j),ai(i,j),asum,acnt
      end if
    enddo
    if(acnt .ne. 0.0)amean = asum/acnt

   end subroutine meanvar
end module stats
