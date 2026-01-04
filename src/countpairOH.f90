subroutine countpairoh(ng, nm, np, maxoh, genofr, pairids, ohrf)
  implicit none
  
  integer, intent(IN) :: ng, nm, np, maxoh
  integer, intent(IN) :: pairids(2*np), genofr(ng*nm)
  integer, intent(INOUT) :: ohrf(np)
  integer :: i,j, l,x, pairs(np, 2), Genos(nm, ng), IsOppHom(-1:2,-1:2)
  
  pairs(:,1) = pairids(1:np)
  pairs(:,2) = pairids((np+1):(2*np))
  
  ! genofr vector to matrix
  Genos = -1
  j = 0
  do l=1,nm      
    do i=1, ng
    j = j+1
    if (GenoFR(j)>=0) then
      Genos(l,i) = GenoFR(j)
    endif
    enddo    
  enddo   
  
  ! count OH
  IsOppHom = 0
  IsOppHom(0, 2) = 1
  IsOppHom(2, 0) = 1
  
  ohrf = -9
  do x=1, np
    if (pairs(x,1) <= 0 .or. pairs(x,2) <= 0)  cycle   ! not genotyped
    ohrf(x) = calcOH(pairs(x,1), pairs(x,2))
  enddo
  

contains
  integer function calcOH(i,j)  result(OH_ij)
    integer, intent(IN) :: i,j
    integer :: l
    
    OH_ij = 0
    do l=1, nm
      OH_ij = OH_ij + IsOppHom(Genos(l,i), Genos(l,j))  
      if (OH_ij > maxOH) exit          
    enddo
  end function calcOH
  
  
end subroutine countpairoh


