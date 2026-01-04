! Author: Jisca Huisman,  jisca.huisman@gmail.com
! Most of this code was written as a post doc in Evolutionary biology 
! at the University of Edinburgh, UK.
!
! This code is available under GNU General Public License v3
!
! The program is described in the paper
! "Pedigree reconstruction from SNP data: 
! Parentage assignment, sibship clustering, and beyond", 
! in Molecular Ecology Resources, 2017
!
! beta-versions are available at  https://github.com/JiscaH , 
! as well as a non-R version, reading & writing to text files.
!
!
! ####################################################################
! @@@@   MODULES   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! ####################################################################
module Global 
implicit none

integer :: nInd, nSnp, nIndLH, maxSibSize, MaxOppHom, MaxMendelE, Hermaphrodites, &
  nC(2), nYears, maxAgePO, nPairs, Complx, quiet, AgePhase, BYzero, chunk_printdot_i(10)
integer, parameter :: mxA=2**6, & ! max no. ancestors considered when testing for pedigree loop
  mxCP = 50, &  ! max no. candidate parents per sex
  XP = 5  ! multiplier nInd --> max no. candidate sib pairs
logical :: DoMtDif                  
logical, allocatable, dimension(:) :: ToCheck, SelfedIndiv
logical, allocatable, dimension(:,:) :: IsNewSibship, mtDif
integer, allocatable, dimension(:) :: Sex, BY, nFS, Mate, YearLast
integer,allocatable,dimension(:,:) :: Genos, Parent, nS, FSID, DumMate, DumClone
integer, allocatable, dimension(:,:,:) :: SibID, GpID
double precision :: TF, TA, OcA(3,-1:2), AKA2P(3,3,3), OKA2P(-1:2,3,3), zero
double precision, parameter ::  missing = 999D0, impossible=777D0, &
  NotImplemented = 444D0, MaybeOtherParent = 222D0, AlreadyAss = 888D0, neg_missing=-999D0
double precision, parameter :: ALR_tiny = -3.001  ! for consistency between R and SA; <0.001 (<log10(3)) rounded to 0
double precision, parameter :: quite_tiny = 0.0001d0  ! arises in ALR calcs due to rounding errors
double precision, allocatable, dimension(:) ::  Lind
double precision, allocatable, dimension(:,:) :: AHWE, OHWE, CLL
double precision, allocatable, dimension(:,:,:) :: AKAP, OKAP, OKOP, PPO, &
  PHS, PFS, LindX, IndBY, AgePriorA
double precision, allocatable, dimension(:,:,:,:) :: DumP, DumBY
double precision, allocatable, dimension(:,:,:,:,:) :: XPr
 
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                 
  contains
pure function MaxLL(V)
double precision, intent(IN) :: V(:)
double precision :: MaxLL
double precision, parameter :: tiny_LL = -0.0001d0        
MaxLL = missing
if (ANY(V < tiny_LL .and. V>-HUGE(0d0))) then  
  MaxLL = MAXVAL(V, mask = V<tiny_LL .and. V>-HUGE(0d0), DIM=1)      
else
  MaxLL = MINVAL(V, mask = V>ABS(tiny_LL), DIM=1) 
  ! 222  could be (related via) other parent
  ! 444  not (yet) implemented
  ! 777  impossible
  ! 888  already assigned
  ! 999  missing/not calculated
endif
end function MaxLL

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pure function addALR(LLg, ALR)   
double precision, intent(IN) :: LLg, ALR
double precision :: addALR

if (LLg > 0) then
  addALR = LLg
else if (ALR == impossible) then
  addALR = impossible
else
  addALR = LLg + ALR
endif

end function addALR

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pure function getPar(A, kA)
integer, intent(IN) :: A, kA
integer :: getPar(2)

if (A > 0) then
  getPar = Parent(A,:)
else if (A < 0) then
  getPar = GpID(:, -A, kA)
else
  getPar = 0
endif

end function getPar

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pure function getMate(A, kA)
integer, intent(IN) :: A, kA
integer :: getMate

if (A==0 .or. Complx/=0) then
  getMate = 0
else if (A > 0) then
  getMate = Mate(A)
else !if (A < 0) then
  getMate = DumMate(-A, kA)
endif

end function getMate

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine setMate(A,kA, M)   ! only called by SelectParent()
integer, intent(IN) :: A, kA, M

if (A==0 .or. Complx/=0)  return

if (A > 0) then
  Mate(A) = M
else !if (A < 0) then
  DumMate(-A, kA) = M
endif

end subroutine setMate

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pure function getNFS(A)
  integer, intent(IN) :: A
  integer :: getNFS, Ai

  getNFS = 0
  if (A<=0) return
  if (nFS(A)>0) then
    getNFS = nFS(A)
  else
    Ai = FSID(maxSibSize+1, A)
    getNFS = nFS(Ai)
  endif
end function getNFS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function getAP(AgeD, Rel, k, m, noGo)  
! Rel: 1=PO, 2=FS, 3=HS, 4=GP, 5=FA, 6=HA
! k: sex of 2nd indiv (except for HA)
! m: related via mat/pat                                        

integer, intent(IN) :: AgeD, Rel, k, m
double precision :: getAP
double precision :: AM(5,5), noGo   
integer :: D2, D3

getAP = zero
if (AgeD == 999) return
if (Rel < 1 .or. Rel > 6)  call Erstop('getAP: illegal Rel', .TRUE.)   
if ((Rel == 1 .and. AgeD <=0) .or. &   ! PO must have ageD >=1
    (Rel == 1 .and. AgeD > MaxAgePO) .or. &  ! PO must have AgeD <= MaxAgePO
    (Rel == 4 .and. AgeD <=1) .or. &     ! GP must have ageD >=2 
    (AgeD <= -MaxAgePO) .or. &   ! neg. AgeD only defined for HA/FA, <-MaxAgePO impossible
    (AgeD > 2*MaxAgePO)) then     ! GGP etc todo separately elsewhere
  getAP = noGo
  return
endif   

if (((m<1 .or.m>4) .and. Rel>2) .or. &
  ((k<1 .or.k>4) .and. (Rel==1 .or. Rel==4 .or. Rel==6)))  then
  call Erstop("getAP: illegal k or m!", .TRUE.)
endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! AgePriorA D2 + D3:
!  1    2     3
!  M   MGM   PGM
!  P   MGP   PGP
! FS   MFA   PFA
! MS  MMHA  PMHA
! PS  MPHA  PPHA
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (Rel < 4) then
  D3 = 1
else
  D3 = 1+m  ! via mother/father
endif
if (Rel==1 .or. Rel==4) then 
  D2 = k  ! sex of parent/GP
else if (Rel==2 .or. Rel==5) then  ! FS/FA
  D2 = 3
else if (Rel==3) then
  D2 = 3+m       
else if (Rel==6) then   ! HA (k HS of parent m)
  D2 = 3+k
else
  D2 = 0
  call ErStop("getAP: illegal Rel", .TRUE.)
endif

AM(:,1:3) = AgePriorA(AgeD,:,:)
AM(:,4) = (AM(:,2) + AM(:,3)) / 2.0  ! m=3: unknown via which parent
AM(:,5) = AM(:,4)    ! m=4: hermaphrodite = unknown sex

if ((Rel==1 .or. Rel==4) .and. k>2) then  ! unknown (grand)parent sex
  getAP = SUM(AM(1:2, D3)) / 2.0
else if ((Rel==3 .and. m>2) .or. (Rel==6 .and. k>2)) then
  getAP = SUM(AM(4:5, D3)) / 2.0
else 
  getAP = AM(D2, D3)
endif

if (getAP/=getAP .or. getAP==zero) then
  getAP = noGo
else
  getAP = LOG10(getAP)
  if (getAP < ALR_tiny)  getAP = noGo                                   
endif

end function getAP

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
pure function AgeDiff(i,j)   

integer, intent(IN) :: i,j
integer :: AgeDiff

if (BY(i)>=0 .and. BY(j)>=0) then
  AgeDiff = BY(i) - BY(j)
else
  AgeDiff = 999
endif

end function AgeDiff

!~~~~~~~~~~~~~~~~~~~
function get_agedif(i,j)  ! returns a double vector of size 0:MaxAgePO*2
  integer, intent(IN) :: i,j
  double precision :: get_agedif(0:MaxAgePO*2)
  double precision :: BYLR_i(nYears), BYLR_j(nYears), prob_agedif(0:nYears*2)
  integer :: x,y,a
  
  if (i<0 .or. j<0)  call ErStop('get_agedif only for >0', .TRUE.)
  
  prob_agedif = 0d0
  if (BY(i)>=0 .and. BY(j)>=0) then
    a = BY(i) - BY(j)
    prob_agedif(a) = 1d0
  
  else
    call getEstBY(i, 3, 5, BYLR_i)  ! k=3, 5=all contributions from all relatives
    call getEstBY(j, 3, 5, BYLR_j)
    
    if (all((BYLR_i - LOG10(1.0D0/nYears)) < quite_tiny) .and. &
     all((BYLR_i - LOG10(1.0D0/nYears)) < quite_tiny)) then   ! no age info for either
      get_agedif = 1.0D0/(1+nYears*2)
      return
   endif   
    
    do x=1, nYears  ! i
      if (BYLR_i(x) < -HUGE(0.0D0)) cycle 
      do y=1,nYears  ! j
        if (BYLR_j(y) < -HUGE(0.0D0)) cycle 
        a = ABS(x-y)
        prob_agedif(a) = prob_agedif(a) + 10**(BYLR_i(x) + BYLR_j(y))
      enddo
    enddo
  endif
  prob_agedif = prob_agedif / sum(prob_agedif)   ! scale to sum to 1
  
  get_agedif = 0d0
  if (nYears < MaxAgePO) then  ! not sure if this can happen??
    get_agedif(0:nYears*2) = prob_agedif(:)
  else
    get_agedif(:) = prob_agedif(0:MaxAgePO*2)
  endif

end function get_agedif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                 
subroutine Rprint(message, IntData, DblData, DataType)            

character(len=*), intent(IN) :: message
integer, intent(IN) :: IntData(:)
double precision, intent(IN) :: DblData(:)
character(3), intent(IN) :: DataType
integer :: nchar, ndata, IntDummy(0)
double precision, allocatable, dimension(:) :: DblDataRounded

nchar = LEN(trim(message))

if (DataType == "DBL") then
  ndata = SIZE(DblData)
  allocate(DblDataRounded(ndata))
  DblDataRounded = DNINT(DblData * 100.0) / 100.0
  call dblepr(trim(message), nchar, DblDataRounded, ndata)
  deallocate(DblDataRounded)
else if (DataType == "INT") then
  ndata = SIZE(IntData)
  call intpr(trim(message), nchar, IntData, ndata)
else if (DataType == "NON") then
  call intpr(trim(message), nchar, IntDummy, 0) 
!  call labelpr(trim(message), nchar)    R >= 4.0.0
else
  call ErStop("invalid DataType for Rprint", .TRUE.)
endif

end subroutine Rprint

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
! print progress dots
subroutine print_dot(ndots)
  integer, intent(IN), optional:: ndots
  integer :: d,n

  if (present(ndots)) then
    n = ndots
  else
    n=1
  endif
  
  do d=1,n
    call rprint_status_tbl_dot()
  enddo
end subroutine print_dot

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
! create a sequence from 1 to n, with length equal to n_steps
function mk_seq(n, n_steps) result(seq) 
  integer, intent(IN) :: n, n_steps
  integer, allocatable :: seq(:)
  integer :: i
  real, allocatable :: probs(:)
  
  allocate(probs(n_steps))
  allocate(seq(n_steps))
  
  probs = (/ (i,i=1,n_steps) /) / real(n_steps)
  seq = NINT( probs * n )  ! round to nearest whole number
  WHERE (seq == 0)  seq = 1
  WHERE (seq >  n)  seq = n  
  deallocate(probs)
end function mk_seq
 

end module Global

!===============================================================================

module OHfun
  use Global, ONLY : nSnp, Genos, maxOppHom, nInd, neg_Missing, OHWE, OKOP
  implicit none
  
  integer, private :: IsOppHom(-1:2,-1:2), Ecnt(-1:2,-1:2, -1:2)
  double precision, private, allocatable :: PPO(:,:,:), QPOsparse_value(:)
  ! using Compressed Sparse Row (CSR) for QPO matrix
  integer, private, allocatable :: QPOsparse_row(:), QPOsparse_col(:)  

  contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine init_OH()  
    integer :: l,i,j
  
    IsOppHom = 0
    IsOppHom(0, 2) = 1
    IsOppHom(2, 0) = 1
    
    ! offspring - dam - sire
    Ecnt = 0
    Ecnt(0,2,:) = 1
    Ecnt(2,0,:) = 1
    Ecnt(0,:,2) = Ecnt(0,:,2) +1
    Ecnt(2,:,0) = Ecnt(2,:,0) +1
    Ecnt(1,0,0) = 1
    Ecnt(1,2,2) = 1   
    
    allocate(PPO(-1:2,-1:2,nSnp))
    PPO = 1D0
    do l=1,nSnp
      do i=0,2  ! obs offspring 
        do j=0,2    ! obs parent
          PPO(i,j,l) = OKOP(i,j,l) / OHWE(i,l)
        enddo
      enddo
    enddo 
    
    call mk_QPO()

  end subroutine init_OH
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine mk_QPO()   
    ! make Sparse Matrix with LR(PO/U) for pairs with OH <= maxOppHom; not conditional on any parents
    ! CSR = Compressed sparse row
    integer :: i,j,l,n, OH_ij
    integer, allocatable :: col_index(:), row_index(:), col_tmp(:)
    double precision, allocatable :: V(:), V_tmp(:)
    
    ! initialise temporary arrays  
    allocate(V(5*nInd), col_index(5*nInd), row_index(nInd+1))    
    V = neg_Missing
    col_index = 0
    row_index = 0
    row_index(1) = 0
    n = 0
    do i=1, nInd-1       
      do j=i+1, nInd
        OH_ij = 0
        do l=1, nSnp
          OH_ij = OH_ij + IsOppHom(Genos(l,i), Genos(l,j))  
          if (OH_ij > maxOppHom) exit          
        enddo
        if (OH_ij <= maxOppHom) then
          n = n+1  
          
          if (n > SIZE(col_index)) then  ! increase size of temporary arrays
            allocate(V_tmp(SIZE(col_index) + nInd), col_tmp(SIZE(col_index) + nInd))
            col_tmp = 0
            col_tmp(1:(n-1)) = col_index(1:(n-1))
            V_tmp = nSnp
            V_tmp(1:(n-1)) = V(1:(n-1))
!            call move_alloc(col_tmp, col_index)   ! Fortran2003
!            call move_alloc(V_tmp, V)   ! Fortran2003
            deallocate(V, col_index)
            allocate(V(SIZE(V_tmp)), col_index(SIZE(col_tmp)))
            V = V_tmp
            col_index = col_tmp
            deallocate(V_tmp, col_tmp)
          endif
          
          V(n) = calc_QLR_PO(i,j)
          col_index(n) = j 
        endif
      enddo
      row_index(i+1) = n
    enddo
    
    ! create output arrays of correct size   
    allocate(QPOsparse_value(n), QPOsparse_col(n), QPOsparse_row(nInd+1))
    QPOsparse_value = V(1:n)
    QPOsparse_col = col_index(1:n)
    QPOsparse_row = row_index
  end subroutine mk_QPO
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure function calc_QLR_PO(i,j) result(QLR_PO)  !quick check, not conditioning on parents.
    integer, intent(IN) :: i,j
    double precision :: QLR_PO
    integer :: l
    double precision :: PrL(nSnp)

    PrL = 0D0
    do l=1,nSnp
      PrL(l) = LOG10(PPO(Genos(l,i), Genos(l,j), l))  
    enddo
    QLR_PO = SUM(PrL)
  end function calc_QLR_PO
    
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function QLR_PO(i,j)   ! Lookup in sparse matrix
    integer, intent(IN) :: i,j
    double precision :: QLR_PO
    integer :: ii, jj, x
       
    if (i == j) then
      QLR_PO = neg_Missing
      return
    else if (i < j) then
      ii = i
      jj = j
    else
      ii = j
      jj = i
    endif
    
    QLR_PO = neg_Missing ! if i+j not found in sparse matrix, OH > maxOppHom, and QLR(PO/U) not calculated
    do x = QPOsparse_row(ii)+1, QPOsparse_row(ii+1)    ! +1 because Fortran indexes starting from 1, not from 0
      if (QPOsparse_col(x) == jj) then
        QLR_PO = QPOsparse_value(x)
        exit
      endif
    enddo 
  end function QLR_PO

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  integer function calcOH(i,j, maxOH)  result(OH_ij)
    integer, intent(IN) :: i,j
    integer, intent(IN), optional :: maxOH
    integer :: l, maxOH_
    
    maxOH_ = maxOppHom
    if (present(maxOH))  maxOH_ = maxOH
    
    OH_ij = 0
    do l=1, nSnp
      OH_ij = OH_ij + IsOppHom(Genos(l,i), Genos(l,j))  
      if (OH_ij > maxOH_) exit          
    enddo
  end function calcOH
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  integer function CalcTrioErr(A,par)
    integer, intent(IN) :: A, Par(2)
    integer :: l, ME

    ME = 0
    do l=1,nSnp
      ME = ME + Ecnt(Genos(l,A), Genos(l, Par(1)), Genos(l, Par(2)))
    enddo
    CalcTrioErr = ME
  end function CalcTrioErr
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  subroutine dealloc_QPOsparse()
    if (allocated(QPOsparse_value))  deallocate(QPOsparse_value)
    if (allocated(QPOsparse_row))  deallocate(QPOsparse_row)
    if (allocated(QPOsparse_col))  deallocate(QPOsparse_col)
    if (allocated(PPO))  deallocate(PPO)   
  end subroutine dealloc_QPOsparse

end module OHfun

!===============================================================================

module CalcLik
  use Global, ONLY : Genos, Parent, nFS, FSID, maxSibSize, hermaphrodites, DumClone,&
    AKA2P, OKA2P
  implicit none

  contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine CalcLind(i)
    use Global
    implicit none

    integer, intent(IN) :: i
    integer :: nFSi, FSi(MaxSibSize), Inbr(2), Anc(2,mxA), PIK, l, x, y, k, z, j
    double precision :: PrL(nSnp), Px(3,2), PrG(3), PrXYZ(3,3,3), PrX(3)
    logical :: Selfed

    !if (real(COUNT(Genos(:,i)/=-1)) < nSnp/50.0) then
    !  return 
    !endif

    nFSi = 1 
    FSi = 0                
    if (Parent(i,1)<0 .or. Parent(i,2)<0) then
      nFSi = nFS(FSID(maxSibSize+1,i))
      FSi = FSID(1:maxSibSize, FSID(maxSibSize+1,i))  
    endif

    ! PO- and GP-mating
    Inbr = 0
    Anc = 0
    PIK = 0
    if (Parent(i,1)/=0 .and. Parent(i,2)/=0) then
      call getAncest(i,1,Anc)
      do k=1,2
        if (Anc(3-k,5-k) == Parent(i,3-k)) then
          Inbr(k) = 1
          if (Parent(i,3-k)>0) then
            PIK = Parent(i,3-k)
          endif
        endif
      enddo
    endif

    Selfed = .FALSE.
    if (hermaphrodites/=0) then
      if (Parent(i,1)==Parent(i,2) .and. Parent(i,1)>0) then
        Selfed = .TRUE.
      else if (all(Parent(i,:) < 0)) then
        if (DumClone(-Parent(i,1),1) == -parent(i,2)) then  
          Selfed = .TRUE.
        endif
      endif
    endif

    PrL = 0D0
    do l=1,nSnp
      do k=1,2
        call ParProb(l, Parent(i,k), k, i,-1, Px(:,k))
        if (Inbr(k)==1) then
          call ParProb(l, Anc(k,k+2), k, PIK,0, PrG)
        endif
      enddo
      PrXYZ = 0D0
      do y=1,3
        do z=1,3
          if (Selfed .and. y/=z)  cycle
          if (ANY(Inbr==1)) then
            do k=1,2
              if (Inbr(k)==-1) then
                PrXYZ(:,y,z) =AKA2P(:,y,z) * SUM(AKA2P(y,z,:)*PrG) *&
                 Px(y,k) * Px(z,3-k)
              endif
            enddo
          else if (Selfed) then   
            PrXYZ(:,y,y) = AKA2P(:, y, y) * Px(y,1)
          else
            PrXYZ(:,y,z) = AKA2P(:, y, z) * Px(y,1) * Px(z,2)
          endif 
          if (nFSi > 1) then
            do j=1, nFSi
              if (FSi(j) == i) cycle
              PrXYZ(:,y,z) = PrXYZ(:,y,z) * OKA2P(Genos(l,FSi(j)),y, z)
            enddo
          endif
        enddo
      enddo

      PrX = 0D0
      if(SUM(PrXYZ) > 0D0) then  
        do x=1,3
          PrX(x) = SUM(PrXYZ(x,:,:))/SUM(PrXYZ)
        enddo
        PrX = PrX / SUM(PrX)
      endif
      PrX = OcA(:,Genos(l,i)) * PrX
      PrL(l) = LOG10(SUM(PrX))
      LindX(:,l, i) = PrX
    enddo
    Lind(i) = SUM(PrL)

    if (Lind(i)> 0D0 .or. Lind(i)/=Lind(i)  .or. Lind(i) < -HUGE(1D0) .or. &   
      any(LindX(:,:,i)/=LindX(:,:,i))) then  
      call Rprint("",(/i,Parent(i,:)/), (/0D0/), "INT")
      if (any(LindX(:,:,i)/=LindX(:,:,i)))  call Rprint("LindX NA", (/0/), (/0D0/), "NON")
      call Erstop("Invalid individual LL - try increasing Err", .FALSE.)
    endif

  end subroutine CalcLind

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  function FSLik(l,i)    ! LL of FS set

    integer, intent(IN) :: l,i
    double precision :: FSLik(3,3), FSLikTMP(3,3)
    integer :: j, x, y

    FSLik = 1D0
    if (nFS(i)==0) then
      FSLik = 1D0
    else
      FSLikTMP = 1D0
      do j=1, nFS(i)
        do y=1,3
          do x=1,3
            FSLikTMP(x,y) = FSLikTMP(x,y) * OKA2P(Genos(l,FSID(j,i)), x, y)
          enddo
        enddo
      enddo
      if (any(FSLikTMP/=FSLikTMP) .or. any(FSLikTMP>1.0D0))  then 
        call Erstop("Invalid FS LL", .TRUE.)
      endif
      FSLik = FSLikTMP
    endif
    
  end function FSLik

end module CalcLik


!===============================================================================

module reshapers
implicit none
! reshape matrices/3D arrays into vectors, before passing to R

contains 
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine MtoVi(M, d1, d2, V)  ! IN, OUT: integer
    integer, intent(IN) :: d1, d2
    integer, intent(IN) :: M(d1, d2)
    integer, intent(OUT) :: V(d1*d2)
    integer :: i, x

    V = 0
    do i= 1, d1
      do x=1, d2
        V(i + (x-1)*d1) = M(i, x)
      enddo
    enddo
  end subroutine MtoVi

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  subroutine MtoVd(M, d1, d2, V)   ! IN, OUT: double
    integer, intent(IN) :: d1, d2
    double precision, intent(IN) :: M(d1,d2)
    double precision, intent(OUT) :: V(d1*d2)
    integer :: i, j

    V = 0D0
    do i=1,d1   
      do j=1,d2
        V((j-1)*d1 + i) = M(i,j)
      enddo
    enddo
  end subroutine MtoVd

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine AtoVi(A, d1, d2, x, V)   ! IN, OUT: integer
    integer, intent(IN) :: d1, d2, x(2)
    integer, intent(IN) :: A(d1,d2,2)
    integer, intent(OUT) :: V(d1*d2*2)
    integer :: i, j, k

    V = 0
    do j=1,d1
      do k=1,2  ! works for d(3)=2
        do i=1, x(k)
          V((j-1)*2*d2 + (k-1)*x(1) + i) = A(j, i, k)
        enddo
      enddo
    enddo
  end subroutine AtoVi

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine AtoVd(A, d1, d2, x, V)  ! IN, OUT: double
    integer, intent(IN) ::  d1, d2, x(2)
    double precision, intent(IN) :: A(d1,d2,2)
    double precision, intent(OUT) :: V(d1*d2*2)
    integer :: i, j, k

    V = 0D0
    do j=1,d1
      do k=1,2  ! works for d(3)=2
        do i=1, x(k)
          V((j-1)*2*d2 + (k-1)*x(1) + i) = A(j, i, k)
        enddo
      enddo
    enddo
  end subroutine AtoVd

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine AAtoVd(A, d1, d2, d3, V)  ! IN, OUT: double ; R-friendly RESHAPE
    integer, intent(IN) ::  d1, d2, d3
    double precision, intent(IN) :: A(d1,d2,d3)
    double precision, intent(OUT) :: V(d1*d2*d3)
    integer :: i, j, k, z

    V = 0D0
    do k=1, d3
      do j=1, d2
        do i=1, d1
          z = ((k-1)*d2 + (j-1))*d1 + i
          V(z) = A(i,j,k)
        enddo
      enddo
    enddo
  end subroutine AAtoVd

end module reshapers


! ####################################################################
! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,1997
! Made F conformant by Walt Brainerd

! Adapted by J Huisman (jisca.huisman@gmail.com) to output rank, to
! enable sorting of parallel vectors, and changed to decreasing rather 
! than increasing order

module sort_module
implicit none
private :: Partition

 contains
recursive subroutine QsortC(A, Rank)
  double precision, intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: Rank
  integer :: iq

  if(size(A) > 1) then
   call Partition(A, iq, Rank)
   call QsortC(A(:iq-1), Rank(:iq-1))
   call QsortC(A(iq:), Rank(iq:))
  endif
end subroutine QsortC

subroutine Partition(A, marker, Rank)
  double precision, intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: Rank
  integer, intent(out) :: marker
  integer :: i, j, TmpI
  double precision :: temp
  double precision :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1
  do
   j = j-1
   do
    if (j < 1) exit
    if (A(j) <= x) exit
    j = j-1
   end do
   i = i+1
   do
    if (i >= size(A)) exit
    if (A(i) >= x) exit
    i = i+1
   end do
   if (i < j) then
    ! exchange A(i) and A(j)
    temp = A(i)
    A(i) = A(j)
    A(j) = temp 
    
    TmpI = Rank(i) 
    Rank(i) = Rank(j)
    Rank(j) = TmpI
   elseif (i == j) then
    marker = i+1
    return
   else
    marker = i
    return
   endif
  end do

end subroutine Partition

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine getRank_i(BYrank)   ! get ranking based on generation number + birthyear (increasing)
use Global
implicit none

integer, intent(OUT) :: BYrank(nInd)
double precision :: SortBY(nInd)
integer :: i, BYx(nInd), Gen(nInd)

Gen = 0
call getGenerations(Gen) 

BYRank = (/ (i, i=1, nInd, 1) /)
SortBY = REAL(Gen, 8)
call QsortC(SortBy, BYRank)

BYx = BY
do i=1, nInd
  if (BY(i) < 0) then
    BYx(i) = MAXLOC(IndBY(:,i,5), DIM=1)  ! incl parents & offspring w unknown BY
  endif
enddo

SortBY = REAL(BYx, 8)
call QsortC(SortBy, BYRank)  ! do earlier birth years before later birth years

end subroutine getRank_i

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine sort_sibships(ss, kk)  ! sort sibships by birthyear (increasing) & size (decreasing)
use Global, ONLY: nC, nS, maxAgePO, DumBY, neg_missing

integer, allocatable, intent(OUT) :: ss(:), kk(:)
double precision, allocatable :: SortBY(:,:), SortBY_tmp(:)
integer, allocatable :: BYrankk(:,:), BYrank_tmp(:)
integer :: k, s, nCT,BYx, nsx, maxns, x1, x2, i

nCT = nC(1) + nC(2)
maxns = MAXVAL(ns)
allocate(SortBY(MAXVAL(nC),2), BYrankk(MAXVAL(nC),2))
BYrankk = 0
SortBY = 0.d0

do k=1,2 
  if (nc(k)==0)  cycle
  allocate(SortBY_tmp(nC(k)), BYrank_tmp(nC(k)))
  BYrank_tmp = (/ (s, s=1, nC(k))/)
  
  ! combine standardised est. birthyear of dummy & standardised sibship size
  do s=1, nC(k)
    BYx = MAXLOC(DumBY(:, s, k, 5), DIM=1)  ! with dummy parents & offspring
    nsx = ns(s,k) - maxns
    SortBY_tmp(s) = (REAL(BYx,8)/maxAgePO) * (REAL(nsx,8)/maxns)
  enddo
  
  ! sort within-sex  
  call QsortC(SortBY_tmp, BYRank_tmp)
  SortBY(1:nC(k), k)  = SortBY_tmp
  BYrankk(1:nC(k), k) = BYRank_tmp
  deallocate(SortBY_tmp, BYrank_tmp)
enddo

! combine
allocate(ss(nCT), kk(nCT))
if (all(nC > 0)) then
  x1 = 1  ! counter within BYrankk(:,1)
  x2 = 1  ! counter within BYrankk(:,2)
  do i=1,nCT
    if (SortbY(x1,1) > SortBY(x2,2)) then
      ss(i) = BYrankk(x1,1)
      kk(i) = 1
      if (x1 < nC(1)) then
        x1 = x1 +1
      else
        SortBY(nC(1),1) = neg_missing  ! don't do this x1 again
      endif
    else
      ss(i) = BYrankk(x2,2)
      kk(i) = 2
      if (x2 < nC(2)) then
        x2 = x2 +1
      else
        SortBY(nC(2),2) = neg_missing
      endif
    endif
  enddo
else if (nc(1) > 0) then
  ss = BYrankk(1:nc(1),1)
  kk = 1
else 
  ss = BYrankk(1:nc(2),2)
  kk = 2
endif


end subroutine sort_sibships

end module sort_module

! #####################################################################
 
subroutine Erstop(message, bug)
use Global
implicit none

character(len=*), intent(IN) :: message
logical, intent(IN) :: bug
 
call DeAllocAll
call intpr(" ", 1, (/0/), 0) 
!call labelpr(" ", 1)
if (bug) then
  call rexit("  ERROR! *** "//message//" *** Please report this bug.")
else
  call rexit("  ERROR! *** "//message//" *** ")
endif

end subroutine Erstop

! ####################################################################

! @@@@   PROGRAMS   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! ####################################################################

subroutine makeped(ng, nm, ny, specsintglb, specsintmkped, specsdbl, errv, genofr, &
  sexrf,  byrf, lyrf, aprf, mtdif_rf, parentsrf, lrrf, ohrf, & 
  nd, dumparrf, dumlrrf, dumbyrf, totll, apout)
use Global
use OHfun
use reshapers
implicit none

integer, intent(IN) :: ng, nm, ny, specsintglb(7), specsintmkped(3)
double precision, intent(IN) :: specsdbl(2), errv(9), aprf(5*ny)
integer, intent(IN) :: genofr(ng*nm), mtdif_rf(ng*ng)
integer, intent(INOUT) :: parentsrf(2*ng), ohrf(3*ng), sexrf(ng), byrf(3*ng), & 
  lyrf(ng), nd(2), dumparrf(2*ng), dumbyrf(3*ng)
double precision, intent(INOUT) :: lrrf(3*ng), dumlrrf(3*ng), totll(42), apout((3*ny+1)*5*3)
integer :: ParSib, i, CalcLLR, AgeEffect, k, s, x, m, j, IndBYmm(3), DumBYmm(3, ng/2, 2)
integer, parameter :: million = 1E6
double precision, allocatable :: AP_TMP(:,:,:), LLR_parent(:,:), LLR_GP(:,:,:)

allocate(LLR_parent(ng,3))
allocate(LLR_GP(3,ng/2,2))

ParSib = SpecsIntMkPed(1) 
AgeEffect = SpecsIntMkPed(2)
CalcLLR = SpecsIntMkPed(3)

call Initiate(Ng,Nm,Ny, SpecsIntGlb, SpecsDbl, ErrV, GenoFR, &
  SexRF, BYRF, LYRF, APRF, parentsRF, DumParRF)

if (quiet<1)  call rprint_status_tbl_header()

if(quiet==-1)  call rprint_tbl_update_a(0,90)  ! "count OH  ")
if(quiet==-1)  call rprint_status_tbl_no_dots() 
!call CalcOppHomAll()   ! also calls CalcPO
if (quiet==-1)  call rprint_status_tbl_eol()

if (any(mtdif_rf == 1)) then
  DoMtDif = .TRUE.
  allocate(mtDif(0:nInd,0:nInd))  ! different mt haplotypes
  mtDif = .FALSE.
  x = 0
  do i=1,nInd
    do j=1,nInd
      x = x+1
      mtDif(i,j) = mtdif_rf(x) == 1  
    enddo
  enddo
else
  DoMtDif = .FALSE.
endif

!=========================
! if (ParSib == 0) then
! ! call duplicates(nDupGenos,DupGenosFR, nMisMFR, DupGenoLR)  ! called directly from R
!=========================
if (ParSib == 1) then
  call parents(TotLL) 
!=========================  
else if (ParSib == 2) then
  call sibships(AgeEffect, TotLL)  
endif

!========================= 
! calculate estimated birth years + 95% CR
if (quiet==-1)  call rprint_tbl_update_a(99,91) ! "est byears")

BYRF = -9
do i=1, nInd
  if (quiet==-1 .and. any(chunk_printdot_i==i)) call print_dot()
  call EstBYrange(i, Sex(i), IndBYmm)
  BYRF(i) = IndBYmm(1)
  BYRF(i +nInd) = IndBYmm(2)
  BYRF(i +2*nInd) = IndBYmm(3)
enddo

SexRF = Sex   ! may include updated sex, if assigned as parent

DumBYmm = -9
if (ParSib == 2) then
  do k=1,2
    if (nC(k)==0)  cycle
    do s=1,nC(k)
      call EstBYrange(-s, k, DumBYmm(:, s,k))
    enddo
  enddo
  call AtoVi(DumBYmm, 3,nInd/2, nC, DumBYRF)
endif
if (quiet==-1)  call rprint_status_tbl_eol()
if (quiet<1)    call rprint_eol()  ! blank line

!========================= 
OhRF = -9
do i=1,nInd
  if (Parent(i,1)>0) OhRF(i) = calcOH(i, Parent(i,1))
  if (Parent(i,2)>0) OhRF(Ng+i) = calcOH(i, Parent(i,2))
  if (Parent(i,1)>0 .and. Parent(i,2)>0) then
    OhRF(2*Ng +i) = CalcTrioErr(i, Parent(i,:))
  endif
enddo

LLR_parent = missing
LLR_GP = missing
if(CalcLLR==1) then
  call rchkusr()
  call UpdateAllProbs() 
  if (quiet<1)  call Rprint("Calculating parent LLR ... ", (/0/), (/0.0D0/), "NON")
  call CalcParentLLR(LLR_parent, LLR_GP)
  if (quiet<1)  call rprint_eol()
  if (quiet<1)  call rprint_eol()  ! blank line
endif


if (any(DumClone /= 0)) then   ! renumber hermaphrodite sibshibs
  do k=1,2
    do x=1,nC(k)
      if (DumClone(x,k) == 0)  cycle
      !call Rprint('Dummy Clone: ', (/k, x, DumClone(x,k)/), (/0d0/), 'INT')
      do i=1, nInd
        if (Parent(i,k) == -x) then
          if (k==1)  Parent(i,k) = -x - million
          if (k==2)  Parent(i,k) = -DumClone(x,k) - million  ! use number of female clone
        endif
      enddo
      do m=1,2
        do s=1,nC(m)
          if (GpID(k,s,m) == -x) then
            if (k==1)  GpID(k,s,m) = -x - million
            if (k==2)  GpID(k,s,m) = -DumClone(x,k) - million  
          endif
        enddo
      enddo
      ! TODO: GPID of hermaphrodite dummies - drop? How to indicate?
    enddo
  enddo
endif

parentsRF = 0
LrRF = missing
do i=1,nInd 
  parentsRF(i) = Parent(i,1)
  parentsRF(nInd+i) = Parent(i,2)
  LrRF(i) = LLR_parent(i, 1)
  LrRF(nInd+i) = LLR_parent(i, 2)
  LrRF(2*nInd+i) = LLR_parent(i, 3)
enddo
if (ParSib == 2) then
  Nd = nC  ! no. dummies
  call AtoVi(GpID, 2,nInd/2, nC, DumParRF)
  call AtoVd(LLR_GP, 3,nInd/2, nC, DumLrRF) 
endif

allocate(AP_TMP(3*ny+1,5,3))   
AP_TMP = 0D0
! MaxAgePO <= ny by definition, as MaxAgePO is derived from input ageprior matrix
AP_TMP(1:(3*MaxAgePO+1), :, :) = AgePriorA(-MaxAgePO : (2*MaxAgePO), :, :)
call AAtoVd(AP_TMP, 3*ny+1, 5, 3, APout)

deallocate(AP_TMP)
deallocate(LLR_parent)
deallocate(LLR_GP)

call DeAllocAll

end subroutine makeped

! ######################################################################

subroutine Initiate(Ng, Nm, Ny, SpecsInt, SpecsDbl, ErrV, GenoFR, SexRF, BYRF, &
   LYRF, APRF, parentsRF, DumParRF)
use Global
use OHfun, only: init_OH
implicit none

integer, intent(IN) :: Ng, Nm, Ny, SpecsInt(7)
double precision, intent(IN) :: SpecsDbl(2), ErrV(9), APRF(5*Ny)
integer, intent(IN) :: GenoFR(Ng*Nm), SexRF(Ng), BYRF(3*Ng), LYRF(Ng), &
  parentsRF(2*Ng), DumParRF(2*Ng)
integer :: i,j,l,k,s, BYrange(Ng,2)
double precision :: AP_IN(Ny,5)  

! set global parameters 
zero = 0.0D0
nInd = Ng
nSnp = Nm
Complx = SpecsInt(1)
Hermaphrodites = SpecsInt(2)
quiet = SpecsInt(3)
maxSibSize = SpecsInt(4)
! MaxMisMatch = SpecsInt(5)  ! only used by duplicates()
MaxOppHom = SpecsInt(6)
MaxMendelE = SpecsInt(7)

TF = SpecsDbl(1)
TA = SpecsDbl(2)

allocate(Genos(nSnp, nInd))
Genos = -1
j = 0
do l=1,nSnp      
  do i=1, nInd
    j = j+1
    if (GenoFR(j)>=0) then
      Genos(l,i) = GenoFR(j)
    endif
  enddo    
enddo  

!=================
! allocate arrays  & prep
call AllocArrays()
call rchkusr()

Sex = SexRF
BY = BYRF(1:nInd)
do i=1,nInd
  BYrange(i,1) = BYRF(i + nInd)
  BYrange(i,2) = BYRF(i + 2*nInd)
enddo 
YearLast = LYRF
!WHERE(YearLast < 0)  YearLast = 999   done in PrepAgeData

AP_IN = 0.0D0
k = 0
do i=1,5
  do j=1, Ny
    k = k+1
    AP_IN(j,i) =  APRF(k)
  enddo
enddo

call PrepAgeData(Ny, AP_IN, BYrange)
call rchkusr()

call PrecalcProbs(ErrV)
call rchkusr()

call init_OH()
call rchkusr()

!=========================
! parents (pedigree prior or prev. assigned parents)
call ReadInputPed(ParentsRF, DumParRF)
call rchkusr()


if (Complx == 0 .and. any(Parent /=0))   call CheckMono()

if (hermaphrodites/=0) then
  do k=1,2
    if (nC(k)==0)  cycle
    do s=1, nC(k)
      call CheckSelfed(-s, k)    
    enddo
  enddo
  do i=1,nInd
    call CheckSelfed(i, Sex(i))
  enddo
endif

call UpdateAllProbs()


! helper array to print progress every 5%, 10%, ... of individuals
chunk_printdot_i =  mk_seq(nInd-1, 10)   ! several subroutines doing double loop, outer over i=1,nInd-1

contains
  subroutine PrecalcProbs(ErrV)

  double precision, intent(IN) :: ErrV(9)
  integer :: h,i,j,k,l,m
  double precision :: OjA(-1:2,3,nSnp), Tmp1(3), Tmp2(3,3), AF(nSnp)

  ! allele frequencies
  AF = 1D0            
  do l=1,nSnp
    if (ALL(Genos(l,:)==-1)) cycle
    AF(l)=dble(SUM(Genos(l,:), MASK=Genos(l,:)/=-1))/(COUNT(Genos(l,:)/=-1)*2)
  enddo


  !###################
  ! Prob. observed (rows) conditional on actual (columns)
  OcA(:, 0) = ErrV(1:3)  ! obs=0
  OcA(:, 1) = ErrV(4:6)  ! obs=1
  OcA(:, 2) = ErrV(7:9)  ! obs=2
  OcA(:,-1) = 1.0D0      ! missing  

  ! OcA(0:2, 1) = (/ (1-Er/2)**2, Er*(1-Er/2), (Er/2)**2 /)  ! act=0
  ! OcA(0:2, 2) = (/ Er/2, 1-Er, Er/2 /)  ! act=1
  ! OcA(0:2, 3) = (/ (Er/2)**2, Er*(1-Er/2),  (1-Er/2)**2 /)  ! act=2

  ! probabilities actual genotypes under HWE
  allocate(AHWE(3,nSnp))
  do l=1,nSnp
    AHWE(1,l)=(1 - AF(l))**2 
    AHWE(2,l)=2*AF(l)*(1-AF(l)) 
    AHWE(3,l)=AF(l)**2 
  enddo

  ! joined probabilities actual & observed under HWE
  do l=1,nSnp
    do i=-1,2    ! obs
      do j=1,3    ! act
        OjA(i, j, l) = OcA(j,i) * AHWE(j, l)
      enddo
    enddo
  enddo

  ! marginal prob. observed genotypes
  allocate(OHWE(-1:2,nSnp))
  do l=1,nSnp
    do i=-1,2
      OHWE(i, l) = SUM(OjA(i, :, l))
    enddo
  enddo

  ! ########################
  ! inheritance conditional on 1 parent
  allocate(AKAP(3,3,nSnp))
  allocate(OKAP(-1:2,3,nSnp))
  allocate(OKOP(-1:2,-1:2,nSnp))

  do l=1,nSnp
    AKAP(1, :, l) = (/ 1-AF(l), (1-AF(l))/2, 0.0D0 /)
    AKAP(2, :, l) = (/ AF(l), 0.5D0, 1-AF(l) /)
    AKAP(3, :, l) = (/ 0D0, AF(l)/2, AF(l) /)
  enddo

  OKAP = 1D0          
  do l=1,nSnp
    do i=0,2  ! obs offspring
      do j=1,3    ! act parent
        Tmp1=0D0
        do k=1,3    ! act offspring
          Tmp1(k) = OcA(k,i) * AKAP(k,j,l)
        enddo
        OKAP(i,j,l) = SUM(Tmp1)
      enddo
    enddo
  enddo

  OKOP = 1D0
  do l=1,nSnp
    do i=0,2  ! obs offspring
      do j=0,2    ! obs parent
        Tmp2=0D0
        do k=1,3    ! act offspring
          do m=1,3    ! act parent
            Tmp2(k,m) = OcA(k,i) * OcA(m,j) * AKAP(k,m,l)
          enddo
        enddo
        OKOP(i,j,l) = SUM(Tmp2) 
      enddo
    enddo
  enddo

  ! #########################
  ! inheritance conditional on both parents

  AKA2P(1,1,:) = dble((/ 1.0, 0.5, 0.0 /))
  AKA2P(1,2,:) = dble((/ 0.5, 0.25, 0.0 /))
  AKA2P(1,3,:) = dble((/ 0.0, 0.0, 0.0 /))

  AKA2P(2,1,:) = dble((/ 0.0, 0.5, 1.0 /))
  AKA2P(2,2,:) = dble((/ 0.5, 0.5, 0.5 /))
  AKA2P(2,3,:) = dble((/ 1.0, 0.5, 0.0 /))

  AKA2P(3,1,:) = dble((/ 0.0, 0.0, 0.0 /))
  AKA2P(3,2,:) = dble((/ 0.0, 0.25, 0.5 /))
  AKA2P(3,3,:) = dble((/ 0.0, 0.5, 1.0 /))

  OKA2P = 1D0
  do i=0,2  ! obs offspring
    do j=1,3    ! act parent 1
      do h=1,3    !act parent 2
        Tmp1=0D0
        do k=1,3    ! act offspring
          Tmp1(k) = OcA(k,i) * AKA2P(k,j,h) 
        enddo
        OKA2P(i,j,h) = SUM(Tmp1)
      enddo
    enddo
  enddo
  
  allocate(PPO(-1:2,-1:2,nSnp))
  PPO = 1D0
  do l=1,nSnp
    do i=0,2  ! obs offspring 
      do j=0,2    ! obs parent
        PPO(i,j,l) = OKOP(i,j,l) / OHWE(i,l)
      enddo
    enddo
  enddo  

  allocate(PHS(-1:2,-1:2,nSnp))
  allocate(PFS(-1:2,-1:2,nSnp))
  PHS = 1D0
  PFS = 1D0
  do l=1,nSnp
    do i=0,2  ! obs offspring 1
      do j=0,2    ! obs offspring 2
        Tmp1=0D0
        Tmp2=0D0
        do m=1,3    !act shared parent 
          Tmp1(m) = OKAP(i,m,l) * OKAP(j,m,l) * AHWE(m,l)
          do h=1,3
            Tmp2(m,h) = OKA2P(i,m,h) * OKA2P(j,m,h) * AHWE(m,l) *AHWE(h,l)
          enddo
        enddo
        PHS(i,j,l) = SUM(Tmp1) / (OHWE(i,l) * OHWE(j,l))
        PFS(i,j,l) = SUM(Tmp2) / (OHWE(i,l) * OHWE(j,l))
      enddo
    enddo
  enddo

  allocate(DumP(3,nSnp, nInd/2,2))
  allocate(XPr(3,3,nSNP, nInd/2,2))
  allocate(LindX(3,nSnp,nInd))  ! used when missing genotype & at start
  XPr = 1D0
  do l=1,nSnp
    do i=1,3  
      LindX(i,l,:) = AHWE(i,l)  
      DumP(i,l,:,:) = AHWE(i,l)
      XPr(2,i,l,:,:) = AHWE(i,l)  ! GP contribution
    enddo
  enddo

  end subroutine PrecalcProbs


end subroutine Initiate

! ####################################################################

subroutine AllocArrays()
use Global
implicit none

integer :: i

allocate(Parent(nInd,2))
Parent = 0
allocate(Lind(nInd))
Lind = 0D0

nC = 0  
allocate(nS(nInd/2,2))
ns = 0
allocate(SibID(maxSibSize, nInd/2, 2))
SibID = 0 
allocate(GpID(2, nInd/2,2))
GpID = 0
allocate(CLL(nInd/2,2))
CLL = missing
allocate(IsNewSibship(nInd/2, 2))
IsNewSibship = .TRUE.
allocate(ToCheck(nInd))
ToCheck = .FALSE.
allocate(SelfedIndiv(nInd))
SelfedIndiv = .FALSE.

allocate(NFS(nInd))
NFS = 1
allocate(FSID(MaxSibSize+1, nInd))
FSID = 0
FSID(1, :) = (/ (i, i=1, nInd) /)
FSID(MaxSibSize+1, :) = (/ (i, i=1, nInd) /)    ! 'primary' sib

allocate(Mate(nInd))
Mate = 0  
allocate(DumMate(nInd/2,2))
DumMate = 0  
allocate(DumClone(nInd/2,2))
DumClone = 0  

allocate(Sex(nInd))
Sex = 3
allocate(BY(nInd))
BY = -999
allocate(YearLast(nInd))
YearLast = +999                   

end subroutine AllocArrays

! ####################################################################

subroutine ReadInputPed(parentV, dumParV)
use Global
implicit none

integer, intent(IN) :: parentV(2*nInd), dumParV(2*nInd)  
integer :: i, ParTmp(2), k, s, j, x

! reset all parents & sibships
Parent = 0  
nC = 0
ns = 0
SibID = 0
GpID = 0
NFS = 1
FSID = 0
FSID(1, :) = (/ (i, i=1, nInd) /)
FSID(MaxSibSize+1, :) = (/ (i, i=1, nInd) /) 

if (all(ParentV==0))  return

if(quiet<1)  call Rprint("Transferring input pedigree ...", (/0/), (/0.0D0/), "NON")

do i=1,nInd
  ParTmp(1) = parentV(i)
  ParTmp(2) = parentV(nInd + i)
  do k=1,2
    if (ParTmp(k) > 0) then
      Parent(i,k) = ParTmp(k)
    else if (ParTmp(k) < 0) then
      s = -ParTmp(k)
      if (nC(k) < s)  nC(k) = s   
        Parent(i,k) = ParTmp(k)
        nS(s,k) = ns(s,k) +1
        SibID(ns(s,k), s, k) = i
    endif
  enddo
enddo

do k=1,2
  DumBY(:,1:nC(k),k,:) = LOG10(1d0/nYears)
enddo

do i=1, nInd
  if (Sex(i)==3) then
    if (ANY(Parent(:,1) == i)) then
      Sex(i) = 1
    else if (ANY(Parent(:,2) == i)) then
      Sex(i) = 2
    endif                                     
  endif
enddo    

! find current FS 
do i=1,nInd-1
  do j=i+1,nInd
    if (Parent(i,1)==Parent(j,1) .and. Parent(i,1)/=0 .and. &
      Parent(i,2)==Parent(j,2) .and. Parent(i,2)/=0) then
      call MakeFS(i, j)
    endif
  enddo
enddo

if (all(DumParV == 0))  return
do k=1,2
  if (nC(k)==0)  cycle
  do s=1, nC(k)
    do i=1,2
      x = (k-1)*2*INT(nInd/2) + (s-1)*2 + i  ! reverse  AtoVi(GpID, 2,nInd/2, nC, DumParRF)
      GpID(i, s, k) = DumParV(x)
    enddo
  enddo
enddo

call rchkusr()

end subroutine ReadInputPed

! ####################################################################

subroutine CheckMono
use Global
implicit none

integer :: i, j, k, nOff, Off(maxSibSize), sxOff(maxSibSize), s, par
logical :: ParOK

nOff = 0
Off = 0
sxOff = 3
par = 0
ParOK = .TRUE.

do i=1, nInd
  do k=1,2
    call getOff(i, k, .FALSE., nOff, Off, sxOff)
    if (nOff == 0)  cycle
    do j=1, nOff
      if (Parent(Off(j),3-k) /= 0) then
        Mate(i) = Parent(Off(j),3-k)
        exit
      endif
    enddo
    if (Mate(i)==0) then
      call NewSibship(Off(1), 0, 3-k)   ! sets Mate(i) = -nC(3-k)
    endif
    do j=1, nOff
      if (Parent(Off(j),3-k) /= Mate(i) .and. Parent(Off(j), 3-k)/=0) then
        call Erstop("Please change to Complex='simp', assigned parents suggest non-monogamy", .FALSE.)
      else if (Parent(Off(j),3-k) == 0) then
        call ChkValidPar(Off(j), sxOff(j), Mate(i), 3-k, ParOK)
        if (ParOK) then
          call SetPar(Off(j), sxOff(j), Mate(i), 3-k)  ! make all half-sibs full sibs
        else 
          call SetPar(Off(j), sxOff(j), 0, k)
        endif
      endif
    enddo
    ! dummy offspring
    call getOff(i, k, .TRUE., nOff, Off, sxOff)
    do j=1, nOff
      if (Off(j) >0) cycle
      if (GpID(3-k, -Off(j), sxOff(j)) /= Mate(i) .and. GpID(3-k, -Off(j), sxOff(j)) /= 0) then
        call Erstop("Please change to Complex='simp', assigned (grand)parents suggest non-monogamy", .FALSE.)
      endif
      call ChkValidPar(Off(j), sxOff(j), Mate(i), 3-k, ParOK)
      if (ParOK) then
        call setPar(Off(j), sxOff(j), Mate(i), 3-k) 
      else
        call setPar(Off(j), sxOff(j), 0, k)   
      endif
    enddo
  enddo
enddo

if (quiet==-1 .and. any(Mate/=0)) then
  call Rprint("Added dummy parents to ensure monogamy...", (/0/), (/0.0D0/), "NON")
endif

do k=1,2
  if (nC(k) == 0)  cycle
  do s=1, nC(k)
    call getFSpar(s, k, .TRUE., par)
    if (par /= 0) then
      DumMate(s,k) = par
      if (any(Parent(SibID(1:ns(s,k),s,k), 3-k) == 0)) then
        do j=1, ns(s,k)
          i = SibID(j,s,k)
          call ChkValidPar(i, 3, Par, 3-k, ParOK)
          if (ParOK) then
            call setPar(i, 3, Par, 3-k)
          else
            call setPar(i, 3, 0, k)  ! remove from sibship
          endif
        enddo
      endif
    else if (any(Parent(SibID(1:ns(s,k),s,k), 3-k) /= 0)) then
      call Erstop("Please change to Complex='simp', assigned (dummy) parents suggest non-monogamy", .FALSE.)
    else
      call NewSibship(SibID(1,s,k), SibID(2,s,k), 3-k)   ! also works if ns=1, and SibID(2)=0
      if (ns(s,k) > 2) then
        do j=3, ns(s,k)
          call setPar(SibID(j,s,k), 3, -nC(3-k), 3-k)         
        enddo
      endif
      ToCheck(SibID(1:ns(s,k),s,k)) = .TRUE.
      DumMate(s,k) = Parent(SibID(1,s,k), 3-k)
    endif
    ! dummy offspring
    call getOff(i, k, .TRUE., nOff, Off, sxOff)
    do j=1, nOff
      if (Off(j) >0) cycle
      if (GpID(3-k, -Off(j), sxOff(j)) /= DumMate(s,k) .and. GpID(3-k, -Off(j), sxOff(j)) /= 0) then
        call Erstop("Please change to Complex='simp', assigned (grand)parents suggest non-monogamy", .FALSE.)
      endif
      call ChkValidPar(Off(j), sxOff(j), DumMate(s,k), 3-k, ParOK)
      if (ParOK) then
        call setPar(Off(j), sxOff(j), DumMate(s,k), 3-k) 
      else
        call setPar(Off(j), sxOff(j), 0, k)   
      endif
    enddo  
  enddo
enddo

end subroutine CheckMono

! ####################################################################

subroutine duplicates(ng,nm,ny, specsint, specsdbl, errv, dupratio, genofr, &  ! in
  sexrf, byrf, aprf, & ! fake, empty in
  ndupgenos, dupgenos, nmismatch, snpdboth, duplr)  ! OUT
use Global
implicit none

integer, intent(IN) :: ng,nm,ny, specsint(7), dupratio
double precision, intent(IN) :: specsdbl(2), errv(9), aprf(5*ny)
integer, intent(IN) :: genofr(ng*nm), sexrf(ng), byrf(3*ng) 
integer, intent(INOUT) :: ndupgenos, dupgenos(2*dupratio*ng), nmismatch(dupratio*ng), snpdboth(dupratio*ng)
double precision, intent(INOUT) :: duplr(dupratio*ng)
integer :: i, j, l, CountMismatch, MaxMisMatchDup, &
   parentsRF(2*Ng), DumParRF(2*Ng), LYRF(Ng)  ! fake & empty
integer :: IsBothScored(-1:2,-1:2), IsDifferent(-1:2,-1:2), SnpdBoth_ij, i_quadraginta(40)
double precision :: LL(7), LLtmp(2), LLX(7)

parentsRF = 0
DumParRF = 0 
LYRF = -999

call Initiate(Ng,Nm,Ny, SpecsInt, SpecsDbl, ErrV, GenoFR, &
  SexRF, BYRF, LYRF, APRF, parentsRF, DumParRF) 
  
MaxMismatchDup = SpecsInt(5)

!====================
! (nearly) identical genotypes?                               
nDupGenos = 0
DupGenos = -9
nMismatch = -9
SnpdBoth = -9
DupLR = missing

IsBothScored = 1
IsBothScored(-1,:) = 0
IsBothScored(:,-1) = 0
IsDifferent = 0
IsDifferent(0, 1:2) = 1
IsDifferent(1, (/0,2/)) = 1
IsDifferent(2, 0:1) = 1

i_quadraginta = mk_seq(nInd-1, 40)
if (quiet==-1)  call rprint_progbar_header()

LLtmp = missing
LL = missing
LLX = missing   
do i=1,nInd-1
  if (MODULO(i,100)==0) call rchkusr()
  if (quiet==-1 .and. any(i_quadraginta==i)) call rprint_progbar_dot()
  do j=i+1, nInd
    CountMismatch = 0
    SnpdBoth_ij = 0               
    do l=1, nSnp
      SnpdBoth_ij = SnpdBoth_ij + IsBothScored(Genos(l,i), Genos(l,j))
      CountMismatch = CountMismatch + IsDifferent(Genos(l,i), Genos(l,j))
      if (CountMismatch > MaxMismatchDup)  exit
    enddo
    if (CountMismatch > MaxMismatchDup)  cycle
    LLtmp = missing               
    call PairSelf(i, j, LLtmp(1))
    call CheckPair(i, j,3,7,LL, LLX)   
    LLtmp(2) = MaxLL(LL)
    if ((LLtmp(1) - LLtmp(2)) > TF)  then   ! what threshold?
      nDupGenos = nDupGenos + 1
      DupGenos(nDupGenos) = i
      DupGenos(nDupGenos + dupratio*ng) = j
      nMisMatch(nDupGenos) = CountMismatch
      SnpdBoth(nDupGenos) = SnpdBoth_ij
      DupLR(nDupGenos) = LLtmp(1) - LLtmp(2)
    endif
    if (nDupGenos==dupratio*ng) then
      if(quiet<1) call rwarn("reached max for duplicates") 
      exit
    endif
    if (nDupGenos==dupratio*ng) exit
  enddo
  if (nDupGenos==dupratio*ng) exit
enddo

if (quiet==-1)  call rprint_eol()
call DeAllocAll    

 end subroutine duplicates
 
! ###################################################################### 

subroutine findambig(ng, nm, ny, specsint, specsintamb, specsdbl, errv, genofr, &
    sexrf, byrf, aprf, parentsrf, dumparrf, &
    namb, ambigid, ambigrel, ambiglr, ambigoh, &
    ntrio, trioids, triolr, triooh)
use Global
use OHfun
implicit none

integer, intent(IN) :: ng, nm, ny, specsint(7), specsintamb(2)
double precision, intent(IN) :: specsdbl(2), errv(9), aprf(5*ny)
integer, intent(IN) :: sexrf(ng), byrf(3*ng), genofr(ng*nm)
integer, intent(INOUT) :: parentsrf(2*ng), dumparrf(2*ng), &
  namb, ambigid(2*specsintamb(2)), ambigrel(2*specsintamb(2)), &
  ambigoh(specsintamb(2)), ntrio, trioids(3*ng), triooh(3*ng)
double precision, intent(INOUT) :: ambiglr(2*specsintamb(2)), triolr(3*ng) 
integer :: ParSib, nAmbMax, i, j, k, x, topX, Anci(2,mxA), Ancj(2,mxA), BYtmp(2), maybe, &
   Lboth, LYRF(Ng), chunk_printdot(40)
double precision :: LL(7), LLtmp(7,3), dLL, LRR(3), LLX(7)
logical :: ParOK(2)

LYRF = -999

call Initiate(Ng,Nm,Ny, SpecsInt, SpecsDbl, ErrV, GenoFR, SexRF,  BYRF, LYRF, &  ! IN
  APRF, parentsRF, DumParRF) 

if(quiet<1)  call Rprint("Counting opposing homozygous loci between all individuals ... ", &
  (/0/), (/0.0D0/), "NON")
!call CalcOppHomAll()   ! also calls CalcPO

ParSib = SpecsIntAmb(1)  ! TODO: pass separately?  
nAmbMax = SpecsIntAmb(2)
  
nAmb = 0
AmbigID = 0
AmbigLR = missing
AmbigOH = -9
AmbigRel = 0   

Anci = 0
Ancj = 0

if(quiet<1) then
  if (ParSib ==1) then
    call Rprint("Checking for non-assigned Parent-Offspring pairs ... ",(/0/), (/0.0D0/), "NON")
  else if (ParSib ==2) then
    call Rprint("Checking for non-assigned relatives ... ",(/0/), (/0.0D0/), "NON")
  endif
endif

if (quiet<1)  call rprint_progbar_header()
chunk_printdot = mk_seq(nInd-1,40)    ! NOTE: does not work well for nInd<41

do i=1,nInd-1
  if (nAmb==nAmbMax)  exit   
  if (MODULO(i,200)==0)   call rchkusr()
   if (quiet<1 .and. any(chunk_printdot==i)) call rprint_progbar_dot()   ! asterisks
  call GetAncest(i,1,Anci)
  
  do j=i+1,nInd    
    Lboth = COUNT(Genos(:,i)/=-1 .and. Genos(:,j)/=-1)  
    if (Lboth < nSnp/2.0)   cycle   ! >1/2th of markers missing
    if (ANY(Parent(i,:)==j) .or. ANY(Parent(j,:)==i)) cycle  ! PO
    if (ALL(Parent(i,:)/=0)) then
      if (Parent(i,1)==Parent(j,1) .and. Parent(i,2)==Parent(j,2)) cycle  ! FS
      if (ANY(Anci(:,3)==j) .and. ANY(Anci(:,4)==j)) cycle  ! double GP
    endif
    if (Parent(j,1)/=0 .and. Parent(j,2)/=0) then
      call GetAncest(j,1,Ancj)
      if (ANY(Ancj(:,3)==i) .and. ANY(Ancj(:,4)==i)) cycle  ! double GP
    endif
    
    LLtmp = missing
    LL = missing
    topX = 0
    dLL = missing
    ParOK = .FALSE.         
    if (ParSib==1) then   ! check if they're not PO only
      if (QLR_PO(i,j) < TA)   cycle   
    else if (All(Parent(i,:)/=0) .and. ALL(Parent(j,:)/=0)) then
      if (QLR_PO(i,j) < 2*TA)   cycle  
    endif
    
    if (ParSib==1) then
      BYtmp(1:2) = BY((/i,j/))
      BY((/i,j/)) = -9
      call ChkValidPar(i,Sex(i), j,Sex(j), ParOK(1))
      call ChkValidPar(j,Sex(j), i,Sex(i), ParOK(2))
      if (.not. ANY(ParOK))  cycle     
      if (ParOK(1))  call CheckPair(i, j, Sex(j), 1, LLtmp(:,1), LLX)
      if (ParOK(2))  call CheckPair(j, i, Sex(i), 1, LLtmp(:,2), LLX)
      do k=1,7  
        LL(k) = MaxLL(LLtmp(k,1:2)) 
      enddo
      call BestRel2(LL, topX, dLL)
      BY((/i,j/)) = BYtmp(1:2)
      if (topX==6 .or. topX==7) cycle   ! conditionally unrelated
      if (COUNT(Parent == 0) > 0.95*nInd .and. topX>2)  cycle  ! else will exceed nAmbMax
    
    else if (ParSib == 2) then 
      maybe = 0
      LRR = missing
      LRR(1) = QLR_PO(i,j)
      topX = 0
      do k=1,2 
        if (Parent(i,k)/=0 .and. Parent(i,k)==Parent(j,k)) cycle
        if (Parent(i,k)>0) then
            if (ANY(Parent(Parent(i,k), :)==j)) cycle
        else if (Parent(i,k)<0) then
            if (ANY(GpID(:, -Parent(i,k), k)==j)) cycle
        endif
        if (Parent(j,k)>0) then
            if (ANY(Parent(Parent(j,k), :)==i)) cycle
        else if (Parent(j,k)<0) then
            if (ANY(GpID(:, -Parent(j,k), k)==i)) cycle
        endif
        call PairQFS(i, j, LRR(2)) 
        call PairQHS(i, j, LRR(3))       
        maybe = 0
        do x=1,3
          if (LRR(x) > 2*TA .and. LRR(x) < missing)  maybe=1  
        enddo
        if (maybe==0)  cycle
        if (BY(i) > BY(j)) then
          call CheckPair(i, j, k, 7, LL, LLX)  
        else
          call CheckPair(j, i, k, 7, LL, LLX)
        endif
        call BestRel2(LL, topX, dLL) 
        if (COUNT(Parent == 0) > 0.95*nInd .and. topX>2) then
          cycle  ! else will exceed nAmbMax)
        else if (topX==6 .or. topX==7) then  ! .or. dLL(2)<TA   .or. topX==8
          maybe = 0
          cycle
        else
          exit
        endif
      enddo
      if (maybe==0) cycle
    endif
    
    nAmb = nAmb + 1
    AmbigID(nAmb) = i
    AmbigID(nAmbMax + nAmb) = j
    AmbigOH(nAmb) = calcOH(i,j, nSnp)
    if (ParSib==1) then
      AmbigLR(nAmb) = QLR_PO(i,j)
      AmbigRel(nAmb) = 1
    else if (ParSib==2) then
      AmbigLR(nAmb) = MAXVAL(LRR, MASK=LRR<MaybeOtherParent)
      AmbigRel(nAmb) = MAXLOC(LRR, MASK=LRR<MaybeOtherParent, DIM=1)
    endif
    AmbigRel(nAmbMax + nAmb) = TopX
    AmbigLR(nAmbMax + nAmb) = dLL 
    if (nAmb==nAmbMax) then
      if(quiet<1) then
!        call rwarn("reached max for maybe-rel")  ! doesn't return any output
        if (ParSib == 1) then
          call Rprint("WARNING - reached max for maybe-par, truncated!",(/0/), (/0.0D0/), "NON")
        else
          call Rprint("WARNING - reached max for maybe-rel, truncated!",(/0/), (/0.0D0/), "NON")
        endif
      endif
      exit
    endif
  enddo
enddo

if (quiet<1)  call rprint_eol()

! triads
ntrio = 0
trioLR = missing 
trioOH = -9  

if (nAmb>1) then   
  if (COUNT(AmbigRel((nAmbMax+1) : (2*nAmbMax)) == 1) > 1) then
    if(quiet<1) then
      call Rprint("Checking for Parent-Parent-Offspring trios ... ",(/0/), (/0.0D0/), "NON")
    endif
    call triads(nAmbMax, AmbigID, AmbigRel((nAmbMax+1) : (2*nAmbMax)),&
     ntrio, trioIDs, trioLR, trioOH)
  endif
endif

call DeAllocAll

end subroutine findambig

! #####################################################################

subroutine triads(nAmbMax, AmbigIDIN, topRel, ntrio, trioIDsOUT, trioLROUT, trioOHOUT)
use Global
use OHfun
use reshapers
implicit none

integer, intent(IN) :: nAmbMax, AmbigIDIN(2*nAmbMax), topRel(nAmbMax)
integer, intent(OUT) :: ntrio, trioIDsOUT(3*nInd), trioOHOUT(3*nInd)
double precision, intent(OUT) :: trioLROUT(3*nInd)
integer :: i, j, u, v, ncp, CandPar(mxCP), m, nAmb, &
 AmbigID(nAmbMax,2), trioIDs(nInd, 3), trioOH(nInd, 3) 
double precision :: trioLR(nInd, 3), LLP(3)
logical :: AncOK

AmbigID(:,1) = AmbigIDIN(1:nAmbMax)
AmbigID(:,2) = AmbigIDIN((nAmbMax+1) : (2*nAmbMax))
nAmb = COUNT(AmbigID(:,1)>0)      
ntrio = 0
trioIDs = 0
trioLR = missing
trioOH = -9
LLP = missing
AncOK = .TRUE.

do i=1, nInd
  if (MODULO(i,500)==0)   call rchkusr()
  if (ntrio == nInd) exit
  if (ANY(Parent(i,:)/=0) .and. Hermaphrodites==0) cycle
  if (ANY(Parent(i,:)>0) .and. Hermaphrodites>0) cycle  
  ncp = 0
  CandPar = 0  
  if ((COUNT(AmbigID(:,1) == i .and. topRel <3) + &   ! PO or FS
    COUNT(AmbigID(:,2) == i .and. topRel <3)) < 2) cycle
  do j=1, nAmb
    if (topRel(j) >2)  cycle 
    if (.not. ANY(AmbigID(j,:) == i))  cycle
    if (ncp == mxCP) exit
    do m=1,2
      if (AmbigID(j,m) == i) then
        if (AgeDiff(i, AmbigID(j,3-m))<=0)  cycle  ! unknown=999. 
        call ChkAncest(AmbigID(j,3-m), sex(AmbigID(j,3-m)), i, sex(i), AncOK)
        if (.not. AncOK)  cycle
        ncp = ncp + 1
        CandPar(ncp) = AmbigID(j,3-m)
      endif
    enddo
  enddo
  
  if (ncp > 1) then   
    do u=1, ncp-1
      do v=u+1, ncp
        if (Sex(CandPar(u))<3 .and. Sex(CandPar(u))==Sex(CandPar(v))) cycle
        if (Sex(CandPar(u))==1 .or. Sex(CandPar(v))==2) then
          call CheckParentPair(i, Sex(i), CandPar( (/u,v/) ), LLP)
        else
          call CheckParentPair(i, Sex(i), CandPar( (/v,u/) ), LLP)
        endif
        if (LLP(3) < -TA .or. LLP(3)==missing)   cycle
        ntrio = ntrio +1
        trioIDs(ntrio,1) = i
        trioIDs(ntrio, 2:3) = CandPar( (/u,v/) )
        trioLR(ntrio,:) = LLP
        
        trioOH(ntrio,1) = calcOH(i, CandPar(u))
        trioOH(ntrio,2) = calcOH(i, CandPar(v))
        trioOH(ntrio,3) = CalcTrioErr(i, CandPar( (/u,v/) ) )
        if (ntrio == nInd) exit        
      enddo
      if (ntrio == nInd) exit
    enddo
  endif
enddo

! fold into vectors
call MtoVi(trioIDs, nInd, 3, trioIDsOUT)
call MtoVi(trioOH, nInd, 3, trioOHOUT)
call MtoVd(trioLR, nInd, 3, trioLROUT)

end subroutine triads

! ######################################################################

subroutine getpedllr(ng, nm, ny, specsint, specsintmkped, specsdbl, errv, genofr, &
 sexrf, byrf, aprf, parentsrf, ohrf, lrrf, snpdboth, dumparrf, dumlrrf, dumbyrf)
use Global
use OHfun
use reshapers
implicit none

integer, intent(IN) :: ng, nm, ny, specsint(7), specsintmkped(2)
double precision, intent(IN) :: specsdbl(2), errv(9), aprf(5*ny)
integer, intent(IN) :: genofr(ng*nm), parentsrf(2*ng), dumparrf(2*ng) 
integer, intent(INOUT) :: ohrf(3*ng), snpdboth(2*ng), sexrf(ng), byrf(3*ng), dumbyrf(3*ng)
double precision, intent(INOUT) :: lrrf(3*ng), dumlrrf(3*ng)
integer :: i,j,l, IndBYmm(3), k, s, CalcLLR, DumBYmm(3, ng/2, 2), LYRF(ng)
double precision, allocatable :: LLR_Parent(:,:), LLR_GP(:,:,:)

allocate(LLR_Parent(ng,3))
allocate(LLR_GP(3, ng/2, 2))
 
!AgeEffect = SpecsIntMkPed(1)
CalcLLR = SpecsIntMkPed(2)
LYRF = -999

if (CalcLLR == 1) then
  call Initiate(Ng,Nm,Ny, SpecsInt, SpecsDbl, ErrV, GenoFR, SexRF,  BYRF, LYRF, &  
    APRF, parentsRF, DumParRF)

else   ! only initiate minimally necessary
  nInd = Ng
  nSnp = Nm
  quiet = SpecsInt(3)

  allocate(Genos(nSnp, nInd))
  Genos = -1
  j = 0
  do l=1,nSnp      
    do i=1, nInd
      j = j+1
      if (GenoFR(j)/=-9) then
        Genos(l,i) = GenoFR(j)                           
      endif
    enddo    
  enddo  
  
  allocate(Parent(nInd,2))
  Parent = 0
  do i=1,nInd
    Parent(i,1) = parentsRF(i)
    Parent(i,2) = parentsRF(nInd + i)
  enddo
endif

!~~ calculate OH & ME ~~
SnpdBoth = -9
OhRF = -9
!MaxOppHom = nSnp
if (quiet<1)  call Rprint("Counting Mendelian errors ... ",(/0/), (/0.0D0/), "NON")
do i=1,nInd
  if (MODULO(i,500)==0)  call rchkusr()
  do k=1,2
    j = i + (k-1)*Ng
    if (Parent(i,k)>0) then
      SnpdBoth(j) = COUNT(Genos(:,i)/=-1 .and. Genos(:,Parent(i,k))/=-1)
      OhRF(j) = calcOH(i, Parent(i,k))
    else if (Parent(i,k) < 0) then
      SnpdBoth(j) = 0
    endif
  enddo
  if (all(Parent(i,:)>0)) then
    OhRF(2*Ng +i) = CalcTrioErr(i, Parent(i,:))
  endif
enddo


!~~ calculate LLR Parent/Not ~~~
LrRF = missing
LLR_Parent = missing
LLR_GP = missing
BYRF = -9
DumBYRF = -9
  
if (CalcLLR == 1) then
  if(quiet<1)  call Rprint("Counting opposing homozygous loci between all individuals ... ", &
    (/0/), (/0.0D0/), "NON")
!  call CalcOppHomAll()   ! also calls CalcPO
  
  if (quiet<1)  call Rprint("Calculating parent LLR ... ",(/0/), (/0.0D0/), "NON")
  call UpdateAllProbs()
  call rchkusr()
  call CalcParentLLR(LLR_parent, LLR_GP)
  if (quiet<1)  call rprint_eol()

  do i=1,nInd 
    LrRF(i) = LLR_parent(i, 1)
    LrRF(nInd+i) = LLR_parent(i, 2)
    LrRF(2*nInd+i) = LLR_parent(i, 3)
  enddo
  
  call AtoVd(LLR_GP, 3,nInd/2, nC, DumLrRF) 
  
  SexRF = Sex   ! may include updated sex, if assigned as parent
  
  ! estimated birth years + 95% CR
  DumBYmm = -9
  do i=1, nInd
    call EstBYrange(i, Sex(i), IndBYmm)
    BYRF(i) = IndBYmm(1)
    BYRF(i +nInd) = IndBYmm(2)
    BYRF(i +2*nInd) = IndBYmm(3)
  enddo
  if (any(nC > 0)) then
    do k=1,2
      if (nC(k)==0)  cycle
      do s=1,nC(k)
        call EstBYrange(-s, k, DumBYmm(:, s,k))
      enddo
    enddo
    call AtoVi(DumBYmm, 3,nInd/2, nC, DumBYRF)
  endif
endif

deallocate(LLR_Parent)
deallocate(LLR_GP)
call DeAllocAll

end subroutine getpedllr

! ######################################################################

subroutine getpairll(ng, nm, ny, np, specsint, specsdbl, errv, nrels, genofr, &
  byrf, aprf, &
  pairids, pairsex, pairagediff, pairfocal, pairk, dropp, parentsrf, dumparrf, llrf)
use Global
use reshapers
!use OHfun, only: CalcOppHomAll
implicit none

integer, intent(IN) :: ng, nm, ny, np, specsint(7), nrels
double precision, intent(IN) :: specsdbl(2), errv(9), aprf(5*ny)
integer, intent(IN) :: genofr(ng*nm), parentsrf(2*ng), dumparrf(2*ng), &
  pairids(2*np), pairsex(2*np), pairagediff(np), pairfocal(np), pairk(np), dropp(2*np)
integer, intent(INOUT) :: byrf(3*ng)
double precision, intent(INOUT) :: llrf(nrels*np)
integer :: x, ij(2), kij(2), a, m, SexRF(ng), LYRF(ng), curPar(2,2), &
  Sex_ij(2), psex(np, 2), pdrop(np, 2), BYtmp(2), chunk_printdot(40), mates_tmp(2,2)
double precision :: LLpair(np, nrels), LLg(7), LLa(7), LLdup
double precision, allocatable, dimension(:,:,:) :: IndBYtmp
logical :: no_can_do, tmp_drop_mates(2)

SexRF = 3
LYRF = -999

call Initiate(Ng,Nm,Ny, SpecsInt, SpecsDbl, ErrV, GenoFR, SexRF,  BYRF, LYRF, &  
  APRF, parentsRF, DumParRF)       

allocate(IndBYtmp(nYears,2,5))   

! fold vectors into 2-column matrices
psex(:,1) = pairsex(1:np)
psex(:,2) = pairsex((np+1) : (2*np))
pdrop(:,1) = dropp(1:np)
pdrop(:,2) = dropp((np+1) : (2*np))

if (nP > 50 .and. quiet < 1) then
  call rprint_progbar_header()  
  chunk_printdot =  mk_seq(nP, 40)
else
  chunk_printdot = 0
endif

LLpair = Missing
do x = 1, nP
  if (MODULO(x,500)==0) call rchkusr()
  if (quiet<1 .and. nP>50 .and. any(chunk_printdot == x))   call rprint_progbar_dot() 

  ij = pairids((/x, nP+x/))
  if (any(ij == 0))  cycle   ! one or both not genotyped & not dummy

  ! check if sex & age difference specified in pairs file do not conflict with pedigree
  no_can_do = .FALSE.
  do a=1,2
    if (ij(a)<0) then
      if (psex(x,a)==3)  no_can_do = .TRUE.     ! unclear if maternal or paternal dummy
    else if (psex(x,a) /= sex(ij(a)) .and. any(parent == ij(a))) then
     ! cannot change sex if assigned as parent in pedigree
      no_can_do = .TRUE.    
!      call Rprint('no can do: sex ', (/x, a, sex(ij(a)), psex(x,a)/), (/0d0/), 'INT')
    endif
  enddo
  if (all(ij>0)) then  ! agedif specifying column in pairs file
    if (pairAgeDiff(x) /= AgeDiff(ij(1),ij(2)) .and. .not. (pairAgeDiff(x)==999 .or. AgeDiff(ij(1),ij(2))==999)) then
      do a=1,2
        if (any(parent(ij(a),:)/=0) .or. any(parent == ij(a)))  no_can_do = .TRUE.
        ! new agedif might make some pedigree links impossible; tricky to check. 
      enddo
 !      call Rprint('no can do: age ', (/x, a, AgeDiff(ij(1),ij(2)), pairAgeDiff(x)/), (/0d0/), 'INT')
    endif
  endif
  if (no_can_do) then
    LLpair(x,:) = NotImplemented  ! 444
    cycle
  endif
  
  ! backup sex and birth years from lifehistory data and pedigree
  ! set sex & age difference specified in pairs DF
  Sex_ij = 0
  BYtmp = -9
  IndBYtmp = 0d0
  do a=1,2
    if (ij(a)>0) then
      Sex_ij(a) = sex(ij(a))       ! backup
      Sex(ij(a)) = psex(x,a)     ! change sex to that specified in pairs file         
    endif
    if (all(ij>0)) then
      BYtmp(a) = BY(ij(a))   ! backup
      IndBYtmp(:,a,:) = IndBY(:,ij(a),:)
    endif
  enddo
  if (all(ij>0)) then      
    if (pairAgeDiff(x) == 999 .and. AgeDiff(ij(1),ij(2))/=999) then
      BY(ij) = -9   ! doesn't matter which one is set to missing? set both to missing to be sure
      IndBY(:,ij,:) = LOG10(1.0D0/nYears)  
    else if (pairAgeDiff(x) /= AgeDiff(ij(1),ij(2))) then
      BY(ij(1)) = pairAgeDiff(x) +1
      BY(ij(2)) = 1    ! AgeDiff = BY(i) - BY(j)
    endif
  endif
      
  ! temp. drop parents
  curPar = 0
  tmp_drop_mates=.FALSE.
  mates_tmp = 0  ! TODO: really necessary? it's just the co-parent isn't it?
  do a=1,2
    curPar(:,a) = getPar(ij(a), psex(x,a))
    if (Complx==0 .and. getNFS(ij(a))==1 .and. pdrop(x,a)>0)  tmp_drop_mates(a)=.TRUE.  ! no FS; do not condition on mate
    do m=1,2  
      if (tmp_drop_mates(a))  mates_tmp(m,a) = getMate(curPar(m,a),m)
      if (pdrop(x,a)==m .or. pdrop(x,a)==3) then 
        call setParTmp(ij(a), psex(x,a), 0, m)
        if (tmp_drop_mates(a))  call setMate(curPar(m,a),m, 0)
      endif
    enddo
  enddo 
  
  ! calc likelihoods
  if (all(ij>0)) then   ! .and. any(pairk(1:np)/=1) ?
    kij = pairk(x)
  else
    kij = psex(x, :)
  endif
  if (pairAgeDiff(x) < 0 .and. pairfocal(x)==7) then  ! swap, else LL_PO & LL_GP not calculated
    call CheckRel(ij(2), kij(2), ij(1), kij(1), pairfocal(x), LLg, LLa)
  else
  !  fclx=pairfocal(x)
 !   if (ij(1)<0 .and. pairfocal(x)==1)  fclx=4  ! GP of sibship. Reverse taken care of by CheckRel()
    call CheckRel(ij(1), kij(1), ij(2), kij(2), pairfocal(x), LLg, LLa)
 !   if (ij(1)<0 .and. pairfocal(x)==1)  LLpair(x,:) = reorderadd2(LLpair(x,:))  ! Done by CalcParentProbs() in R. 
  endif  

  LLdup = missing
  if (nrels==8) then   ! also calc LL to be duplicates
    if (all(ij>0)) then
      call PairSelf(ij(1), ij(2), LLdup) 
    else if (ij(1)>0) then
      call AddParent(ij(1), -ij(2), kij(2), LLdup)   ! does NOT consider age (?)
    else if (ij(2)>0) then
      call AddParent(ij(2), -ij(1), kij(1), LLdup)
    else if (kij(1)==kij(2)) then
      call MergeSibs(-ij(1), -ij(2), kij(1), LLdup)
    else
      LLdup = impossible
    endif
    LLpair(x,1) = LLdup
    LLpair(x,2:8) = LLg
  else
    LLpair(x,:) = LLg
  endif 
  
  ! restore sex, BY, IndBY, Parents    
  do a=1,2
    if (ij(a) > 0)  sex(ij(a)) = sex_ij(a) 
    do m=1,2
      call setParTmp(ij(a), psex(x,a), curPar(m,a), m)
      if (tmp_drop_mates(a))  call setMate(curPar(m,a),m, mates_tmp(m,a))
    enddo
  enddo
  if (all(ij>0)) then 
    BY(ij) = BYtmp(1:2)
    IndBY(:,ij,:) = IndBYtmp
  endif
enddo 

call MtoVd(LLpair, np, nrels, llrf)

if (nP > 50 .and. quiet < 1)  call rprint_eol() 
deallocate(IndBYtmp)
call DeAllocAll


contains
  ! function reorderadd2(LL)  result(LLtmp)
    ! ! reorder output from CheckAdd for compatibility with CheckPair
    ! double precision, intent(IN) :: LL(7)
    ! double precision :: LLtmp(7)

    ! LLtmp = missing
    ! LLtmp(1) = LL(4)
    ! LLtmp(2:3) = LL(5:6)
    ! LLtmp(5) = LL(2)
    ! LLtmp(6) =  MaxLL(LL((/1,3/)))   
    ! LLtmp(7) = LL(7) 
  ! end function reorderadd2
end subroutine getpairll

! ######################################################################

subroutine getbyprobs(ng, nx, ny, nYearsIn, byrf, lyrf, aprf, parentsrf, byprobv)
use Global
use reshapers
use sort_module, ONLY: getRank_i
implicit none

integer, intent(IN) :: ng, nx, ny, nYearsIn
integer, intent(IN) :: byrf(3*ng), lyrf(ng), parentsrf(2*ng)
double precision, intent(IN) :: aprf(5*ny)
double precision, intent(INOUT) :: byprobv(nx*nYearsIn)
integer :: i, j, k, BYrange(Ng,2), BYrankI(Ng), dumParV(2*Ng)
double precision :: AP_IN(ny, 5), BYLR(nYearsIN)
double precision, allocatable :: BYprobM(:,:)

allocate(BYprobM(nx, nYearsIn))

nInd = Ng
maxSibSize = 500   ! TODO? make user-setable?

call AllocArrays()

do i=1,nInd
  BY(i) = BYRF(i)
  BYrange(i,1) = BYRF(i + nInd)
  BYrange(i,2) = BYRF(i + 2*nInd)
enddo  
YearLast = LYRF

AP_IN = 0.0D0
k = 0
do i=1,5
  do j=1, ny
    k = k+1
    AP_IN(j,i) =  APRF(k)
  enddo
enddo

! grandparent & avuncular ageprior columns
call PrepAgeData(Ny, AP_IN, BYrange)
if (nYears /= nYearsIN)  call ErStop("nYears differ", .TRUE.)

! process pedigree
dumParV = 0
call ReadInputPed(ParentsRF, dumParV)
call rchkusr()

! estimate birth years; do oldest known-BY individuals first
call getRank_i(BYrankI)
do k=1,7   ! iterate a few times to ensure convergence
  do j=1, nInd
    i = BYRankI(j)
    call setEstBY(i, Sex(i))                              
  enddo
enddo

! retrieve estimates for those with unknown BY
j = 0
BYprobM = 0D0
BYLR = missing
do i=1, nInd
  if (BY(i) > 0)  cycle
  j = j+1
  if (j > nx)  call ErStop("number w/o BY does not match", .TRUE.)
  call getEstBY(i, 0, 5, BYLR)
  BYprobM(j,:) = 10**BYLR
enddo

! matrix to vector
call MtoVd(BYprobM, nx, nYearsIN, byprobv)

deallocate(BYprobM)
call DeAllocAll

end subroutine getBYprobs

! ######################################################################

subroutine parents(TotLL)
use sort_module
use Global
implicit none

double precision, intent(INOUT) :: TotLL(42)
integer :: i, j, k, Round, isP(2), BYRank(nInd), BYRank_Selfed(nInd)
double precision :: SortBY(nInd)
logical :: withPrior

AgePhase = 0                                       
 
call rchkusr()     
call UpdateAllProbs()

if (quiet<1) then
  call rprint_tbl_update_a(0,100) ! 'initial   ')
  call rprint_status_tbl_no_dots()
  call rprint_tbl_update_b()
endif

if (ANY(Parent /= 0)) then
  if(quiet<1)  call rprint_tbl_update_a(0,101) ! 'ped check ') 
  call CheckPedigree(.TRUE.)
  call UpdateAllProbs()
  if(quiet<1)  call rprint_tbl_update_b()
  withPrior = .TRUE.
else
  withPrior = .FALSE.
endif

if(quiet==0)  call rprint_tbl_update_a(0,102) ! 'parents   ')

!============================
! get birthyear ranking (increasing)
call getRank_i(BYrank)  ! do earlier birth years before later birth years

if (hermaphrodites/=0) then  ! do selfed individuals first to reduce false positives
  BYRank_Selfed = (/ (i, i=1, nInd, 1) /)
  do i=1, nInd
    call IsSelfed(i, .FALSE., SortBY(i))
    SortBY(i) = -SortBY(i)
  enddo
  call QsortC(SortBy, BYRank_Selfed)
endif

TotLL = 0D0
TotLL(1) = SUM(Lind)
do Round=1,41
  call rchkusr()
  if(quiet==-1)  call rprint_tbl_update_a(Round,102) ! 'parents   ')
  if (hermaphrodites/=0) then
    call Parentage(BYRank_Selfed, withPrior .and. Round==1)
  else
    call Parentage(BYrank, withPrior .and. Round==1)       
  endif
  call UpdateAllProbs()
  if(quiet==-1)  call rprint_tbl_update_b()  
 
  if (any (BY < 0))  call getRank_i(BYrank)
  
  do i=1,nInd
    if (Sex(i)==3) then
      isP = 0
      do k=1,2
        do j=1,nInd
          if (Parent(j,k) == i) then
            isP(k) = isP(k) + 1
          endif
        enddo
      enddo
      if (isP(1)>0 .and. isP(2)>0) then
        call rwarn("individual assigned as both dam & sire")
      else
        do k=1,2
          if (isP(k)>1) then
            Sex(i) = k
          endif
        enddo
      endif
    endif     
  enddo

  TotLL(Round + 1) = SUM(Lind)
  if (TotLL(Round + 1) - TotLL(Round) < ABS(TF)) exit ! convergence
  if (Round==41) then
    call Erstop("parentage not converging - need better SNP data", .FALSE.)
  endif
enddo

if(quiet==0) call rprint_status_tbl_no_dots()
if(quiet==0) call rprint_tbl_update_b()  

end subroutine parents

! ####################################################################

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! ####################################################################

subroutine sibships(AgeEffect, TotLL)
use Global
implicit none

integer, intent(INOUT) :: AgeEffect   ! IN
double precision, intent(INOUT) :: TotLL(42)
integer :: Round, RX, k, PairID(XP*nInd,2), PairType(XP*nInd), MaxRounds                                 
! character(len=10), allocatable :: stepname(:)   ! done in pretty_R_print.c

RX = 1  ! no. of initial rounds, pairs-cluster-merge only    
MaxRounds = 42   ! failsafe
! 0: no ageprior; 1: with + w/o; 2: extra age in last rounds
if (AgeEffect==0) then  !  .or. nYears==1
  AgePhase = 0  ! do not use ageprior
else
  AgePhase = 1  ! use age prior
endif     

if (Complx == 0 .and. any(Parent /=0))   call CheckMono()
call UpdateAllProbs()     

! allocate(stepname(12))
! stepname(1) = 'find pairs'
! stepname(2) = 'clustering'
! stepname(3) = 'GP pairs'
! stepname(4) = 'merging'
! stepname(5) = 'P of sibs'
! stepname(6) = 'GP Hsibs'
! stepname(7) = 'GP Fsibs'
! stepname(8) = 'find/check'
! stepname(9) = '(all)'
! stepname(10) = 'initial'
! stepname(11) = 'ped check'
! stepname(12) = 'end'

if (quiet<1) then
  call rprint_tbl_update_a(0,100)  ! initial
  call rprint_status_tbl_no_dots()
  call rprint_tbl_update_b()
endif

TotLL = 0D0
TotLL(1) = SUM(Lind)

if (any(Parent /=0)) then
  if(quiet==-1)  call rprint_tbl_update_a(0,101)  ! ped check 
  call CheckPedigree(.FALSE.)  ! double check parents, using updated ageprior    
  call UpdateAllProbs()
  if(quiet==-1)  call rprint_tbl_update_b()
endif     

do Round=1, MaxRounds                        
  call rchkusr()
  if(quiet==0)  call rprint_tbl_update_a(Round,300)  ! all
  
!  if(quiet==-1)  call Rprint("--- Round "//RoundChars(Round)//" start ---", (/0/), (/0.0D0/), "NON")
  
  ! == find pairs of potential siblings ==
  if(quiet==-1)  call rprint_tbl_update_a(Round,1)  ! find pairs
  call FindPairs(PairID, PairType)
  if(quiet==-1)  call rprint_tbl_update_b()
  call rchkusr()
  
  ! == cluster pairs into sibships ==
  if(quiet==-1)  call rprint_tbl_update_a(Round,2)
  call Clustering(PairID, PairType)  
  call UpdateAllProbs() 
  if(quiet==-1)  call rprint_tbl_update_b()  
  call rchkusr() 
  
  ! == grandparent - grandoffspring pairs ==
  if ((Round > RX+1) .or. (sum(nC)==0)) then
    if(quiet==-1)  call rprint_tbl_update_a(Round,3)
    call GGpairs()                  
    call UpdateAllProbs()
    if(quiet==-1)  call rprint_tbl_update_b() 
    call rchkusr()
  endif       
  
  if (sum(nC)==0)   RX=-1   ! no. of initial rounds, pairs-cluster-merge only   
  
  ! == merge sibship clusters ==
  if(quiet==-1)  call rprint_tbl_update_a(Round,4)
  call Merging()
  call UpdateAllProbs()  
  if(quiet==-1)  call rprint_tbl_update_b()   
  call rchkusr()
  
  ! == replace sibship dummy parents by genotyped individuals ==
  if(quiet==-1)  call rprint_tbl_update_a(Round,5)
  call SibParent()   
  call UpdateAllProbs()
  if(quiet==-1)  call rprint_tbl_update_b()   
  call rchkusr()
  
  ! == sibship grandparents ==
  if (Round > RX .or. Round==MaxRounds) then  
    if(quiet==-1)  call rprint_tbl_update_a(Round,6)
    call SibGrandparents()         
    call UpdateAllProbs() 
    if(quiet==-1)  call rprint_tbl_update_b()   
    call rchkusr()  
    
    if (Round > RX+1) then
      if(quiet==-1)  call rprint_tbl_update_a(Round,7)
      call FsibsGPs()
      call UpdateAllProbs()
      if(quiet==-1)  call rprint_tbl_update_b()   
      call rchkusr()      
    endif                    
    do k=1,2
      IsNewSibship(1:nC(k), k) = .FALSE.
    enddo 
  endif   
  
  ! == assign additional parents & double check new assignments ==
  if(quiet==-1)  call rprint_tbl_update_a(Round,8)
  call MoreParent() 
  call UpdateAllProbs()
  if(quiet==-1)  call rprint_tbl_update_b()
  call rchkusr()
  
  if(quiet==0)  call rprint_status_tbl_no_dots()
  if(quiet==0)  call rprint_tbl_update_b()
  if(quiet==-1)  call Rprint('', (/0/), (/0d0/), 'NON')  ! blank line
  
  TotLL(Round +1) = SUM(Lind)
  if (Round == MaxRounds) then
    exit
  else if (TotLL(Round +1) - TotLL(Round) < 2D0*ABS(TF) .and. Round > RX+1) then  
    if (AgeEffect==2 .and. AgePhase==1) then
      AgePhase = 2
    else             
      exit
    endif
  endif  
enddo  

end subroutine Sibships

subroutine rprint_tbl_update_a(round, step)
  implicit none
  integer, intent(IN) :: round, step
  integer :: date_time_values(8)

  call date_and_time(VALUES=date_time_values)
  call rprint_status_tbl_a(date_time_values(5:7), Round, step)  ! pretty_R_print.c
end subroutine rprint_tbl_update_a

subroutine rprint_tbl_update_b()
  use Global
  implicit none
  integer :: nparents(2), ngp
  double precision :: totlik
   
  nparents = count(Parent/=0, DIM=1)
  ngp = count(GpID/=0)
  totlik = SUM(Lind)
  call rprint_status_tbl_b(nparents, ngp, totlik)  ! defined in pretty_R_print.c   
end subroutine rprint_tbl_update_b

  
! #####################################################################

! @@@@   SUBROUTINES   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! #####################################################################

subroutine CalcP2(A,kA, B, C, kB, LLR)   ! LR of A; B & C both parent vs both U  
use Global
use OHfun
implicit none

integer, intent(IN) :: A, kA, B, C, kB
double precision, intent(OUT) :: LLR
integer :: l, x, y, z, PP(2), k
double precision :: PrL(nSnp,2), PrXYZ(3,3,3,2), PrP(3,2), PrA(3)
logical :: Selfed                 

if (A==0 .or. (B==0 .and. C==0)) then
  LLR = 0D0
  return
endif

if (kB==1) then
  PP = (/B,C/)
else
  PP = (/C,B/)
endif

LLR = missing
do k=1,2
  if (A>0 .and. PP(k)>0) then
    if (QLR_PO(A,PP(k)) < 5.0*TF) then
      LLR = impossible 
      return
    endif
  endif
enddo
  
if (A>0 .and. all(PP>0))  then 
  if (CalcTrioErr(A, PP) > MaxMendelE) then
    LLR = impossible
    return
  endif
endif

if (all(PP<0)) then   ! if many sibs shared, this approx is not valid
  if (any(Parent(SibID(1:ns(-PP(1),1),-PP(1),1), 2) == PP(2))) then
    LLR = 0d0
    return
  endif
endif

Selfed = .FALSE.
if (hermaphrodites/=0 .and. A>0) then
  if (SelfedIndiv(A))  Selfed = .TRUE.
endif

PrL = 0D0
do l=1,nSnp
  call ParProb(l, A, kA, -4, 0, PrA)
  do k=1,2
    call ParProb(l, PP(k), k, 0,0, PrP(:,k))
  enddo
  PrXYZ = 0D0
  do x=1,3
    do y=1,3
      if (Selfed) then
        PrXYZ(x,y,y,1) = PrA(x) * AKA2P(x,y,y) * PrP(y,1)
        PrXYZ(x,y,y,2) = PrA(x) * AKA2P(x,y,y) * AHWE(y,l)
      else
        do z=1,3
          PrXYZ(x,y,z,1) = PrA(x) * AKA2P(x,y,z) * PrP(y,1) * PrP(z,2)
          PrXYZ(x,y,z,2) = PrA(x) * AKA2P(x,y,z) * AHWE(y,l) * AHWE(z,l)
        enddo
      endif
    enddo
  enddo
  PrL(l,1) = LOG10(sum(PrXYZ(:,:,:,1)))
  PrL(l,2) = LOG10(sum(PrXYZ(:,:,:,2)))
enddo

if (SUM(PrL(:,1)) < -HUGE(0D0)) then
  LLR = impossible
else
  LLR = SUM(PrL(:,1)) - SUM(PrL(:,2)) 
endif

end subroutine CalcP2

! ######################################################################

subroutine ChkValidPar(A, kA, P, kP, OK)  ! age, ancestor, & OH check
use Global
implicit none

integer, intent(IN) :: A, kA, P, kP
logical, intent(OUT) :: OK
logical :: AncOK, MoreChk
double precision :: ALR, LRQ
integer :: ParA(2), i, mateP, mateParA, ParP(2)

if (A==0 .or. P==0) then
  OK = .TRUE.
  return
endif

if (A==P .and. (A>0 .or. kA==kP)) then
  OK = .FALSE.
  return
endif

ParA = getpar(A, kA)
if (kP < 3) then
  if (ParA(kP) == P) then
    OK = .TRUE.
    return
  endif
endif

if (Complx==0) then    ! TODO: test effect on assignment. 
  mateP = getMate(P,kP)
  if (mateP/=0 .and. mateP/=ParA(3-kP) .and. ParA(3-kP)/=0) then
    OK = .FALSE.
    return
  endif
  mateParA = getMate(ParA(3-kP),3-kP)
  if (mateParA/=0 .and. mateParA/=P) then
    OK = .FALSE.
    return
  endif
endif 

OK = .FALSE.
AncOK = .FALSE.
ALR = missing
LRQ = missing

if (A>0 .and. P>0) then
  if (AgeDiff(A,P) <=0)  return   ! P younger than A  (unknown coded as >0)
endif

call ChkAncest(P, kP, A, kA, AncOK)  ! check that A is not an ancestor of P
if (.not. AncOK)  return


if (Complx < 2) then   ! creating inbred configs invalid  
  ParP = getpar(P, kP)
  if (A>0) then
    if (ParP(3-kP)==ParA(3-kP) .and. ParA(3-kP)/=0)  return 
  else
    if (any(Parent(SibID(1:ns(-A,kA),-A,kA), kP) == P)) return
    if (ParP(3-kA)/=0) then
      if (any(Parent(SibID(1:ns(-A,kA),-A,kA), 3-kA) == ParP(3-kA)))  return
    endif
  endif
endif  


call CalcAgeLR(A,kA, P, kP, 0,1, .TRUE., ALR)
if (ALR == impossible)  return

if (kP < 3) then
  call CalcP2(A, kA, P, ParA(3-kP), kP, LRQ)  
else
  call CalcP2(A, kA, P, 0, kP, LRQ)  
endif
if (LRQ == impossible)  return

if (A>0 .and. P<0 .and. ALR < -TA) then   ! check age diff with future sibs
  do i=1, ns(-P,kP)
    call CalcAgeLR(A, kA, SibID(i,-P,kP), 3, kP, 3, .TRUE., ALR)
    if (ALR == impossible)  return
  enddo
endif

OK = .TRUE.

MoreChk = .FALSE.
if (any(OKA2P < TINY(0D0))) then  ! some combi of obs. offspring genotype + act parents is impossible
  if (kP < 3) then
    if (ParA(kP)==0)  MoreChk = .TRUE.
  endif
endif

if (MoreChk) then
  if (A > 0  .and. P < 0 .and. ParA(3-kP)<0) then
    Parent(A, kP) = P
    nS(-P,kP) = nS(-P,kP) +1
    SibID(nS(-P,kP), -P,kP) = A  
    call CalcCLL(-P,kP)
    if (CLL(-P,kP) < -HUGE(0D0))  OK = .FALSE.    ! inbreeding / FindEE loop
    SibID(nS(-P,kP), -P,kP) = 0  
    nS(-P,kP) = nS(-P,kP) -1
    Parent(A,kP) = 0
    call CalcCLL(-P,kP)

  else if (A < 0) then
    GpID(kP, -A, kA) = P
    call CalcCLL(-A,kA)
    if (CLL(-A,kA) < -HUGE(0D0))  OK = .FALSE.    ! e.g. inbreeding
    GpID(kP, -A, kA) = 0
    call CalcCLL(-A,kA)
  endif
endif

end subroutine ChkValidPar

! ######################################################################

subroutine ChkValidGP(A, kA, B, kB, m, OK, ALR)   ! age & ancestor check
use Global
implicit none

integer, intent(IN) :: A, kA, B, kB, m
logical, intent(OUT) :: OK
double precision, intent(OUT) :: ALR
logical :: AncOK
integer :: AncA(2,mxA), ParA(2), kBx(2), mateB, ParB(2)
                  
ALR = 0d0
if (A==0 .or. B==0) then
  OK = .TRUE.
  return
endif

OK = .FALSE.

if (A==B .and. (A>0 .or. kA==kB))   return

if (B<0 .and. .not. (kB==1 .or. kB==2))  call ErStop("ChkValidGP: kB must be 1 or 2 if B<0", .TRUE.)

if (A>0 .and. B>0) then
  if (AgeDiff(A, B) <= 0)  return   ! B younger than A  (unknown coded as >0)
endif

parA = getPar(A,kA)
if (parA(m)/=0) then
  call ChkAncest(B, kB, ParA(m),m, AncOK)  
else
  call ChkAncest(B, kB, A, kA, AncOK)  ! check that A is not an ancestor of B
endif   
if (.not. AncOK)  return

if (Complx==0) then  ! monogamous
  mateB = getMate(B,kB)
  if (mateB/=0 .and. ParA(3-kB)==mateB)  return
  !mateParA = getMate(ParA(3-m),3-m)   ! TODO?
endif

if (Complx < 2 .and. ParA(3-m)/=0) then   ! creating inbred configs invalid
  if (ParA(3-m)==B .and. 3-m==kB)  return 
  ParB = getPar(B,kB)
  if (ParA(3-m) == ParB(3-m))  return
endif  

call CalcAgeLR(A,kA, B, kB, m,4, .TRUE., ALR) 
if (ALR == impossible)  return

call GetAncest(A, kA, AncA)
if (kB==1 .or. kB==2) then
  kBx = (/kB,kB/)
else
  kBx = (/1,2/)   ! also covers case of hermaprodites
endif
if (any(AncA(kBx, 3:4) == B)) then  ! is already GP.
  OK = .TRUE.
  return
else if (all(AncA(kBx, 3:4) /= 0)) then
  OK = .FALSE.  ! both grandfathers/mothers already assigned, and neither is B
  return
endif

OK = .TRUE.

end subroutine ChkValidGP

! ######################################################################

subroutine ChkDoQuick(s,k, DoQuick)
use Global
implicit none

integer, intent(IN) :: s,k
integer, intent(OUT) :: DoQuick
integer :: nOff, Offspr(maxSibSize), sxOff(maxSibSize), i, OpPar, &
  UseEE(ns(s,k)), Sibs(ns(s,k)), MatePar(ns(s,k))

!  1: all parents 3-k >=0, no inbreeding
!  0: some parents 3-k <0, no inbreeding
! -1: inbreeding: Parent(Offspr(i),3-k) == GpID(3-k,s,k)
! -2: all are FS; ns(s,k) = ns(..,3-k)
! -3: s has a dummy clone
!  2: inbreeding: Parent(Bj,3-k) = Offspr(i)
!  3: Parent(Bi,3-k) close relative of Parent(Bj,3-k)

DoQuick = 1
if (ns(s,k)==0)  return   ! may happen temporarily
!if (Complx==1)  return                      

if (any(Parent(SibID(1:ns(s,k),s,k),3-k) < 0)) then 
  DoQuick = 0
else
  DoQuick = 1
endif

OpPar = 0
if (all(Parent(SibID(1:ns(s,k),s,k),3-k) < 0)) then
  call getFSpar(s, k, .TRUE., OpPar)
  if (OpPar < 0) then
    if (ns(-opPar, 3-k) == ns(s,k))  DoQuick = -2
  endif
endif

if (Hermaphrodites/=0) then
  if (DumClone(s,k)/=0)  DoQuick = -3   
endif

nOff = 0
Offspr = 0
sxOff = 3
call getOff(-s,k, .TRUE., nOff, Offspr, sxOff)  
do i=1, nOff
  if (Offspr(i) < 0 .and. sxOff(i)==3-k) then
    if (any(Parent(SibID(1:ns(s,k),s,k),3-k) == Offspr(i))) then
      DoQuick = 2
    endif
  else if (Offspr(i) > 0) then
    if (Parent(Offspr(i),3-k) == GpID(3-k,s,k) .and. GpID(3-k,s,k)/=0) then
      DoQuick = -1
      exit
    endif
  endif
enddo

UseEE = 0
if (DoQuick==0 .and. Complx/=1) then
  Sibs = SibID(1:ns(s,k), s, k)
  call FindEE(Sibs, ns(s,k), 0, k, UseEE, MatePar)
  if (any(UseEE/=0))  DoQuick = 3
endif

end subroutine ChkDoQuick

! ######################################################################

subroutine CalcPX2(A, kA, P1, P2, LLR)   ! joined LR of A+B+C; B & C as parent vs either or both U  
use Global
implicit none

integer, intent(IN) :: A, kA, P1, P2
double precision, intent(OUT) :: LLR
integer :: x, curPar(2)   ! May or may not be P1, P2
double precision :: LLY(2,2), LLU(4), LLcor(3,2)

curPar = getPar(A, kA)
do x=1,2
  call setParTmp(A, kA, 0, x)
enddo

LLY = missing
call Calc4U((/P1, P2/), 0,0, A,kA, LLU, LLcor)
LLY(1,1) = LLcor(3,1) + LLU(4)  ! no parents. LLU(4) = CLL(-A,kA)

call setParTmp(A,kA,P1,1)
if (Complx/=0) then
  call CalcU(A,kA, P1,1, LLY(2,1))
  LLY(2,1) = LLY(2,1) + LLcor(1,1)   ! only dam
endif

call setParTmp(A,kA,P2,2)        
call CalcU(A,kA, P1,1, LLY(2,2))
LLY(2,2) = LLY(2,2) + LLcor(1,1)   ! dam + sire
 
if (Complx/=0) then    
  call setParTmp(A,kA,0,1)
  call CalcU(A,kA, P2,2, LLY(1,2))
  LLY(1,2) = LLY(1,2) + LLcor(2,2)   ! only sire
endif

do x=1,2
  call setParTmp(A, kA, curPar(x), x)
enddo

LLR = LLY(2,2) - MaxLL((/LLY(1,:), LLY(2,1)/))

if (hermaphrodites/=0) then   ! A>0 & A<0 (?)
  if (P1<0 .and. P2<0) then
    if (DumClone(-P1,1) == -P2) then
      LLR = LLY(2,2) - LLY(1,1)
    endif
  else if (P1>0 .and. P2>0) then
    if (P1 == P2) then
      LLR = LLY(2,2) - LLY(1,1)
    endif
  endif
endif

end subroutine CalcPX2

! ######################################################################

subroutine CheckPedigree(ParOnly)
use Global
use sort_module, ONLY: getRank_i                                
implicit none

logical, intent(IN) :: ParOnly  ! T:only parents / F:also dummy parents
integer :: i, k, x, curPar(2), BYrank(nInd)
logical :: parOK(2), dropS
double precision :: LLRpair!, LL(7,2)

call getRank_i(BYrank)

do x=1, nInd
  if (MODULO(x,100)==0)  call rchkusr()
  if (quiet==-1 .and. any(chunk_printdot_i==x)) call print_dot()
  
  i = BYRank(x)
  curPar = Parent(i,:)

  do k=1,2
    call setParTmp(i,Sex(i), 0,k)  ! also drops i from sibship
  enddo
  
  parOK = .TRUE.
  do k=1,2
    if (curPar(k)>0) then
      if (Sex(curPar(k)) /= k .and. Sex(curPar(k)) < 3) then
        ParOK(k) = .FALSE.    ! male assigned as dam or female as sire.
      endif
    endif
    if (ParOK(k)) then
      call ChkValidPar(i, Sex(i), curPar(k), k, parOK(k))
    endif
  enddo
  
  ! check parent-pair
  if (Parent(i,1)/=0 .and. Parent(i,2)/=0 .and. all(ParOK)) then
    call CalcPX2(i, Sex(i), curPar(1), curPar(2), LLRpair)
    if (LLRpair < TA) ParOK = .FALSE.
  endif
  
  ! do k=1,2
    ! if (ParOK(k) .and. curPar(k)/=0) then
      ! LL = missing
      ! call CheckRel(i,Sex(i),curPar(k),k, 1, LL(:,1), LL(:,2))  
      ! if (LL(1,2) - MaxLL(LL(:,2)) < TF)  ParOK(k) = .FALSE.
    ! endif
  ! enddo 
    
  do k=1,2
    if (ParOK(k)) then
      call setPar(i,Sex(i), curPar(k),k)
    else if (curPar(k) < 0) then
      dropS = .FALSE.
      call CheckDropSibship(-curPar(k), k, DropS)
!      if (hermaphrodites/=0 .and. .not. DropS) then
!        call CheckSelfed(curPar(k),k)
!      endif
    endif
  enddo
    
  if (Parent(i,1)/=curPar(1) .or. Parent(i,2)/=curPar(2)) then
    if (Complx == 0)  call UpdateMate(i, Sex(i), curPar, ParOnly) 
    if (hermaphrodites /= 0 .and. .not. ParOnly) then
      call CheckSelfed(i,Sex(i))
      if (SelfedIndiv(i) .and. all(Parent(i,:)==0)) then
        do k=1,2
          call NewSibship(i, 0, k)
        enddo
        do k=1,2
          DumClone(-Parent(i,k), k) = -Parent(i,3-k)
        enddo
      endif
    endif
  endif
enddo

end subroutine CheckPedigree

! ######################################################################

subroutine Parentage(BYrank, withPrior)
use Global
use sort_module 
use OHfun, ONLY: QLR_PO           
implicit none

integer, intent(IN) :: BYrank(nInd)
logical, intent(IN) :: withPrior
integer, parameter :: mxxCP = 5*mxCP
integer :: i, j, x, y, k, CandPar(mxxCP, 2), nCP(2), curPar(2), SexTmp(2), &
  CP_rank(mxxCP), CP_tmp(mxxCP)
double precision :: ALR, LLR_CP(mxxCP)
logical :: AncOK

AncOK = .FALSE.
ALR = missing    
CP_tmp = 0    

do x=1, nInd
  if (MOD(x,200)==0) call rchkusr() 
  if (quiet==-1 .and. any(chunk_printdot_i==x)) call print_dot()
  i = BYRank(x)
  if (Parent(i,1)>0 .and. Parent(i,2)>0)  cycle
  curPar = Parent(i,:)
  nCP = 0
  CandPar = 0
  SexTmp = 3 
  do k=1,2
    if(Parent(i,k)/=0) then
      nCP(k) = nCP(k) + 1
      CandPar(nCP(k), k) = Parent(i,k)
    endif
  enddo  
  
  do y=1,nInd 
    j = BYRank(y)
    if (i==j) cycle
    if (sex(j)<3) then
      if (nCP(sex(j))==mxxCP)  cycle
    endif
    if (ANY(Parent(j,:)==i)) cycle
    if (ANY(CandPar==j) .and. (Sex(j)<3 .or. withPrior)) cycle  ! already included
    if (QLR_PO(i,j) < TF) cycle
    if (AgeDiff(i,j) <= 0)  cycle                                                           
    call ChkAncest(j,sex(j), i,sex(i), AncOK)  ! check that i is not an ancestor of j
    if (.not. AncOK)  cycle      
    call CalcAgeLR(i,sex(i), j,sex(j), 0,1, .TRUE., ALR) 
    if (ALR == impossible)  cycle
    
    do k=1,2
      if (Sex(j) < 3 .and. Sex(j)/= k) cycle
      if (nCP(k)==mxxCP) cycle
      if (ANY(CandPar(:,k) == j))  cycle    
      if (DoMtDif) then
        if (k==1 .and. Sex(j)>2 .and. mtDif(i,j))  cycle    
      endif  
      nCP(k) = nCP(k) + 1
      CandPar(nCP(k), k) = j
      if (Complx==0 .and. Mate(j)/=0) then
        if ((.not. any(CandPar == Mate(j))) .and. nCP(3-k)<mxxCP .and. Mate(j)>0) then
          nCP(3-k) = nCP(3-k) +1
          CandPar(nCP(3-k), 3-k) = Mate(j)
        endif
      endif
    enddo
  enddo

  if (ALL(nCP <=1) .and. ALL(candPar(1,:) == Parent(i,:)))  cycle
  ! highly unlikely different sire will be assigned in absence of any cand. dam and v.v.
  if (nCP(1)==0 .and. curPar(2)/=0)  cycle  
  if (nCP(2)==0 .and. curPar(1)/=0)  cycle                                
  
  if (any(nCP > mxCP)) then ! sort by QLR_PO
    do k=1,2
      if (nCP(k) <= mxCP)  cycle
      LLR_CP = missing
      do y = 1, nCP(k)
        LLR_CP(y) = QLR_PO(CandPar(y,k), i)
      enddo
      CP_rank = (/ (y, y=1, mxxCP, 1) /)
      LLR_CP = -LLR_CP   ! to sort decreasing
      call QsortC(LLR_CP(1:nCP(k)), CP_rank(1:nCP(k)))
      CP_tmp(1:nCP(k)) = CandPar(1:nCP(k),k)
      CandPar(1:nCP(k),k) = CP_tmp(CP_rank(1:nCP(k)))
      nCP(k) = mxCP
    enddo
  endif

  call SelectParent(i, Sex(i), nCP, CandPar(1:mxCP,:), .TRUE.)  ! does actual assignment
  
enddo

end subroutine Parentage

! ######################################################################

subroutine CheckParentPair(A, kA, Par, dLL)
use Global
implicit none

integer, intent(IN) :: A, kA, Par(2)
double precision, intent(OUT) :: dLL(3)  ! 1:dam, 2:sire, 3: both (vs none or either)
integer :: NowPar(2), m
double precision :: LLRP(2), gLL(4,2)

NowPar = getPar(A, kA)
do m=1,2
  call setParTmp(A, kA, 0, m)
enddo

dLL = missing
LLRP = missing
gLL = missing

call CalcP2(A, kA, Par(1), Par(2), 1, LLRP(1))

if (LLRP(1) > 2*TF .and. LLRP(1) /= impossible)  then
  call CalcPX2(A, Sex(A), Par(1), Par(2), LLRP(2))
  
  if (LLRP(2) > TA)  then
    call setParTmp(A, kA, Par(1), 1)
    call CalcPOGPZ(A, kA, Par(2), 2, gLL)
    dLL = gLL(1:3,1)  ! no ageeffect
  endif
endif

do m=1,2
  call setParTmp(A, kA, NowPar(m), m)
enddo

end subroutine CheckParentPair

! ######################################################################

subroutine SelectParent(A, kAIN, nCP, CandPar, ParOnly)   
! assigns parent / grandparent as side effect
use Global
implicit none

integer, intent(IN) :: A, kAIN, nCP(2), CandPar(mxCP, 2)
logical, intent(IN) :: ParOnly
integer :: m, u, v, best(2), par, AG(2), kA, curpar(2), uv(2), mate_v,mate_u, &
   i,j, n, fcl, nSingle, ToDoSingleCheck(MAXVAL(nCP),2), swapMate(MAXVAL(nCP),2)
double precision :: LRS, LLRX(nCP(1),nCP(2)), LLRY(nCP(1),nCP(2)), ALR, & 
  LLRZpair(nCP(1),nCP(2),2), LLRZsingle(2,MAXVAL(nCP),2), gLL(4,2), TAx(2), dLLrev(2), &
   LLA(7,7,2,2), LRStmp(2,2)
logical :: AgeAUnk, SexUnk(MAXVAL(nCP), 2), MonoPair(nCP(1),nCP(2)), emptySibship(MAXVAL(nCP), 2), &
 PairCand(nCP(1),nCP(2)), AgeUnk(2), maybeRev(2), DoSingle, SingleCand(MAXVAL(nCP),2), &
  DoneCombo(SUM(nCP), SUM(nCP)), DoubtAssign
  
if (ALL(nCP==0))  return

DoubtAssign = .FALSE.   ! TODO: implement 

if (A > 0) then
  fcl = 1
else !if (A < 0) then
  fcl = 4
endif
if (kAIN >2 .or. kAIN==0) then    ! (only?) necessary for checkMaybeRev
  if (Sex(A)<3) then
    kA = Sex(A)
  else
    kA = 1
  endif
else
  kA = kAIN
endif

LRS = missing
if (hermaphrodites/=0 .and. A>0) then
  call IsSelfed(A, .FALSE., LRS)
endif

curpar = getPar(A,kA)
do m=1,2
  call setParTmp(A, kA, 0, m)
  call SetEstBY(curPar(m), m)  
enddo
call SetEstBY(A, kA)   ! conditional on no parents

AgeAUnk = .FALSE.                  
if (A < 0) then
  AgeAUnk = .TRUE.
else if (BY(A) < 0) then
  AgeAUnk = .TRUE.
endif

if (AgePhase == 0) then     ! only [genetics] must pass TA
  TAx = (/TA, LOG(zero)/) 
  AG = (/1,1/)
else if (AgePhase==1) then  ! [genetics] & [genetics + age] must pass TA
  AG = (/1,2/)       ! changed 2024-09-23 from (/2/) to (/1,2/)
  TAx = (/TA, TA/)  
else                        ! only [genetics + age] must pass TA,
  AG = (/2,2/)             
  TAx = (/-TA, TA/)        ! added 2024-12-24: [genetics] must pass -TA 
endif

if (A < 0) then   ! new 2025-08-09: also 2*TAx for grandparent-pairs of singletons
  if (ns(-A, kA) == 1) TAx = 2*TAx   
endif
 

SexUnk = .FALSE.
MonoPair = .FALSE.                  
EmptySibship = .FALSE.   
par = 0           
do m=1,2
  do u = 1, nCP(m)
    if (CandPar(u,m) > 0) then
      if (Sex(CandPar(u,m)) > 2) then
        if (A>0 .and. .not. any(CandPar(:,3-m) == CandPar(u,m))) then  ! pedigree prior
          SexUnk(u,m) = .FALSE.  
        else
          SexUnk(u,m) = .TRUE.
        endif
        if (DoMtDif) then
          if (A>0) then
            if (m==2 .and. mtDif(A, CandPar(u,m)))  SexUnk(u,m) = .FALSE.   ! cannot be maternal relative
          else if (A < 0) then
            if (m==2 .and. mtDif(SibID(1,-A,kA), CandPar(u,m)))  SexUnk(u,m) = .FALSE. 
          endif
        endif
      endif
    else if (CandPar(u,m) < 0) then
      v = -CandPar(u,m)
      if (ns(v,m) == 0 .and. all(GpID(:,v,m)==0)) then
        EmptySibship(u,m) = .TRUE.
      else if (ALL(Parent(SibID(1:nS(v,m), v, m), 3-m) < 0)) then
        call getFSpar(v, m, .TRUE.,par)
        if (par < 0) then
          if (nS(-par, 3-m) == nS(v,m)) then  ! cannot tell if mat or pat
            SexUnk(u,m) = .TRUE.   ! TODO?
            if (DoMtDif .and. A>0) then
              if (m==2 .and. mtDif(A, SibID(1,v,m)))  SexUnk(u,m) = .FALSE. 
            endif                 
            if (ANY(CandPar(:,3-m) == par)) then
              do v=1, nCP(3-m)
                if (CandPar(v,3-m) == par) then
                  if (m==1)  MonoPair(u,v) = .TRUE.  ! monogamous parent pair
                  if (m==2)  MonoPair(v,u) = .TRUE.
                endif
              enddo
            endif
          endif
        endif          
      endif
    endif
  enddo
enddo

PairCand = .FALSE.   ! check if candidates form an eligible parent-pair
mate_u = 0                                                                      
if (ALL(nCP>0)) then   
  do u=1, nCP(1)
    if (complx==0)  mate_u = getMate(CandPar(u,1),1)
    do v=1, nCP(2)
      PairCand(u,v) = .TRUE.
      
      if (Complx==0) then       
        mate_v = getMate(CandPar(v,2),2)       
        if (mate_u == candpar(v,2))  MonoPair(u,v) = .TRUE.                      
        if (mate_u/=0 .and. mate_u/=candpar(v,2))  PairCand(u,v) = .FALSE.  
        if (mate_v/=0 .and. mate_v/=candpar(u,1))  PairCand(u,v) = .FALSE.    
      endif
      
      if (SexUnk(u,1) .and. SexUnk(v,2)) then
        if (hermaphrodites==0 .and. .not. MonoPair(u,v))  PairCand(u,v) = .FALSE.
        if (hermaphrodites==1) then  ! only if selfed
          if (CandPar(u,1) > 0) then
            if (CandPar(u,1)/=CandPar(v,2))  PairCand(u,v) = .FALSE.
          else
            if (DumClone(-CandPar(u,1),1) /= -CandPar(v,2))  PairCand(u,v) = .FALSE.
          endif
        endif
      endif
      
      if (hermaphrodites/=0 .and. A>0) then
        if (LRS > TA .and. CandPar(u,1)/=CandPar(v,2) .and. (CandPar(u,1)>0 .or. CandPar(v,2)>0)) then
          PairCand(u,v) = .FALSE.
        endif
        if (LRS > TA .and. CandPar(u,1)>0) then
          if (Sex(CandPar(u,1)) /= 4) PairCand(u,v) = .FALSE.
        endif
        if (LRS < 5*TF .and. CandPar(u,1)==CandPar(v,2) .and. CandPar(u,1)>0) then
          PairCand(u,v) = .FALSE.
        endif
      endif

      if (hermaphrodites/=0 .and. CandPar(u,1)>0 .and. CandPar(v,2)>0) then    
        if (Hermaphrodites == 2 .and. any(CandPar(1:(u-1),1) ==CandPar(v,2)) .and. &
          any(CandPar(:,2) == CandPar(u,1)))   PairCand(u,v) = .FALSE.
          ! don't care if dam or sire; drop if pair already considered.
        if (Hermaphrodites == 1 .and. any(CandPar(:,1) == CandPar(v,2)) .and. &
          any(CandPar(:,2) == CandPar(u,1)) .and. &   ! can't distinguish dam-sire vs sire-dam
          .not. candPar(u,1) == CandPar(v,2))   PairCand(u,v) = .FALSE.   ! exception: selfing 
      endif
     enddo
  enddo
endif

swapMate = 0
if (Hermaphrodites==2) then
! monogamous dummy partners are both included with multiple partners: 
! keep only one from each set of defacto identical pairs
  do u=1, nCP(1)
    do v=1, nCP(2)
      if ((candpar(u,1)>0 .and. candPar(u,1) == candPar(v,2)) .or. &
        (candpar(u,1)<0 .and. candpar(v,2)<0 .and. monopair(u,v))) then
          swapMate(u,1) = v
          swapMate(v,2) = u
      endif
    enddo
  enddo
  
  do u=1, nCP(1)
    if (swapMate(u,1)==0)  cycle
    j = swapMate(u,1)
    do v=1,nCP(2)
      if (swapMate(v,2)==0)  cycle
      if (.not. PairCand(u,v))  cycle  
      i = swapMate(v,2)
      if (u==i .and. v==j)  cycle
      PairCand(i,j) = .FALSE.
    enddo
  enddo
endif
 

LLRX = missing  ! LR(P/U), w/o parents, ignoring inbreeding etc
LLRY = missing  ! LR(P/U), w parents + their relatedness
if (ANY(PairCand)) then   ! find plausible parent-pairs
  do u=1, nCP(1)
    do v=1, nCP(2)
      if (.not. PairCand(u,v))  cycle
      if (candPar(u,1) < 0) then
        if (ns(-CandPar(u,1),1) <= 1)  cycle   ! CalcPX2 unreliable (?)
      endif
      if (candPar(v,2) < 0) then
        if (ns(-CandPar(v,2),2) <= 1)  cycle  
      endif
      ! calc LLR parent pair / both unrelated --> LLRX(u,v)
      call CalcP2(A, kA, candPar(u,1), CandPar(v,2), 1, LLRX(u,v))
      if (LLRX(u,v)==Impossible .or. LLRX(u,v) < 2*TF) then
        PairCand(u,v) = .FALSE.
      else       
        call CalcPX2(A, kA, candPar(u,1), candPar(v,2), LLRY(u,v))
        if (A <0 .and. candPar(u,1) < 0 .and. candPar(v,2) < 0) then
          if (LLRY(u,v) < 2*TF)  PairCand(u,v) = .FALSE.  
        else if (candPar(u,1) < 0 .and. candPar(v,2) < 0) then
          if (LLRY(u,v) < -TA)  PairCand(u,v) = .FALSE.   
        else if (LLRY(u,v) < TA) then      
          PairCand(u,v) = .FALSE.
        endif 
      endif
    enddo
  enddo
endif
     
! age compatibility check
if (AgeAunk .and. ANY(PairCand)) then 
  ALR = missing
  do u=1, nCP(1)
    call setParTmp(A, kA, CandPar(u,1), 1)
    call SetEstBY(A, kA)
    do v=1, nCP(2)
      if (.not. PairCand(u,v))  cycle
      call CalcAgeLR(A, kA, CandPar(v,2), 2, 0,1, .TRUE., ALR)
      if (ALR == impossible)  PairCand(u,v) = .FALSE.
    enddo
  enddo
  call setParTmp(A, kA, 0, 1)
  call SetEstBY(A, kA)
endif


LLRZpair = missing  ! LR(P/next-most-likely)
LLRZsingle = missing
gLL = missing
DoneCombo = .FALSE.                   
if (ANY(PairCand)) then  ! possibly parent-pair
  do u=1, nCP(1)
    do v=1, nCP(2)
      if (.not. PairCand(u,v))  cycle
      if (Complx==0 .and. A<0 .and. MonoPair(u,v)) then
        mate_v = getMate(CandPar(v,2), 2)
        if (mate_v==0)  call setParTmp(A,kA, CandPar(u,1),1)
        if (mate_v/=0)  call setParTmp(A,kA, 0,1)  ! addGP will consider CandPar(v,2) + its mate
        call CalcPOGPZ(A, kA, CandPar(v,2), 2, gLL) 
        if (mate_v==0)  LLRZpair(u,v,:) = gLL(3,:)
        if (mate_v/=0)  LLRZpair(u,v,:) = gLL(2,:)
        cycle
      endif
      if (Complx==0 .and. A>0) then
        if (EmptySibship(v,2)) then
          call CalcPOGPZ(A, kA, CandPar(u,1), 1, gLL)
        else if (EmptySibship(u,1) .or. (CandPar(u,1) < 0 .and. CandPar(v,2) < 0)) then
          call setParTmp(A, kA, 0, 1)
          call CalcPOGPZ(A, kA, CandPar(v,2), 2, gLL)
        else
          call setParTmp(A, kA, CandPar(u,1), 1)
          call CalcPOGPZ(A, kA, CandPar(v,2), 2, gLL)
          DoneCombo(u, nCP(1)+v) = .TRUE. 
!          DoneCombo(nCP(1)+v, u) = .TRUE.          
        endif  
      else
        call setParTmp(A, kA, CandPar(u,1), 1)
        call CalcPOGPZ(A, kA, CandPar(v,2), 2, gLL)
        DoneCombo(u, nCP(1)+v) = .TRUE.   
!        DoneCombo(nCP(1)+v, u) = .TRUE.        
      endif
      LLRZpair(u,v,:) = gLL(3,:)  ! vs one or no parents
      if (any(gLL(1,AG) < LLRZsingle(AG,u,1))) then
        LLRZsingle(:,u,1) = gLL(1,:)   
      endif
      if (any(gLL(2,AG) < LLRZsingle(AG,v,2))) then
        LLRZsingle(:,v,2) = gLL(2,:)
      endif
    enddo
  enddo
endif
call setParTmp(A, kA, 0, 1) 

do u=1, nCP(1)
  do v=1, nCP(2)
    if (.not. PairCand(u,v))  cycle
    if (any(LLRZpair(u,v,AG) >= MaybeOtherParent)) then
      PairCand(u,v) = .FALSE.     
    else if (DoubtAssign .and. ABS(LLRZpair(u,v,1))<0.01 .and. &
     (COUNT(PairCand(u,:))>1 .or. COUNT(PairCand(:,v))>1)) then
      PairCand(u,v) = .TRUE.   ! not-yet-merged sibship / duplicated sample involved
      ! TODO: make more efficient. 
    else if (any(LLRZpair(u,v,:) < TAx(:))) then
      PairCand(u,v) = .FALSE.  
    endif
  enddo
enddo


if (ANY(PairCand)) then    ! check if parent-offspring not reversed
  do u=1, nCP(1)
    do v=1, nCP(2)
      if (.not. PairCand(u,v))  cycle
      !if (MonoPair(u,v) .and. hermaphrodites==2 .and. count(PairCand==1))  cycle   
 !     if (ParOnly)  cycle    ! reverse check by CheckPair()                
      AgeUnk = .FALSE.
      maybeRev = .FALSE.
      dLLrev = missing
      uv = (/u, v/)
      
      do m=1,2
        if (CandPar(uv(m),m) < 0) then
          AgeUnk(m) = .TRUE.
        else if (BY(CandPar(uv(m),m)) < 0) then
          AgeUnk(m) = .TRUE.
        endif
      enddo

      if (AgeAUnk .or. any(AgeUnk)) then                                  
        do m=1,2
          if (CandPar(uv(m),m) < 0) then             
            if (ns(-CandPar(uv(m),m),m) == 0)  cycle   ! A is/was only member of sibship
          endif
          call setParTmp(A, kA, CandPar(uv(3-m), 3-m), 3-m)
          call SetEstBY(A, kA)
          call CheckMaybeRev(A, kA, CandPar(uv(m),m), m, maybeRev(m), dLLrev)
          call setParTmp(A, kA, 0, 3-m)
          call SetEstBY(A, kA)  
          if (maybeRev(m) .and. (all(dLLrev==missing) .or. &
           any(LLRZpair(u,v,:) - dLLrev(:) < 2*TAx(:)))) then
            PairCand(u,v) = .FALSE.
            exit
          endif
        enddo
      endif
    enddo
  enddo
endif

best = 0
if (ANY(PairCand)) then
  best = MAXLOC(LLRZpair(:,:,AG(2)), MASK=PairCand)      
  do m=1,2
    call setParTmp(A, kA, CandPar(best(m), m), m)
  enddo 
endif  

  ! ~~~~  check when >1 eligible pair  ~~~~
if (COUNT(PairCand) > 1) then              
  do v=1, nCP(2)
    do u=1, nCP(1)   
      !if (DoneCombo(best(1), nCP(1)+v) .and. DoneCombo(u, nCP(1)+best(2)))  cycle
    
      mate_u = -999999
      if ((MonoPair(u,v) .or. Complx==0) .and. PairCand(u,v) .and. u/=best(1)) then   ! MonoPair(best(1), best(2)) .and. 
        if (Complx==0) then  ! TODO: .and. .not. DoneCombo(...)
          ! assign temporary mate if currently none
          mate_u = getMate(CandPar(u,1), 1)
          if (mate_u==0)  call setMate(CandPar(u,1), 1, CandPar(v,2))
        endif
        call CalcPOGPZ(A, kA, CandPar(u,1), 1, gLL)
        if (Complx==0 .and. mate_u==0)  call setMate(CandPar(u,1), 1, 0)   ! restore   
        if (any(gLL(4,AG) < LLRZpair(best(1),best(2),AG))) then
          LLRZpair(best(1),best(2),:) = gLL(4,:)
        endif  
        if (any(gLL(4,:) < TAx(:))) then  ! best(1) + (2) not most likely
          call setParTmp(A, kA, CandPar(u,1), 1)
          call setParTmp(A, kA, CandPar(v,2), 2)
          call CalcPOGPZ(A, kA, CandPar(best(1),1), 1, gLL)
          DoneCombo(u, nCP(1)+v) = .TRUE.  ! used for singles check only
          DoneCombo(u, best(1)) = .TRUE.
          DoneCombo(best(1), u) = .TRUE.
!          DoneCombo(nCP(1)+v, best(1)) = .TRUE.     
          if (any(gLL(4,AG) < LLRZpair(u,v,AG))) then
            LLRZpair(u,v,:) = gLL(4,:)
          endif              
          if (all(gLL(4,:) > TAx(:) .and. gLL(4,:)<555d0)) then      ! TODO: why 555?     
            best(1) = u
            best(2) = v
          else  ! restore
            call setParTmp(A, kA, best(1), 1)          
            call setParTmp(A, kA, best(2), 2) 
          endif
        endif
        cycle
      endif

      if (PairCand(best(1), v) .and. v/=best(2) .and. &
        (.not. (DoneCombo(best(1), nCP(1)+v) .and. DoneCombo(nCP(1)+best(2), nCP(1)+v)))) then
        call CalcPOGPZ(A, kA, CandPar(v,2), 2, gLL)
        DoneCombo(best(1), nCP(1)+v) = .TRUE.
        DoneCombo(nCP(1)+best(2), nCP(1)+v) = .TRUE.
        DoneCombo(nCP(1)+v, nCP(1)+best(2)) = .TRUE.            
        if (any(gLL(4,AG) < LLRZpair(best(1),best(2),AG))) then
          LLRZpair(best(1),best(2),:) = gLL(4,:)
        endif
        if (any(gLL(3,AG) < LLRZpair(best(1),v,AG))) then
          LLRZpair(best(1),v,:) = gLL(3,:)
        endif
        if (gLL(3,AG(1)) > gLL(4,AG(1))) then   ! best(1) + CandPar(v,2) more likely than best(1) + (2)
          call setParTmp(A, kA, CandPar(v,2), 2)
          best(2) = v           
        endif  
      endif
      
      if (PairCand(u, best(2)) .and. u/=best(1) .and. &
        (.not. (DoneCombo(u, nCP(1)+best(2))) .and. DoneCombo(u, best(1))))  then
        call CalcPOGPZ(A, kA, CandPar(u,1), 1, gLL) 
        DoneCombo(u, nCP(1)+best(2)) = .TRUE. 
        DoneCombo(u, best(1)) = .TRUE. 
        DoneCombo(best(1),u) = .TRUE.        
        if (any(gLL(4,AG) < LLRZpair(best(1),best(2),AG))) then
          LLRZpair(best(1),best(2),:) = gLL(4,:)
        endif
        if (any(gLL(3,AG) < LLRZpair(u, best(2),AG))) then
          LLRZpair(u, best(2),:) = gLL(3,:)
        endif
        if (gLL(3,AG(2)) > gLL(4,AG(2))) then
          call setParTmp(A, kA, CandPar(u,1), 1)
          best(1) = u
        endif  
      endif
      
    enddo
  enddo

  best = MAXLOC(LLRZpair(:,:,AG(2)), MASK=PairCand)
  do m=1,2
    call setParTmp(A, kA, CandPar(best(m), m), m)
  enddo
  ! if multiple parent pairs tested, tiny difference likely due to not-yet-merged sibships -- keep assignment
endif     
if (all(best>0) .and. MAXVAL(LLRZpair(:,:,AG(2)), MASK=PairCand) < TA) then
  do m=1,2
    call setParTmp(A, kA, 0, m)
  enddo  
endif

curPar = getPar(A,kA)
if (ALL(curPar/=0)) then  ! make it official
  do m=1,2
    call setPar(A,kA, curPar(m),m)
  enddo
  return
else
  do m=1,2
    call setParTmp(A,kA, 0,m)
  enddo
endif


DoSingle = .TRUE.
if (hermaphrodites/=0 .and. LRS > TA) then  ! selfed
  DoSingle = .FALSE.
else if (Complx==0 .and. A<0) then
  DoSingle = .FALSE.
  do m=1,2
    do u = 1, nCP(m)
      if (CandPar(u,m) > 0) then  
        if (Mate(CandPar(u,m)) == 0) then   ! exception: no genotyped offspring --> no mate yet
          DoSingle = .TRUE.
          exit
        endif
      endif  ! CandPar(u,m) < 0 w/o DumMate shouldn't exist
    enddo
    if (DoSingle)  exit
  enddo
endif

if (.not. DoSingle) then
  do m=1,2
    call setPar(A,kA, 0,m)   
  enddo  
  return
endif


! ~~~~  single parent  ~~~~
SingleCand = .FALSE.
do m=1,2
  do u=1, nCP(m)
    if (SexUnk(u,m) .and. hermaphrodites/=2)  cycle
    if (hermaphrodites==2 .and. m==2) then
      if (CandPar(u,m)>0) then
        if (Sex(CandPar(u,m))==4)  cycle  ! assign as dam
      else if (any(CandPar(1:nCP(1),1) == -DumClone(-CandPar(u,2),2))) then
        cycle ! assign as dam
      endif
    endif
    if (ANY(LLRZsingle(:,u,m) < TAx(:)))  cycle   

    if (Complx==0) then
      if (CandPar(u,m) < 0) then
        if (DumMate(-candpar(u,m),m)/=0)  cycle
      else
        if (Mate(CandPar(u,m))/=0)  cycle
      endif
    endif

    SingleCand(u,m) = .TRUE.
  enddo
enddo     
  
nSingle = COUNT(SingleCand)  ! SingleCand edited inside loop

if (nSingle >= 1) then
  ToDoSingleCheck = 0
  do m=1,2
    do u=1,nCP(m)
      if (SingleCand(u,m))  ToDoSingleCheck(u,m) = MAX(3, nSingle-1)  ! check each valid candidate with at least 3 others
    enddo
  enddo

  if (SUM(ToDoSingleCheck)>0 .and. (SUM(nCP)==1 .or. Complx==0)) then                                                                   
    do m=1,2
      do u=1,nCP(m)
        if (ToDoSingleCheck(u,m)==0)  cycle  
        LRStmp = Missing      
        call CalcCandParLL(A, kA, CandPar(u,m), m, LLA)
        call LLA2LLR(LLA, fcl, LRStmp)                             
        if (any(LRStmp(AG,1) < LLRZsingle(AG,u,m)))  LLRZsingle(:,u,m) = LRStmp(:,1)
        if (any(LLRZsingle(:,u,m) < TAx(:))) then
          SingleCand(u,m) = .FALSE.   
          ToDoSingleCheck(u,m) = 0
        else 
          ToDoSingleCheck(u,m) = ToDoSingleCheck(u,m) -1
        endif
      enddo
    enddo
  endif
   if (SUM(ToDoSingleCheck)> 0 .and. sum(nCP)>1 .and. Complx > 0) then  ! for mono, do not also condition on other CP being parent 
    do m=1,2
      do u=1,nCP(m)
        if (all(ToDoSingleCheck == 0))  exit                                       
        if (ToDoSingleCheck(u,m)==0)  cycle
        i = (m-1)*nCP(1) +u
        LRStmp = Missing 
        
        do n=1,2 
          do v=1,nCP(n)                          
            if (n==m .and. v==u)  cycle
            if (CandPar(u,m)>0 .and. CandPar(u,m) == CandPar(v,n))  cycle  ! when sex=4 
            j = (n-1)*nCP(1) + v
            if (DoneCombo(i,j) .or. DoneCombo(j,i))  cycle  
            if (ToDoSingleCheck(v,n)==0 .and. COUNT(ToDoSingleCheck>0)>1) cycle  ! at least 1 other to pair i with   
         
            call setParTmp(A, kA, CandPar(v,n), n) 
            call CalcCandParLL(A, kA, CandPar(u,m), m, LLA)
            call LLA2LLR(LLA, fcl, LRStmp) 
            DoneCombo(i,j) = .TRUE.
  
            if (any(LRStmp(AG,1) < LLRZsingle(AG,u,m)))  LLRZsingle(:,u,m) = LRStmp(:,1)
            if (any(LRStmp(AG,2) < LLRZsingle(AG,v,n)))  LLRZsingle(:,v,n) = LRStmp(:,2)
            if (ANY(LLRZsingle(:,u,m) < TAx(:))) then   ! changed 2024-09-23 from all to any
              SingleCand(u,m) = .FALSE. 
              ToDoSingleCheck(u,m) = 0
            else 
              ToDoSingleCheck(u,m) = ToDoSingleCheck(u,m) -1
            endif
            if (ANY(LLRZsingle(:,v,n) < TAx(:))) then
              SingleCand(v,n) = .FALSE.
              ToDoSingleCheck(v,n) = 0
            else 
              ToDoSingleCheck(v,n) = ToDoSingleCheck(v,n) -1
            endif
            if (all(ToDoSingleCheck == 0))  exit                                       
          enddo  ! v       
          call setParTmp(A, kA, 0, n)
        enddo  ! n
      
      enddo  ! u
    enddo  ! m
  endif
endif

do m=1,2
  do u=1,nCP(m)
    if (ANY(LLRZsingle(:,u,m) < TAx(:)) .or. ANY(LLRZsingle(:,u,m) >= 222D0)) then
      SingleCand(u,m) = .FALSE.
    else if (A<0 .and. CandPar(u,m)<0) then
      if ((ns(-A,kA) == 1 .or. ns(-CandPar(u,m),m) == 1) .and. ANY(LLRZsingle(:,u,m) < 2*TA)) then  
        SingleCand(u,m) = .FALSE.  ! GGpairs sensitive to false pos. added 2024-12-24; Maybe too conservative?
      endif
    endif
  enddo
enddo
           
! check if PO pair could be flipped, or relies strongly on (uncertain) age estimate
dLLrev = missing 
do m=1,2
  do u=1,nCP(m)
    if (.not. SingleCand(u,m))  cycle
    call CheckMaybeRev(A, kA,CandPar(u,m), m, maybeRev(1), dLLrev)  ! includes age check  
    if (maybeRev(1)) then
      if (all(dLLrev==missing) .or. any(LLRZsingle(AG,u,m) - dLLrev(AG) < TAx(AG)) .or. &
        ALL(dLLrev > TAx(AG))) then  
        SingleCand(u,m) = .FALSE.
      else   ! sometimes LLRZsingle for B is from LL(CP|B), which is not directly comparable
        call CalcPOGPZ(A, kA, CandPar(u,m), m, gLL)
        if (any(gLL(1,AG) - dLLrev(AG) < TAx(AG)))  SingleCand(u,m) = .FALSE.
      endif
    endif
  enddo
enddo 

best=0
best = MAXLOC(LLRZsingle(AG(2),:,:), MASK=SingleCand)  ! u,m

if (COUNT(SingleCand) == 1) then
  call setPar(A, kA, CandPar(best(1), best(2)), best(2))
  if (Complx == 0 .and. A>0 .and. .not. ParOnly) then  ! create dummy mate if assigned parent has no mate yet
    call NewSibship(A, 0, 3-best(2))
  endif
else
  do m=1,2
    call setParTmp(A,kA, 0,m)   ! do no funny bussiness w temp assigned ex-cand-parents
  enddo
  do m=1,2
    call setPar(A,kA, 0,m)
  enddo
endif

end subroutine SelectParent

! ######################################################################

subroutine CalcPOGPZ(A, kA, B, kB, pLLR)
use Global
implicit none

integer, intent(IN) :: A, kA, B, kB
double precision, intent(OUT) :: pLLR(4,2)  ! dam - sire - new pair - old pair , w/o - w age
integer :: AG, fcl, notfcl(6), x, curpar(2), mateB
double precision :: LLA(7,7, 2,2), LLS, LLP(4,2,2) 
logical :: ParClone               

pLLR = missing
curpar = getPar(A, kA)
if (curPar(kB) == B)  return

ParClone = .FALSE.
if (Hermaphrodites/=0 .and. A>0 .and. curPar(kB)==0 .and. &
 B>0 .and. curpar(3-kB) == B) then  
  ParClone = .TRUE.
  call setParTmp(A, kA, 0, 3-kB)
endif

LLA = missing
call CalcCandParLL(A, kA, B, kB, LLA)  

if (A > 0) then
  fcl = 1
  notfcl = (/2,3,4,5,6,7/)
else 
  fcl = 4
  notfcl = (/1,2,3,5,6,7/)
endif
 
if (ParClone) then  ! A>0 .and. B >0
  call setParTmp(A, kA, curPar(3-kB), 3-kB)
  call PairPO(A, B, kB, 1, LLS)   
  LLA(fcl,fcl,3-kB,:) = LLS   ! TODO: ALR?
endif

LLP = missing
do AG=1,2
if (curPar(kB)==0) then   ! --> LLA(:,:,kB,:) empty
   ! B + curPar(3-kB)
    LLP(3   ,1,AG) = LLA(fcl,fcl,3-kB,AG)
    LLA(fcl,fcl,3-kB,AG) = 555D0
    LLP(3   ,2,AG) = MaxLL(RESHAPE(LLA(:,:,3-kB,AG), (/7*7/)))
    ! only B
    LLP(kB  ,1,AG) = MaxLL(LLA(fcl,notfcl,3-kB,AG)) 
    LLP(kB  ,2,AG) = MaxLL(RESHAPE(LLA(notfcl,:,3-kB,AG), (/6*7/))) 
    ! only curPar(3-kB) 
    LLP(3-kB,1,AG) = MaxLL(LLA(notfcl,fcl,3-kB,AG))
    LLP(3-kB,2,AG) = MaxLL(RESHAPE(LLA(:,notfcl,3-kB,AG), (/6*7/))) 
  else if (all(curPar /= 0)) then                                      
    ! B + curPar(3-kB)
    LLP(3   ,1,AG) = MaxLL((/LLA(fcl,fcl,3-kB,AG), LLA(fcl,:,kB,AG)/))
    LLP(3   ,2,AG) = MaxLL((/RESHAPE(LLA(notfcl,:,:,AG), (/6*7*2/)), &
      LLA(fcl,notfcl,3-kB,AG)/))
    ! curPar pair                             
    LLP(4   ,1,AG) = MaxLL(RESHAPE(LLA(notfcl,fcl,:,AG), (/6*2/)))
    LLP(4   ,2,AG) = MaxLL((/RESHAPE(LLA(:,notfcl,:,AG),(/7*6*2/)), &
      LLA(fcl,fcl,3-kB,AG)/))
    LLA(fcl,fcl,:,AG) = 555D0     
    ! only B
    LLP(kB  ,1,AG) = MaxLL(LLA(fcl,notfcl,3-kB,AG))
    LLP(kB  ,2,AG) = MaxLL((/RESHAPE(LLA(notfcl,:,:,AG),(/2*6*7/)),&
      LLA(fcl,:,kB,AG), LLA(fcl,fcl,3-kB,AG)/)) 
    ! only curPar(3-kB)
    LLP(3-kB,1,AG) = MaxLL(RESHAPE(LLA(notfcl,notfcl,kB,AG), (/6*6/)))
    LLP(3-kB,2,AG) = MaxLL((/RESHAPE(LLA(fcl,:,:,AG), (/2*7/)), &
      RESHAPE(LLA(:,fcl,:,AG), (/2*7/))/))                        
  else 
    LLP(kB  ,1,AG) = LLA(fcl,7,3-kB,AG)
    LLP(kB  ,2,AG) = MaxLL(LLA(notfcl,7,3-kB,AG))     
  endif
enddo
            
pLLR = 555D0
do x=1,4
  if (all(LLP(x,:,:) < 0))  pLLR(x,:) = LLP(x,1,:) - LLP(x,2,:)
enddo  

if (Complx==0) then  ! monogamous  
  mateB = getMate(B,kB)
  if (mateB/=0 .and. LLP(3,1,1)>0d0)  pLLR(3,:) = pLLR(kB,:)    ! newpair = B + B's mate
  if (all(curPar/=0) .and. LLP(4,1,1)>0d0)   pLLR(4,:) = pLLR(3-kB,:)  ! oldpair = curpar(3-kB) + its mate
endif

end subroutine CalcPOGPZ

! ######################################################################

subroutine LLA2LLR(LLA, fcl, LLR)   ! output array from CalcCandParLL -> LR(PO/not)
use Global
implicit none

double precision, intent(IN) :: LLA(7,7,2,2)
integer, intent(IN) :: fcl
double precision, intent(OUT) :: LLR(2,2)  ! UseAge; vs not-PO / vs U  ;B, curPar (kB or 3-kB)
integer :: notfcl(6), a
double precision :: LLPO

if (fcl==1) then   ! A>0
  notfcl = (/2,3,4,5,6,7/)
else !if (fcl==4) then  ! A<0
  notfcl = (/1,2,3,5,6,7/)
endif

LLR = Missing
do a=1,2
  if (ALL(LLA(fcl,7,:,a) > 0)) then   ! maybeOtherParent / Impossible w/o co-parent
    LLR(a,1) = MINVAL(LLA(fcl,7,:,a))
  else
    LLPO = MaxLL( RESHAPE(LLA(fcl,:,:,a), (/2*7/)) )   ! B parent
    LLR(a,1) = LLPO - MaxLL( RESHAPE(LLA(notfcl,:,:,a), (/2*6*7/)) )
  endif

  if (ALL(LLA(7,fcl,:,a) > 0)) then  
      LLR(a,2) = MINVAL(LLA(7,fcl,:,a))
  else
    LLPO = MaxLL( RESHAPE(LLA(:,fcl,:,a), (/2*7/)) ) ! curPar parent
    LLR(a,2) = LLPO - MaxLL( RESHAPE(LLA(:,notfcl,:,a), (/2*6*7/)) )
  endif
enddo

end subroutine LLA2LLR

! ######################################################################

subroutine CalcCandParLL(A, kA, B, kB, LLA)
use Global
implicit none

! calc LL over A, candidate parent B, + current parent of A, under a set of diff relationships
! combination of prev. CalcPOZ & CalcGPZ

integer, intent(IN) :: A, kA, B, kB
double precision, intent(OUT) :: LLA(7,7,2,2)  ! B, curPar(3-n), curPar(n), Age
integer :: focal, m, curPar(2), mid(5), parB(2), r, fclx(3), z, curGP(2,2), BC(2), kBC(2)
double precision :: LLcp(3,2), LLU(4), U, ALR(4,2), LLtmp, LLtrio(6), &
 LLX(2,3), ALRtrio(2,6), ALRHS(2,2), ALRGP(2,3,3), LLcorner(2,2), ALRrev(2)
logical :: ParOK, GpOK(2), AreConnected(3), swapcorner(2), parrevOK(2), do_trios

LLA = missing
curPar = getPar(A, kA)      
if (curPar(kB) == B)  return
curGP = 0    
parB = getPar(B, kB) 

if (A > 0) then
  focal = 1
  mid = (/2,3,4,5,6/)
else if (A < 0) then
  focal = 4
  mid = (/1,2,3,5,6/)
else
  return
endif               
fclx = focal
if (A>0) then
  do m=1,2
    if (curPar(m)<0)  fclx(m) = 6
  enddo
  if (B<0)  fclx(3) = 6
endif

call CalcU(A, kA, 0, 0, U)  ! ensure LL(A) is up to date

if (ALL(curPar==0)) then
  call CheckRel(A, kA, B, kB, fclx(3), LLA(:,7,3-kB,1), LLA(:,7,3-kB,2))

else
  ParOK = .TRUE.
  do m=1,2
    curGP(:,m) = getPar(curPar(m), m)
  enddo   
  
  LLCP = 0D0
  LLU = missing
  call Calc4U(curPar, B,kB, A,kA, LLU, LLCP)
  ! LLU: curPar(1), curPar(2), B, A  ; 1-3 w/o A if <0
  ! LLCP: correction factors        
  
  ALR = 0D0
  do m=1,2
    call CalcAgeLR(A,kA, curPar(m), m, 0,1, .FALSE., ALR(m,2))  ! D2: without/with age LR
  enddo 
  call CalcAgeLR(A,kA, B, kB, 0,1, .FALSE., ALR(3,2))  
  ! above is conditional on curpar, in case A becomes FS with B's & different ALR for FS than HS
  ! TODO: check if becomes FS. if so, no ALR correction needed - already included in checkrel() (RIGHT ???)
  ! if so: ALR(3,2)=0
  
   
  do m=1,2 ! sex currently assigned par
    if (curPar(m)==0) cycle 
    
    ! switch regarding which LL to pick for corner where 1 is parent
    ! (most inclusive consideration of all possibilities (POHA,..) and all linked sibships)
    if (A>0) then
      swapcorner(1) = .FALSE.
      swapcorner(2) = .TRUE.  ! A+B considers POHA; A+curpar | B does not  
    else 
      swapcorner(1) = .TRUE.   ! LLA(7,focal); only curpar(m) PO; default from A + curpar(m)
      swapcorner(2) = .FALSE.  ! LLA(focal,7); only B PO; default from A + curpar(m)
      if (curpar(m)<0) then
        if (m/=kB .and. ANY(Parent(SibID(1:nS(-curPar(m),m), -curPar(m),m), kB)==B)) then
          swapcorner(1) = .FALSE.  ! LLA(7,focal)
          !swapcorner(2) = .TRUE.  ?
        endif
      else if (B<0) then
        if (m/=kB .and. ANY(Parent(SibID(1:nS(-B,kB), -B,kB),m)==curPar(m))) then
          swapcorner(2) = .TRUE.  ! LLA(focal,7)
          ! swapcorner(1) = .FALSE. ?
        endif
      endif
    endif      
   
    call checkRel(A, kA, B, kB, fclx(3), LLA(:,focal,m,1), LLA(:,focal,m,2))   ! curPar(m)=GP + A_7 
    LLcorner(:,1) = LLA(7,focal,m,:)  
    call setParTmp(A, kA, 0, m)
    call checkRel(A, kA, B, kB, fclx(3), LLA(:,7,m,1), LLA(:,7,m,2))   ! A_7
    LLcorner(:,2) = LLA(focal,7,m,:)
    call checkRel(A, kA, curPar(m), m, fclx(m), LLA(7,:,m,1), LLA(7,:,m,2))  ! curPar(m)_7
    if (swapcorner(1))  LLA(7,focal,m,:) = LLcorner(:,1)  !!  .or. LLA(7,focal,m,:)>0d0 ?
    call CalcAgeLR(A,kA, B, kB, 0,1, .FALSE., ALR(4,2))   ! not conditional on curpar(m)
    
    call ChkValidPar(A, kA, B, kB, ParOK)
    if (ParOK) then 
      call setParTmp(A, kA, B, kB)
      call checkRel(A, kA, curPar(m),m, fclx(m), LLA(focal,:,m,1), LLA(focal,:,m,2)) 
      if (swapcorner(2))  LLA(focal,7,m,:) = LLcorner(:,2)  
    endif

    call setParTmp(A, kA, curPar(kB), kB)
    call setParTmp(A, kA, curPar(3-kB), 3-kB)
    
    ! fix: B or curpar parent as side-effect of curpar resp B as HS    
    if (A>0) then
      if (curPar(m) == parB(m))  LLA(2:3,7,m,:) = impossible  ! HS implies curPar = Par
      if (curGP(kB,m) == B)      LLA(7,2:3,m,:) = impossible  ! HS implies B = Par
      if (curPar(3-m)/=0) then     
        if (curPar(3-m) == parB(3-m))  LLA(2,7,m,:) = impossible
        if (curGP(kB,3-m) == B)        LLA(7,2,m,:) = impossible
      endif
    endif
    if (Complx==0) then
      LLA(focal,7,m,:) = impossible
      LLA(7,focal,m,:) = impossible
    endif

    do z=1,2                        
      WHERE (LLA(mid,focal,m,z)<0) LLA(mid,focal,m,z) = LLA(mid,focal,m,z) + LLcp(3,m) + ALR(m,z)   ! ALR(:,1)=0
      WHERE (LLA(mid,7,m,z)<0)     LLA(mid,7,m,z)     = LLA(mid,7,m,z)     + LLcp(3,m)
      IF  (LLA(focal,focal,m,z)<0)  LLA(focal,focal,m,z) = LLA(focal,focal,m,z) + LLcp(m,m) + ALR(3,z)
      WHERE (LLA(focal,mid,m,z)<0)  LLA(focal,mid,m,z) = LLA(focal,mid,m,z) + LLcp(m,m) + ALR(4,z)
      WHERE (LLA(7,(/mid,7/),m,z)<0)      LLA(7,(/mid,7/),m,z)     = LLA(7,(/mid,7/),m,z)     + LLcp(m,m) + ALR(3-m,z)
      if (swapcorner(1)) then
        if (LLA(7,focal,m,z)<0)  LLA(7,focal,m,z) = LLA(7,focal,m,z) + LLcp(3,m) + ALR(m,z)  
      else
        if (LLA(7,focal,m,z)<0)  LLA(7,focal,m,z) = LLA(7,focal,m,z) + LLcp(m,m) + ALR(3-m,z)
      endif
      if (swapcorner(2)) then
        if (LLA(focal,7,m,z)<0)  LLA(focal,7,m,z) = LLA(focal,7,m,z) + LLcp(3,m)
      else
        if (LLA(focal,7,m,z)<0)  LLA(focal,7,m,z) = LLA(focal,7,m,z) + LLcp(m,m) + ALR(4,z)
      endif
!      LLcorner(z,1) = LLcorner(z,1) + LLcp(3,m)              ! for debug only
!      LLcorner(z,2) = LLcorner(z,2) + LLcp(3,m) + ALR(m,z)   ! for debug only
    enddo
    
    if (m == kB) then
      WHERE (LLA(focal,mid,m,2)<0) LLA(focal,mid,m,2) = LLA(focal,mid,m,2) + ALR(3-m,2)
      LLA(focal,focal,m,:) = impossible ! cannot have 2 same-sex parents
    endif
    
     if (A<0 .and. ParB(m)==curPar(m)) then  ! then 'B FS of A' identical with & without curpar
      LLA(5,7,m,:) = impossible
    endif

    call setParTmp(A, kA, 0, m)  ! for ageLR & connect check
    ! trio likelihoods. Simplified, assuming A,B,C are unconnected
    if (parB(m)==curPar(m) .or. curPar(3-m)/=0 .or. &
      (kB==m .and. (B==curPar(m) .or. B==curGP(m,m)))) then
      do_trios = .FALSE.
    else
      call connected(A,kA, B,kB, AreConnected(1))
      call connected(A,kA, curPar(m),m, AreConnected(2))
      call connected(B,kB, curPar(m),m, AreConnected(3))
      if (any(AreConnected)) then
        do_trios = .FALSE.
      else
        do_trios = .TRUE.
      endif
    endif
    LLtrio = Missing
    ALRtrio = 0D0        
    if (do_trios) then                       
      ! check if B & curPar(m) both FS / both GP of A, + HS & FA if all >0
      do r=2,6
        if (r/=5 .and. any(ParB == curGP(:,m) .and. parB/=0))  cycle           
        if (r==2 .and. A<0 .and. (any(parB/=0) .or. any(curGP(:,m)/=0)))  cycle
        if (r==3 .and. (Complx==0 .or. any(parB/=0) .or. any(curGP(:,m)/=0)))  cycle     
        if (r==4 .and. A<0)  cycle
        if (r==5)  cycle           
        if (r/=6) then
          call CalcAgeLR(A,kA, B,kB, kB,r, .FALSE., ALRtrio(1,r))
          call CalcAgeLR(A,kA, CurPar(m),m, kB,r, .FALSE., ALRtrio(2,r))
          if (any(ALRtrio(:,r)==impossible)) then
            if (r==4 .and. (B<0 .or. curpar(m)<0)) then  ! proxy for other rels
              if (ALRtrio(1,r)==impossible)  ALRtrio(1,r) = 0d0
              if (ALRtrio(2,r)==impossible)  ALRtrio(2,r) = 0d0
            else
              cycle
            endif
          endif                          
        endif
        if (r==2)  call trioFS(A,kA, B,kB, curPar(m),m, LLtrio(2))
        if (r==3)  call trioHS(A,kA, B,kB, curPar(m),m, LLtrio(3))  ! B & curPar HS to eachother or not.
        if (r==4)  call trioGP(A,kA, B,kB, curPar(m),m, LLtrio(4))  ! A>0 only
        if (r==5)  call trioFA(A,kA, B,kB, curPar(m),m, LLtrio(5)) 
        if (r==6) then
          if (A<0) then
            call CalcAgeLR(A,kA, B,kB, kB,4, .FALSE., ALRtrio(1,r))
            call CalcAgeLR(A,kA, CurPar(m),m, kB,4, .FALSE., ALRtrio(2,r))
            if (any(ALRtrio(:,r)==impossible)) then
              ALRtrio(:,r) = 0D0   ! use as proxy for other 3rd degree rels
            endif
            call trioGP(A,kA, B,kB, curPar(m),m, LLtrio(6))           
          else if (B>0 .and. curPar(m)>0) then  
            if (all(ParB/=0) .or. all(curGP(:,m)/=0)) then
              call trioGGP(A,kA, B,curPar(m), LLtrio(6))
              call CalcAgeLR(A,kA, B,kB, kB,4, .FALSE., ALRtrio(1,r))
              call CalcAgeLR(A,kA, CurPar(m),m, kB,4, .FALSE., ALRtrio(2,r))
              do z=1,2  ! TODO: ALR for GGP 
                if (ALRtrio(z,r)==impossible)  ALRtrio(z,r)=0.0D0
              enddo
            else
              call CalcAgeLR(A,kA, B,kB, kB,5, .FALSE., ALRtrio(1,r))
              call CalcAgeLR(A,kA, CurPar(m),m, kB,5, .FALSE., ALRtrio(2,r)) 
              if (any(ALRtrio(:,r)==impossible)) cycle              
              call trioHA(A,kA, B,curPar(m), LLtrio(6))                         
            endif 
          else
            LLtrio(6) = NotImplemented   ! (yet)
            cycle
          endif
        endif        
        if (LLtrio(r) < -HUGE(0D0) .or. LLtrio(r)>=NotImplemented)  cycle
        LLA(r,r,m,:) = LLtrio(r) + LLU(3-m)! &  ! curPar(3-m)             
        LLA(r,r,m,2) = addALR( addALR(LLA(r,r,m,2), ALRtrio(1,r)), ALRtrio(2,r) )
      enddo
    endif
      ! HS + GP
      LLX = Missing
      GpOK = .FALSE. 
      ALRHS = Missing             
      ALRGP = neg_missing           
      if (curPar(3-m)==0 .and. Complx>0) then   ! A<0  ! TODO? curPar(3-m)/=0
        do z=1,2
          call ChkValidGP(A,kA,CurPar(m),m, z, GpOK(1), ALRGP(z,1,2))      
          call ChkValidGP(B,kB,CurPar(m),m, z, GpOK(2), ALRGP(z,3,2))
          if (.not. all(GpOK))  cycle          
          call CalcAgeLR(A,kA, B,kB, z,3, .FALSE., ALRHS(1,z)) 
          if (ALRHS(1,z)==impossible .or.ALRHS(1,z) < 5*TF)  cycle
          call trioHSGP(A,kA, B,kB, CurPar(m),m, z, LLX(1,z))
        enddo
        do z=1,2
          call ChkValidGP(A,kA,B,kB, z, GpOK(1), ALRGP(z,1,3))
          call ChkValidGP(CurPar(m),m, B,kB, z, GpOK(2), ALRGP(z,2,3))
          if (.not. all(GpOK))  cycle  
          call CalcAgeLR(A,kA, CurPar(m),m, z,3, .FALSE.,ALRHS(2,z))
          if (ALRHS(2,z)==impossible .or. ALRHS(2,z) < 5*TF)  cycle
          call trioHSGP(A,kA, CurPar(m),m, B,kB, z, LLX(2,z))         
        enddo
        if (any(LLX < 0d0)) then
          LLA(6,5,m,1) = MaxLL((/LLX(:,1), LLX(:,2)/))
          do z=1,2
            if (ABS(LLX(1,z) - LLA(6,5,m,1)) < TINY(0.0d0)) then
              LLA(6,5,m,2) = LLA(6,5,m,1) + ALRHS(1,z) + ALRGP(z,1,2) + ALRGP(z,3,2)
            else if (ABS(LLX(2,z) - LLA(6,5,m,1)) < TINY(0.0d0)) then
              LLA(6,5,m,2) = LLA(6,5,m,1) + ALRHS(2,z) + ALRGP(z,1,3) + ALRGP(z,2,3)
            endif
          enddo
        endif
      endif  

      if (Complx/=1 .and. curpar(3-m)==0) then
        if (all(ParB == 0))       call trioFSFA(A,kA, B,kB, curpar(m),m, LLX(1,3))
        if (all(curGP(:,m) == 0)) call trioFSFA(A,kA, curpar(m),m, B,kB, LLX(2,3))
        LLA(2,5,m,1) = LLX(1,3)
        LLA(2,5,m,2) = addALR( addALR(LLX(1,3), ALRtrio(1,2)), ALRtrio(2,5) )
        LLA(5,2,m,1) = LLX(2,3)
        LLA(5,2,m,2) = addALR( addALR(LLX(2,3), ALRtrio(2,2)), ALRtrio(1,5) )
      endif      
               
    call setParTmp(A, kA, curPar(m), m)  
        
    if (kA>2)  cycle
    ParrevOK = .FALSE.
    if (curpar(3-m)==0 .and. (A>0 .or. B>0 .or. curpar(m)>0)) then   ! m/=kB .and. 
      ! check if B is offspring & curpar(m) parent, or vv, or both offspring
      call setParTmp(A, kA, 0, m)    
      if (parB(kA)==0)  call ChkValidPar(B, kB, A, kA, ParrevOK(1))
      if (curGP(kA,m)==0) call ChkValidPar(curPar(m), (m), A, kA, ParrevOK(2))
      BC = (/B, curpar(m)/)
      kBC = (/kB, m/)
      do z=1,2    ! 1: B offspring of A; 2: curpar(m) offspring of A
        if (.not. parrevOK(z))  cycle
        do r=1,2  ! 1: other is parent of A; 2: other is offspring of B/C 
          if (r==2) then
            if (z==1 .and. curGP(kB,m)/=0)  cycle
            if (z==2 .and. parB(m)/=0)  cycle
          endif
          if (A>0 .or. r==2) then
            call CalcU(B,kB, curPar(m),m, LLX(1,1)) 
          else
            call CalcU(A,kA, BC(3-z), kBC(3-z), LLX(1,1))
          endif
          
          call setParTmp(BC(z), kBC(z), A, kA)  
          ParOK = .FALSE.
          if (r==1)  call ChkValidPar(A, kA, BC(3-z), kBC(3-z), ParOK)
          if (r==2)  call ChkValidPar(BC(3-z), kBC(3-z), BC(z), kBC(z), ParOK)
          if (.not. ParOK) then
            call setParTmp(BC(z), kBC(z), 0, kA)  
            cycle
          endif          
          if (r==1)  call setParTmp(A, kA, BC(3-z), kBC(3-z))
          if (r==2)  call setParTmp(BC(3-z), kBC(3-z), BC(z), kBC(z))   
          
          if (A>0 .or. r==2) then
            call CalcU(B,kB, curPar(m),m, LLX(2,1))       
          else
            call CalcU(A,kA, BC(3-z), kBC(3-z), LLX(2,1))    
          endif
          
          call setParTmp(BC(z), kBC(z), 0, kA)
          if (r==1)  call setParTmp(A, kA, 0, kBC(3-z))
          if (r==2)  call setParTmp(BC(3-z), kBC(3-z), 0, kBC(z)) 

          LLtmp = LLX(2,1) - LLX(1,1) + LLA(7,7,m,1)   ! change in LL relative to both unrelated    
          if (LLtmp > LLA(6,6,m,1) .or. LLA(6,6,m,1)>0) then
            LLA(6,6,m,1) = LLtmp  
            call CalcAgeLR(BC(z), kBC(z), A, kA, 0,1, .FALSE., ALRrev(1))
            if (r==1)  call CalcAgeLR(A, kA, BC(3-z), kBC(3-z), 0,1, .FALSE., ALRrev(2))
            if (r==2)  call CalcAgeLR(BC(3-z), kBC(3-z), BC(z), kBC(z), 0,1, .FALSE., ALRrev(2))
            LLA(6,6,m,2) = addALR( addALR(LLA(6,6,m,1), ALRrev(1)), ALRrev(2))
          endif
        enddo
      enddo
      call setParTmp(A, kA, curPar(m), m)  
    endif

  enddo  ! m
endif

end subroutine CalcCandParLL

! #####################################################################

subroutine CheckMaybeRev(A, kA, candP, kP, maybe, dLL)   ! check if A could be parent of candP instead.
use Global
implicit none

integer, intent(IN) :: A, kA, candP, kP
logical, intent(OUT) :: maybe
double precision, intent(OUT) :: dLL(2)
integer :: ParCP(2), ParA(2), fcl, n, notfcl(6), fclx
double precision :: ALR(2), LLrev(7,2), LLtmp(3), LLrevX(7,2)
logical :: SexUnk, ParOK                   

dLL = missing
maybe = .TRUE.

if (CandP == 0)  return

SexUnk = .FALSE.
if (A > 0) then
  if (Sex(A)>3)  SexUnk = .TRUE.
endif                              
ParCP = getPar(candP, kP)
ParA = getPar(A,kA)    
 
if (ALL(ParCP > 0)) then
  if (hermaphrodites/=0 .and. candP>0 .and. A>0) then
    if (SelfedIndiv(CandP)) then
      maybe = .TRUE.
    else
      maybe = .FALSE.
    endif
  else
    maybe = .FALSE.
  endif
  if (.not. maybe)  return
else if (Complx==0 .and. candP < 0 .and. ParA(3-kP)/=0) then
  if (DumMate(-CandP, kP) == ParA(3-kP)) then
    maybe = .FALSE.
    return
  endif
else if (A>0 .and. .not. SexUnk) then
  if (ParCP(kA) > 0 .and. Sex(A)<3) then
    maybe = .FALSE.
    return
  endif
endif

ParOK = .TRUE.
call ChkValidPar(candP, kP, A, kA, ParOK)
if (.not. ParOK) then
  maybe = .FALSE.
  return
endif

ALR = missing
call CalcAgeLR(A, kA, candP, kP, 0,1, .TRUE., ALR(1))
call CalcAgeLR(candP, kP, A, kA, 0,1, .TRUE., ALR(2))

if (ALR(2)==impossible .or. ALR(1)-ALR(2) > 3.0*ABS(TF)) then
  maybe = .FALSE.
  return
else if (ALL(ABS(ALR) < 0.01)) then
  if (all(ParA==0) .and. all(ParCP==0))  return   ! No way to tell which is parent & which offspring
endif

LLrev = missing
if (hermaphrodites/=0 .and. ALL(ParCP > 0) .and. ParCP(1)==ParCP(2)) then
  call ChkValidPar(ParCP(1), 1, A, kA, ParOK)
  if (ParOK) then
    call CheckRel(ParCP(1), 1, A, kA, 1, LLrev(:,1), LLrev(:,2)) 
    notfcl = (/ (n, n = 2, 7) /)    
    do n=1,2
      if (LLrev(1,n) < 0) then
        maybe = .TRUE.
        dLL(n) = LLrev(1, n) - MaxLL(LLrev(notfcl, n))
      else
        maybe = .FALSE.
      endif
    enddo
  else
    maybe = .FALSE.
  endif
  return
endif

fclx = 7
if (candP>0) then
  fcl = 1
  notfcl = (/ (n, n = 2, 7) /)
  if (A>0)  fclx = 1
else ! if (candP<0) then
  fcl = 4
  notfcl = (/1,2,3,5,6,7/)
endif

LLrev = missing
LLtmp = missing
call setParTmp(candP, kP, 0, kA)
call ChkValidPar(candP, kP, A, kA, ParOK)
if (ParOK) then  
  if (ParCP(3-kA) < 0) then  ! include changes in CLL & validPar chk with other sibs
    fcl = 1
    call CheckRel(ParCP(3-kA), 3-kA, A, kA, fcl, LLrev(:,1), LLrev(:,2))   
  else
    call CheckRel(candP, kP, A, kA, fclx, LLrev(:,1), LLrev(:,2))
  endif
endif                                   
call setParTmp(candP, kP, parCP(kA), kA)

LLrevX = missing                
if (SexUnk .and. ParOK .and. Hermaphrodites/=2) then  ! implies A>0; also consider as parent of other sex
  call ChkValidPar(candP, kP, A, 3-kA, ParOK)
  if (ParOK) then
    if (candP > 0) then
      call CheckRel(candP, 3-kA, A, 3-kA, fclx, LLrevX(:,1), LLrevX(:,2))
    else
      call CheckRel(candP, kP, A, 3-kA, fclx, LLrevX(:,1), LLrevX(:,2))
    endif
    if (ParCP(kA) < 0) then  ! include changes in CLL
      call CalcU(candP, kP, parCP(kA),kA, LLtmp(1))
      call setParTmp(candP, kP, A, 3-kA)
      call CalcU(candP, kP, parCP(kA),kA, LLtmp(2))
      LLrevX(fcl,1) = LLrevX(fcl,1) + (LLtmp(2) - LLtmp(1))
      call setParTmp(candP, kP, parCP(3-kA), 3-kA)
    endif
  endif
endif     

do n=1,2
  if (LLrev(fcl,n) < 0) then
    maybe = .TRUE.
    dLL(n) = LLrev(fcl, n) - MaxLL(LLrev(notfcl, n))
    if (LLrevX(fcl,n) < 0) then
      dLL(n) = MAX(dLL(n), LLrevX(fcl, n) - MaxLL(LLrevX(notfcl, n)))
    endif
  else
    maybe = .FALSE.
  endif
enddo

end subroutine CheckMaybeRev

! #####################################################################

subroutine FindPairs(PairID, PairType)
use Global
use sort_module
implicit none

integer, intent(OUT) :: PairID(XP*nInd,2), PairType(XP*nInd)        
logical :: UseAge, cPair, matpat(2)
integer :: k, i, j, top, PairTypeTmp(XP*nInd), PairIDtmp(XP*nInd,2), x
double precision :: dLL, PairLLRtmp(XP*nInd), LL(7), LLg(7), LRS(2), PairLLR(XP*nInd), ALR  
integer, allocatable, dimension(:) :: Rank
double precision, allocatable, dimension(:) :: SortBy                                     

nPairs = 0
PairID = -9
PairLLR = missing
PairType = 0
UseAge = AgePhase > 0

do i=1,  nInd-1  
  if (MODULO(i,100)==0) call rchkusr()
  if (quiet==-1 .and. any(chunk_printdot_i==i)) call print_dot()
  if (ALL(Parent(i,:)/=0)) cycle
  do j=i+1,nInd
    if (hermaphrodites==1 .and. ((ANY(parent(i,:)/=0) .and. ALL(parent(j,:)==0)) .or. &
     (ALL(parent(i,:)==0) .and. ANY(parent(j,:)/=0))))  cycle  
    LRS = 0D0
    matpat = .FALSE.
    do k=1,2
      if (Parent(i,k)/=0 .or. Parent(j,k)/=0) cycle
      if (Parent(i,k)==j .or. Parent(j,k)==i) cycle
      if (UseAge) then
        ALR = getAP(AgeDiff(i,j), 3, 0, k, Impossible)
        if (ALR == Impossible)  cycle
      endif 
      matpat(k) = .TRUE.
    enddo
    if (.not. any(matpat))  cycle                             
    if (DoMtDif) then
      if (mtDif(i,j))  matpat(1) = .FALSE.   ! not sharing same mt haplotype --> not mat sibs    
      if (.not. any(matpat))  cycle
    endif   
    if (Complx==0) then
      call PairQFS(i, j, LRS(2))  ! quick check
      if (LRS(2) < 2*TF) cycle  
    else
      call PairQHS(i, j, LRS(1))
      if (LRS(1) < 2*TF) cycle  
    endif
    if (Complx>0 .and. ((ALL(Parent(i,:)==0) .and. ALL(Parent(j,:)==0) .and. &
      UseAge .and. ALL(matpat)) .or. &
       (Hermaphrodites/=0 .and. (ALL(Parent(i,:)==0) .or. ALL(Parent(j,:)==0))))) then
      call PairQFS(i, j, LRS(2)) 
    endif
    
    cPair = .FALSE.    
    do k=1,2
      if ((.not. matpat(k)) .or. cPair)  cycle
      if (Complx==0 .and. k==2)  cycle
      if (k==2 .and. matpat(1)) then
        if (.not. UseAge) cycle
        if (LRS(2) < TF) cycle  
      endif    
      if (hermaphrodites==1 .and. (ALL(Parent(i,:)==0) .or. ALL(Parent(j,:)==0))) then
        if (LRS(2) < TF) cycle 
      endif  
      if (Complx==0) then
        x = 2
      else
        x = 3
      endif
      LLg = missing
      LL = missing
      if (AgeDiff(i,j)>=0) then
        call CheckPair(i, j, k, x, LLg, LL)
      else
        call CheckPair(j, i, k, x, LLg, LL)
      endif
      top = 0
      dLL = missing
      if (UseAge)  call BestRel(LL, 3, top, dLL)
      if (.not. UseAge)  call BestRel(LLg, 3, top, dLL)
      if (hermaphrodites==1 .and. (ALL(Parent(i,:)==0) .or. &
        ALL(Parent(j,:)==0)) .and. top/=2)  cycle
      if (top==2 .or. top==3) then  
        if (nPairs >= XP*nInd) cycle  ! do in next round
        nPairs = nPairs+1
        PairID(nPairs, :) = (/ i, j /)
        PairLLR(nPairs) = dLL
        if (k==1 .and. matpat(2)) then
          pairType(nPairs) = 3
          cPair = .TRUE.
        else
          PairType(nPairs) = k
        endif
      endif
    enddo
  enddo
enddo

! sort by decreasing dLL
PairIDtmp = 0
PairLLRtmp = 0D0
allocate(Rank(nPairs))
allocate(SortBy(nPairs))
Rank = (/ (i, i=1, nPairs, 1) /)
SortBy = PairLLR(1:nPairs)
 
 call QsortC(SortBy, Rank(1:nPairs))
do i=1,nPairs
  PairTypeTmp(i) = PairType(Rank(nPairs-i+1))  ! decreasing order
  PairIDtmp(i,1:2) = PairID(Rank(nPairs-i+1), 1:2)  
  PairLLRtmp(i) = PairLLR(Rank(nPairs-i+1)) 
enddo 

PairType = PairTypeTmp
PairID = PairIDtmp 
PairLLR = PairLLRtmp
deallocate(Rank)
deallocate(SortBy)
                                                              
end subroutine FindPairs

! #####################################################################

subroutine PairQFS(A, B, LR)  !quick check, not conditioning on parents.
use Global
implicit none

integer, intent(IN) :: A, B
double precision, intent(OUT) :: LR
integer :: l
double precision :: PrL(nSnp)

PrL = 0D0
do l=1,nSnp
  PrL(l) = LOG10(PFS(Genos(l,A), Genos(l,B), l))  
enddo
LR = SUM(PrL)

end subroutine PairQFS

! #####################################################################

subroutine PairQHS(A, B, LR)  !quick check, not conditioning on parents.
use Global
implicit none

integer, intent(IN) :: A, B
double precision, intent(OUT) :: LR
integer :: l
double precision :: PrL(nSnp)

PrL = 0D0
do l=1,nSnp
  PrL(l) = LOG10(PHS(Genos(l,A), Genos(l,B), l)) ! note: >0 for FS 
enddo
LR = SUM(PrL)

end subroutine PairQHS

! #####################################################################

subroutine IsSelfed(A, withFS, LR)  ! A is product of selfing, not conditioning on parents.
use Global
use CalcLik           
implicit none

integer, intent(IN) :: A
logical, intent(IN) :: withFS
double precision, intent(OUT) :: LR
integer :: l, x,u, z, v,y
double precision :: PrX(3), PrXY(3,3,3), PrL(nSnp,4), LLtmp(4), PrZV(3,3), PrA(3,3) 

PrL = 0D0
do l=1,nSnp
  do x=1,3
    if (withFS) then
      PrA = FSLik(l,A)
    else
      PrA = OKA2P(Genos(l,A), :, :)
    endif
    PrX(x) = PrA(x, x) * AHWE(x,l)  ! selfed
    PrXY(x,:,1) = PrA(x, :) * AHWE(x,l) * AHWE(:,l)   ! parents U
    PrXY(x,:,2) = PrA(x, :) * AKAP(x,:,l) * AHWE(:,l)   ! parents PO
    do y=1,3
      do z=1,3
        do v=1,3
          PrZV(z,v) = PrA(x,y) * AKA2P(x,z,v) * AKA2P(y,z,v) * &
            AHWE(z,l) * AHWE(v,l)    ! parents FS. identical to parents PO?
        enddo
      enddo
      PrXY(x,y,3) = SUM(PrZV)
    enddo
  enddo
  PrL(l,1) = LOG10(SUM(PrX))
  do u=1,3
    PrL(l,u+1) = LOG10(SUM(PrXY(:,:,u)))
  enddo
enddo
LLtmp = SUM(PrL,DIM=1)
LR = LLtmp(1) - MAXVAL(LLtmp(2:4))

end subroutine IsSelfed

! #####################################################################

subroutine CheckPair(A, B, kIN, focal, LLg, LL) 
! joined LL A,B under each hypothesis
use Global
implicit none

integer, intent(IN) :: A,B,kIN, focal
double precision, intent(OUT) :: Llg(7), LL(7)  ! PO,FS,HS,GG,FAU,HAU,U
integer :: x, cgp(2), k, AB(2), i, z
double precision :: ALR(7), LLGR(3), ALRgr(3), LLGGP(5), ALRgg(5), &
  LLtmpAU(2,3), ALRAU(2,3), LLCC, LLFC, ALRp(2,2), LLX(4), LLPA(3), LLFAx(2), &
  LLZ(8), LLC(7,2), LLHH(2,2),  LLP(2,2), LLFA(2), LLPS(2,2), LLU, ALRf(2), ALRs, LLS 
logical :: fclsib, hasparent(2,2)

LLg = missing
ALR = missing
LL = missing   

if (kIN>2) then
  k = 1
else
  k = kIN
endif     
AB = (/A, B/) 

HasParent = .FALSE.
do i=1,2
  do z=1,2
    if (Parent(AB(i),z) > 0) then
      HasParent(z,i) = .TRUE.
    else if (Parent(AB(i),z) < 0) then
      if (ns(-Parent(AB(i),z),z) >1 .or. any(GpID(:,-Parent(AB(i),z),z)/=0)) then
        HasParent(z,i) = .TRUE.   ! singleton sibship w/o GPs does not count
      endif
    endif
  enddo
enddo

if (focal==1 .and. Sex(B)/=k .and. Sex(B)<3) then 
  LLg(1) = impossible
  LL(1) = impossible
  return
endif

call CalcU(A,k,B,k, LLg(7))
LL(7) = LLg(7)

if (hermaphrodites >0) then
  if (any(Parent(B,:) == A)) then
    LLg(2) = impossible
    if (focal==1 .or. focal==4) then  ! reverse (B-A) is possible
      LLg(1) = impossible
      LLg(4) = impossible
    endif
  else if (any(Parent(A,:) == B)) then
    LLg(2) = impossible
  endif
  LL = LLg
  if (LLg(focal) == impossible)  return
endif

do x=2,4   ! not x=1: done by ChkValidPar()
  if (focal /= x)  cycle
  call CalcAgeLR(A,Sex(A), B, Sex(B), k, x, .TRUE., ALR(x)) 
  if (ALR(x)==impossible) then
    LLg(x) = impossible
    LL = LLg
    return
  endif
enddo

if (focal == 2) then
  call PairFullSib(A, B, LLg(2))
else if (focal == 3) then
  call PairHalfSib(A, B, k, LLg(3)) 
  if (any(parent(A,:)==Parent(B,:) .and. Parent(A,:)/=0)) then
    call PairFullSib(A, B, LLg(2)) 
  endif
else if (focal == 4) then
  call PairGP(A, B, k, focal, LLg(4))
endif

if (focal<5 .and. (LLg(focal)==impossible .or. ((LLg(focal)-LLg(7)) < TA .and. .not. &
  (focal==3 .and. LLg(2)<impossible .and. (LLg(2)-LLg(7)) >= TA)))) then   
  LL = LLg
  return
endif

do x=1,4
  if (ALR(x) /= missing) cycle
  call CalcAgeLR(A,Sex(A), B, Sex(B), k, x, .TRUE., ALR(x)) 
  if (ALR(x) == impossible) then
    LLg(x) = impossible
    LL(x) = impossible
  endif
enddo

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 1-4 non-focals

if (LLg(1)==missing)  call PairPO(A, B, k, focal, LLg(1)) 
if (LLg(2)==missing)  call PairFullSib(A, B, LLg(2)) 
if (LLg(3)==missing .and. Complx>0)   call PairHalfSib(A, B, k, LLg(3))    
if (LLg(4)==missing)  call PairGP(A, B, k, focal, LLg(4)) 
      
do x=1,4
  LL(x) = addALR(LLg(x), ALR(x)) 
enddo
!~~~  GGP  ~~~~~~~
LLGGP = missing
ALRgg = 0D0        
if (ALR(4)/=impossible) then   ! no ageprior for GGP
  if (AgeDiff(A,B)>=3) then
    call PairGGP(A, B, k, focal, LLGGP(1))
    call PairGGP(A,B,3-k, focal, LLGGP(2))    ! TODO drop?
  endif
  do x=1,3  ! hf
    if (Complx==0 .and. x<3)  cycle
    call PairGA(A, B, k, x, LLGGP(2+x))   ! HS of GP 4th degree rel, but giving false pos. 
  enddo
endif
if (LLGGP(1) <0D0) then
  if (Parent(A,k)/=0) then
    call CalcAgeLR(Parent(A,k),k, B,Sex(B), k, 4, .TRUE., ALRgg(1))
    ALRgg(2) = ALRgg(1)
    call CalcAgeLR(Parent(A,k),k, B,Sex(B), k, 6, .TRUE., ALRgg(3))
    ALRgg(4) = ALRgg(3)
    call CalcAgeLR(Parent(A,k),k, B,Sex(B), k, 5, .TRUE., ALRgg(5))
  else
    ! TODO ageprior
  endif
endif

!~~~  FA/HA  ~~~~~~~
LLtmpAU = missing
ALRAU = missing
LLFAx = missing
! call CalcAgeLR(A,k, B,k, 0, 6, .TRUE., ALR(6))
do i=1,2   ! A, B
  do x=1,3  ! mat, pat, FS
    if (x<3 .and. Complx>0) then
      call CalcAgeLR(AB(i),k, AB(3-i),x, k, 6, .TRUE., ALRAU(i,x))
    else
      call CalcAgeLR(AB(i),k, AB(3-i),3, k, 5, .TRUE., ALRAU(i,x))
    endif
    if (ALRAU(i,x)/=impossible .and. .not. (Complx==0 .and. x<3)) then
      call PairUA(AB(i), AB(3-i),k, x, LLtmpAU(i,x))
      if (Complx/=1)  call pairFAx(AB(i), AB(3-i), LLFAx(i))  ! A inbred   
    endif
  enddo  
enddo

LLg(5) = MaxLL((/LLtmpAU(:,3), LLFAx/))
if (Complx>0)  LLg(6) = MaxLL(RESHAPE(LLtmpAU(:,1:2), (/2*2/) ))
do i=1,2   
  do x=1,3
    LLtmpAU(i,x) = addALR(LLtmpAU(i,x), ALRAU(i,x))
  enddo
  LLFAx(i) = addALR(LLFAx(i), ALRAU(i,3))
enddo
LL(5) =  MaxLL((/LLtmpAU(:,3), LLFAx/))
if (Complx>0)  LL(6) = MaxLL(RESHAPE(LLtmpAU(:,1:2), (/2*2/) ))
ALR(6) = MAXVAL((/ALRAU(1,1:2), ALRAU(2,1:2)/))

LLCC = missing
call PairCC(A, B, k, LLCC)   ! full 1st cousins

LLg(6) = MaxLL( (/LLg(6), LLGGP, LLCC/) )  ! most likely 3rd degree
do x=1,5
  LLGGP(x) = addALR(LLGGP(x), ALRgg(x))
enddo
LL(6) = MaxLL( (/LL(6), LLGGP, LLCC/) )

!~~~ FS + HC ~~~
LLFC = missing
if (Complx==2 .and. LLg(2)<0D0 .and. Parent(A,3-k)==Parent(B,3-k) .and. &
  (MaxLL(LLtmpAU(:,3)) - MaxLL(LLg(2:3)) > -TA)) then
  call FSHC(A, B, k, LLFC)
  if (LLFC > LLg(2) .and. LLFC<0D0) then
    LLg(2) = LLFC   ! note: no ageprior for FC
    LL(2) = addALR(LLg(2), ALR(2))
  endif
endif  

!~~~  Parent 3-k / B-A  ~~~~~~~  
LLP = missing
LLPS = missing
ALRp = missing 
LLU = missing 
if (all(Parent(A,:)/=B) .and. focal/=7 .and. .not. &
  (focal==1 .and. .not. any(HasParent(:,2)))) then ! if B no parents, LL indistinguishable.
  ! NOT: .and. all(Parent(A,:)==0) --> inconsistency when called by CalcCandPar                             
  do x=1,2
    do i=1,2
      if (i==1 .and. x==k .and. hermaphrodites==0)  cycle  ! 'default'  
      if (AgeDiff(AB(i),AB(3-i)) <=0 .or. Sex(AB(3-i))==3-x) then
        LLP(i,x) = impossible  ! known, invalid agedif / sex
        cycle
      endif 
      call CalcAgeLR(AB(i),Sex(AB(i)), AB(3-i),Sex(AB(3-i)), x, 1, .TRUE., ALRp(i,x))  
      if (ALRp(i,x) == impossible) then
        LLP(i,x) = impossible
        cycle
      endif
      if (hermaphrodites/=0) then
        call PairPOX(AB(i), AB(3-i), x, focal, LLPS(i,x))    ! AB(3-i) result of selfing
        if (i==1 .and. (x==k .or. hermaphrodites==2))  cycle
      endif  
      if (Parent(AB(i),x) == 0) then
        call PairPO(AB(i), AB(3-i), x, 0, LLP(i,x))   ! WAS: FOCAL=0 ?  
      else if (Parent(AB(i),x) < 0) then
        call AddParent(AB(3-i), -Parent(AB(i),x), x, LLP(i,x))
        call CalcU(AB(3-i), Sex(AB(3-i)), Parent(AB(i),x), x, LLU)
        if (LLP(i,x) < 0)  LLP(i,x) = LLP(i,x) - LLU + LLg(7)
      endif
    enddo
  enddo   
  
  if (focal==1) then
    z=6   ! obv. not 3rd degree rel, but most convenient place to store w/o having to call CheckPair again
  else
    z=1
  endif
  LLg(z) = MaxLL((/LLg(z), LLP(:,1), LLP(:,2), LLPS(:,1), LLPS(:,2)/))
  do x=1,2
    do i=1,2
      LLP(i,x) = addALR(LLP(i,x), ALRp(i,x))
      LLPS(i,x) = addALR(LLPS(i,x), ALRp(i,x))
    enddo
  enddo
  LL(z) = MaxLL((/LL(z), LLP(:,1), LLP(:,2), LLPS(:,1), LLPS(:,2)/))
endif

!~~~  Sibs 3-k  ~~~~~~~
LLS = missing
ALRs = missing
if ((ANY(Parent(A,:)/=0) .or. ANY(Parent(B,:)/=0)) &  ! else: could be both. 
  .and. .not. (Parent(A,k) == Parent(B,k))) then                                              
  call CalcAgeLR(A,Sex(A), B, Sex(B), 3-k, 3, .TRUE., ALRs)
  if (ALRs /= impossible) then
    call PairHalfSib(A, B, 3-k, LLS)
    if (LLS<0D0 .and. (LLS - LLg(3) > TA .or. LLg(3)>0D0)) then
      if (focal == 3) then
        LLg(3) = MaybeOtherParent
        LL(3) = MaybeOtherParent
      else
        LLg(3) = LLS
        LL(3) = addALR(LLS, ALRs)
      endif
    endif
  endif
endif

!~~~  GP 3-k  ~~~~~~~
LLGR = missing 
ALRgr = missing    
if (((focal==3 .and. (MaxLL(LLg(2:3))>LLg(4) .or. LLg(4)>0D0)) .or. &
   (focal==4 .and. (LLg(4) > MaxLL(LLg(2:3))) .or. MaxLL(LLg(2:3))>0D0) .or. &
   (focal==1 .and. (LLg(1) > LLg(4) .or. LLg(4)>0))) .and. &
  Sex(B)<3 .and. .not. (Parent(A,3-k)==Parent(B,3-k) .and. Parent(A,3-k)/=0)) then
  cgp = getPar(Parent(A,3-k), 3-k)
  if (cgp(Sex(B)) == 0) then
    call CalcAgeLR(A,Sex(A),B,sex(B),3-k,4,.TRUE.,ALRgr(1))
    if (ALRgr(1)/=impossible) then
      call PairGP(A, B, 3-k, focal, LLGR(1))
    endif   
    if (focal==3 .or. focal==1) then
      LLg(4) = MaxLL((/ LLg(4), LLGR(1) /)) 
      LL(4) = MaxLL((/LL(4), addALR(LLGR(1), ALRgr(1))/))
    else if (focal==4 .and. hermaphrodites/=2) then
      if (LLGR(1)<0D0 .and. (LLg(4) - LLGR(1)) < TA) then
        LLg(4) = MaybeOtherParent
        LL(4) = MaybeOtherParent
      endif
    endif
  endif
endif

!~~~  various double rel  ~~~~~~~
LLPA = missing 
LLZ = missing
LLC = missing   
LLHH = missing 
LLX = missing 
if (focal==2 .or. focal==3) then   
  fclsib = .TRUE.
else
  fclsib = .FALSE.
endif  

if (Complx==2 .and. .not. (SelfedIndiv(A) .or. SelfedIndiv(B))) then     
  if (AgeDiff(A,B)>0 .and. LLg(1)/=impossible) then
    do x=1,3
      if (ALRAU(1,x)==impossible)  cycle
      call PairPOHA(A, B, k, x, LLPA(x))
    enddo
    if (any(LLPA > LLg(1) .and. LLPA < 0d0)) then
      LLg(1) = MaxLL(LLPA)
      do x=1,3
        LLPA(x) = addALR( addALR(LLPA(x), ALR(1)), ALRAU(1,x) )
      enddo
      LL(1) = MaxLL(LLPA)
    endif
  endif
  
  if (.not. any(HasParent)) then
    if (AgeDiff(A,B)>0) call PairHSPO(A,B,LLX(1))   ! HS via k, PO via 3-k  !  .and. Sex(B)/=k
    if (AgeDiff(B,A)>0) call PairHSPO(B,A,LLX(2))
    if (any(LLX(1:2)<0D0)) then
      if (fclsib) then
        LLg(3) = MaxLL((/LLg(3), LLX(1:2)/))
      else
        LLg(1) = MaxLL((/LLg(1), LLX(1:2)/))
      endif
      if (ALRp(1,3-k)==missing)  call CalcAgeLR(A,sex(A),B,sex(B),3-k,1,.TRUE.,ALRp(1,3-k))
      if (ALRp(2,3-k)==missing)  call CalcAgeLR(B,sex(B),A,sex(A),3-k,1,.TRUE.,ALRp(2,3-k))
      do x=1,2
        LLX(x) = addALR( addALR(LLX(x), ALR(3)), ALRp(x,3-k))
      enddo
      if ((LLg(2) - MaxLL(LLX(1:2))) < TA)  LLg(2) = MaybeOtherParent
    endif
  endif

  if (LLg(2)<0 .and. LLg(2) - LL(7) > TA) then  ! check if inbred FS (& boost LL) 
    call PairFSHA(A, B, k, LLZ(1))
    call PairFSHA(A, B, 3-k, LLZ(2))
    if (any(LLZ(1:2) > LLg(2) .and. LLZ(1:2) < 0d0)) then
      LLg(2) = MaxLL(LLZ(1:2))
      LL(2) = addALR(LLg(2), ALR(2))
    endif
  endif
    
  if (fclsib .and. (ABS(MaxLL(LLg) - LLg(2))<0.01 .or. LLg(3) - MaxLL(LLg(5:6)) < TA) .and. &
  (Parent(A,3-k)<=0 .or. Parent(B,3-k)<=0) .and. ALR(3)/=impossible .and. ALR(6)/=impossible) then
    do i=1,2
      if (i==2 .and. .not. any(HasParent)) then  ! symmetrical
        LLHH(:,2) = LLHH(:,1)
      else
        call PairHSHA(A, B, i, LLHH(1,i), (/.FALSE.,.FALSE./))    
        call PairHSHA(B, A, i, LLHH(2,i), (/.FALSE.,.FALSE./))
      endif
    enddo
    if (fclsib .and. any(LLHH(:,k)<0d0) .and. (MaxLL(LLHH(:,k)) - MaxLL(LLHH(:,3-k))) > TA .and. &
     MaxLL(LLHH(:,k)) > LLg(3)) then  ! .and. fclsib ?
      LLg(3) = MaxLL(LLHH(:,k))
      LL(3) = addALR( addALR(LLg(3), ALR(3)), ALR(6))   ! TODO: fix ageprior
    endif 
    if (any(LLHH(:,3-k)<0d0)) then  ! focal/=2 .and. 
      LLg(6) = MaxLL((/LLg(6), LLHH(:,3-k)/))
      LL(6) = addALR( addALR(LLg(6), ALR(3)), ALR(6))
    endif
  endif

  if (LLg(3)/=impossible .and. focal/=1 .and. LLGR(1)/=impossible .and. &
   (.not. all(HasParent(3-k,:)))) then 
   if (ALRgr(1)==Missing)  call CalcAgeLR(A,Sex(A),B,sex(B),3-k,4,.TRUE.,ALRgr(1)) 
   if (ALRgr(1)/=impossible)  call PairHSGP(A, B,k, LLX(3))  ! HS via k, GP via 3-k
   call CalcAgeLR(B,Sex(B),A,Sex(A),3-k,4,.TRUE.,ALRgr(2))
   if (ALRgr(2)/=impossible)  call PairHSGP(B, A,k, LLX(4))
   if (any(LLX(3:4)<0D0)) then
      if ((LLg(2) -MaxLL(LLX(3:4)))<TA  .and. fclsib) then
        LLg(2) = MaybeOtherParent
      endif
      if (MaxLL(LLX(3:4)) > LLg(3)) then       
        LLg(3) = MaxLL(LLX(3:4))
        LL(3) = addALR( addALR(LLg(3), ALR(3)), ALRgr(1)) 
      endif
    endif
  endif
  
  if ((LLg(2)/=MaybeOtherParent .and. .not. any(HasParent)) .or. &
   (focal==1 .and. LLg(1) - MaxLL(LLg) > -TA .and. ANY(ALRAU /= impossible)))  then
    call pairFAHA(A, B, .FALSE., LLZ(5))
    call pairFAHA(B, A, .FALSE., LLZ(6))
    if (ANY(LLZ(5:6) < 0D0)) then
      if ((LLg(2) - MaxLL(LLZ(5:6))) < TA  .and. fclsib) then
        LLg(2) = MaybeOtherParent
      endif
      if (MaxLL(LLZ(5:6)) > LLg(5)) then
        LLg(5) = MaxLL(LLZ(5:6))
        LL(5) = addALR(LLg(5), ALR(6))
      endif
    endif
  endif
    
  if (LLg(2)/=MaybeOtherParent .and. fclsib) then  ! check if GG in any way. can't be FS and GP
    if(LLGR(1)==missing)  call PairGP(A, B, 3-k, focal, LLGR(1))
    if (AgeDiff(A,B)==missing) then
      call PairGP(B, A, k, focal, LLGR(2))
      call PairGP(B, A, 3-k, focal, LLGR(3))
    endif
    if (MaxLL(LLGR)<0D0 .and. (LLg(4) - MaxLL(LLGR)) <TA) then
      LLg(4) = MaxLL(LLGR)
      if (ALRgr(1)==missing) then
        call CalcAgeLR(A,Sex(A),B,sex(B),3-k,4,.TRUE.,ALRgr(1))
        LLGR(1) = addALR(LLGR(1), ALRgr(1))
      endif
      do x=1,2
        call CalcAgeLR(B,Sex(B),A,Sex(A),x,4,.TRUE.,ALRgr(x+1))
         LLGR(x+1) = addALR(LLGR(x+1), ALRgr(x+1))
      enddo
      LL(4) = MaxLL(LLGR) 
    endif
  endif

  if (LLg(2)==MaybeOtherParent)  LL(2) = MaybeOtherParent
endif

if (hermaphrodites/=0) then
  call PairFSSelfed(A, B, LLZ(2))
  if (LLZ(2) > LLg(2) .and. LLZ(2) < 0d0) then
    LLg(2) = LLZ(2)
    LL(2) = addALR(LLZ(2), ALR(2))
  endif
endif

if (Complx==2 .and. (MaxLL(LL)>=LL(3) .or. MaxLL(LL)==LL(2))) then  
  call pairHSHAI(A, B, k, LLZ(3)) ! HS + inbr HA
  call pairHSHAI(B, A, k, LLZ(4))
  if (MaxLL(LLZ(3:4)) < 0D0 .and. MaxLL(LLZ(3:4)) > LLg(3)) then
    LLg(3) = MaxLL(LLZ(3:4))
    LL(3) = addALR(MaxLL(LLZ(3:4)), ALR(3))
  endif
    
  call PairHSCC(A,B,k, LLZ(7))
  if (LLZ(7) < 0 .and. LLZ(7) > LLg(3)) then
    LLg(3) = LLZ(7)
    LL(3) = addALR(LLZ(7), ALR(3))
  endif  
endif 

if (complx/=1 .and. LLg(2)<0 .and. LLg(2)>LLg(7)) then
  call pairFSFC(A,B, LLZ(8))
endif
LLFA = missing
ALRf = missing            
if (focal==1) then  ! check if FS of other parent  
  call CalcAgeLR(A,3-k, B,0, 3-k, 5, .TRUE., ALRf(1))
  if (ALRf(1) /= impossible) then  
    call PairUA(A, B, 3-k, 3, LLFA(1))
  endif
  ! check if FA of parent
  cgp = getPar(Parent(A,3-k), 3-k)
  if (Parent(A,3-k)>0 .and. any(cgp==0)) then
    if (cgp(1)==0) then
      z = 1
    else 
      z = 2
    endif
    call CalcAgeLR(Parent(A,3-k),3-k, B,0, z, 5, .TRUE., ALRf(2))
    if (ALRf(2) /= impossible) then  
      call PairUA(Parent(A,3-k), B, z, 3, LLFA(2))
      if (LLFA(2) < 0) then
        call CalcU(Parent(A,3-k),3-k, B,sex(B), LLU)
        LLFA(2) = LLFA(2) - LLU + LLg(7)
      endif
    endif
  endif
  
  if (LL(5)>0 .or. ((LL(5) - addALR(LLFA(1), ALRf(1))) < TA .and. LLFA(1)<0) .or. &
   ((LL(5) - addALR(LLFA(2), ALRf(2))) < TA .and. LLFA(2)<0)) then
    LLg(5) = MaxLL(LLFA)
    LL(5) = MaxLL( (/ addALR(LLFA(1), ALRf(1)), addALR(LLFA(2), ALRf(2)) /) )
  endif
endif

! check if add to opp. sibship  
if (fclsib .and. ANY(LLg(2:3)<0d0) .and. (MaxLL(LLg(2:3)) - MaxLL(LLg((/1,4,5,6,7/))))>TA .and. &
  (.not. all(HasParent(3-k,:))) .and. ANY(Parent(AB,3-k)<0)) then  
  if (Parent(A,3-k) < 0 .and. .not. HasParent(3-k,2)) then  
    call CheckAdd(B, -Parent(A,3-k), 3-k, 3, LLC(:,1), LLC(:,2))  !! DANGER
  else if (Parent(B,3-k)<0 .and. .not. HasParent(3-k,1)) then
    call CheckAdd(A, -Parent(B,3-k), 3-k, 3, LLC(:,1), LLC(:,2))  !! DANGER 
  endif
  if (LLC(2,1)<0D0 .and. (LLC(2,1) - MaxLL(LLC((/1,3,4,5,6,7/),1))) < TA) then     ! TODO: use as check for GP ?
    LLg(2) = MaybeOtherParent
    LL(2) = MaybeOtherParent
  endif
  if (LLC(3,1)<0D0 .and. (MaxLL(LLC((/1,2,4,5,6,7/),1)) - LLC(3,1)) > TA) then  
    LLg(3) = MaybeOtherParent
    LL(3) = MaybeOtherParent
  endif
endif  

end subroutine CheckPair

! #####################################################################

subroutine PairSelf(A, B, LL)  ! A==B; currently only called w/o parents
use Global
implicit none

integer, intent(IN) :: A, B
double precision, intent(OUT) :: LL
integer :: l
double precision :: PrL(nSnp)

PrL = 0D0
forall (l=1:nSnp)  PrL(l) = LOG10(SUM( LindX(:,l,A) * OcA(:,Genos(l,B)) ))
LL = SUM(PrL)

end subroutine PairSelf

! #####################################################################

subroutine GetPOconfigs(A, B, k, focal, Maybe)  ! TODO: check for PairPOX
use Global
implicit none

integer, intent(IN) :: A, B, k, focal
logical, intent(OUT) :: Maybe(5)
integer :: GA(2), PAB, AncB(2,mxA)

Maybe = .FALSE.  ! 1: non-inbred, 2: B PO & GP; 3: B PO & HS, 
                 ! 4: Parent(A,3-k) ancestor of B, 5: B selfing
! B HS of Parent(A,3-k) : see PairPOHA                                      
Maybe(1) = .TRUE. 
if (Complx==0)  return 
GA = getPar(Parent(A,3-k), 3-k)   
if (GA(k)==B .or. (Complx==2 .and. GA(k) == 0)) then
  Maybe(2) = .TRUE.
endif

if ((Complx==2 .and. (Parent(A,3-k)==Parent(B,3-k) .or. Parent(A,3-k)==0 &
  .or. Parent(B,3-k)==0)) .or. (Parent(A,3-k)==Parent(B,3-k) .and. focal==7)) then
    if ((focal==1 .or. focal==7) .and. Parent(A,3-k)==0 .and. Parent(B,3-k)/=0) then
      Maybe(3) = .FALSE.
    else if (Parent(A,3-k)/=0) then
      if (ANY(GA == B)) then
      Maybe(3) = .FALSE.
    else
      Maybe(3) = .TRUE. 
    endif
  else
    Maybe(3) = .TRUE.
  endif
endif

if (Parent(A,3-k)==Parent(B,3-k) .and. Parent(A,3-k)/=0) then
  Maybe(1) = .FALSE.  ! becomes config (3)
else if (Parent(A,3-k) > 0) then
  if (any(GA == B)) then
    Maybe(1) = .FALSE.  ! becomes config (2).
  endif
else if (Parent(A,3-k) < 0) then
  if (any (GpID(:, -Parent(A,3-k), 3-k) == B)) then
    Maybe(1) = .FALSE.
  endif
endif

if (Parent(A,3-k)/=0) then
  AncB = 0
  call getAncest(B, k, AncB)   
  if (ANY(AncB(3-k, 3:4) == Parent(A,3-k))) then
    Maybe(1) = .FALSE. 
    Maybe(4) = .TRUE.  
  endif
endif

if (hermaphrodites>0) then
  if (Parent(A,3-k) == B) then
    Maybe(1:4) = .FALSE.
    Maybe(5) = .TRUE.
  else if (focal/=1 .and. focal/=7 .and. Parent(A,3-k)<=0 .and. Sex(B)==4) then
    Maybe(5) = .TRUE.  
  else
    Maybe(5) = .FALSE.     ! else single/double parent same LL   
  endif
  if (SelfedIndiv(B)) then
    Maybe(3) = .FALSE.
  endif
endif

PAB = Parent(A,3-k)
if(Maybe(3) .and. Parent(A,3-k)==0) then  ! B PO & HS
  PAB = Parent(B,3-k)
  if (PAB < 0) then
    if (any(parent(SibID(1:ns(-PAB,3-k),-PAB,3-k),k) == A)) then
      Maybe(3) = .FALSE.  ! not implemented
    endif
  endif
endif       

end subroutine GetPOconfigs

! #####################################################################

subroutine PairPO(A, B, k, focal, LL)
use Global
use CalcLik           
implicit none

integer, intent(IN) :: A, B, k, focal
double precision, intent(OUT) :: LL
integer :: l, x, y,m, curPar(2), AncB(2,mxA), PAB, GA(2), Gtmp(2), Btmp(2), Bx
double precision :: PrL(nSnp,5), PrX(3,3,5), PrPA(3), PrB(3),PrPB(3,2),&
  LLtmp(5), PrPAB(3), PrPAX(3), PrG(3), LLX(5)
logical :: Maybe(5), ParOK                 

LL = missing
if (Sex(B)<3 .and. Sex(B)/=k) then
  LL = impossible
else if(Parent(A,k)>0) then  ! allow dummy to be replaced (need for AddFS)
  if (Parent(A,k)==B) then
    LL = AlreadyAss
  else if (focal==1) then    ! else do consider (in case current parent wrong)
    LL = impossible
  endif
else if (Parent(A,k)<0) then
  if (any(SibID(:,-parent(A,k),k) == B)) then
    LL = impossible
  endif
endif
if (LL/=missing) return

ParOK = .TRUE.
call ChkValidPar(A, Sex(A), B, k, ParOK)
if (.not. ParOK .and. focal/=7) then 
  LL = impossible
  return
endif

Maybe = .TRUE.
call GetPOconfigs(A, B, k, focal, Maybe)
! 1: non-inbred, 2: B PO & GP; 3: B PO & HS, 
! 4: Parent(A,3-k) ancestor of B, 5: B selfing
GA = getPar(Parent(A,3-k), 3-k)  
PAB = Parent(A,3-k)
if(Maybe(3) .and. Parent(A,3-k)==0) then  ! B PO & HS
  PAB = Parent(B,3-k)
endif   


LLtmp = missing     
LLX = missing     
if ((Parent(A,3-k)<0 .or. (Maybe(4) .and. .not. Maybe(1))) .and. ParOK) then                      
  curPar = Parent(A, :)
  Maybe(1) = .TRUE.
  if (Parent(B,k)<0) then
    Bx = Parent(B,k)
  else
    Bx = B
  endif
  
  call setParTmp(A, Sex(A), 0, k)     ! B vs none.    
  call CalcU(A,3-k, B,k, LLX(1))
  call CalcU(Parent(A,3-k),3-k, Bx,k, LLX(2))

  call setParTmp(A, Sex(A), B, k)  
  call CalcU(Parent(A,3-k),3-k, Bx,k, LLX(3))
  LLtmp(1) = LLX(1) + (LLX(3) - LLX(2)) 
  
  if (Maybe(2) .and. .not. ANY(GA == B) .and. parent(A,3-k)/=0) then    ! B PO & GP ? (if GA==B, covered by LLtmp(1))
    AncB = 0                                                                                
    call GetAncest(B, Sex(B), AncB)
    call ChkValidPar(Parent(A,3-k), 3-k, B, k, ParOK)
    if (ParOK .and. .not. ANY(AncB(3-k,:)==Parent(A,3-k))) then
      Gtmp = getPar(Parent(A,3-k), 3-k)
      call setParTmp(Parent(A,3-k), 3-k, B, k)
      call CalcU(Parent(A,3-k),3-k, Bx,k, LLX(4))
      LLtmp(2) = LLX(1) + (LLX(4) - LLX(2))
      call setParTmp(Parent(A,3-k), 3-k, Gtmp(k), k)
    endif
  endif
  
  ! 3-PO & HS 
  if (Maybe(3) .and. Parent(B,3-k)==0 .and. parent(A,3-k)/=0) then
    call ChkValidPar(B,k, Parent(A,3-k), 3-k, ParOK)
    if (ParOK) then
      Btmp = getPar(B,sex(B))
      call setParTmp(B,sex(B), Parent(A,3-k),3-k)
      call CalcU(Parent(A,3-k),3-k, Bx,k, LLX(5))
      LLtmp(3) = LLX(1) + (LLX(5) - LLX(2))
      call setParTmp(B,sex(B), Btmp(3-k),3-k)   
    endif
  endif
  
  ! TODO: 5-selfing
  
  call setParTmp(A, Sex(A), curPar(k), k)  ! restore

else
  PrL = 0D0                        
  do l=1,nSnp    
    if (Complx==0 .and. Mate(B)/=0) then
      call ParProb(l, Mate(B), 3-k, A, 0, PrPA)
    else
      call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPA)
    endif
    PrB = LindX(:,l,B)  ! unscaled
    if (Maybe(2)) then
      call ParProb(l, Parent(A,3-k), 3-k, A, -4, PrPAX)
      if (GA(3-k)==B) then  ! hermaphrodite  
        m = k
      else
        m = 3-k
      endif
      if (Parent(A,3-k)>0) then
        call ParProb(l, GA(m), m, Parent(A,3-k), 0, PrG)  
      else 
        call ParProb(l, GA(m), m, 0, 0, PrG)       
      endif
    endif
    if (Maybe(3)) then
      do m=1,2  
        call ParProb(l, Parent(B,m), m, B, 0, PrPB(:,m))
      enddo
      call ParProb(l, PAB, 3-k, A, B, PrPAB)
    endif

    do x=1,3  ! B
      do y=1,3  ! parent(A,3-k)
        PrX(x,y,:) = OKA2P(Genos(l, A), x, y)
        if (Maybe(1))  PrX(x,y,1) = PrX(x,y,1) * PrB(x) * PrPA(y)
        if (Maybe(2)) then  ! B PO & GP
          PrX(x,y,2) = PrX(x,y,2) * PrB(x)* PrPAX(y) * SUM(AKA2P(y,x,:) * PrG)
        endif
        if (Maybe(3)) then  ! B PO & HS
          PrX(x,y,3) = PrX(x,y,3) * PrPAB(y) * SUM(AKA2P(x,y,:) * PrPB(:,k)) * OcA(x,Genos(l,B))
        endif
        if (Maybe(5)) then  ! B selfing
          if (x/=y)  PrX(x,y,5) = 0D0
          if (x==y)  PrX(x,y,5) = PrX(x,y,5) * PrB(x)
        endif
      enddo
    enddo
    do x=1,5
      if (Maybe(x))  PrL(l,x) = LOG10(SUM(PrX(:,:,x)))
    enddo
  enddo
  LLtmp = SUM(PrL, DIM=1)
  if (Parent(A,3-k) > 0  .and. Maybe(2)) then
    LLtmp(2) = LLtmp(2) - Lind(Parent(A,3-k))
  endif
endif

do x=1,5
  if (.not. Maybe(x) .or. LLtmp(x)>=0)  LLtmp(x) = impossible
enddo

LL = MaxLL(LLtmp)

end subroutine PairPO

! #####################################################################

subroutine PairPOX(A, B, k, focal, LL)    ! B parent of A; B result of selfing
use Global
implicit none

integer, intent(IN) :: A, B, k, focal
double precision, intent(OUT) :: LL
integer :: l, x, y, z, GA(2), PAB, m
double precision :: PrL(nSnp,5), PrPB(3), PrPA(3), PrXYZ(3,3,3,5), LLtmp(5), &
  PrG(3), PrPAX(3)
logical :: ParOK, Maybe(5)

LL = missing
if (hermaphrodites==0) then
  LL = impossible
else if (Parent(A,k)/=0) then
  LL = impossible
else if (Parent(B,1) >0) then
  if (Parent(B,1)/=Parent(B,2) .and. Parent(B,2)/=0) then
    LL = impossible
  else if (Sex(Parent(B,1))/=4) then
    LL = impossible
  endif
endif
if (LL == impossible)  return

ParOK = .TRUE.
call ChkValidPar(A, Sex(A), B, k, ParOK)
if (.not. ParOK .and. focal/=7) then
  LL = impossible
  return
endif

Maybe = .TRUE.
call GetPOconfigs(A, B, k, focal, Maybe)
Maybe(4) = .FALSE.   ! not implemented with B selfed
GA = getPar(Parent(A,3-k), 3-k)  
PAB = Parent(A,3-k)
if(Maybe(3) .and. Parent(A,3-k)==0) then  ! B PO & HS
  PAB = Parent(B,3-k)
endif

PrL = 0D0 
PrPAX = missing
do l=1,nSnp
  call ParProb(l, Parent(B,1), 1, B, 0, PrPB)
  call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPA)
  if (Maybe(2)) then
    call ParProb(l, Parent(A,3-k), 3-k, A, -4, PrPAX)
    if (GA(3-k)==B) then  ! hermaphrodite    
      m = k
    else
      m = 3-k
    endif
    if (Parent(A,3-k)>0) then
      call ParProb(l, GA(m), m, Parent(A,3-k), 0, PrG)  
    else 
      call ParProb(l, GA(m), m, 0, 0, PrG)    
    endif
  endif
  
  do x=1,3  ! B
    do y=1,3  ! parent A 3-k
      do z=1,3  ! parent B
        PrXYZ(x,y,z,:) = OKA2P(Genos(l, A), x, y) * OcA(x,Genos(l,B)) * AKA2P(x,z,z) * PrPB(z)
        if (Maybe(1))  PrXYZ(x,y,z,1) = PrXYZ(x,y,z,1) * PrPA(y)
        if (Maybe(2))  PrXYZ(x,y,z,2) = PrXYZ(x,y,z,2) * SUM(PrPAX(y) * AKA2P(y,x,:) * PrG)
        if (Maybe(3) .and. y/=z)  PrXYZ(x,y,z,3) = 0D0   ! else as is
        if (Maybe(5) .and. x/=y)  PrXYZ(x,y,z,5) = 0D0   ! else as is
      enddo
    enddo
  enddo
  do x=1,5
    if (Maybe(x))  PrL(l,x) = LOG10(SUM(PrXYZ(:,:,:,x)))
  enddo 
enddo
LLtmp = SUM(PrL, DIM=1)
if (Parent(A,3-k) > 0 .and. Maybe(2)) then
  LLtmp(2) = LLtmp(2) - Lind(Parent(A,3-k))
endif
do x=1,5
  if (.not. Maybe(x) .or. LLtmp(x)>=0 .or. LLtmp(x) < -HUGE(0D0))  LLtmp(x) = impossible
enddo

LL = MaxLL(LLtmp)

end subroutine PairPOX

! #####################################################################

subroutine PairFullSib(A, B, LL)
use Global           
use CalcLik           
implicit none

integer, intent(IN) :: A, B
double precision, intent(OUT) :: LL
integer :: x, y, l,k, Par(2), ix, i, AlreadyHS, Ei, DoQuick, AB(2), curPar(2,2) 
double precision :: PrL(nSnp), PrXY(3,3), Px(3,2), LUX(2), LLtmp, &
  dx(maxSibSize), LRQ(2), PrE(3), PrXYb(3,3,2), PrFS(3,3)
logical :: AncOK(2)

LL = missing
Par = 0  ! joined parents of A & B
if (Parent(A,1)==Parent(B,1) .and. Parent(A,1)/=0 .and. &
 Parent(A,2)==Parent(B,2) .and. Parent(A,2)/=0) then ! already FS
  LL = AlreadyAss
  return
else 
  do k=1,2
    if (Parent(A,k) == B .or. Parent(B,k) == A) then
      LL = impossible
      return
    else if (Parent(A,k)/=Parent(B,k) .and. .not. (Parent(A,k)==0 .or. &
      Parent(B,k)==0)) then
      LL = impossible
      return
    else if (Parent(A,k)/=0) then
      Par(k) = Parent(A,k)
    else
      Par(k) = Parent(B,k)
    endif
    if (Complx==0 .and. Par(k)/=0) then
      Par(3-k) = getMate(Par(k),k) 
    endif
  enddo
endif  

if (Complx==0) then
  do k=1,2  
    if (Par(k)==A .or. Par(k)==B) then
      LL = impossible
      return
    endif
  enddo
endif

LRQ = missing
call CalcP2(A, Sex(A), Par(1), Par(2), 1, LRQ(1))
call CalcP2(B, Sex(B), Par(1), Par(2), 1, LRQ(2))
if (any(LRQ == impossible)) then 
  LL = impossible
  return
endif

AncOK = .TRUE.
call ChkAncest(A,0,B,0, AncOK(1))
call ChkAncest(B,0,A,0, AncOK(2))
if (any(.not. AncOK)) then  ! cannot be both ancestors of eachother & full siblings
  LL = impossible
  return
endif

AlreadyHS = 0
do k=1,2
  if (Parent(A,k) == Parent(B,k) .and. Parent(A,k)<0) then
    AlreadyHS = k
  endif
enddo

PrL = 0D0 
LUX = 0D0
LLtmp = missing
dx = missing           
curPar=0              
if ((Par(1) < 0 .or. Par(2)<0) .and. (AlreadyHS==0 .or. (Par(1)/=0 .and. Par(2)/=0))) then
  AB = (/A, B/)
  forall (x=1:2)  curpar(:,x) = Parent(AB(x),:)
  do k=1,2 
    if (Par(k) >= 0)  cycle
    ! set parents 3-k to 0 so that LUX(1) is consistent with addFS
    do x=1,2
      call setParTmp(AB(x),1,0,3-k)
      call CalcU(A, 1, B, 1, LUX(x))
      call setParTmp(AB(x),1,curPar(3-k,x),3-k) 
    enddo
    do x=1,2
      if (Parent(AB(x),k)==0 .and. Parent(AB(3-x),k)==Par(k)) then
        call addFS(AB(x), -Par(k), k, 0, k, LLtmp, ix, dx)   ! call AddFS
        do i=1, nS(-Par(k),k)
          if (SibID(i,-Par(k),k) == AB(3-x)) then
            if (dx(i) < impossible) then
              LL = LUX(x) + dx(i)
            else
              LL = dx(i)
            endif 
            return
          endif
        enddo
      endif
    enddo      
  enddo   
  
else if (AlreadyHS/=0 .and. (Par(1) < 0 .or. Par(2)<0)) then
  k = AlreadyHS
  DoQuick = 1
  call ChkDoQuick(-Parent(A,k),k, DoQuick) 
  if (DoQuick>1 .or. DoQuick==-1 .or. DoQuick==-3) then
    LL = NotImplemented
    Return
  endif
   
  do l=1, nSnp   
    call ParProb(l, Par(k), k, -1, 0, Px(:,k))
    call ParProb(l, Par(3-k), 3-k, 0, 0, Px(:,3-k))   ! other parent ==0   
    do x=1,3  ! already shared parent
      do y=1,3  ! other parent
        PrXYb(x,y,:) = Px(x,k) * Px(y,3-k)
        do i=1, ns(-Par(k),k)
          Ei = SibID(i,-Par(k), k)
          if (nFS(Ei)==0 .or. Ei==A .or. Ei==B)  cycle
          PrFS = FSLik(l,Ei)
          call ParProb(l, Parent(Ei, 3-k), 3-k, Ei,-1, PrE)              
          PrE = PrE * PrFS(:,x)
          if (.not. ALL(PrE==1D0))  PrXYb(x,y,:) = PrXYb(x,y,:) * SUM(PrE)
        enddo
        PrXYb(x,y,2) = PrXYb(x,y,2) * OKA2P(Genos(l,A), x, y) * OKA2P(Genos(l,B), x, y)
      enddo
    enddo
    PrL(l) = LOG10(SUM(PrXYb(:,:,2))) - LOG10(SUM(PrXYb(:,:,1)))
  enddo
  LL = SUM(PrL) 

else  
  do l=1, nSnp
    do k=1,2
      if (Parent(A,k)==Parent(B,k)) then  
        call ParProb(l, Par(k), k, A, B, Px(:,k))
      else if (Parent(A,k)==Par(k)) then
        call ParProb(l, Par(k), k, A, 0, Px(:,k))
      else if (Parent(B,k)==Par(k)) then
        call ParProb(l, Par(k), k, B, 0, Px(:,k))
      else
        call ParProb(l, Par(k), k, 0, 0, Px(:,k))
      endif       
    enddo 
  
    do x=1,3
      do y=1,3
        PrXY(x,y) = Px(x,1) * Px(y,2) * OKA2P(Genos(l,A), x, y) * OKA2P(Genos(l,B), x, y)
      enddo
    enddo
    PrL(l) = LOG10(SUM(PrXY))
  enddo
  LL = SUM(PrL) 
endif

end subroutine PairFullSib

! #####################################################################

subroutine PairHalfSib(A, B, k, LL)
use Global
implicit none

integer, intent(IN) :: A, B, k
double precision, intent(OUT) :: LL
integer :: x,y, l, Par, Inbr, AB(2), i, PP, ZZ
double precision :: PrL(nSnp), PrX(3), PrPx(3,2), PrXY(3,3), LLtmp(2), LLU
logical :: ParOK

LL = missing
Par = 0  ! parent K
if (Parent(A,k)/=0) then
  if (Parent(A,k)/=Parent(B,k) .and. Parent(B,k)/=0) then
    LL = impossible ! mismatch
  else if (Parent(A,k)==Parent(B,k)) then
    LL = AlreadyAss ! already captured under H0
  else
    Par = Parent(A,k)
    call ChkValidPar(B,Sex(B), Par,k, ParOK)  ! includes age & ancestor check
    if (.not. ParOK)  LL = impossible
  endif
else if (Parent(B,k)/=0) then
  Par = Parent(B,k)
  call ChkValidPar(A, Sex(A), Par,k, ParOK)                 
  if (.not. ParOK)  LL = impossible
endif
if (LL/=missing) return

if (Parent(A,3-k)==Parent(B,3-k) .and. Parent(A,3-k)/=0) then
  LL = impossible  ! would become full sibs
  return
endif

AB = (/ A, B /)
LLtmp = missing
LLU = missing
if (Par < 0 .and. (all(Parent(A,:)>=0) .or. all(Parent(B,:)>=0))) then
  do i=1,2
    if (Parent(AB(i),k) == Par) then
      call addSib(AB(3-i), -Par, k, LLtmp(2))
      if (LLtmp(2) < 0) then
        call CalcU(AB(3-i), 3, Par, k, LLtmp(1))        
        call CalcU(AB(1),3, AB(2),3, LLU)
        LL = LLU + (LLtmp(2) - LLtmp(1))
      else
        LL = LLtmp(2)
      endif
      return
    endif
  enddo

else if (Par < 0 .and. any(Parent(AB,3-k) < 0)) then   
  PP = 0
  ZZ = Par
  do i=1,2
    if (Parent(AB(i),3-k) < 0) then
      PP = Parent(AB(i),3-k)
      if (Par >= 0)  ZZ = AB(3-i)
    endif
  enddo
  do i=1,2
    if (Parent(AB(i),k) == Par .and. parent(AB(3-i),k)==0) then
      call CalcU(ZZ, k, PP, 3-k, LLtmp(1))
      call SetParTmp(AB(3-i),3, Par, k)
      call CalcU(ZZ, k, PP, 3-k, LLtmp(2))
      call SetParTmp(AB(3-i),3, 0, k)
      call CalcCLL(-Par,k)
      call CalcCLL(-PP, 3-k)
      call CalcU(AB(1),3, AB(2),3, LLU)
      LL = LLU + (LLtmp(2) - LLtmp(1))
      return
    endif
  enddo
endif

Inbr = 0
do i=1,2
  if (Parent(AB(i),3-k) == AB(3-i)) Inbr = i
enddo

PrL = 0D0
do l=1,nSnp
  if (Par==Parent(A,k) .and. Par/=0) then
    call ParProb(l, Par, k, A, 0, PrX)
  else if (Par==Parent(B,k) .and. Par/=0) then
    call ParProb(l, Par, k, B, 0, PrX)
  else
    call ParProb(l, Par, k, 0, 0, PrX)    
  endif
  call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPx(:,1))
  call ParProb(l, Parent(B,3-k), 3-k, B, 0, PrPx(:,2))
  if (Inbr==0) then
    do x=1,3
      do y=1,3 
        PrXY(x,y) = PrX(x) * PrPX(y,1) * OKA2P(Genos(l,A),x,y) * &
          SUM(OKA2P(Genos(l,B),x,:) * PrPX(:,2))
      enddo
    enddo
    PrL(l) = LOG10(SUM(PrXY))
  else 
    do x=1,3
      do y=1,3
        PrXY(x,y)=PrX(x) * SUM(AKA2P(y, x, :) * PrPx(:,3-Inbr)) * &
          OKA2P(Genos(l,AB(Inbr)), x, y) * OcA(y,Genos(l, AB(3-Inbr)))
      enddo
    enddo
    PrL(l) = LOG10(SUM(PrXY))
  endif
enddo
LL = SUM(PrL)

end subroutine PairHalfSib

! #####################################################################

subroutine pairHSHA(A, B, k, LL, withFS)  !HS via k, & parent A is HS of B via 3-k
! returns max of hf 1: HSHA, 2: HSFA, 
! NOT 3: A inbred -- see pairHSHAI                              
use Global
use CalcLik  
implicit none

integer, intent(IN) :: A,B, k
logical, intent(IN) :: withFS(2)  ! A, B                  
double precision, intent(OUT) :: LL
integer :: l, x, y, z, PAB, Ai, Bj, exclFS(2), h, parA(2), parB(2), GA(2)
double precision :: PrL(nSnp,2), PrXYZ(3,3,3,3), PrY(3), PrZ(3), PrX(3), &
  PrTmp, PrW(3), LLtmp(2), PrA(3,3), PrB(3,3)
logical :: ParOK                

LL = missing
PAB = 0
if (Parent(A,3-k)/=0 .and. Parent(B,3-k)/=0) then
  call ChkValidPar(Parent(A,3-k),3-k, Parent(B,3-k),3-k, ParOK)
  if (.not. ParOK) then
    LL = Impossible
    return
  endif
endif

if (Parent(A,k)/=Parent(B,k)) then                                  
  if (Parent(A,k)/=0) then
    if(Parent(B,k)/=0) then
      LL = impossible
     else
      PAB = Parent(A,k)
      call chkvalidPar(B,sex(B), PAB,k, ParOK)
      if (.not. ParOK)  LL = impossible
    endif
  else
    PAB = Parent(B,k)
    call chkvalidPar(A,sex(A), PAB,k, ParOK)
    if (.not. ParOK)  LL = impossible
  endif
endif 
if (LL == impossible)  return   

if (Parent(A,3-k)<0 .or. (PAB < 0 .and. Parent(B,3-k)<0)) then  ! TODO: .or. 
  parA = getPar(A,sex(A))
  parB = getPar(B,sex(B))
  if (Parent(A,k)/=PAB)  call setParTmp(A,sex(A), PAB,k)
  if (Parent(B,k)/=PAB)  call setParTmp(B,sex(B), PAB,k)
  Ai = A
  Bj = B
  if (withFS(1) .and. Parent(A,3-k)<0)  Ai = Parent(A,3-k)
  if (withFS(2) .and. Parent(B,3-k)<0)  Bj = Parent(B,3-k)  
  call pairUA(Bj,Ai,3-k,3-k,LL)
  call setParTmp(A,sex(A), ParA(k),k)
  call setParTmp(B,sex(B), ParB(k),k)   
  
  return
endif  

if (PAB < 0 .and. hermaphrodites/=0) then
  if (DumClone(-PAB,k)/=0) then
    LL = NotImplemented
    Return
  endif
endif  

Ai = FSID(maxSibSize+1, A)
Bj = FSID(maxSibSize+1, B)
exclFS = 0
do x=1,2
  if (withFS(x))  exclFS(x) = -1
enddo

GA = getPar(Parent(A,3-k),3-k)

PrL = 0D0
PrXYZ = 0d0
do l=1, nSnp
  call ParProb(l, Parent(A,3-k), 3-k, A,-4, PrX)
  call ParProb(l, Parent(B,3-k), 3-k, B, exclFS(2), PrY)
  if (Parent(A,3-k)>0) then
    call ParProb(l, GA(k), k, Parent(A,3-k), 0, PrW)
  else
    call ParProb(l, GA(k), k, 0, 0, PrW)
  endif
  if (Parent(A,k)==Parent(B,k)) then
    call ParProb(l, PAB, k, A, B, PrZ)
  else if (Parent(A,k) == PAB) then
    call ParProb(l, Parent(A,k), k, A, exclFS(1), PrZ)
  else if (Parent(B,k) == PAB) then
    call ParProb(l, Parent(B,k), k, B, exclFS(2), PrZ)
  endif
  
  if (withFS(1)) then
    PrA = FSlik(l,Ai)
  else
    PrA = OKA2P(Genos(l,A),:,:)
  endif
  if (withFS(2)) then
    PrB = FSlik(l,Bj)
  else
    PrB = OKA2P(Genos(l,B),:,:)
  endif

  do x=1,3
    do y=1,3    
      do z=1,3
        PrTmp = PrY(y) * PrZ(z) * PrX(x) * PrA(x,z) * PrB(y,z)
        PrXYZ(x,y,z,1) = PrTmp * SUM(AKA2P(x,y,:) * PrW)  ! HSHA
        if (GA(k)==0 .or. GA(k)==PAB) then
          PrXYZ(x,y,z,2) = PrTmp * AKA2P(x,y,z)  ! HSFA
        endif
!        PrXYZ(x,y,z,3) = PrTmp *  AKAP(x,z,l)  ! A inbred  -- see HSHAI
      enddo
    enddo
  enddo
  do h=1,2
    PrL(l,h) = LOG10(SUM(PrXYZ(:,:,:,h)))
  enddo
enddo

LLtmp = SUM(PrL,DIM=1)
if (GA(k)/=0 .and. GA(k)/=PAB)   LLtmp(2) = Impossible
LL = MaxLL(LLtmp)
if (Parent(A,3-k)>0)  LL = LL - Lind(Parent(A,3-k))
if (LL < -HUGE(0D0))  LL = impossible

end subroutine pairHSHA

! #####################################################################

subroutine pairHSHAI(A, B, k, LL)  !HS via k, & A inbred
use Global
implicit none

integer, intent(IN) :: A, B, k
double precision, intent(OUT) :: LL
integer :: l, x, y, f,g
double precision :: PrL(nSnp,2,2), PrXY(3,3,2,2), PrPA(3), PrPAX(3), PrPB(3), LLU, LLtmp(2)
logical :: AncOK

LL = Missing
if (Parent(A,k)>0 .or. Parent(B,k)>0) then
  LL = NotImplemented
endif 
if (Parent(A,k) <0) then
  if (ns(-parent(A,k),k)/=1 .or. any(GpID(:,-parent(A,k),k)/=0))  LL = NotImplemented
endif
if (Parent(B,k) <0) then
  if (ns(-parent(B,k),k)/=1 .or. any(GpID(:,-parent(B,k),k)/=0))  LL = NotImplemented
endif
if (Parent(A, 3-k)==B) LL = impossible
if (LL/=Missing) return

call ChkAncest(A,0,B,0, AncOK)
if (.not. AncOK)  then
  LL = impossible
  return
endif
if (Parent(A,3-k)==Parent(B,3-k) .and. Parent(A,3-k)/=0) then
    LL = NotImplemented  ! likely picked up elsewhere
    return
endif

PrL = 0D0
do l=1, nSnp
  call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPA)
  call ParProb(l, Parent(A,3-k), 3-k, A, -4, PrPAX)  ! OcA if Parent(A,3-k)>0
  call ParProb(l, Parent(B,3-k), 3-k, B, 0, PrPB)
  do x=1,3
    do y=1,3    
      PrXY(x,y,1,:) = AKAP(x,y,l) * PrPAX(x) * AHWE(y,l) 
      PrXY(x,y,2,:) = AKAP(y,x,l) * PrPA(x)
      PrXY(x,y,:,:) = PrXY(x,y,:,:) * OKA2P(Genos(l,A), x,y)
      PrXY(x,y,:, 1) = PrXY(x,y,:, 1) * SUM(OKA2P(Genos(l,B), y,:)*PrPB)
      PrXY(x,y,:, 2) = PrXY(x,y,:, 2) * SUM(OKAP(Genos(l,B),:,l)*PrPB)
    enddo
  enddo
  do f=1,2
    do g=1,2
      PrL(l,f,g) = LOG10(SUM(PrXY(:,:,f,g)))
    enddo
  enddo
enddo

LLU = missing
LLtmp = missing  
call CalcU(A,k,B,k, LLU)
do f=1,2
  if (SUM(PrL(:,f,2)) > LLU) then
    LLtmp(f) = SUM(PrL(:,f,1)) - SUM(PrL(:,f,2)) + LLU
  else
    LLtmp(f) = SUM(PrL(:,f,1))
  endif
enddo
LL = MaxLL(LLtmp)

end subroutine pairHSHAI

! #####################################################################

subroutine pairFAHA(A, B, withFS, LL)  !B FA via k & HA via 3-k ; A inbred. 
use Global
implicit none

integer, intent(IN) :: A, B
logical, intent(IN) :: withFS
double precision, intent(OUT) :: LL
integer :: l, k, x, y, z, v, i, AA(maxSibSize), BB(maxSibSize), nA, nB
double precision :: PrL(nSnp, 4), PrXY(3,3,3,3, 4), PrPB(3,2), LLU, LLtmp(4)
logical :: ParOK

LL = missing
if (ANY(Parent(A,:)>0)) then  
  LL = NotImplemented   
  return
else
  do k=1,2
    if (Parent(A,k)<0) then
      if (any(GpID(:,-Parent(A,k),k) /=0))  LL = NotImplemented
      if (Parent(A,k) == Parent(B,k))  LL = Impossible
    endif
  enddo
  if (LL /= Missing)  return
endif

do k=1,2
  call ChkValidPar(A,Sex(A), Parent(B,k),k, ParOK)
  if (.not. ParOK) then
    LL = impossible
    return
  endif
enddo

AA = 0
BB = 0
if (withFS) then
  nA = nFS(FSID(maxSibSize+1, A))  
  AA(1:nA) = FSID(1:nA, FSID(maxSibSize+1, A))
  nB = nFS(FSID(maxSibSize+1, B))
  BB(1:nB) = FSID(1:nB, FSID(maxSibSize+1, B))
else
  nA = 1
  AA(1) = A
  nB = 1
  BB(1) = B
endif

PrL = 0D0
do l=1, nSnp
  do k=1,2
    if (withFS) then
      call ParProb(l, Parent(B,k), k, B, -1, PrPB(:,k))
    else
      call ParProb(l, Parent(B,k), k, B, 0, PrPB(:,k))
    endif
  enddo
  do x=1,3  ! Par A, FS of B
    do y=1,3   ! par B, double GP of A 
      do z=1,3  ! par B, GP of A
        do v=1,3  ! Par A, HS of B
          PrXY(x,y,z,v,1) = PrPB(y,1) * PrPB(z,2) * AKAP(v,y,l) * AKA2P(x,y,z)
          PrXY(x,y,z,v,2) = PrPB(y,2) * PrPB(z,1) * AKAP(v,y,l) * AKA2P(x,y,z)
          PrXY(x,y,z,v,3) = PrPB(y,1) * PrPB(z,2) * AHWE(v,l) * AHWE(x,l)  ! A, B unrelated
          PrXY(x,y,z,v,4) = PrPB(y,1) * PrPB(z,2) * SUM(AKAP(v,:,l) * AKAP(x,:,l) * AHWE(:,l))  
          ! A, B unrelated; A inbred
          do i=1, nA
            PrXY(x,y,z,v,:) = PrXY(x,y,z,v,:) * OKA2P(Genos(l,AA(i)), x, v)
          enddo
          do i=1, nB
            PrXY(x,y,z,v,:) = PrXY(x,y,z,v,:) * OKA2P(Genos(l,BB(i)), y, z)
          enddo
        enddo
      enddo
    enddo                         
  enddo
  do i=1,4
    PrL(l,i) = LOG10(SUM(PrXY(:,:,:,:,i)))
  enddo
enddo
LLtmp = SUM(PrL, DIM=1)

if (.not. withFS) then
  LL = MAXVAL(LLtmp(1:2))
else
  LLU = missing
  call CalcU(A, 0, B, 0, LLU)
  LL = MAXVAL(LLtmp(1:2)) - LLtmp(3) + LLU   !MAXVAL(LLtmp(3:4)) + LLU
endif

end subroutine pairFAHA

! #####################################################################

subroutine pairHSPO(A, B, LL)   ! HS via k, & PO via 3-k
use Global
implicit none

integer, intent(IN) :: A,B
double precision, intent(OUT) :: LL
integer :: l, x, y,m,k
double precision :: PrL(nSnp), PrXY(3,3), PrPB(3,2)
logical :: ParOK

if (ANY(Parent(A,:)/=0)) then   !  .or. ANY(Parent(B,:)/=0)
  LL = NotImplemented
  return   ! else not necessary.
endif  

call ChkValidPar(A, Sex(A), B,Sex(B), ParOK)
if (.not. ParOK) then
  LL = impossible
  return
endif

k=1
if (Sex(B)<3)  k=3-Sex(B)

PrL = 0D0
do l=1, nSnp
  do m=1,2
    call ParProb(l, Parent(B,m), m, B, 0, PrPB(:,m))
  enddo
  do x=1,3 
    do y=1,3    ! B
      PrXY(x,y) = OKA2P(Genos(l,A), x,y) * OcA(y,Genos(l,B)) * PrPB(x,k) * SUM(AKAP(y,x,:) * PrPB(:,3-k)) 
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXY))
enddo
LL = SUM(PrL)

end subroutine pairHSPO

! #####################################################################

subroutine PairPOHA(A, B, k, hf, LL)  ! B parent of A via k, and 'hf' sib of A's 3-k parent
use Global
implicit none

integer, intent(IN) :: A, B, k, hf
double precision, intent(OUT) :: LL
integer :: l, x, y,z,m, ParB(2), GA(2), GG(2)
double precision :: PrL(nSnp), PrXYZ(3,3,3), PrPAX(3), PrGG(3,2), PrPB(3), PrGA(3), LLX(3)

ParB = Parent(B,:)
GA = getPar(Parent(A,3-k), 3-k)
GG = 0
LL = missing

if (Parent(A,k)/=0 .or. (Parent(A,3-k) == ParB(3-k) .and. ParB(3-k)/=0)) then
  LL = impossible
  return
endif

do m=1,2
  if (m/=hf .and. hf/=3)  cycle
  if (GA(m) == B) then
    LL = impossible
  else if (ParB(m) == GA(m) .or. GA(m) == 0) then
    GG(m) = ParB(m)
  else if (ParB(m) == 0) then
    GG(m) = GA(m)
  else
    LL = impossible
  endif
  if (GG(m) < 0 .and. hermaphrodites/=0) then
    if (DumClone(-GG(m),m)/=0) then
      LL = NotImplemented
    endif
  endif
enddo
if (LL /= Missing)  return                          

if (hf<3) then
  if (GG(hf)<0) then
    call CalcU(A,sex(A),B,sex(B), LLX(1))
    call CalcU(A,sex(A),GG(hf),hf, LLX(2))
    call setParTmp(A,sex(A), B,k)  ! validity check done by checkpair
    call PairUA(A,GG(hf),3-k,hf,LLX(3))
    LL = LLX(3) - LLX(2) + LLX(1)
    call setParTmp(A,sex(A), 0,k)
    return
  endif
endif    

PrL = 0D0  
PrPB = missing
PrGA = missing
do l=1,nSnp       
  call ParProb(l, Parent(A,3-k), 3-k, A, -4, PrPAX)  ! offspring only, except A 
  do m=1,2
    if (m/=hf .and. hf/=3)  cycle
    if (Parent(A,3-k)>0) then
      call ParProb(l, GG(m), m, B, Parent(A,3-k), PrGG(:,m))
    else
      call ParProb(l, GG(m), m, B, 0, PrGG(:,m))  
    endif
  enddo
  if (hf < 3) then
    call ParProb(l, Parent(B,3-hf), 3-hf, B, 0, PrPB)
    if (Parent(A,3-k)>0) then
      call ParProb(l, GA(3-hf), 3-hf, Parent(A,3-k), 0, PrGA)
    else
      call ParProb(l, GA(3-hf), 3-hf, 0, 0, PrGA)
    endif
  endif
  
  PrXYZ = 0D0
  do x=1,3  ! B
    do y=1,3  ! parent(A,3-k)
      do z=1,3  ! double grandparent (dam if hf=3)
        if (hf < 3) then
          PrXYZ(x,y,z) = SUM(AKA2P(x,z,:) * PrPB) * &
            PrPAX(y) * SUM(AKA2P(y,z,:) * PrGA) * PrGG(z, hf)
        else
          PrXYZ(x,y,z) = SUM(AKA2P(x,z,:) * PrPAX(y) * AKA2P(y,z,:) * &
            PrGG(z, 1) * PrGG(:,2))
        endif
        PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l, A), x, y) * OcA(x,Genos(l, B))
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXYZ))
enddo

LL = SUM(PrL)
if (LL < -HUGE(0D0))  LL = impossible

if (Parent(A,3-k)>0 .and. LL/=impossible) then
  LL = LL - Lind(Parent(A,3-k))
endif
end subroutine PairPOHA

! #####################################################################

subroutine clustHSHA(SA, SB, k, LL)   ! HS via 3-k, & SB parent of SA; SA,SB FS
use Global
implicit none

integer, intent(IN) :: SA,SB, k
double precision, intent(OUT) :: LL
integer :: l, x, y, z,i, Par(2), GC(2), u
double precision :: PrL(nSnp), PrXYZ(3,3,3), PrGA(3),PrGC(3,2),PrUZ(3,3)

! all checks done by CheckMerge.

! grandparents of opp. parent
 call getFSpar(SA, k, .TRUE., Par(1))  ! TODO: or strict=.FALSE.?
 call getFSpar(SB, k, .TRUE., Par(2))
GC = 0
do i=1,2
  if (Par(1)<0) then
    GC(i) = GpID(i, -Par(1),3-k)
    if (GpID(i, -Par(1),3-k)/=0) then
      if (Par(2) < 0) then
       if (GpID(i,-Par(2),3-k)/=GC(i) .and. GpID(i,-Par(2),3-k)/=0) then
          GC(i) = 0   ! shouldn't happen
        else if (GC(i)==0 .and. GpID(i, -Par(2),3-k)/=0) then
          GC(i) = GpID(i, -Par(2),3-k)
        endif
      endif
    endif
  endif
enddo

PrL = 0D0
do l=1, nSnp
  call ParProb(l, GpID(3-k,SA,k), 3-k, 0, 0, PrGA)
  do i=1,2
    call ParProb(l, GC(i), i, 0, 0, PrGC(:,i))
  enddo
  do z=1,3
    do u=1,3
        PrUZ(u,z) = SUM(AKA2P(z,u,:) * PrGC(u,1) * PrGC(:,2))
    enddo
    do x=1,3    
      do y=1,3
        PrXYZ(x,y,z) = SUM(AKA2P(x,y,:) * PrGA) * XPr(2,y,l,SB,k) *&
          SUM(PrUZ(:,z))
        do i=1,nS(SA,k)
          PrXYZ(x,y,z) = PrXYZ(x,y,z) *OKA2P(Genos(l,SibID(i,SA,k)),x,z)
        enddo 
        do i=1,nS(SB,k)
          PrXYZ(:,y,z) =PrXYZ(:,y,z) *OKA2P(Genos(l, SibID(i,SB,k)),y,z)
        enddo
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXYZ))
enddo
LL = SUM(PrL)

end subroutine clustHSHA

! #####################################################################

subroutine FSHC(A, B, k, LL)  ! FS + parents are HS; B may be neg
use Global
implicit none

integer, intent(IN) :: A,B, k
double precision, intent(OUT) :: LL
integer :: l, x, y, z, Par(2), m, GG(2,2), kG, i, PM(2)
double precision :: PrL(nSnp), PrG(3,2), PrXYZ(3,3,3), PrZ(3)

LL = missing
Par = 0
if (B < 0 .and. A>0) then
  Par(k) = B
  if (Parent(A,k)/=0) then
    LL = impossible
    return
  endif
  if (ALL(Parent(SibID(1:nS(-B,k), -B, k), 3-k) == 0)) then
    Par(3-k) = Parent(A, 3-k)
  else
    call getFSpar(-B, k, .TRUE., Par(3-k))
    if (Par(3-k)==0 .or. (Parent(A,3-k)/=Par(3-k) .and. Parent(A, 3-k)/=0)) then
      LL = impossible
      return
    endif
  endif
else if (B > 0 .and. A>0) then
    do m=1,2
        if (Parent(B,m)==0) then
            Par(m) = Parent(A,m)
        else if (Parent(B,m) /= Parent(A,m) .and. Parent(A,m)/=0) then
            LL = impossible
            return
        else
            Par(m) = Parent(B,m)
        endif
    enddo
else if (B<0 .and. A<0) then
    if (ANY(GpID(:,-B,k)/=0) .or. ANY(GpID(:,-A,k)/=0)) then
        LL = NotImplemented
        return
    else
        Par = 0
    endif
endif

if (hermaphrodites/=0) then
  LL = NotImplemented
  return
endif

GG = 0
kG = 0
PM = 0
do m=1,2
  GG(:, m) = getPar(Par(m), m)
enddo
do m=1,2
    if (GG(m,1)==0 .or. GG(m,2)==0) then  ! GG(m,1)==GG(m,2) not needed
      kG = m
    endif
    if (Par(m)>0) then
      PM(m) = Par(m)
    endif
enddo
if (kG==0) then
    LL = AlreadyAss
    return
endif

PrL = 0D0
do l=1, nSnp
  do m=1,2
    call ParProb(l, GG(3-kG, m), 3-kG, PM(m),0, PrG(:,m))
  enddo
  if (GG(kG,1)/=0) then
    call ParProb(l, GG(kG,1), kG, PM(1),0, PrZ)
  else
    call ParProb(l, GG(kG,2), kG, PM(2),0, PrZ)
  endif
  do x=1,3  ! Par(1)
    do y=1,3  !Par(2)
      do z=1,3  ! GG(kG)
        PrXYZ(x,y,z) = PrZ(z) * SUM(AKA2P(x, z, :) * PrG(:,1)) * &
          SUM(AKA2P(y, z, :) * PrG(:,2))
        if (A>0) then
          PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,A), x, y)
        else
          do i=1, nS(-A,k)
            PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,SibID(i,-A,k)), x, y)
          enddo
        endif
        if (B>0) then
          PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,B), x, y)
        else
          do i=1, nS(-B,k)
            PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,SibID(i,-B,k)), x, y)
          enddo
        endif
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXYZ))
enddo
LL = SUM(PrL)

end subroutine FSHC

! #####################################################################

subroutine pairFSHA(A, B, k, LL) !inbred FS: par k offspring of par 3-k
use Global
implicit none

integer, intent(IN) :: A,B,k
double precision, intent(OUT) :: LL
integer :: l, x, y
double precision :: PrL(nSnp), PrXY(3,3), PrY(3)

if (Parent(A,k)/=0 .or. Parent(B,k)/=0) then
  LL = NotImplemented 
  return
endif  

PrL = 0D0
do l=1, nSnp
  call ParProb(l, Parent(A,3-k), 3-k, -1,0, PrY) 
  do x=1,3
    do y=1,3    
      PrXY(x,y) = PrY(y) * AKAP(x, y, l) * OKA2P(Genos(l,B), x, y) * OKA2P(Genos(l,A), x, y)
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXY))
enddo
LL = SUM(PrL)

end subroutine pairFSHA

! #####################################################################

subroutine pairFSFC(A, B, LL) !inbred FS: FS + parents are FS (or FA)
use Global
implicit none

integer, intent(IN) :: A,B
double precision, intent(OUT) :: LL
integer :: l, x, y,z,v,w
double precision :: PrL(nSnp,2), PrXY(3,3,3,3), PrXYZ(3,3,3,3,3)

if (any(Parent(A,:)/=0) .or. any(Parent(B,:)/=0)) then
  LL = NotImplemented 
  return
endif  

PrL = 0D0
do l=1, nSnp
  do x=1,3
    do y=1,3
      do v=1,3
        do w=1,3
          PrXY(w,v,y,x) = AHWE(v,l) * AHWE(w,l) * AKA2P(x,v,w) * &
           OKA2P(Genos(l,B), x, y) * OKA2P(Genos(l,A), x, y)
          do z=1,3
            PrXYZ(z,w,v,y,x) = PrXY(w,v,y,x) * AKAP(y,z,l) * AKA2P(z,v,w)
          enddo
          PrXY(w,v,y,x) = PrXY(w,v,y,x) * AKA2P(y,v,w)
        enddo
      enddo
    enddo
  enddo
  PrL(l,1) = LOG10(SUM(PrXY))
  PrL(l,2) = LOG10(SUM(PrXYZ))
enddo

LL = MAXVAL( SUM(PrL, DIM=1) )

end subroutine pairFSFC

! #####################################################################

subroutine pairFSselfed(A, B, LL) !A & B both product of selfing by same parent
use Global
implicit none

integer, intent(IN) :: A,B
double precision, intent(OUT) :: LL
integer :: l, x
double precision :: PrL(nSnp), PrX(3)

if (any(Parent(A,:)/=0) .or. any(Parent(B,:)/=0)) then
  LL = NotImplemented 
  return
endif  

PrL = 0D0
do l=1, nSnp
  do x=1,3
    PrX(x) = OKA2P(Genos(l,B), x, x) * OKA2P(Genos(l,A), x, x) * AHWE(x,l)
  enddo
  PrL(l) = LOG10(SUM(PrX))
enddo
LL = SUM(PrL)

end subroutine pairFSselfed

! #####################################################################

subroutine trioGP(A,kA, B, kB, C, kC, LL)  ! B & C both GP of A, A>0 or A<0
use Global
implicit none

integer, intent(IN) :: A,B,C, kA,kB, kC
double precision, intent(OUT) :: LL
integer :: l, x, y, v, BC(2), kBC(2), oppar_BC(2), w, i
double precision :: PrL(nSnp,3), PrXY(3,3), PrXYV(3,3,3), PrB(3), PrC(3),PrA(3), &
  LLBC, PrMB(3), PrMC(3), PrX(3,3), PrY(3,3), LLtmp(3)
logical :: withMate(2)

BC = (/B,C/)
kBC = (/kB, kC/)
opPar_BC = 0
withMate = .FALSE.
if (Complx==0) then
 do x=1,2
   if (BC(x)>0)  opPar_BC(x) = Mate(BC(x))
   if (BC(x)<0)  opPar_BC(x) = DumMate(-BC(x), kBC(x))
   if (opPar_BC(x)/=0) then
    LL = NotImplemented  ! something breaks below. TODO. 
    return
   endif
 enddo
endif
if (B<0 .and. C<0) then  !  with mate of B & C
  do x=1,2
    !if (BC(x)>0)  cycle        
    call getFSpar(-BC(x),kBC(x), .TRUE., opPar_BC(x))
    if (opPar_BC(x)/=0 .or. all(parent(SibID(1:ns(-BC(x),kBC(x)),-BC(x),kBC(x)), 3-kBC(x)) == 0)) then
      withMate(x) = .TRUE.
    endif
  enddo
endif

! parents are -not- ignored, but assumes Par(A,:)=0 
PrL = 0D0
do l=1, nSnp
  PrXY = 0D0
  call OffProb(l, A, kA, PrA)  ! OcA if >0, XPr(1,:,) if <0
  call ParProb(l, B, kB, 0, 0, PrB)
  call ParProb(l, C, kC, 0, 0, PrC)
  do x=1,3
    do y=1,3
      ! one mat GP, other pat GP 
      PrXY(x,y) = SUM(PrA * AKA2P(:, x, y)) * SUM(AKAP(x,:,l) * PrB) * &
        SUM(AKAP(y,:,l) * PrC)
      ! both mat or both pat GP
      do v=1,3
        PrXYV(x,y,v) = SUM(PrA * AKA2P(:, x, y)) * AHWE(y,l) * SUM(AKA2P(x,v,:) * PrB(v) * PrC)
      enddo
    enddo
  enddo
  PrL(l,1) = LOG10(SUM(PrXY))
  PrL(l,2) = LOG10(SUM(PrXYV))
enddo

if (all(withMate)) then   ! TODO?: any   ! B<0 & C<0
  do l=1, nSnp
    PrXY = 0D0
    call OffProb(l, A, kA, PrA) 
    call ParProb(l, B, kB, -1, 0, PrB)
    call ParProb(l, opPar_BC(1), 3-kB, -1, 0, PrMB)
    call ParProb(l, C, kC, -1, 0, PrC)
    call ParProb(l, opPar_BC(2), 3-kC, -1, 0, PrMC)
    do x=1,3
      do y=1,3
        PrXY(x,y) = SUM(PrA * AKA2P(:, x, y))
        do v=1,3
          do w=1,3
            PrX(v,w) = AKA2P(x,v,w) * PrB(v) * PrMB(w)
            do i=1,ns(-B,kB)
              PrX(v,w) = PrX(v,w) * OKA2P(Genos(l,SibID(i,-B,kB)),v,w)
            enddo
          enddo
        enddo
        do v=1,3
          do w=1,3
            PrY(v,w) = AKA2P(y,v,w) * PrC(v) * PrMC(w)
            do i=1,ns(-C,kC)
              PrX(v,w) = PrX(v,w) * OKA2P(Genos(l,SibID(i,-C,kC)),v,w)
            enddo
          enddo
        enddo
        PrXY(x,y) = PrXY(x,y) * SUM(PrX) * SUM(PrY)
      enddo
    enddo
    PrL(l,3) = LOG10(SUM(PrXY))
  enddo
endif

LLtmp = SUM(PrL,DIM=1)
if (.not. all(withMate))  LLtmp(3) = NotImplemented
LLBC = missing
call CalcU(B, kB, C, kC, LLBC)
LLtmp(1:2) = LLtmp(1:2) + LLBC

LL = MAXLL(LLtmp)
if (LL < -HUGE(0D0))  LL = impossible

end subroutine trioGP

! #####################################################################

subroutine trioFS(A,kA, B,kB, C,kC, LL) 
use Global
implicit none

integer, intent(IN) :: A,kA, B,kB, C,kC
double precision, intent(OUT) :: LL
integer :: l, x, y, iX(3), kX(3), s, ParT(2), kP, ParX(2), i_excl(2)
double precision :: PrL(nSnp), PrXY(3,3), PrS(3,3), PrP(3,2)

iX = (/A, B, C/)
kX = (/kA, kB, kC/)

if (Hermaphrodites/=0) then
  do s=1,3
    if (iX(s) >=0)  cycle
    if (DumClone(-iX(s), kX(s))/=0) then
      LL = NotImplemented
      Return
    endif
  enddo
endif

ParT = 0  ! trio parents
i_excl = 0  ! sibship members to exclude when calling ParProb()                                                               
do s=1,3  
  ParX = getPar(iX(s), kX(s))
  do kP=1,2
    if (ParX(kP)/= ParT(kP) .and. ParT(kP)/=0) then
      LL = Impossible
      Return
    else
      ParT(kP) = ParX(kP)
      if (iX(s) > 0 .and. ParT(kP) < 0) then
        i_excl(kP) = iX(s)   ! TODO?: 2 already have ParT as parent
      endif
    endif
  enddo
enddo

do kP=1,2
  do s=1,3  
    if (kX(s)/=kP)  cycle
    if (ParT(kP) == iX(s)) then
      LL = Impossible
      Return 
    endif
  enddo
enddo

PrL = 0D0
do l=1, nSnp
  do kP=1,2
    call ParProb(l,ParT(kP),kP, i_excl(kP),0, PrP(:,kP))
  enddo
  do s=1,3
    call OffProb(l,iX(s), kX(s), PrS(:,s))
  enddo
  do x=1,3
    do y=1,3
      PrXY(x,y) = PrP(x,1) * PrP(y,2) 
      do s=1,3       
        PrXY(x,y) = PrXY(x,y) * SUM(PrS(:,s) * AKA2P(:,x,y))
      enddo   
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXY))
enddo

LL = SUM(PrL)

end subroutine trioFS

! #####################################################################

subroutine trioHS(A,kA, B,kB, C,kC, LL) 
use Global
implicit none

integer, intent(IN) :: A,kA, B,kB, C,kC
double precision, intent(OUT) :: LL
integer :: l, x, y, iX(3), kX(3), s
double precision :: PrL(nSnp, 2), PrX(3), PrXY(3,3), PrS(3,3)

iX = (/A, B, C/)
kX = (/kA, kB, kC/)

PrL = 0D0
do l=1, nSnp
  do s=1,3
    call OffProb(l,iX(s), kX(s), PrS(:,s))    ! OcA if >0, XPr(1,:,) if <0
  enddo

  ! all 3 are halfsibs
  do x=1,3
    PrX(x) = AHWE(x,l)
    do s=1,3
      PrX(x) = PrX(x) * SUM(PrS(:,s) * AKAP(:,x,l))
    enddo
  enddo
  PrL(l,1) = LOG10(SUM(PrX))
  
  ! A-B HS and A-C HS, but not B-C
  do x=1,3
    do y=1,3
      PrXY(x,y) = SUM(PrS(:,1) * AKA2P(:,x,y)) * SUM(PrS(:,2) * AKAP(:,x,l)) * &
        SUM(PrS(:,3) * AKAP(:,y,l)) * AHWE(x,l) * AHWE(y,l)
    enddo
  enddo
  PrL(l,2) = LOG10(SUM(PrXY))      
enddo

LL = MAXVAL(SUM(PrL, DIM=1))

end subroutine trioHS

! #####################################################################

subroutine trioFA(A,kA, B,kB, C,kC, LL)   
! if B<0, then *-B* is considered as FA, NOT the siblings in SB (contrast with pairUA)
! so members of SB become full cousins of A 
use Global
implicit none

integer, intent(IN) :: A,kA, B,kB, C,kC
double precision, intent(OUT) :: LL
integer :: l, k, x, y, z, v,u,w, parB(2), parC(2) !, i, pari
double precision :: PrL(nSnp,2), PrXY(3,3,2), PrPB(3,2), PrPC(3,2), PrUV(3,3), &
   PrWZ(3,3), PrA(3), PrB(3), PrC(3)
logical :: do_indep, do_BC_FS, AncOK

parB = getPar(B,kB)
parC = getPar(C,kC)

do_indep = .TRUE.
do_BC_FS = .TRUE.
do k=1,2
   if (parB(k)/=0 .and. parC(k)/=0) then
    if (parB(k)==parC(k))  do_indep = .FALSE.
    if (parB(k)/=parC(k))  do_BC_FS = .FALSE.
  endif
enddo
LL = Missing
if (.not. do_indep .and. .not. do_BC_FS)  return        

if (do_BC_FS) then
  call ChkAncest(B,kB, C,kC, AncOK)
  if (.not. AncOK)  do_BC_FS = .FALSE.   ! would become own ancestor
  call ChkAncest(C,kC, B,kB, AncOK)
  if (.not. AncOK)  do_BC_FS = .FALSE. 
  if (.not. do_indep .and. .not. do_BC_FS)  return
endif

PrL = 0D0
do l=1, nSnp
  call OffProb(l, A, kA, PrA)  ! OcA if >0, XPr(1,:,) if <0
  call OffProb(l, B, kB, PrB)
  call OffProb(l, C, kC, PrC)
  do k=1,2
    call ParProb(l,parB(k),k,-1,0,PrPB(:,k))
    call ParProb(l,parC(k),k,-1,0,PrPC(:,k))
  enddo

  do x=1,3
    do y=1,3
      do u=1,3
        do v=1,3
          PrUV(u,v) = PrPB(u,1) * PrPB(v,2) * AKA2P(x,u,v) * SUM(PrB * AKA2P(:,u,v))
        enddo
      enddo
      
      if (do_indep) then
        do w=1,3
          do z=1,3
            PrWZ(w,z) = PrPC(w,1) * PrPC(z,2) * AKA2P(y,w,z) * SUM(PrC * AKA2P(:,w,z))
          enddo
        enddo
        PrXY(x,y,1) = SUM(PrA * AKA2P(:,x,y)) * SUM(PrUV) * SUM(PrWZ)
      endif
      
      if (do_BC_FS) then
        do u=1,3
          do v=1,3
            PrUV(u,v) = PrUV(u,v) * SUM(PrC * AKA2P(:,u,v))
          enddo
        enddo
        PrXY(x,y,2) = SUM(PrA * AKA2P(:,x,y)) * SUM(PrUV) * AHWE(y,l)
      endif                                     

    enddo
  enddo 
  if (do_indep)  PrL(l,1) = LOG10(SUM(PrXY(:,:,1)))
  if (do_BC_FS)  PrL(l,2) = LOG10(SUM(PrXY(:,:,2)))
enddo

LL = MaxLL( SUM(PrL, DIM=1) )

end subroutine trioFA

! #####################################################################

subroutine trioFSFA(A,kA, B,kB, C,kC, LL)   ! A+B FS, C FA of both
use Global
implicit none

integer, intent(IN) :: A,kA, B,kB, C,kC
double precision, intent(OUT) :: LL
integer :: l, k, x, y, u, v, parC(2), parB(2)
double precision :: PrL(nSnp), PrXY(3,3,3,3), PrPC(3,2), PrA(3), PrB(3), PrC(3), ALRGP
logical :: GPOK(2)

parB = getPar(B,kB)
if (any(parB/=0)) then  ! else indirect parent assignment to A, messes up CalcCandPar
  LL = NotImplemented
  return
endif

parC = getPar(C,kC)
do k=1,2
  GPOK = .TRUE.
  if (ParC(k)==0)  cycle
  call chkvalidGP(A,kA, parC(k),k, k, GPOK(1), ALRGP)   ! age & ancestor check
  call chkvalidGP(B,kB, parC(k),k, k, GPOK(2), ALRGP)
  if (any(.not. GPOK)) then
    LL = Impossible
    return
  endif
enddo

PrL = 0D0
do l=1, nSnp
  call OffProb(l, A, kA, PrA)  ! OcA if >0, XPr(1,:,) if <0
  call OffProb(l, B, kB, PrB)
  call OffProb(l, C, kC, PrC)
  do k=1,2
    call ParProb(l,parC(k),k,-1,0,PrPC(:,k))
  enddo
  do x=1,3
    do y=1,3     
      do u=1,3
        do v=1,3
          PrXY(v,u,y,x) = SUM(PrA * AKA2P(:,x,y)) * SUM(PrB * AKA2P(:,x,y)) * &
            AHWE(x,l) * AKA2P(y,v,u) * SUM(PrC * AKA2P(:,v,u)) * PrPC(u,1) * PrPC(v,2) 
        enddo
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXY))
enddo

LL = SUM(PrL)

end subroutine trioFSFA

! #####################################################################

subroutine trioHA(A,kA, B,C, LL)  ! B & C > 0
use Global
implicit none

integer, intent(IN) :: A,kA, B, C
double precision, intent(OUT) :: LL
integer :: l, x, y, z, v, PB, PC, kB, kC
double precision :: PrL(nSnp,2), PrXYZ(3,3,3,2), PrV(3), PrA(3), PrPB(3), PrPC(3)

PB = 0
PC = 0
if (all(Parent(B,:)/=0) .or. all(Parent(C,:)/=0)) then
  LL = NotImplemented   ! TODO? 
  return
else
  if (Parent(B,1)/=0) then
    PB = Parent(B,1)
    kB=1
  else
    PB = Parent(B,2)
    kB = 2
  endif
  if (Parent(C,1)/=0) then
    PC = Parent(C,1)
    kC=1
  else
    PC = Parent(C,2)
    kC = 2
  endif
endif

PrL = 0D0
do l=1, nSnp
  call OffProb(l, A, kA, PrA)  ! OcA if >0, XPr(1,:,) if <0
  call ParProb(l, PB,kB, B, 0, PrPB)
  call ParProb(l, PC,kC, C, 0, PrPC)
  do x=1,3
    do y=1,3
      do z=1,3
        PrXYZ(x,y,z,:) = SUM(PrA * AKA2P(:,x,y)) * AKAP(x,z,l) * AHWE(z,l) * &
          SUM(OKA2P(Genos(l,B),z,:) * PrPB)
        ! B+C HS
        PrXYZ(x,y,z,1) = PrXYZ(x,y,z,1) * AHWE(y,l) * SUM(OKA2P(Genos(l,C),z,:) * PrPC)         
        ! B+C U
        do v=1,3
          PrV(v) = AKAP(y,v,l) * AHWE(v,l) * SUM(OKA2P(Genos(l,C),v,:) * PrPC) 
        enddo
        PrXYZ(x,y,z,2) = PrXYZ(x,y,z,2) * SUM(PrV)
      enddo
    enddo
  enddo
  do v=1,2
    PrL(l,v) = LOG10(SUM(PrXYZ(:,:,:,v)))
  enddo
enddo

LL = MAXVAL(SUM(PrL, DIM=1))

end subroutine trioHA

! #####################################################################

subroutine trioGGP(A,kA, B,C, LL)  ! B & C > 0
use Global
implicit none

integer, intent(IN) :: A,kA, B, C
double precision, intent(OUT) :: LL
integer :: l, x, y, z, v
double precision :: PrL(nSnp), PrXV(3,3,3,3), PrA(3), PrB(3), PrC(3)

PrL = 0D0
do l=1, nSnp
  call OffProb(l, A, kA, PrA)  ! OcA if >0, XPr(1,:,) if <0
  call ParProb(l, B,1, 0,0, PrB)
  call ParProb(l, C,2, 0,0, PrC)
  do x=1,3
    do y=1,3
      do z=1,3
        do v=1,3
          PrXV(x,y,z,v) = SUM(PrA * AKA2P(:,x,y)) * AKAP(x,z,l) * AKAP(y,v,l) * &
            SUM(AKAP(z,:,l) * PrB) * SUM(AKAP(v,:,l) * PrC)
        enddo
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXV))
enddo

LL = SUM(PrL, DIM=1) + Lind(B) + Lind(C)

end subroutine trioGGP

! #####################################################################

subroutine trioHSGP(A,kA, B,kB, C,kC, k, LL)  ! A & B HS, C GP of both, via k
use Global
implicit none

integer, intent(IN) :: A,kA, B,kB, C,kC, k
double precision, intent(OUT) :: LL
integer :: l, x, PB(2), GB(2), BB, y,w, PC(2)
double precision :: PrL(nSnp), PrX(3), PrA(3), PrPB(3,2), PrC(3), PrMC(3), PrY(3),&
  PrB(3), PrW(3)

PB = getPar(B,kB)
GB = getPar(PB(k),k)
if (GB(kC)/=0 .and. GB(kC)/=C) then
  LL = impossible
  return
endif

PC = getPar(C, kC)
if (PC(kB) == B .or. PB(kC)==C) then
  LL = NotImplemented
  return
endif

PrL = 0D0
if (B>0) then
  BB = B
else
  BB = 0
endif
do l=1, nSnp
  call OffProb(l, A, kA, PrA)  ! OcA if >0, XPr(1,:,) if <0
  call OffProb(l, B, kB, PrB)
  call ParProb(l, C, kC, 0, 0, PrC)
  call ParProb(l, PB(3-k),3-k,BB,0,PrPB(:,3-k))
  call ParProb(l, PB(k),k, BB, -4, PrPB(:,k))
  if (PB(k)>0) then
    call ParProb(l,GB(3-kC),3-kC,PB(k),0,PrMC)
  else
    call ParProb(l,GB(3-kC),3-kC,0,0,PrMC)
  endif
  do x=1,3
    do y=1,3
      PrY(y) = SUM(AKA2P(x,y,:) * PrC(y) * PrMC)
    enddo
    do w=1,3
      PrW(w) = SUM(PrB(w) * AKA2P(w,x,:) * PrPB(:,3-k))                         
    enddo
    PrX(x) = SUM(PrA * AKAP(:,x,l)) * SUM(PrW) * PrPB(x,k) * SUM(PrY)
  enddo
  PrL(l) = LOG10(SUM(PrX))
enddo

if (C > 0) then
  LL = SUM(PrL) + Lind(C)
else
  LL = SUM(PrL) + CLL(-C,kC)
endif
if (PB(k)>0)   LL = LL - Lind(PB(k))

end subroutine trioHSGP

! #####################################################################

subroutine pairHSGP(A, B,k, LL)   ! HS via k, B is GP of A via 3-k
use Global
implicit none

integer, intent(IN) :: A,B,k
double precision, intent(OUT) :: LL
integer :: l, x, y, z
double precision :: PrL(nSnp), PrXYZ(3,3,3), PrPB(3)

if (Parent(A,3-k)/=0) then
  LL = NotImplemented
  return
endif 

PrL = 0D0
do l=1, nSnp
  call ParProb(l, Parent(B,3-k), 3-k, B, 0, PrPB)
  do x=1,3  ! parent 3-k of A, offspring of B
    do y=1,3  ! shared parent k 
      do z=1,3  ! B
        PrXYZ(x,y,z) =AKAP(x,z,l)*AHWE(y,l)*SUM(AKA2P(z,y,:)*PrPB) * &
          OKA2P(Genos(l,A),x, y) * OcA(z,Genos(l,B))
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXYZ))
enddo
LL = SUM(PrL)

end subroutine pairHSGP

! #####################################################################

subroutine PairGP(A, B, k, focal, LL)  
! calculates LL that B is maternal(k=1) or paternal(k=2) gp of A
use Global
use CalcLik           
implicit none

integer, intent(IN) :: A,B,k, focal
double precision, intent(OUT) :: LL
integer :: l, x, y, curGP(2,2), m, z, v, i, ParAB, Ei,j
double precision :: PrL(nSnp,8), PrPA(3,2), PrG(3), LLtmp(8),&
   PrXZ(3,3,3,8), PrB(3), PrGx(3), PrPB(3,2), PrV(3), PrPAX(3,2), ALR, LLU(2), PrE(3)
logical :: cat(8), GPOK, Aselfed, AncOK
 
LL = missing
call chkValidGP(A,sex(A), B,sex(B), k, GPOK, ALR)  ! age & ancestor check
if (.not. GPOK) then
  LL = impossible
  return                                 
endif 

curGP = 0 
do z=1,2
  curGP(:,z) = getPar(Parent(A,z), z)
enddo

if (Sex(B)<3) then
  m = Sex(B)
else if (curGP(1,k) == 0) then
  m=1    ! doesn't really matter.
else
  m = 2
endif

if (ANY(curGP(:,k) == B)) then
  LL = AlreadyAss 
else if (focal == 4) then
  if (curGP(m,k) /=0)  LL = impossible
else
  if (curGP(m,k) > 0)  LL = impossible
endif
if (LL/=missing) return

if (Complx==0 .and. curGP(3-m,k)==0 .and. Parent(A,3-m)==0) then
  curGP(3-m,k) = Mate(B)
endif

LLtmp = missing
if (Parent(A,k) /= 0) then
  call CalcAgeLR(Parent(A,k),k, B, Sex(B), m, 1, .FALSE., ALR) 
  if (ALR == impossible) then   ! B too young/old to be parent of Parent(A,k)
    LL = impossible 
  else if (Parent(A,k) > 0) then
    call PairPO(Parent(A,k), B, m, focal, LLtmp(1))
    if (LLtmp(1) > 0) then    ! impossible
      LL = impossible
    else 
      call CalcU(Parent(A,k), k, B,k, LLtmp(2))
      if (LLtmp(1) - LLtmp(2) < TA) then
        LL = impossible
      endif
    endif
  else if (Parent(A,k) < 0) then
    if (any(SibID(:,-parent(A,k),k) == B)) then
      LL = impossible
      return
    endif
    if (ns(-Parent(A,k),k) > 1) then
      call AddGP(B, -Parent(A,k), k, LLtmp(1))
      if (LLtmp(1) < 0) then
        ParAB = Parent(A,k)
        ! most-likely alternative: only ParAB PO of A, or only B GP of sibship ParAB
        call CalcU(ParAB, k, B, sex(B), LLtmp(2))
        call setParTmp(A,sex(A),0,k)
        !call CalcU(ParAB, k, A, sex(A), LLtmp(3))
        LLtmp(3) = LLtmp(3)
        call AddGP(B, -ParAB, k, LLtmp(3))
        call CalcU(ParAB, k, B, sex(B), LLtmp(4))
        LLtmp(5) = LLtmp(1) - LLtmp(2)
        LLtmp(6) = LLtmp(3) - LLtmp(4)
        call setParTmp(A,sex(A),ParAB,k)  ! restore
        call CalcU(A, sex(A), B, sex(B), LLtmp(7))
        LL = LLtmp(7) + MINVAL(LLtmp(5:6))
        ! new 2025-12-29; TODO double check if better or worse?
        !LL = LLtmp(1) - CLL(-Parent(A,k), k) + Lind(A)  OLD
      else
        LL = LLtmp(1)
      endif
      return
    endif
  endif
endif
if (LL/=missing) return

Aselfed = .FALSE.
if (hermaphrodites/=0) then
  if (all(Parent(A,:) > 0)) then
    if (Parent(A,k) == Parent(A,3-k))  Aselfed = .TRUE.
  else if (all(Parent(A,:) < 0)) then
    if (DumClone(-Parent(A,k), k) == -Parent(A,3-k))  Aselfed = .TRUE. 
  endif
endif

! cat: 1: non-inbred, 2: double GP, 3: GP+HS, 4: GP+PO, 5 & 6: P-O mating   
if (complx < 2 .or. Aselfed) then
  cat = .FALSE.
  cat(1) = .TRUE.
else
  cat = .TRUE.   
  if (Parent(B,3-k)==Parent(A,3-k) .and. Parent(A,3-k)/=0) then
    cat(1) = .FALSE.
  endif
  if ((focal==1 .or. focal==4) .and. all(parent(A,:)==0)) then
  ! if no parents assigned yet, double GP indistinguishable from PO
    cat(2) = .FALSE.    ! cat(2:3) ?
  endif
  if (Parent(A,3-k)/=0) then
    if (curGP(m, 3-k) == B) then
      cat = .FALSE.
      cat(2) = .TRUE.    ! assignment would make it double GP
    else if (curGP(m,3-k)/=0) then
      cat(2) = .FALSE.  ! already has another GP                                  
    else  
      call ChkAncest(B, Sex(B), Parent(A,3-k), 3-k, AncOK)
      if (.not. AncOK) then
        cat(2) = .FALSE.
      else
        call CalcAgeLR(Parent(A,3-k), 3-k, B,k, 0, 1, .TRUE., ALR)
        if (ALR == impossible .or. ALR < 3.0*TF)  cat(2) = .FALSE.
      endif 
    endif
  endif
  if (focal==3) then  
    cat(3) = .FALSE.
  else if ((focal==1 .or. focal==7) .and. Parent(A,3-k)==0 .and. Parent(B,3-k)/=0) then
    cat(3) = .FALSE.   ! GP+HS; consistency w getPOconfigs() PO+HS
  else if (Parent(B,3-k)==0 .or. Parent(A,3-k)==0 .or. & 
   Parent(B,3-k)==Parent(A,3-k)) then
    if (Parent(A,3-k)/=0 .and. Parent(B,3-k)==0) then
      call ChkAncest(Parent(A,3-k), 3-k, B, 0, AncOK)
      if (.not. AncOK) then
        cat(3) = .FALSE.
      else
        call CalcAgeLR(B,k, Parent(A,3-k), 3-k, 0, 1, .TRUE., ALR)
        if (ALR == impossible .or. ALR < 3.0*TF)  cat(3) = .FALSE.
      endif
    endif
  else
    cat(3) = .FALSE.
  endif
  if (cat(3)) then
    call CalcAgeLR(A,0, B,0, 3-k, 3, .TRUE., ALR) 
    if (ALR == impossible .or. ALR < 3.0*TF)  cat(3) = .FALSE.
  endif 

  if (Parent(A,3-k)==0 .and. Sex(B)==3-k .and. focal/=1 .and. focal/=7) then 
    cat(4) = .FALSE. !.TRUE.   ! not implemented. is already implemented in PairPO?
  else
    cat(4) = .FALSE.
  endif
  if (any(Parent(A,:) /= 0)) then
    cat(5:6) = .FALSE.   ! possible but not necessary to consider (?)
  endif
endif

ParAB = 0 
if (cat(2) .and. Parent(A,3-k)<0) then
  ParAB = Parent(A,3-k)
else if (cat(3)) then
  if (Parent(A,3-k)==0) then
    ParAB = Parent(B,3-k)
  else
    ParAB = Parent(A,3-k)
  endif
endif

cat(7:8) = .FALSE.
if (Sex(B)==4 .and. focal/=7) then
!  if (Parent(A,3-k)/=B)  cat(7) = .TRUE.   ! in-between parent is selfed; identical to LL(PO)
  if (all(parent(A,:)==0) .or. Aselfed)  cat(8) = .TRUE.   ! A also selfed
endif
  
PrL = 0D0
do l=1,nSnp
  call ParProb(l, curGP(3-m,k), 3-m, 0,0, PrG) 
  call ParProb(l, Parent(A,k), k, A, -4, PrPA(:, k)) 
  if (Aselfed) then
    PrPA(:,3-k) = 1D0
  else if (Parent(A,3-k)==Parent(B,3-k) .and. Parent(A,3-k)/=0) then
    call ParProb(l, Parent(A,3-k), 3-k, A, B, PrPA(:, 3-k))
  else
    call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPA(:, 3-k))
  endif
  call ParProb(l, B, 0, 0, 0, PrB)
  if (cat(3)) then
    call ParProb(l, Parent(B,k), k, B, 0, PrPB(:,k))
    call ParProb(l, ParAB, 3-k, -1, 0, PrPB(:,3-k))   ! grandparent contribution only
  endif
  if (cat(2) .or. cat(6)) then
    if (Parent(A,3-k)>0) then
      call ParProb(l, Parent(A,3-k), 3-k, A, -4, PrPAX(:, 3-k))   ! exclude both GPs & A    
      call ParProb(l, curGP(3-k,3-m), 3-m, Parent(A,3-k), 0, PrGx)      
    else
      PrPAx(:,3-k) = 1d0
      call ParProb(l, curGP(3-k,3-m), 3-m, 0, 0, PrGx)
    endif
  endif

  PrXZ = 0D0
  do x=1,3  ! PA(k)     
    do y=1,3  ! PA(3-k)
      if (Aselfed .and. x/=y)  cycle                              
      do z=1,3  !  PrG(3-m)
        PrXZ(x,y,z,:) = OKA2P(Genos(l,A),x,y) * PrPA(x,k) * PrG(z)
        if (cat(1)) then
          PrXZ(x,y,z,1) = PrXZ(x,y,z,1) * PrPA(y,3-k) *&
           SUM(AKA2P(x, :, z) * PrB)  !non-inbred
        endif
        if (cat(2)) then   !inbreeding loop; B double gp
          do v=1,3
            PrV(v) = SUM(AKA2P(y,v,:) * PrGx) * AKA2P(x,v,z) * PrB(v)
          enddo
          PrXZ(x,y,z,2) =PrXZ(x,y,z,2) * PrPAX(y,3-k) * SUM(PrV)
        endif
        if (cat(3)) then  !B GP and HS of A
          do v=1,3
            PrV(v) = AKA2P(x, v,z) * SUM(AKA2P(v,y,:)*PrPB(:,k))*PrPB(y,3-k) * OcA(v,Genos(l,B))
          enddo
          PrXZ(x,y,z,3) =  PrXZ(x,y,z,3) * SUM(PrV)
        endif
        if (ParAB < 0) then     ! implies cat(2) and/or cat(3)
          do i=1, ns(-ParAB,3-k)
            Ei = SibID(i,-ParAB,3-k)
            if (nFS(Ei)==0)  cycle
            call ParProb(l, Parent(Ei,k), k, Ei, -1, PrE)  ! -1: exclude all FS of A
            do j=1,nFS(Ei)
              if (FSID(j,Ei)==A .or. FSID(j,Ei)==B)  cycle
              PrE = PrE * OKA2P(Genos(l,FSID(j,Ei)), :, y)
            enddo
            if (cat(3) .and. .not. all(PrE==1D0)) PrXZ(x,y,z,3) =  PrXZ(x,y,z,3) * SUM(PrE) 
            if (cat(2) .and. .not. all(PrE==1D0)) PrXZ(x,y,z,2) =  PrXZ(x,y,z,2) * SUM(PrE)    
          enddo
        endif
        if(cat(4)) then
        ! TODO?
        endif
        if (cat(5)) then  ! only when Parent(A,:) == 0
          if (x /= z) then
            PrXZ(x,y,z,5) = 0D0
          else
            PrXZ(x,y,z,5) = SUM(OKA2P(Genos(l,A),x,y) * AKA2P(x,y,:) * PrB * AHWE(y,l))
          endif
        endif
        if (cat(6)) then  ! only when Parent(A,:) == 0
          do v=1,3
            PrV(v) = SUM(AKA2P(y,x,:) * PrGx) * AKA2P(x,v,z) * PrB(v)
          enddo
          PrXZ(x,y,z,6) = PrXZ(x,y,z,6) * SUM(PrV)
        endif
        if (cat(7)) then  ! B double GP, A's parent is selfed
          PrXZ(x,y,z,7) = OKA2P(Genos(l,A),x,y) * PrPA(x,k) * PrPA(y,3-k) * &
            AKA2P(x,z,z) * PrB(z)
        endif
        if (cat(8)) then  ! as 7, A also selfed
          if (x/=y) then
            PrXZ(x,y,:,8) = 0D0
          else
            PrXZ(x,x,z,8) = OKA2P(Genos(l,A),x,x) * AKA2P(x,z,z) * PrB(z)
          endif
        endif
      enddo
    enddo
  enddo
  do i=1,8  ! inbred/non-inbred
    if (cat(i))   PrL(l,i) = LOG10(SUM(PrXZ(:,:,:,i)))
  enddo
enddo

LLtmp = SUM(PrL, DIM=1)
WHERE(LLtmp((/1,2,5,6,7,8/)) <0)  LLtmp((/1,2,5,6,7,8/)) = LLtmp((/1,2,5,6,7,8/)) + Lind(B)
if (Parent(A,k)>0) then
  WHERE(LLtmp <0)  LLtmp = LLtmp - Lind(Parent(A,k))
endif
if (cat(2) .and. Parent(A,3-k)>0) then
    LLtmp(2) = LLtmp(2) - Lind(Parent(A,3-k))
endif

LLU = missing           
if (cat(3) .and. ParAB < 0) then   ! parAB may or may not be Parent(A,3-k)
  call CalcU(A,sex(A), B,sex(B), LLU(1))  ! intended reference LL
  if (Parent(A,3-k)==parAB .and. Parent(B,3-k)==0) then  ! actual reference LL
    call CalcU(Parent(A,3-k), 3-k, B,sex(B), LLU(2))      
  else if (Parent(A,3-k)==0 .and. Parent(B,3-k)==parAB) then
    call CalcU(Parent(B,3-k), 3-k, A,sex(A), LLU(2))  
  else if (Parent(A,3-k)==parAB .and. Parent(B,3-k)==ParAB) then
    call CalcU(parAB, 3-k, 0,0, LLU(2))    ! includes A & B
  endif
  LLtmp(3) = SUM(PrL(:,3)) - LLU(2) + LLU(1)         
endif                          

if (cat(2) .and. Parent(A,3-k)<0) then
  call CalcU(A,sex(A), B,sex(B), LLU(1))  ! intended reference LL
  call CalcU(Parent(A,3-k), 3-k, B,sex(B), LLU(2))
  LLtmp(2) = LLtmp(2) - LLU(2) + LLU(1)   ! LLtmp(2)  includes Lind(B)
endif

LL = MaxLL(LLtmp)
if (LL >= 0) then
  LL = impossible
endif

end subroutine PairGP

! #####################################################################

subroutine LRGG(A,k, B,kB, LR)
use Global
implicit none

integer, intent(IN) :: A,k,B,kB
double precision, intent(OUT) :: LR
integer :: x, y, l
double precision :: PrXY(3,3,2), PrL(nSnp), PrPA(3), PrB(3)

PrL = 0D0
do l=1,nSnp
  call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPA)
  call ParProb(l, B, kB, 0, 0, PrB)
  PrXY = 1D0
  do x=1,3  ! PA(k)
    do y=1,3  ! B
      PrXY(x,y,:) = SUM(OKA2P(Genos(l,A),x,:) * PrPA) * PrB(y)
      PrXY(x,y,1) = PrXY(x,y,1) * AKAP(x,y,l)
      PrXY(x,y,2) = PrXY(x,y,2) * AHWE(x,l)
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXY(:,:,1))) - LOG10(SUM(PrXY(:,:,2)))
enddo
LR = SUM(PrL)

end subroutine LRGG

! #####################################################################

subroutine PairGGP(A, B, k, fcl, LL)   
! calculates LL that B is maternal(k=1) paternal(k=2), or double(k3) ggp
use Global
implicit none

integer, intent(IN) :: A,B,k, fcl
double precision, intent(OUT) :: LL
integer :: l, x, y,z,w, m, n, AncA(2,mxA), MateB
double precision :: PrL(nSnp,4), PrXY(3,3), PrXZ(3,3,3), PrXW(3,3,3,3,2), LLtmp(4),&
  PrPA(3),PrB(3), PrG(3,2), PrPAX(3,2), ALR
logical :: MaybeLoop(2), AncOK(2)
  
LL = missing
LLtmp = missing
AncOK = .TRUE.
if (AgeDiff(A, B) <= 0) then  ! B younger than A
  LL = impossible
else if (B==Parent(A,k)) then
  LL = impossible
else
  call ChkAncest(B,k, A, 0, AncOK(1))
  if (.not. AncOK(1)) then
    LL = impossible
    return
  endif
  if (Parent(A,3-k)/=0 .and. Parent(A,3-k)==Parent(B,3-k)) then
    LL = NotImplemented   ! not impossible, just unlikely.
    return
  endif
  if (Parent(A,k)/=0) then
    call ChkAncest(B,k,Parent(A,k),k, AncOK(2))
    if (.not. AncOK(2))   LL = impossible
  endif
endif
if (LL==impossible .or. LL==NotImplemented) return

if (Complx==0 .and. sex(B)<3) then  ! monogamous
  mateB = getMate(B,sex(B))
  if (mateB/=0 .and. Parent(A,3-sex(B))==mateB) then
    LL=impossible
    return
  endif
endif

if (Parent(A,k)>0) then 
  if (ANY(Parent(Parent(A,k), :)/=0)) then
    LL = NotImplemented    ! should be picked up elsewere
  else
    call PairGP(Parent(A,k), B, k, 4, LLtmp(1))
    if (LLtmp(1) > 0) then    
      LL = LLtmp(1)
    endif
  endif
else if (Parent(A,k)<0) then
  if (ANY(GpID(:,-Parent(A,k),k)/=0)) LL = NotImplemented
endif
if (LL/=missing) return

MaybeLoop = .FALSE.
if (fcl/=4) then  ! double GGP indistinguishable from GP
  if (ALL(Parent(A,:)==0)) then
    MaybeLoop = .TRUE.
  else
    do m=1,2
      if (ALL(getPar(Parent(A,m),m)==0)) then
        call CalcAgeLR(Parent(A,m),m,B,Sex(B),k,4, .TRUE., ALR)
        if (ALR/=impossible .and. ALR>TF) then
          MaybeLoop(m) = .TRUE.
        endif
      endif
    enddo
  endif
endif     

AncA = 0
call getAncest(A, k, AncA) 

PrL = 0D0 
do l=1,nSnp
  call ParProb(l, B, 0, 0, 0, PrB)
  call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPA)  ! with GP con
  PrPAX = 1D0
  do n=1,2
    if (n/=k .and. .not. MaybeLoop(n)) cycle
    if (Parent(A,n)/=0) then
      call ParProb(l, Parent(A,n), n, A, -4, PrPAX(:,n))  !=OcA if >0
    endif
    if (ALL(MaybeLoop)) then
      PrG(:,n) = AHWE(:,l)
      do m=1,2
        if (AncA(m, n+2)/=0) then ! either GP ==0
          if (Parent(A,n)>0) then
            call ParProb(l,AncA(m, n+2), m,Parent(A,n),0, PrG(:,n))   
          else
            call ParProb(l,AncA(m, n+2), m,0,0, PrG(:,n))    
          endif
        endif
      enddo
    endif
  enddo
  
  PrXY = 0D0
  PrXZ = 0D0
  PrXW = 0D0
  do x=1,3  
    do y=1,3 
      PrXY(x,y) = SUM(OKA2P(Genos(l,A),x,:) * PrPA) * PrPAX(x,k) * &
       AKAP(x,y,l) * SUM(AKAP(y, :, l) * PrB)
      if (ANY(MaybeLoop)) then
        do z=1,3  !consider double GGP (k & 3-k, or 2x k)
         if (ALL(MaybeLoop)) then
          do w=1,3
            PrXW(x,y,z,w,:) = OKA2P(Genos(l,A),x,z) *PrPAX(x,k) *&
             PrPAX(z,3-k) * SUM(AKA2P(z,w,:)* PrG(:,3-k))* &
             SUM(AKA2P(x,y,:)*PrG(:,k))
            PrXW(x,y,z,w,1) = PrXW(x,y,z,w,1) * SUM(AKAP(y,:,l) * AKAP(w,:,l) * PrB)
            PrXW(x,y,z,w,2) = PrXW(x,y,z,w,2) * SUM(AKAP(y,:,l) * AKAP(w,:,l) * AHWE(:,l))
          enddo
          endif
          if (MaybeLoop(k)) then
            PrXZ(x,y,z) = SUM(OKA2P(Genos(l,A),x,:) *PrPA) * PrPAX(x,k) *&
             AKA2P(x,y,z) * SUM(AKAP(y,:,l) * AKAP(z,:,l) * PrB)
          endif
        enddo
      endif
    enddo
  enddo
  PrL(l,1) = LOG10(SUM(PrXY))
  if (ALL(MaybeLoop)) then
    PrL(l,2) = LOG10(SUM(PrXW(:,:,:,:,1)))
    PrL(l,4) = LOG10(SUM(PrXW(:,:,:,:,2)))
  endif
  if (MaybeLoop(k))    PrL(l,3) = LOG10(SUM(PrXZ))
enddo

LLtmp = SUM(PrL, DIM=1)   
do n=1,2
  if (Parent(A,n)>0) then
    if (n==k) then
      LLtmp(1) = LLtmp(1) - Lind(Parent(A,n))
      if (MaybeLoop(k)) then
        LLtmp(3) = LLtmp(3) - Lind(Parent(A,n))
      endif                             
    endif
    if (ALL(MaybeLoop)) then
      LLtmp(2) = LLtmp(2) - Lind(Parent(A,n))
      LLtmp(4) = LLtmp(4) - Lind(Parent(A,n))
    endif
  endif
enddo

if (ALL(MaybeLoop)) then  ! compare to inbreeding loop w/o B at the top. 
  if (LLtmp(4) > Lind(A)) then
    LLtmp(2) = LLtmp(2) - LLtmp(4) + Lind(A)
  endif
endif

LL = MaxLL(LLtmp(1:3)) + Lind(B)

end subroutine PairGGP

! #####################################################################

 subroutine PairGA(A, B, k, hf, LL)   ! B FS/HS of GP
use Global
implicit none

integer, intent(IN) :: A,B,k, hf
double precision, intent(OUT) :: LL
integer :: l, x, y,v,w, m, n
double precision :: PrL(nSnp), PrX(3,3,3,3), PrPA(3), PrGG(3, 2) 
logical :: AncOK

LL = missing
if (AgeDiff(A, B) <= 0) then  ! B younger than A
  LL = impossible
else if (ANY(Parent(A,:)==B)) then
  LL = NotImplemented
else
  AncOK = .TRUE.
  call ChkAncest(B,0,A,0,AncOK)
  if (.not. AncOK) then
    LL = impossible
  else if (Parent(A,3-k)/=0 .and. Parent(A,3-k)==Parent(B,3-k)) then
    LL = NotImplemented   ! not impossible, just unlikely.
  endif
endif
if (LL/=missing) return

m = 3-k  ! most neutral, doesn't matter in most cases
if (Parent(A,k)/=0) then
  LL = NotImplemented
  return
endif

if (Complx==0 .and. sex(B)<3) then
  if (Mate(B)/=0 .and. Parent(A,3-sex(B))==Mate(B)) then
    LL = NotImplemented
    return
  endif
endif

PrL = 0D0    
do l=1,nSnp
  PrX = 0D0
  call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPA) 
  do n=1,2
    call ParProb(l, Parent(B, n), n, B, -1, PrGG(:,n))
  enddo

  do v=1,3
    do w=1,3
      do x=1,3  
        do y=1,3
          PrX(x,y,v,w) = AKAP(x,y,l) * PrGG(v,1) * PrGG(w,2)
          if (hf==3) then
            PrX(x,y,v,w) = PrX(x,y,v,w) * AKA2P(y, v, w)
          else if (hf==1) then
            PrX(x,y,v,w) = PrX(x,y,v,w) * AKAP(y, v, l)
          else if (hf==2) then
            PrX(x,y,v,w) = PrX(x,y,v,w) * AKAP(y, w, l)
          endif
          PrX(x,y,v,w) = PrX(x,y,v,w) * SUM(OKA2P(Genos(l,A),x,:) * PrPA)
        enddo
      enddo
      PrX(:,:,v,w) = PrX(:,:,v,w) * OKA2P(Genos(l,B), v, w)
    enddo
  enddo         
  PrL(l) = LOG10(SUM(PrX))
enddo

LL = SUM(PrL)                            

end subroutine PairGA

! #####################################################################

subroutine PairUA(A, B, kA, kB, LL)
! B half sib or full sib (kB=3) of parent kA of A?
use Global
use CalcLik          
use sort_module
implicit none

integer, intent(IN) :: A,B,kA, kB  ! kB=3 : full sibs
double precision, intent(OUT) :: LL
integer :: AA(maxSibSize), PA, nA, GA(2), BB(maxSibSize), nB, PB(2), BBx(maxSibSize, 2), &
  nBx(2), Mates(maxSibSize, 2), GG(2), AncG(2, 2,mxA), AncA(2,mxA), AB(2*maxSibSize), &
  catA(maxSibSize), catG(2), catB(maxSibSize), GGP, GGG(2), doneB(maxSibSize), &
  Bj, w, l, x, g, y, z, i, r,u,j,e,Ei,m, PAx, DoQuickB
double precision :: PrL(nSnp), PrG(3,2), PrXYZ(3,3,3), PrPA(3), PrA(3), &
  PrPB(3), PrGA(3), PrAB(3,3,3,2), PrE(3), PrH(3), PrGG(3), PrEW(3,3), PrW(3), PrXY(3,3), PrFS(3,3)     
integer, allocatable, dimension(:) :: UseEE, MateABpar, TypeEE
double precision, allocatable, dimension(:,:,:) :: PrEE
logical :: PAselfed, SIMPL, AncOK, MateLoop(maxSibSize,2), DoAZ

AA = 0  
if (A>0) then  
  PA = Parent(A, kA)
  if (Complx==0 .and. PA==0 .and. Parent(A,3-kA)/=0) then
    PA = getMate(Parent(A,3-kA),3-kA)
    if (PA==B .and. (B>0 .or. kA==kB)) then
      LL = impossible
      return
    endif
  endif
  if (PA<0) then
    nA = ns(-PA,kA)
    AA(1:nA) = SibID(1:nA, -PA, kA)
  else 
    nA = 1
    AA(1) = A
  endif
else
  nA = nS(-A, kA)
  AA(1:nA) = SibID(1:nA, -A, kA)
  PA = A
endif
GA = getPar(PA, kA)                   
  
if (nA==0) then
  LL = NotImplemented
  return
endif

BB = 0
nB = 0
PB = 0
BBx = 0
nBx = 0
Mates = 0
if (B > 0) then
  nB = 1
  BB(1) = B  ! for cat checks
  PB = Parent(B,:)
  do m=1,2
    if (kB<3 .and. m/=kB)  cycle
    if (Parent(B, m) >=0) then
      nBx(m) = 1
      BBx(1, m) = B
    else 
      nBx(m) = nS(-Parent(B, m), m)
      BBx(1:nBx(m), m) = SibID(1:nBx(m), -Parent(B, m), m)  ! half sibs
    endif
    do j=1,nBx(m)
      Mates(j,m) = Parent(BBx(j, m), 3-m)
    enddo
    if (Complx==0 .and. PB(m)==0 .and. PB(3-m)/=0) then  ! TODO: for dummies?
      PB(m) = getMate(PB(3-m),3-m)
      if (PB(m)==A .and. (kA==m .or. A>0)) then
        LL = impossible
        return
      endif
    endif
  enddo
else ! if (B < 0) then
  nB = nS(-B, kB)
  BB(1:nB) = SibID(1:nB, -B, kB)
  PB(kB) = B
  PB(3-kB) = 0
  do j=1,nB
    Mates(j,kB) = Parent(BB(j), 3-kB)
  enddo
endif

if (PB(kA) == PA .and. PA/=0 .and. (kB==kA .or. kB==3)) then
  LL = impossible
  return
endif

if (Hermaphrodites/=0 .and. PA<0) then
  if (DumClone(-PA,kA)/=0) then
    if (DumClone(-PA,kA) == -PB(3-kA))  then  
      if (kB==3 .or. kB==3-kA) then
        LL = impossible   ! would become own parent
      else
        LL = NotImplemented  ! TODO?
      endif
      return
    else if (kA==kB .and. any(DumClone(-PA,kA) == -Mates(:,kA))) then
      LL = NotImplemented
      return
    endif
  else 
    do m=1,2
      if (PB(m)>=0)  cycle
      if (DumClone(-PB(m),m)==0)  cycle
      if (any(DumClone(-PB(m),m) == -Parent(AA(1:nA), 3-m))) then
        LL = NotImplemented
        return  
      endif
    enddo
  endif
endif
if (Hermaphrodites/=0 .and. kB==3) then   ! implemented only for kB<3
  do m=1,2
    if (PB(m)>=0)  cycle 
    if (DumClone(-PB(m),m)/=0) then
      LL = NotImplemented
      return  
    endif  
  enddo
endif


LL = missing
if (kB < 3) then
  if (B < 0 .and. GA(kB) == B) then
    LL = AlreadyAss
  else if (B > 0 .and. GA(kB) == PB(kB) .and. PB(kB)/=0)  then
    LL = AlreadyAss
  else if (GA(3-kB)==PB(3-kB) .and. PB(3-kB)/=0) then ! B>0; FA not HA
    LL = impossible
  endif
else if (GA(1)==PB(1) .and. GA(2)==PB(2) .and. GA(1)/=0 .and. &
  GA(2)/=0) then  ! kB==3
  LL = AlreadyAss
endif
if (LL /= missing) return

GG = 0  ! parent of B, GP of A
AncG = 0
do x=1,2
  if (x/=kB .and. kB/=3) cycle
  if (GA(x)==0) then
    GG(x) = PB(x)
  else if (GA(x)/=PB(x) .and. PB(x)/=0) then
    LL = impossible
  else
    GG(x) = GA(x)
  endif
  if (ANY(AA(1:nA)==GG(x))) then
    LL = impossible
  endif
  do j=1,nBx(x)
    if (PA<0 .and. Parent(BBx(j,x), 3-x)<0) then
      if (GpID(kA, -Parent(BBx(j,x), 3-x), 3-x) == PA) then
        LL = NotImplemented
      endif
    endif
  enddo
enddo
if (LL /= missing) return

do x=1,2
  call GetAncest(GG(x), x, AncG(x, :, :))
  if (PA/=0 .and. (x==kB .or. kB==3)) then
     if (ANY(AncG(x,kA,:) == PA)) then
      LL = impossible
      return
    endif
  endif
enddo

AncA = 0
if (A > 0) then
  if (ANY(AncG == A)) then
    LL = impossible
  endif
else if (A < 0) then
  if (ANY(AncG(:, kA, 2:mxA) == A)) then
    LL = impossible
  endif
endif
if (B > 0) then
  if (ANY(AncG == B)) then
    LL = impossible
  endif
else if (B < 0) then
  if (ANY(AncG(:, kB, 3:mxA) == B)) then
    LL = impossible
  endif
endif
if (kB<3 .and. any(GA < 0)) then
  call GetAncest(A, kA, AncA)   
  if (B<0) then
    if (ANY(AncA(kB,3:mxA)==B))  LL = NotImplemented  ! B is GGP; 
  else if (B>0) then
    if (ANY(AncA(:,3:mxA)==B))  LL = NotImplemented
  endif
endif
if (LL /= missing) return

do x=2,mxA
  do y=1,2
    do g=1,2
      if (AncG(g,y,x) > 0) then
        if (A > 0) then
          if (AgeDiff(A, AncG(g,y,x)) < 0)  LL = impossible  ! A older than putative ancestor 
        else if (A<0) then
         do i=1, ns(-A,kA)
            if (AgeDiff(SibID(i,-A,kA),AncG(g,y,x)) < 0)  LL = impossible
          enddo
        endif
        if (x==2) cycle 
        if (B > 0) then
          if (AgeDiff(B, AncG(g,y,x)) < 0)  LL = impossible  
        else if (B<0) then
          do i=1, ns(-B,kB)
            if(AgeDiff(SibID(i,-B,kB),AncG(g,y,x)) < 0)  LL = impossible
          enddo 
        endif
      endif
    enddo
  enddo
enddo
if (LL /= missing) return

PAselfed = .FALSE.
if (hermaphrodites/=0) then
  if (kB==3) then
    if (all(PB < 0)) then
      if (DumClone(-PB(1),1) == -PB(2))  PAselfed = .TRUE.
    endif
  else if (GA(3-kB) < 0 .and. PB(kB)<0) then
    if (DumClone(-GA(3-kB),3-kB) == -PB(kB))  PAselfed = .TRUE.
  endif
endif

!==============================================

allocate(UseEE(nA+nB))
allocate(TypeEE(nA+nB))
allocate(PrEE(3,3, nA+nB))
allocate(MateABpar(nA+nB))
UseEE = 0
TypeEE = 0
PrEE = 0D0
MateABpar = 0
AB = 0

if (kA==kB) then
  AB(1:nA) = AA(1:nA)
  AB((nA+1):(nA+nB)) = BB(1:nB)
  call FindEE(AB(1:(nA+nB)), nA, nB, kA, UseEE, MateABpar)  ! may reorder AA, BB
  AA(1:nA) = AB(1:nA)
  BB(1:nB) = AB((nA+1):(nA+nB))
  TypeEE = 3-kA
else if (kA/=kB  .and. kB/=3) then
  call FindEE(AA(1:nA), nA, 0, kA, UseEE(1:nA), MateABpar(1:nA)) 
  call FindEE(BB(1:nB), nB, 0, kB, UseEE((nA+1):(nA+nB)), MateABpar((nA+1):(nA+nB)))
  do i=1, nB
    if (UseEE(nA+i)/=0) then
      UseEE(nA+i) = nA + UseEE(nA+i)
    endif
  enddo
  TypeEE(1:nA) = 3-kA
  TypeEE((nA+1):(nA+nB)) = 3-kB
endif

! TODO: kB==3, BBx i.o. BB
!==============================================

PrL = 0D0
PAx = 0
if (nA>0)  PAx = Parent(AA(1),3-kA)                       
 
if (A>0 .and. B>0 .and. PA>=0 .and. PAx>=0 .and. ALL(PB>=0) .and. all(GG>=0)) then  ! quicker.
  do l=1, nSnp
    call ParProb(l, PAx, 3-kA, A, 0, PrPA)
    if (kB == 3) then
      do g=1,2
        call ParProb(l, GG(g), g, 0, 0, PrG(:,g))  ! >=0
      enddo        
    else
      call ParProb(l, GG(kB), kB, 0, 0, PrG(:,kB))  ! >=0
      if (PA>0) then
        call ParProb(l, GA(3-kB), 3-kB, PA, 0, PrGA)
      else
        PrGA = AHWE(:,l)
      endif
      call ParProb(l, PB(3-kB), 3-kB, B, 0, PrPB)
    endif    
    
    PrXYZ = 0D0
    do z=1,3
      do y=1,3
        do x=1,3
          if (kB == 3) then
            PrXYZ(x,y,z) = AKA2P(x,y,z)*PrG(y,1)*PrG(z,2)
          else
            PrXYZ(x,y,z) = AKA2P(x,y,z)*PrG(y,kB)*PrGA(z)   
          endif           
          PrXYZ(x,y,z) = PrXYZ(x,y,z) *SUM(OKA2P(Genos(l,A), x, :) * PrPA)
          if (kB==3) then
            PrXYZ(x,y,z) = PrXYZ(x,y,z) *OKA2P(Genos(l,B),y,z)
          else
            PrXYZ(x,y,z) = PrXYZ(x,y,z) *SUM(OKA2P(Genos(l,B),y,:) * PrPB)                       
          endif
          if (PA>0) then
            PrXYZ(x,y,z) = PrXYZ(x,y,z) * OcA(x,Genos(l,PA))
          endif
        enddo
      enddo
    enddo
    PrL(l) = LOG10(SUM(PrXYZ))
  enddo
  
  LL = SUM(PrL)
  if (PA>0) then
    LL = LL - Lind(PA)
  endif
  deallocate(UseEE)
  deallocate(PrEE)
  deallocate(MateABpar)
  deallocate(TypeEE)
  return
endif

!==============================================
if (A>0 .and. B<0 .and. kB/=3) then
  SIMPL = .TRUE.
  AncOK = .TRUE.
  if(ANY(Parent(A,:)<0) .or. Parent(A,kA)/=0) then   ! ANY(Parent(A,:)/=0)
    SIMPL = .FALSE.
  endif
  do j=1,nB
    if (nFS(BB(j))==0)  cycle
    Bj = BB(j)
    call ChkAncest(Bj, 0, A, Sex(A), AncOK)
    if (.not. AncOK)  SIMPL = .FALSE.
    call ChkAncest(A, Sex(A), Bj, 0, AncOK)
    if (.not. AncOK)  SIMPL = .FALSE.
    if (.not. SIMPL)  exit
  enddo
  if (SIMPL) then
    call ChkDoQuick(-B, kB, DoQuickB)
    do l=1, nSnp
      call ParProb(l,Parent(A,3-kA),3-kA,0,0,PrE)
      if (DoQuickB == -2) then  ! all are FS; ns(s,k) = ns(..,3-k)
        call ParProb(l, Parent(Bj,3-kB),3-kB,-1,0,PrH)
        PrFS = FSLik(l,Bj)                  
      endif
      do y=1,3
        do x=1,3
          if (DoQuickB == -2) then  
            do z=1,3
              PrW(z) = PrFS(z,y) * XPr(2,y,l,-B,kB) * PrH(z) * AKAP(x,y,l) * &
               SUM(OKA2P(Genos(l,A),x,:) * PrE)
            enddo
            PrXY(x,y) = SUM(PrW)         
          else
            PrXY(x,y) = XPr(3,y,l, -B,kB) * AKAP(x,y,l) * SUM(OKA2P(Genos(l,A),x,:) * PrE)
          endif
        enddo
      enddo
      PrL(l) = LOG10(SUM(PrXY))
    enddo
    LL = SUM(PrL)
    deallocate(UseEE)
    deallocate(PrEE)
    deallocate(MateABpar)
    deallocate(TypeEE)
    return
  endif
endif 

!==============================================

if (A<0 .and. B>0 .and. all(UseEE==0)) then
  SIMPL = .TRUE.
  AncOK = .TRUE.              
  if (any(Parent(B,:) < 0)) then
    SIMPL = .FALSE.
  else if (kB < 3) then
    if (Parent(B,kB)/=0)  SIMPL = .FALSE.
  else
    if (any(Parent(B,:)/=0))  SIMPL = .FALSE.
  endif
  if (any(GpID(:,-A,kA)/=0))  SIMPL = .FALSE.
  if (any(Parent(AA(1:nA), 3-kA) < 0))  SIMPL = .FALSE.
  if (SIMPL) then
    do j=1,nA
      call ChkAncest(AA(j), 0, B, Sex(B), AncOK)
      if (.not. AncOK)  SIMPL = .FALSE.
      call ChkAncest(B, Sex(B), AA(j), 0, AncOK)
      if (.not. AncOK)  SIMPL = .FALSE.
    enddo
  endif
  if (SIMPL) then
    do l=1, nSnp
      if (kB < 3)  call ParProb(l,Parent(B,3-kB),3-kB,0,0,PrE)
      do x=1,3
        do y=1,3
          PrXY(x,y) = XPr(1,x,l, -A,kA) * AHWE(y,l)
          if (kB<3) then
            PrXY(x,y) = PrXY(x,y) * AKAP(x,y,l) * SUM(OKA2P(Genos(l,B),y,:) * PrE)
          else
            PrXY(x,y) = PrXY(x,y) * SUM(AKA2P(x,y,:) * OKA2P(Genos(l,B),y,:) * AHWE(:,l))
          endif
        enddo
      enddo
      PrL(l) = LOG10(SUM(PrXY))
    enddo
    LL = SUM(PrL)
    deallocate(UseEE)
    deallocate(PrEE)
    deallocate(MateABpar)
    deallocate(TypeEE)
    return
  endif
endif 

!============================================

catA=0
catG = 0
catB = 0
GGP = 0
GGG = 0
do i = 1, nA
  if (Parent(AA(i), 3-kA)==0) cycle
  if (hermaphrodites/=0 .and. PA<0) then
    if (Parent(AA(i), 3-kA) == -DumClone(-PA,kA)) then
      catA(i) = 12
      cycle
    endif
  endif
  if (kA/=kB .and. Parent(AA(i), 3-kA) == GG(3-kA)) then  !incl. kB=3
    catA(i) = 1  
  else if (kA==kB .and. Parent(AA(i), 3-kA)==GA(3-kA)) then
    catA(i) = 2
    UseEE(i) = 0
  else if (hermaphrodites/=0 .and. PB(kA)<0) then
    if (Parent(AA(i), 3-kA) == -DumClone(-PB(kA),kA)) then
      LL = NotImplemented  ! TODO
      return
    endif
  else 
    if (Parent(AA(i), 3-kA)<0) then
      if (GpID(kA,-Parent(AA(i), 3-kA),3-kA) == PA .and. PA/=0) then
        catA(i) = 7  ! Ai inbred
      endif
    endif
    if (kA==kB) then
      do j=1, MAX(nB, nBx(kB))
        if (B>0) then
          Bj = BBx(j,kB)
        else
          Bj = BB(j)
        endif
        if (AA(i) == Bj) cycle
        if (Parent(AA(i), 3-kA) == Parent(Bj, 3-kA) .and. Parent(Bj, 3-kA) < 0) then
          if (B<0)  catA(i) = 3
        else if (Parent(AA(i), 3-kA) == Bj) then
          catA(i) = -j
        endif
      enddo
    endif
  endif
  do g=1,2
    if (kB/=g .and. kB/=3) cycle
    if (Parent(AA(i), 3-kA) < 0) then
      if (GpID(g,-Parent(AA(i), 3-kA),3-kA) == GG(g) .and. GG(g)/=0) then
        if (g==kB .or. (kB==3 .and. g==3-kA)) then
          if (ALL(GpID(:,-Parent(AA(i), 3-kA),3-kA) == GG) .and. ALL(GG/=0)) then
            catA(i) = 10
          else
            catA(i) = 8  ! via y
          endif
        else
          catA(i) = 9  ! via z
        endif
      endif
    endif
    GGG = getPar(GG(g), g)
    if (Parent(AA(i),3-kA) == GGG(3-kA)) then
      catA(i) = 5  ! TODO? 4+g when kB==3
      catG(g) = 2
      GGP = GGG(kA)
      UseEE(i) = 0
    ! TODO: parent(parent(A,3-kA),kB) == B), B<0
    endif
  enddo
enddo
if (kB/=3) then   ! TODO: for kB==3
  do j=1, MAX(nB, nBx(1), nBx(2))
    if (B>0) then
      Bj = BBx(j,kB)
    else
      Bj = BB(j)
    endif
    if (Parent(Bj,3-kB)==0) cycle
    if (hermaphrodites/=0 .and. PB(kB) < 0) then
      if (Parent(Bj, 3-kB) == -DumClone(-PB(kB),kB) .and. DumClone(-PB(kB),kB)/=0) then
        catB(j) = 12
        cycle
      endif
    endif
    if (Parent(Bj, 3-kB) == GA(3-kB) .and. GA(3-kB)/=0) then  
      catB(j) = 2
      if (B<0) UseEE(nA+j) = 0
      if (B>0) UseEE(nA+1) = 0
    else if (B<0 .and. Parent(Bj,3-kB)<0) then
      if (GpID(kB, -Parent(Bj,3-kB),3-kB) == GG(kB) .and. GG(kB)/=0) then
        catB(j) = 7
      else if (GpID(kA, -Parent(Bj,3-kB),3-kB) == PA .and. PA/=0) then
        catB(j) = 8
      endif
    endif
    do g=1,2
      GGG = getPar(GG(g), g)
      if (Parent(Bj,3-kB) == GGG(3-kB)) then
        catB(j) = 5  
        if(catG(g)==0)  catG(g) = 3
        GGP = GGG(kB)
        if (B<0)  UseEE(nA+j) = 0  ! ??
        if (B>0) UseEE(nA+1) = 0
      endif
    enddo
    if (ANY(catA == 8) .and. kB/=3 .and. catB(j)==0) then
      do i=1,nA
        if (PA<0 .and. NFS(AA(i))==0) cycle
        if (Parent(AA(i), 3-kA)>=0) cycle
        if (GpID(3-kB,-Parent(AA(i), 3-kA),3-kA) == Parent(Bj,3-kB)) then
          catB(j) = -i
        endif
      enddo
    endif       
  enddo
endif

if (Complx==1 .and.  (ANY(catA==1) .or. ANY(catA==3))) then 
  LL = NotImplemented  ! explicit consideration of close inbreeding
  deallocate(UseEE)    ! TODO: extend to more?
  deallocate(PrEE)
  deallocate(MateABpar)
  deallocate(TypeEE)
  return
endif

if (kB/=3) then
  GGG = getPar(GG(kB), kB)
  if (GGG(3-kB) == GA(3-kB) .and. GA(3-kB)/=0) then
    catG(kB) = 1
    GGP = GGG(kB)
  endif
endif

MateLoop = .FALSE.
if (B>0) then    ! TODO: B<0   superseded by UseEE -- TODO convert mateloop > UseEE
  do m=1,2
    if (kB/=3 .and. m/=kB)  cycle
    do j=1, nBx(m)
      Bj = BBx(j, m)
      if (nFS(Bj)==0) cycle    
      if (kB==3 .and. Parent(Bj,1)==GG(1) .and.  Parent(Bj,2)==GG(2))  cycle
      if (Parent(Bj,m)<0 .and. Parent(Bj,3-m)<0) then
        do g=1, nS(-Parent(Bj, 3-m),3-m)
          Ei = SibID(g,-Parent(Bj,3-m),3-m)
          if (Parent(Ei,m)>=0 .or. Parent(Ei,m)==Parent(Bj,m)) cycle
          if (ANY(Mates(:,3-m) == Parent(Ei, m))) then
            MateLoop(j,m) = .TRUE.
          endif
        enddo
      endif
    enddo
  enddo
endif

!==============================================

PrL = 0D0        
DoneB = 0
SIMPL = ALL(catA==0) .and. ALL(catG==0) .and. ALL(catB==0) .and. .not. ANY(MateLoop) .and. &
  ALL(GG >=0) .and. all(Parent(AA(1:nA),3-kA) >=0) .and. .not. PAselfed .and. &
  .not. (all(PB < 0) .and. kB==3)
if (SIMPL .and. ANY(UseEE /= 0)) then
  if (ALL(PB >= 0)) then
    SIMPL = .TRUE.
  else
    SIMPL = .FALSE.
  endif
endif
if (SIMPL .and. A>0 .and. PA<0)  SIMPL = .FALSE.   
if (SIMPL .and. B>0 .and. any(PB <0))  SIMPL = .FALSE.

DoAZ = any(catA==2 .or. catA==9 .or. catA==10) .or. catG(kA)==2
if (nA==1 .and. catA(1)==8 .and. PA==0)  catA(1) = 0   ! can't find the bug... 

do l=1,nSnp
  do g=1,2
    if (g/=kB .and. kB/=3) cycle
    if (catG(g)==0) then
      if (SIMPL .and. B>0) then
        call ParProb(l, GG(g), g, B, -1, PrG(:,g)) 
      else
        call ParProb(l, GG(g), g, -1, 0, PrG(:,g))
      endif
    else
      if (GG(g) > 0) then 
        PrG(:,g) = OcA(:,Genos(l,GG(g)))
        if (catG(g)==1) then
          call ParProb(l, GGP, kB, GG(g), 0, PrGG) 
        else if (catG(g)==2) then
          call ParProb(l, GGP, kA, GG(g), 0, PrGG)
        else if (catG(g)==3) then
          call ParProb(l, GGP, kB, GG(g), 0, PrGG)
        endif
      else if (GG(g) < 0) then
        PrG(:,g) = 1D0
        if (catG(g)==1) then
          call ParProb(l, GGP, kB, 0, 0, PrGG)
        else if (catG(g)==2) then
          call ParProb(l, GGP, kA, 0, 0, PrGG)
        else if (catG(g)==3) then
          call ParProb(l, GGP, kB, 0, 0, PrGG)
        endif
      else 
        PrG(:,g) = AHWE(:,l)
      endif
    endif
  enddo
  if (kB/=3) then 
    if (ANY(catA==2) .or. ANY(catB==2)) then
      call ParProb(l, GA(3-kB), 3-kB, -1, 0, PrGA)
    else if (catG(kB)==1 .and. GG(kB)>0) then
      call ParProb(l, GA(3-kB), 3-kB, GG(kB), 0, PrGA)
    else if (PA>0) then
      call ParProb(l, GA(3-kB), 3-kB, PA, 0, PrGA)
    else
      call ParProb(l, GA(3-kB), 3-kB, 0, 0, PrGA)
    endif
    if (B>0) then
      call ParProb(l, PB(3-kB), 3-kB, B, 0, PrPB)   ! -1?
    endif
  endif

  ! === 
  
  PrA = 1D0
  if (SIMPL) then
    if (A < 0) then
      PrA = Xpr(1,:,l, -A,kA)
    else if (A>0) then
      call ParProb(l, Parent(A,3-kA), 3-kA, A, 0, PrPA)
      call Parprob(l, Parent(A,kA), kA, A, -4, PrH)  ! exclude both GPs & A
      do x=1,3 
        PrA(x) = SUM(OKA2P(Genos(l,A),x,:) * PrPA * PrH(x))
      enddo
    endif
  
    do x=1,3  ! PA, kA
      do y=1,3  ! PrG, kB
        do z=1,3  ! PrGA, 3-kB / PrG, 3-kB
          if (kB==3 .and. B>0) then  ! SA/PA FS of B; 
            PrXYZ(x,y,z) = PrA(x) * AKA2P(x,y,z)*PrG(y,3-kA) * PrG(z,kA) * &
              OKA2P(Genos(l,B), y, z)            
          else 
            PrXYZ(x,y,z) = PrA(x) * AKA2P(x,y,z) * PrGA(z)
            if (B>0) then
              PrXYZ(x,y,z) = PrXYZ(x,y,z) * PrG(y, kB) *SUM(OKA2P(Genos(l,B),y,:) * PrPB)
            else if (B<0) then
              PrXYZ(x,y,z) =PrXYZ(x,y,z) *XPr(3,y,l,-B,kB)
            endif
          endif     
        enddo
      enddo
    enddo
    PrL(l) = LOG10(SUM(PrXYZ))   

  else 

    PrAB = 0D0 
    do y=1,3  ! PrG, kB
      PrEE = 0D0
      do x=1,3  ! PA, kA
        if (PAselfed) then
          if (kB==3) then
            PrAB(x,y,y,:) = AKA2P(x,y,y) * PrG(y,kA)
          else
            PrAB(x,y,y,:) = AKA2P(x,y,y) * PrG(y,kB)
          endif
        else
          do z=1,3
            if (kB==3) then
              PrAB(x,y,z,:) = AKA2P(x,y,z) *PrG(z,kA) *PrG(y,3-kA)
            else if (catG(kB)==1) then
              PrAB(x,y,z,:) = AKA2P(x,y,z) * PrGA(z) * &
                SUM(AKA2P(y,z,:) * PrGG) * PrG(y, kB)  
            else if (catG(kB)==2) then
              PrAB(x,y,z,:) = AKA2P(x,y,z)*PrGA(z)
            else
              PrAB(x,y,z,:) = AKA2P(x,y,z) *PrGA(z) * PrG(y,kB) 
            endif
          enddo
        endif     
        if (PA>0) then
          PrAB(x,y,:,:) = PrAB(x,y,:,:) * OcA(x,Genos(l,PA)) 
        endif  
      enddo  

      do x=1,3
       DoneB = 0
       do z=1,3
        if (z>1 .and. .not. DoAZ)  cycle
        do r=1, nA
          if (PA<0 .and. NFS(AA(r))==0) cycle
          if (catA(r)==7) then
            call ParProb(l, GpID(3-kA,-Parent(AA(r),3-kA),3-kA), 3-kA, 0,0,PrH)
            do e=1,3
              PrE(e) = SUM(AKA2P(e,x,:) * PrH)
            enddo
          else if (catA(r)>=8 .and. catA(r)<=10) then
            if (kB < 3) then
              if (catA(r)==8) then
                g=kB
              else if (catA(r)==9) then
                g=3-kB
              endif
            else 
              if (catA(r)==8) then
                g=3-kA
              else if (catA(r)==9) then
                g=kA
              endif
            endif
            if (catA(r) < 10) then
              call ParProb(l, GpID(3-g,-Parent(AA(r),3-kA),3-kA), 3-g, 0,0,PrH)
            endif
            do e=1,3
              if (catA(r)==8) then
                PrE(e) = SUM(AKA2P(e,y,:) * PrH)
              else if (catA(r)==9) then
                PrE(e) = SUM(AKA2P(e,:,z) * PrH)
              else if (catA(r)==10) then
                PrE(e) = AKA2P(e,y,z)
              endif
            enddo  
            PrE = PrE/SUM(PrE)
          else if (catA(r)==42) then
            cycle ! do with B  (catB(j) = -i)
          else if (UseEE(r)/=0) then
            call ParProb(l, MateABpar(r), 3-TypeEE(r), 0,0,PrH)
            do e=1,3
              do u=1, 3
                PrW(u) = SUM(AKA2P(e,u,:) * PrEE(u,x,UseEE(r)) * PrH)
              enddo
              PrE(e) = SUM(PrW)
            enddo
            PrE = PrE/SUM(PrE)
          else if (catA(r) < 0) then
            if (kB==3) then
              PrH = PrG(:,kA)  
            else 
              if (B>0) then
                Bj = BBx(-catA(r),kA)
              else
                Bj = BB(-catA(r))
              endif
              call ParProb(l, Parent(Bj,3-kB), 3-kB, Bj,0,PrH)
            endif
            do e=1,3
              PrE(e) = SUM(AKA2P(e,y,:) * PrH)
            enddo
            PrE = PrE * OcA(:,Genos(l,Bj))
            PrE = PrE/SUM(PrE)
          else if (catA(r)==0 .or. (catA(r)>2 .and. catA(r)<7)) then
            call ParProb(l,Parent(AA(r),3-kA),3-kA,-1,0,PrE)
          else
            PrE = 1D0
          endif

          if (Parent(AA(r), 3-kA) < 0 .and. Parent(AA(r), 3-kA)/=GG(3-kA)) then     
            if (A>0) then
              do i=1, MAX(nFS(AA(r)),1)
                if (FSID(i, AA(r))==A) cycle      
                if (ANY(GG == FSID(i, AA(r)))) cycle
                PrE=PrE*OKA2P(Genos(l,FSID(i,AA(r))),x,:)  ! FS of A
              enddo
            endif
            
            do e=1,3
              if (catA(r)==1 .and. e/=y)  cycle
              if (catA(r)==2 .and. e/=z)  cycle
              if (catA(r)==12 .and. e/=x)  cycle                                  
              do g=1, nS(-Parent(AA(r), 3-kA), 3-kA)
                Ei = SibID(g, -Parent(AA(r), 3-kA),3-kA)
                if (nFS(Ei) == 0) cycle 
                if (nFS(Ei) == 1 .and. (Ei==A .or. Ei==B))  cycle
                if (Parent(Ei,kA)==PA .and. PA/=0) cycle
                if (kA==kB .and. Parent(Ei,kA)==GG(kA) .and. GG(kA)/=0)  cycle
                if (kB<3) then
                  if (Parent(Ei, kB)== PB(kB) .and. PB(kB)/=0) cycle
                endif
                call ParProb(l,Parent(Ei,kA),kA,Ei,-1,PrH)
                if (catA(r)==5 .and. Parent(Ei,kA)==GGP .and. GGP/=0) then  ! FS of GG
                  do i=1, nFS(Ei)
                    if (any(GG == FSID(i, Ei))) cycle
                    if (FSID(i, Ei)==A .or. FSID(i, Ei)==B) cycle 
                    PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)),:,e)
                  enddo
                  if (catG(kA)==2 .and. kB==3) then 
                    PrE(e) = PrE(e) * SUM(AKA2P(z,e,:) * PrH)  
                  else if (ANY(catG==2)) then
                    PrE(e) = PrE(e) * SUM(AKA2P(y,e,:) * PrH)
                  endif                
                else
                  do i=1, nFS(Ei)
                    if (FSID(i, Ei)==A .or. FSID(i, Ei)==B) cycle 
                    PrH=PrH*OKA2P(Genos(l,FSID(i,Ei)),:,e)
                  enddo
                  if (.not. all(PrH==1D0))  PrE(e) = PrE(e) * SUM(PrH)
                endif
              enddo  ! g
            enddo  ! e
            if (catA(r)==3 .and. B>0) then  ! Parent(AA(r), 3-kA) == Parent(Bj, 3-kA)
              do j=1,nBx(kA)
                if (Parent(BBx(j,kA), 3-kA) /= Parent(AA(r), 3-kA)) cycle
                do i=1, MAX(nFS(BBx(j,kA)),1)
                  if (FSID(i,BBx(j,kA))==B) cycle   
                  PrE = PrE * OKA2P(Genos(l,FSID(i,BBx(j,kA))), y, :)
                enddo
              enddo
            endif   
          endif
          if (catA(r)==5 .and. (Parent(AA(r), 3-kA)>0 .or. GGP==0)) then
            do e=1,3
              if (catG(kA)==2 .and. kB==3) then 
                PrE(e) = PrE(e) * SUM(AKA2P(z,e,:) * PrGG)
              else
                PrE(e) = PrE(e) * SUM(AKA2P(y,e,:) * PrGG)
              endif
            enddo
          endif
          
          if (catA(r)==1) then 
            if (DoAZ)        PrAB(x,y,z,1) = PrAB(x,y,z,1) * PrE(y)
            if (.not. DoAZ)  PrAB(x,y,:,1) = PrAB(x,y,:,1) * PrE(y)
          else if (catA(r)==2) then 
            PrAB(x,y,z,1) = PrAB(x,y,z,1) * PrE(z)
          else if (catA(r)==12) then 
            if (DoAZ)        PrAB(x,y,z,1) = PrAB(x,y,z,1) * PrE(x)
            if (.not. DoAZ)  PrAB(x,y,:,1) = PrAB(x,y,:,1) * PrE(x)
          else if (.not. all(PrE==1D0)) then  
            if (DoAZ)        PrAB(x,y,z,1) = PrAB(x,y,z,1) * SUM(PrE)
            if (.not. DoAZ)  PrAB(x,y,:,1) = PrAB(x,y,:,1) * SUM(PrE)
          endif       

          do i=1, MAX(nFS(AA(r)),1)
            if (A>0 .and. FSID(i, AA(r))/=A .and. .not. ANY(BB==FSID(i, AA(r)))) cycle
            PrE = PrE * OKA2P(Genos(l,FSID(i,AA(r))), x, :)  ! <- A
          enddo

          if ((catA(r)==3 .or. (catA(r)==5 .and. ANY(catB==5)) .or. &
           (catA(r)==2 .and. ANY(catB==2))) .and. kB<3) then 
            do j=1, MAX(nB, nBx(kB))
              if (B>0) then
                Bj = BBx(j,kB)
              else
                Bj = BB(j)
              endif
              if (Parent(Bj,3-kA) /= Parent(AA(r),3-kA)) cycle
              if (ANY(AA == Bj))  cycle
              PrE = PrE * OKA2P(Genos(l,Bj), y, :)
              DoneB(j) = 1
            enddo 
          endif
          
          if (catA(r)==1) then 
            if (DoAZ)        PrAB(x,y,z,2) = PrAB(x,y,z,2) * PrE(y)
            if (.not. DoAZ)  PrAB(x,y,:,2) = PrAB(x,y,:,2) * PrE(y)
          else if (catA(r)==2) then 
            PrAB(x,y,z,2) = PrAB(x,y,z,2) * PrE(z)
          else if (catA(r)==12) then 
            if (DoAZ)        PrAB(x,y,z,2) = PrAB(x,y,z,2) * PrE(x)
            if (.not. DoAZ)  PrAB(x,y,:,2) = PrAB(x,y,:,2) * PrE(x)
          else if (.not. all(PrE==1D0)) then    
            if (DoAZ)        PrAB(x,y,z,2) = PrAB(x,y,z,2) * SUM(PrE)
            if (.not. DoAZ)  PrAB(x,y,:,2) = PrAB(x,y,:,2) * SUM(PrE)
          endif
          PrEE(:,x,r) = PrE
        enddo  ! r
       enddo  ! z (if DoAZ)
      enddo  ! x 
      
      do x=1,3
        if (x>1 .and. ALL(UseEE==0) .and. all(catB>=0))  cycle  !  .and. catB(j)/=5
       if (B<0) then            
        do j=1,nB
          if (nFS(BB(j))==0) cycle
          if (DoneB(j)==1) cycle
          if (kA/=kB .and. PA<0 .and. Parent(BB(j), 3-kB)==PA) cycle
!          DoneB(j) = 2  ! for output check only
          if (catB(j)==2 .or. catB(j)==12) then
            PrE = 1D0
          else if (catB(j)==7) then
            call ParProb(l, GpID(3-kB,-Parent(BB(j),3-kB),3-kB), 3-kB, 0,0,PrH)
            do e=1,3
              PrE(e) = SUM(AKA2P(e,y,:) * PrH)
            enddo
          else if (catB(j)==8) then
            call ParProb(l, GpID(3-kA,-Parent(BB(j),3-kB),3-kB), 3-kA, 0,0,PrH)
            do e=1,3
              PrE(e) = SUM(AKA2P(e,x,:) * PrH)
            enddo
          else if (UseEE(nA+j)/=0) then
            call ParProb(l, MateABpar(nA+j), 3-TypeEE(nA+j), 0,0,PrH)
            do e=1,3
              do u=1, 3
                PrW(u) = SUM(AKA2P(e,u,:) * PrEE(u,x,UseEE(nA+j)) * PrH)
              enddo
              PrE(e) = SUM(PrW)
            enddo
            PrE = PrE/SUM(PrE)
          else
            call ParProb(l, Parent(BB(j), 3-kB), 3-kB, -1,0,PrE)
          endif

          if (Parent(BB(j), 3-kB) < 0) then  
            do e=1,3
              do g=1, nS(-Parent(BB(j), 3-kB), 3-kB)
                Ei = SibID(g, -Parent(BB(j), 3-kB),3-kB)
                if (nFS(Ei) == 0) cycle
                if (Parent(Ei, kB) == PB(kB) .and. PB(kB)/=0) cycle  
                if (Parent(Ei, kA)== PA .and. PA/=0) cycle
                call ParProb(l,Parent(Ei,kB),kB,Ei,-1,PrH) 
                do i=1, nFS(Ei)
                  if (FSID(i, Ei)==A) cycle  
                  PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e)
                enddo
                if (.not. all(PrH == 1D0))  PrE(e) = PrE(e) * SUM(PrH) 
              enddo 
            enddo                   
          endif
          
          if (catB(j)==5) then
            do e=1,3
              PrE(e) = PrE(e) * SUM(AKA2P(y,e,:) * PrGG)  ! TODO: FS of SB
            enddo
          endif 
                    
          if (ANY(UseEE/=0) .or. ANY(catB<0)) then
            if (catB(j)==2) then
              PrAB(x,y,:,1) = PrAB(x,y,:,1) * PrE
            else if (.not. all(PrE==1D0)) then
              PrAB(x,y,:,1) = PrAB(x,y,:,1) * SUM(PrE)
            endif
          else if (catB(j)==2) then
            do z=1,3
              PrAB(:,y,z,1) = PrAB(:,y,z,1) * PrE(z)
            enddo
          else if (catB(j)==12) then
            PrAB(:,y,:,1) = PrAB(:,y,:,1) * PrE(y)
!          else if (catB(j)==5 .and. SUM(PrE)<3) then
!            PrAB(x,y,:,1) = PrAB(x,y,:,1) * SUM(PrE)
          else if (.not. all(PrE==1D0)) then
            PrAB(:,y,:,1) = PrAB(:,y,:,1) * SUM(PrE)
          endif           

          do u=1, nFS(BB(j))
            PrE = PrE * OKA2P(Genos(l,FSID(u,BB(j))), y, :)
          enddo                 
          
          if (ANY(UseEE/=0) .or. ANY(catB<0)) then
            if (catB(j)==2) then
              PrAB(x,y,:,2) = PrAB(x,y,:,2) * PrE
            else if (.not. all(PrE==1D0)) then
              PrAB(x,y,:,2) = PrAB(x,y,:,2) * SUM(PrE)
            endif
            PrEE(:,x,nA+j) = PrE
          else if (catB(j)==2) then 
            do z=1,3
              PrAB(:,y,z,2) = PrAB(:,y,z,2) * PrE(z)
            enddo                                                   
          else if (catB(j)==12) then
            PrAB(:,y,:,2) = PrAB(:,y,:,2) * PrE(y)
          else if (.not. all(PrE==1D0)) then
            PrAB(:,y,:,2) = PrAB(:,y,:,2) * SUM(PrE) 
          endif  
        enddo   ! j

      else if (B>0) then 
        do z=1,3  
          do m=1,2
            if (kB/=3 .and. m/=kB)  cycle
            do j=1, nBx(m)
              Bj = BBx(j,m)
              if (nFS(Bj)==0 .and. parent(Bj,m)<0) cycle   
              if (ANY(FSID(:,Bj)==B) .and. DoneB(1)==1)  cycle  
              if (kA/=kB .and. PA<0 .and. Parent(Bj, kA)==PA) cycle
              if (kB==3 .and. Parent(Bj, 3-m)==GG(3-m) .and. GG(3-m)/=0) then  ! FS of B
                if (Parent(Bj,1)<0 .and. Parent(Bj,2)<0 .and. m==2) cycle
 !               DoneB(1) = 2  ! for output check only
                do u=1, Max(nFS(Bj), 1)
                  if (FSID(u,Bj)==B) cycle  
                  if (ANY(AA(1:nA)==FSID(u,Bj))) cycle
                  if (ALL(UseEE==0)) then
                    PrAB(:,y,z,:) = PrAB(:,y,z,:) * OKA2P(Genos(l,FSID(u,Bj)),y,z) 
                  else
                    PrAB(x,y,z,:) = PrAB(x,y,z,:) * OKA2P(Genos(l,FSID(u,Bj)),y,z)
                  endif
                enddo
              else if (kB==3 .and. Parent(Bj,1)<0 .and. Parent(Bj,2)<0 &
               .and. MateLoop(j,m)) then  
                if (m==2) cycle
                call ParProb(l,Parent(Bj,m),m,-1,0,PrE)
                call ParProb(l,Parent(Bj,3-m),3-m,-1,0,PrW)
                
                do g=1, nS(-Parent(Bj, 3-m),3-m)
                  Ei = SibID(g,-Parent(Bj,3-m),3-m)
                  if (nFS(Ei) == 0) cycle
                  do e=1,3
                    do w=1,3
                      PrEW(e,w) = PrE(e) * PrW(w)
                      if (Parent(Ei,m)==Parent(Bj,m) .and. &
                       Parent(Ei,3-m)==Parent(Bj,3-m)) then
                        do i=1, nFS(Ei)
                          if (FSID(i, Ei)==B) cycle
                          PrEW(e,w) = PrEW(e,w) * OKA2P(Genos(l,FSID(i,Ei)), e, w)
                        enddo
                      else if (Parent(Ei,m)==Parent(Bj,m)) then
                        call ParProb(l, Parent(Ei,3-m),3-m, Ei, -1, PrH)
                        do i=1, nFS(Ei)
                          if (FSID(i, Ei)==B) cycle  
                          PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e)
                        enddo
                        if (.not. all(PrH==1D0))  PrEW(e,w) = PrEW(e,w) * SUM(PrH)
                      else if (Parent(Ei,3-m)==Parent(Bj,3-m)) then
                        call ParProb(l, Parent(Ei,m),m, Ei, -1, PrH)
                        do i=1, nFS(Ei)
                          if (FSID(i, Ei)==B) cycle  
                          PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, w)
                        enddo
                        if (.not. all(PrH==1D0))  PrEW(e,w) = PrEW(e,w) * SUM(PrH)
                      endif
                    enddo  ! w
                  enddo  ! e
                enddo  ! sib g
                if (ALL(UseEE==0)) then
                  PrAB(:,y,z,:) = PrAB(:,y,z,:) * SUM(PrEW)
                else
                  do e=1,3
                    PrEE(:,x,nA+1) = SUM(PrEW(e,:))
                  enddo
                  PrAB(x,y,z,:) = PrAB(x,y,z,:) * SUM(PrEW)
                endif

              else
                if (kB==3) then
                  call ParProb(l,Parent(Bj,3-m),3-m,-1,0,PrE)
                else if (catB(j)==2 .or. catB(j)==12) then
                  PrE = 1D0                         
                else if (catB(j)==7) then
                  call ParProb(l, GpID(3-kB,-Parent(Bj,3-kB),3-kB), 3-kB, 0,0,PrH)
                  do e=1,3
                    PrE(e) = SUM(AKA2P(e,y,:) * PrH)
                  enddo              
                else if (catB(j)==8) then
                  call ParProb(l, GpID(3-kA,-Parent(Bj,3-kB),3-kB), 3-kA, 0,0,PrH)
                  do e=1,3
                    PrE(e) = SUM(AKA2P(e,x,:) * PrH)
                  enddo
                else if (UseEE(nA+1)/=0 .and. Bj==B) then  
                  call ParProb(l, MateABpar(nA+1), 3-TypeEE(nA+1), 0,0,PrH)
                  do e=1,3
                    do u=1, 3
                      PrW(u) = SUM(AKA2P(e,u,:) * PrEE(u,x,UseEE(nA+1)) * PrH)
                    enddo
                    PrE(e) = SUM(PrW)
                  enddo
                  PrE = PrE/SUM(PrE)
                else
                  call ParProb(l,Parent(Bj,3-m),3-m,-1,0,PrE)
                endif
                
                if (Parent(Bj,3-m)<0 .and. Parent(Bj,3-m)/=GG(3-m)) then
                  do g=1, nS(-Parent(Bj, 3-m),3-m)
                    Ei = SibID(g,-Parent(Bj,3-m),3-m)
                    if (nFS(Ei) == 0) cycle
                    if (ANY(AA(1:nA)==Ei)) cycle
                    if (Parent(Ei,m)==GG(m) .and. GG(m)/=0) cycle
                    do e=1,3
                      call ParProb(l, Parent(Ei, m), m, Ei, -1, PrH)
                      do i=1, nFS(Ei)
                        if (ANY(AA(1:nA)==FSID(i, Ei))) cycle
                        if (FSID(i, Ei)==B) cycle   
                        if (FSID(i, Ei)==PA) cycle
                        PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e)
                      enddo
                      if (.not. all(PrH==1D0))  PrE(e) = PrE(e) * SUM(PrH)
                    enddo
                  enddo 
                endif
                do u=1, MAX(nFS(Bj),1)
                  if (FSID(u,Bj)==B) cycle 
                  if (FSID(u, Bj)==PA) cycle 
                  if (ANY(AA(1:nA)==FSID(u,Bj))) cycle
                  if (kB==3 .and. m==kA) then
                    PrE = PrE * OKA2P(Genos(l,FSID(u,Bj)), z, :) 
                  else
                    PrE = PrE * OKA2P(Genos(l,FSID(u,Bj)), y, :) 
                  endif
                enddo
                if (catB(j)==5 .and. kB/=3) then
                  do e=1,3
                    PrE(e) = PrE(e) * SUM(AKA2P(y,e,:) * PrGG)
                  enddo
                endif
                
                if (kB<3) then
                  if (ANY(UseEE/=0) .or. ANY(catB<0)) then
                    if (catB(j)==2) then 
                      PrAB(x,y,z,1) = PrAB(x,y,z,1) * PrE(z)
                    else if (SUM(PrE)<3.0) then
                      PrAB(x,y,z,1) = PrAB(x,y,z,1) * SUM(PrE)
                    endif
                  else if (catB(j)==2) then 
                    PrAB(:,y,z,1) = PrAB(:,y,z,1) * PrE(z)
                  else if (catB(j)==12) then 
                    PrAB(:,y,z,1) = PrAB(:,y,z,1) * PrE(y)
                  else if (SUM(PrE)<3.0) then
                    PrAB(:,y,z,1) = PrAB(:,y,z,1) * SUM(PrE)
                  endif
                  
                  do u=1, MAX(nFS(Bj),1)
                    if (FSID(u,Bj)==B) then
                      PrE = PrE * OKA2P(Genos(l,B), y, :) 
                    endif
                  enddo
                  
                  if (ANY(UseEE/=0) .or. ANY(catB<0)) then
                    if (catB(j)==2) then 
                      PrAB(x,y,z,2) = PrAB(x,y,z,2) * PrE(z)
                    else if (SUM(PrE)<3.0) then
                      PrAB(x,y,z,2) = PrAB(x,y,z,2) * SUM(PrE)
                    endif
                    if (any(FSID(1:nFS(Bj),Bj)==B))  PrEE(:,x,nA+1) = PrE
                  else if (catB(j)==2) then 
                    PrAB(:,y,z,2) = PrAB(:,y,z,2) * PrE(z)
                  else if (catB(j)==12) then 
                    PrAB(:,y,z,2) = PrAB(:,y,z,2) * PrE(y)
                  else if (SUM(PrE)<3.0) then
                    PrAB(:,y,z,2) = PrAB(:,y,z,2) * SUM(PrE)
                  endif
                else if (SUM(PrE)<3.0) then  ! kB==3
                  PrAB(:,y,z,:) = PrAB(:,y,z,:) * SUM(PrE)
                endif
              endif
            enddo  ! j_m
          enddo  ! m
          if (kB==3) then
            PrAB(:,y,z,2) = PrAB(:,y,z,2) * OKA2P(Genos(l,B), y, z)
          endif
        enddo  ! z
      endif  ! B <0 vs B>0
      enddo  ! x (ANY(UseEE/=0)  only)
    enddo  ! y

    PrL(l) = LOG10(SUM(PrAB(:,:,:,2))) - LOG10(SUM(PrAB(:,:,:,1)))
  endif
enddo
LL = SUM(PrL)

deallocate(UseEE)
deallocate(PrEE)
deallocate(MateABpar)
deallocate(TypeEE)

end subroutine PairUA 

! #####################################################################

subroutine addFA(A, SB, kB, LL)  ! SB & partner-of-SB are GP's of A
use Global
implicit none

integer, intent(IN) :: A,SB,kB
double precision, intent(OUT) :: LL
integer :: x, y, z, mateofSB, i, l, Bi
double precision :: PrL(nSnp), PrXYZ(3,3,3), PrY(3), PrZ(3), PrPA(3), & ! LLU, &
   ALRq, PrPAx(3)
logical :: AncOK, Ainbr    

LL = missing
if (Parent(A,kB) == -SB)  LL = Impossible
if (hermaphrodites/=0)  LL = NotImplemented
if ((any(SibID(:,SB,kB) == Parent(A,1)) .and. Parent(A,1)/=0) .or. &
  (any(SibID(:,SB,kB) == Parent(A,2)) .and. Parent(A,2)/=0))  LL = AlreadyAss   
if (LL /= missing)  return 

AncOK = .TRUE.
call ChkAncest(-SB,kB, A,0, AncOK)
if (.not. AncOK  .or. Parent(A,kB)<0) then 
  LL = NotImplemented
  return
endif

mateofSB = 0
if (ns(SB,kB)>0) then
  call getFSpar(SB, kB, .TRUE., mateofSB)  ! TODO: strict=FALSE ?  needs check if Parent(B1,3-k)==mateofSB
  if (mateofSB == 0 .or. (any(Parent(SibID(1:ns(SB,kB),SB,kB), 3-kB) /= mateofSB .and. &
    Parent(SibID(1:ns(SB,kB),SB,kB), 3-kB) /= 0))) then  ! NOTE: all sibs made FS as side-effect
    LL = impossible
  else if (mateofSB/=0 .and. Parent(A,3-kB)==mateofSB) then  ! can't be both parent & GP
    LL = impossible  
  else
    call CalcAgeLR(A,Sex(A), mateofSB,3-kB, kB,4, .TRUE., ALRq) 
    if (ALRq < 2.0*TF .or. ALRq==impossible)  LL = impossible
    if (LL /= impossible) then
      call ChkAncest(mateofSB, 3-kB, A, 0, AncOK)
      if (.not. AncOK)  LL = impossible
    endif
  endif
  if (LL == impossible)  return
endif

do i=1, ns(SB, kB)
  if (nFS(SibID(i,SB,kB))==0)  cycle                                  
  if (Parent(SibID(i,SB,kB), 3-kB) == mateofSB) then
    Bi = SibID(i,SB,kB)
    exit
  endif
enddo

if (Parent(A,3-kB)/=0 .and. Parent(A,3-kB)== mateofSB) then
  Ainbr = .TRUE.  ! A would become inbred
else
  Ainbr = .FALSE.
endif

PrL = 0D0
do l=1, nSnp
  call ParProb(l, -SB, kB, -1, 0, PrY)
  call ParProb(l, mateofSB, 3-kB, Bi, -1, PrZ)   
  call ParProb(l, Parent(A,3-kB), 3-kB, A, 0, PrPA)
  call ParProb(l, Parent(A,kB), kB, A, -4, PrPAx)       
  do x=1,3
    do y=1,3
      do z=1,3
        PrXYZ(x,y,z) = PrY(y) * PrZ(z) * PrPAx(x) * AKA2P(x,y,z)  
!        PrXYZ(x,y,z,2) = PrY(y) * PrZ(z) * PrPAx(x) * AHWE(x,l)   !! NO: parent(A,k) may already have parents!
        if (Ainbr) then
          PrXYZ(x,y,z) = PrXYZ(x,y,z) *  OKA2P(Genos(l,A), x, z)
        else
          PrXYZ(x,y,z) = PrXYZ(x,y,z) *  SUM(OKA2P(Genos(l,A), x, :) * PrPA)
        endif
        do i=1, ns(SB,kB)
          PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l, SibID(i,SB,kB)), y,z)
        enddo  ! NOT FSLik: some siblings may have parent(3-kB)=0
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXYZ))
enddo

LL = SUM(PrL)
if (Parent(A,kB)>0)  LL = LL - Lind(Parent(A,kB))

! LLU = missing
! if (Parent(A,kB) /= 0 .or. Ainbr) then
  ! call CalcU(-SB, kB, A, Sex(A), LLU)
  ! LL = SUM(PrL(:,1)) - SUM(PrL(:,2)) + LLU 
! else
  ! LL = SUM(PrL(:,1))
! endif                    

end subroutine addFA

! #####################################################################

 subroutine pairFAx(A, B, LL)  
! A result of FS mating, B offspr of A's double GPs  (complex = 'mono')
use Global
implicit none

integer, intent(IN) :: A, B
double precision, intent(OUT) :: LL
integer :: l, x, y, v, w
double precision :: PrL(nSnp), PrXY(3,3,3,3)

LL = Missing
if (any(Parent(A,:)/=0 .or. any(parent(B,:)/=0))) then
  LL = NotImplemented
  return
endif

PrL = 0D0
do l =1, nSnp
  do x=1,3  ! dam of A
    do y=1,3  ! sire of A
      do v=1,3  ! grandmother
        do w=1,3  ! grandfather
          PrXY(x,y,v,w) = OKA2P(Genos(l, A), x, y) * AKA2P(x,v,w) * AKA2P(y,v,w) * &
            OKA2P(Genos(l,B),v,w) * AHWE(w,l) * AHWE(v,l)
        enddo
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXY))
enddo

LL = SUM(PrL)

end subroutine pairFAx

! #####################################################################

subroutine HSmating(A, kA, B, kB, hf, LL)   ! SB/PB and HS-of-PB are parents of SA/PA
use Global
implicit none

integer, intent(IN) :: A, kA, B, kB, hf
double precision, intent(OUT) :: LL
integer :: PA, PB, GA(2), GB(2), GGA(2), l, x, y, z, v, DoQuick, GGP  ! , m
double precision :: PrL(nSnp,2), PrXZ(3,3,3,3,2), PrA(3), PrPA(3), PrGA(3), &
  PrB(3), PrPB(3), PrGB(3), PrGGA(3), PrGGP(3), LLU

if (A > 0) then
  PA = Parent(A, kA)
else
  PA = A
endif
if (B > 0) then
  PB = Parent(B, kB)
else
  PB = B
endif

GA = getPar(PA, kA)
GGA = getPar(GA(3-kB), 3-kB)                            
GB = getPar(PB, kB)

LL = missing
if (PA > 0 .or. PB > 0) then
  LL = NotImplemented
else if (GA(kB) == PB .and. PB/=0) then
  LL = AlreadyAss
else if (GA(kB)/=0) then
  LL = impossible
else if (GB(kA) == PA .and. PA/=0) then
  LL = impossible
endif
if (LL /= Missing)  return

if (GGA(hf)/=0) then
  if (GB(hf)==0 .or. GB(hf)==GGA(hf)) then
    GGP = GGA(hf)
  else
    LL = impossible
  endif
else
  GGP = GB(hf)
endif

if (LL /= Missing)  return

if (PA==0 .or. (PB==0 .and. GGP==0) ) then
  LL = NotImplemented   ! causes too many false neg 
  return                ! (indistinguishable from just PA inbred)
endif

DoQuick = 1
if (PA < 0) then
  call ChkDoQuick(-PA, kA, DoQuick)
  if (DoQuick /=1)  LL = NotImplemented
endif
if (PB < 0) then
  call ChkDoQuick(-PB, kB, DoQuick)
  if (DoQuick /=1)  LL = NotImplemented
endif  
if (LL /= Missing)  return

PrL = 0D0
do l=1, nSnp
  call ParProb(l, GA(3-kB), 3-kB, -4, 0, PrGA)  ! offspring contribution only
  if (GA(3-kB) > 0) then
    if (PB > 0) then
      call ParProb(l, GGP, hf, GA(3-kB), PB, PrGGP) 
    else
      call ParProb(l, GGP, hf, GA(3-kB), 0, PrGGP) 
    endif
  else
    if (PB > 0) then
      call ParProb(l, GGP, hf, 0, PB, PrGGP) 
    else
      call ParProb(l, GGP, hf, 0, 0, PrGGP) 
    endif
  endif
  call ParProb(l, GB(3-hf), 3-hf, 0, 0, PrGB)
  call ParProb(l, GGA(3-hf), 3-hf, 0, 0, PrGGA)
  if (A < 0) then
    PrA = XPr(1,:,l,-A,kA)
  else 
    call ParProb(l, PA, kA, A, -4, PrA) ! exclude both GPs & A
    call ParProb(l, Parent(A,3-kA), 3-kA, A, 0, PrPA)
    do x=1,3
      PrA(x) = PrA(x) * SUM(OKA2P(Genos(l,A), x, :) * PrPA)
    enddo
  endif
  if (B < 0) then
    PrB = XPr(1,:,l,-B,kB)
  else 
    call ParProb(l, PB, kB, B, -4, PrB)
    call ParProb(l, Parent(B,3-kB), 3-kB, B, 0, PrPB)
    do y=1,3
      PrB(y) = PrB(y) * SUM(OKA2P(Genos(l,B), y, :) * PrPB)
    enddo
  endif
  
  PrXZ = 0D0
  do x=1,3  ! SA/PA
    do y=1,3  ! SB
      do z=1,3  ! HS & mate of SB
        do v=1,3   ! parent of SB & mate-of-SB
          PrXZ(x,y,z,v,1) = PrA(x) * PrB(y) * AKA2P(x,y,z) * PrGA(z) * PrGGP(v) * &
            SUM(AKA2P(y,v,:) * PrGB) * SUM(AKA2P(z,v,:) * PrGGA)   
         ! SA similarly inbred, but B unrelated
         if (PB==0) then  
           PrXZ(x,y,z,v,2) = PrA(x) * AKA2P(x,y,z) * PrGA(z) * PrGGP(v) * &
              AKAP(y,v,l) * SUM(AKA2P(z,v,:) * PrGGA)
          endif
         
        enddo
      enddo
    enddo
  enddo
  PrL(l,1) = LOG10(SUM(PrXZ(:,:,:,:,1)))
  PrL(l,2) = LOG10(SUM(PrXZ(:,:,:,:,2)))
enddo

if (PB==0) then  ! B>0
  call CalcU(A, kA, B, kB, LLU)
  LL = SUM(PrL(:,1)) - (SUM(PrL(:,2)) + Lind(B)) + LLU
else
  LL = SUM(PrL(:,1))
endif

end subroutine HSmating

! #####################################################################

subroutine addHAselfed(SA,kA, SB,kB, LL)  ! SA result of selfing by SB
use Global
implicit none

integer, intent(IN) :: SA, kA, SB, kB
double precision, intent(OUT) :: LL
integer :: l, x, y
double precision :: PrL(nSnp), PrXY(3,3)

if (any(GpID(:,SA,kA) /= 0)) then
  LL = impossible
  return
endif

if (kA/=kB) then
  if (any(parent(SibID(1:ns(SA,kA),SA,kA), 3-kA) == -SB)) then
    LL = NotImplemented
    return
  endif
endif
! TODO: check if otherwise connected

PrL = 0D0
do l=1,nSnp
  do x=1,3
    do y=1,3
      PrXY(x,y) = XPr(1,x,l,SA,kA) * AKA2P(x,y,y) * XPr(3,y,l,SB,kB)
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXY))
enddo
LL = SUM(PrL)

end subroutine addHAselfed

! #####################################################################

subroutine pairCC(A,B,k, LL)  ! full 1st cousins
use Global
implicit none

integer, intent(IN) :: A,B,k
double precision, intent(OUT) :: LL
integer :: l, x, y, u, v, z
double precision :: PrL(nSnp, 5), PrXY(3,3), PrPA(3), PrPB(3), PrXYh(3,3), &
  PrC(3,3,5), PrZ(3), PrXYf(3,3,2,2), LLself(2), LLU, LLtmp(3), PrX(3), PrY(3)
logical :: AncOK(2), AreHS, MaybeInbr(2) 
  
LL = missing
if (Parent(A,k)/=0 .and. Parent(B,k)/=0) then
  if (Parent(A,k)==Parent(B,k)) then
    LL = impossible
    return
  endif
endif

if (Complx==0) then
  if (Sex(A)<3) then
    if (Mate(A)/=0 .and. Parent(B,3-sex(A))==Mate(A)) LL=impossible
  endif
  if (Sex(B)<3) then
    if (Mate(B)/=0 .and. Parent(A,3-sex(B))==Mate(B)) LL=impossible
  endif
  if (LL==impossible)  return
endif

if (Parent(A,k)<0 .or. Parent(B,k)<0) then
  LL = NotImplemented  ! covered by ParentHFS
  return
endif

AncOK = .TRUE.
call ChkAncest(A,0,B,0, AncOK(1))
call ChkAncest(B,0,A,0, AncOK(2))
if (any(.not. AncOK)) then
  LL = impossible
  return
endif

AreHS = .FALSE.
MaybeInbr = .TRUE.                 
if (Parent(A,3-k)/=0 .and. Parent(A,3-k)==Parent(B,3-k)) then
  AreHS = .TRUE.
endif
if (hermaphrodites>0) then
  LLself = missing
  if (Parent(A,3-k)>0) then
    call PairSelf(B, Parent(A,3-k), LLself(1))
    call PairFullSib(B, Parent(A,3-k), LLself(2))
    if (LLself(1) > LLself(2)) then
      MaybeInbr(1) = .FALSE.
    endif
  endif
  if (Parent(B,3-k)>0) then
    call PairSelf(A, Parent(B,3-k), LLself(1))
    call PairFullSib(A, Parent(B,3-k), LLself(2))
    if (LLself(1) > LLself(2)) then
      MaybeInbr(2) = .FALSE.
    endif
  endif
endif

PrL = 0D0
PrXYf = 0D0  ! 4D: x, y, A/B inbred, CC/not
PrXYh = 0d0
do l=1, nSnp
  if (.not. AreHS) then
    call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPA)
    call ParProb(l, Parent(B,3-k), 3-k, B, 0, PrPB)
  else
    call ParProb(l, Parent(A,3-k), 3-k, A, B, PrPA)
  endif
  PrX = 1d0
  PrY = 1d0
  if (Parent(A,k)>0)  PrX = OcA(:,Genos(l,Parent(A,k)))
  if (Parent(B,k)>0)  PrX = OcA(:,Genos(l,Parent(B,k)))                              
  
  PrC = 0D0
  do u=1,3  ! GG 3-k
    do v=1,3  ! GG k
      do x=1,3  !PA
        do y=1,3    !PB
          PrXY(x,y) = AKA2P(x,u,v) * AKA2P(y,u,v) * AHWE(u,l) * AHWE(v,l) * PrX(x) * PrY(y)
          if (.not. AreHS) then
            if (any(MaybeInbr)) then
              PrXYf(x,y,1,1) = PrXY(x,y) * PrPA(u) / AHWE(u,l) 
              PrXYf(x,y,2,1) = PrXY(x,y) * PrPB(u) / AHWE(u,l)
              ! A or B inbred, not CC (reference)
              PrXYf(x,y,1,2) =  AKA2P(x,u,v) * AHWE(v,l) * AHWE(y,l) * PrPA(u) 
              PrXYf(x,y,2,2) =  AKA2P(y,u,v) * AHWE(v,l) * AHWE(x,l) * PrPB(u)
            endif
            PrXY(x,y) = PrXY(x,y) * SUM(OKA2P(Genos(l,A), x, :) * PrPA) * &
              SUM(OKA2P(Genos(l,B), y, :) * PrPB)
            if (any(MaybeInbr)) then
              PrXYf(x,y,1,:) = PrXYf(x,y,1,:) * OKA2P(Genos(l,A), x, u) * &
                SUM(OKA2P(Genos(l,B), y, :) * PrPB)  ! A inbred
              PrXYf(x,y,2,:) = PrXYf(x,y,2,:) * SUM(OKA2P(Genos(l,A), x, :) * PrPA) * &
                OKA2P(Genos(l,B), y, u) 
              ! not considered: both A & B inbred. 
            endif
          else if (AreHS) then
            do z=1,3
              PrZ(z) = PrPA(z) * OKA2P(Genos(l,A),x,z) * OKA2P(Genos(l,B),y,z)
            enddo
            PrXY(x,y) = PrXY(x,y) * SUM(PrZ)
            ! HS+HC
            PrXYh(x,y) = AKAP(x,u,l) * AKAP(y,u,l) * AHWE(u,l) * PrX(x) * PrY(y)
            PrXYh(x,y) = PrXYh(x,y) * SUM(PrZ)
          endif
        enddo
      enddo
      PrC(u,v,1) = SUM(PrXY)
      if (.not. AreHS) then
        if (MaybeInbr(1))  PrC(u,v,2) = SUM(PrXYf(:,:,1,1))
        if (MaybeInbr(2))  PrC(u,v,3) = SUM(PrXYf(:,:,2,1))
        if (MaybeInbr(1))  PrC(u,v,4) = SUM(PrXYf(:,:,1,2))
        if (MaybeInbr(2))  PrC(u,v,5) = SUM(PrXYf(:,:,2,2))
      else if (u==v) then
        PrC(u,v,2) = SUM(PrXYh)
      endif
    enddo
  enddo 
  do z=1,5
    PrL(l,z) = LOG10(SUM(PrC(:,:,z)))
  enddo
enddo

LLtmp = 999D0
LLtmp(1) = SUM(PrL(:,1))
if (areHS) then
  LLtmp(2) = SUM(PrL(:,2))
else if (any(MaybeInbr)) then
  call CalcU(A, k, B, k, LLU)
  if (MaybeInbr(1))  LLtmp(2) = SUM(PrL(:,2)) - SUM(PrL(:,4)) + LLU
  if (MaybeInbr(2))  LLtmp(3) = SUM(PrL(:,3)) - SUM(PrL(:,5)) + LLU
endif
LL = MaxLL(LLtmp)

end subroutine pairCC

! #####################################################################

subroutine pairHSCC(A,B,k,LL)  ! HS via k + full 1st cousins via 3-k
use Global
implicit none

integer, intent(IN) :: A, B, k
double precision, intent(OUT) :: LL
integer :: x, y, z, v, w, l, curpar(2,2), curGP(2,2), joinpar, joinGP(2)
double precision :: PrL(nSnp), PrAB(3,3,3,3,3), PrX(3), PrY(3), PrZ(3), PrV(3), PrW(3)
logical :: ParOK

curpar(:,1) = Parent(A,:)
curpar(:,2) = Parent(B,:)
if (curPar(k,1)/=0 .and. curPar(k,2)/=0 .and. curPar(k,1)/=curPar(k,2)) then
  LL = impossible   ! different parents --> cannot be HS
else if (curPar(k,1)/=0) then
  joinpar = curPar(k,1)
  call ChkValidPar(B,sex(B), joinpar,k, ParOK)
  if (.not. ParOK)  LL = impossible
else
  joinpar = curPar(k,2)
  call ChkValidPar(A,sex(A), joinpar,k, ParOK)
  if (.not. ParOK)  LL = impossible
endif
if (LL == impossible)  return

curGP = 0
curGP(:,1) = getPar(curPar(3-k,1), 3-k)
curGP(:,2) = getPar(curPar(3-k,2), 3-k)
do x=1,2
  if (curGP(x,1)/=0 .and. curGP(x,2)/=0 .and. curGP(x,1)/=curGP(x,2)) then
    LL = impossible   ! different GPs at 3-k side: cannot be full cousins
    return
  else if (curGP(x,1)/=0) then
    joinGP(x) = curGP(x,1)
  else
    joinGP(x) = curGP(x,2)
  endif
enddo

WHERE (curPar < 0)  curPar = 0

PrL = 0D0
do l=1, nSnp
  call ParProb(l, Parent(A,3-k), 3-k, A, -4, PrX)  ! excl A & GPs
  call ParProb(l, Parent(B,3-k), 3-k, B, -4, PrY)
  call ParProb(l, joinPar, k, A, B, PrZ) 
  call ParProb(l, joinGP(1), 1, curPar(3-k,1), curPar(3-k,2), PrV)   ! excl if curpar>0  
  call ParProb(l, joinGP(2), 2, curPar(3-k,1), curPar(3-k,2), PrW)
  do x=1,3
    do y=1,3
      do z=1,3
        do v=1,3
          do w=1,3
            PrAB(x,y,z,v,w) = PrX(x) * PrY(y) * AKA2P(x,v,w) * AKA2P(y,v,w) * PrV(v) * PrW(w) * &
              PrZ(z) * OKA2P(Genos(l,A),x,z) * OKA2P(Genos(l,B),y,z)
          enddo
        enddo
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrAB))
enddo
LL = SUM(PrL)

end subroutine pairHSCC

! #####################################################################

subroutine pairDHC(A, kA, B, withFS, LL)   ! B double half cousin of A
use Global
use CalcLik           
implicit none

integer, intent(IN) :: A, kA, B
logical, intent(IN) :: withFS
double precision, intent(OUT) :: LL
integer :: x, y, z, v, w, ParB(2), m, l
double precision :: PrL(nSnp), PrX(3,3,3,3,3), PrPA(3), PrB(3,3)

if (B<0) then
  LL = NotImplemented
  Return
endif

parB = getPar(B,0)
LL = Missing
if (Parent(A,kA)/=0 .or. any(parB>0)) then
  LL = NotImplemented
else  
  do m=1,2
    if (ParB(m) < 0) then
      if (any(GpID(:,-ParB(m),m)/=0)) then
        LL = NotImplemented
      endif
    endif
  enddo
endif
if (LL == NotImplemented)  return
! also assumed that B has no additional siblings, asside from full sibs

PrL = 0D0
do l =1, nSnp
  call ParProb(l, Parent(A,3-kA), 3-kA, A, 0, PrPA)
  if (withFS) then  
    PrB = FSLik(l,B)
  else
    PrB = OKA2P(Genos(l,B),:,:)
  endif                
  do x=1,3  ! parent of A
    do y=1,3  ! grandparent 1
      do z=1,3  ! grandparent 2
        do v=1,3  ! dam of B
          do w=1,3  ! sire of B  
            PrX(x,y,z,v,w) = SUM(OKA2P(Genos(l,A), x, :) * PrPA) * AKA2P(x,y,z) * &
              AKAP(v,y,l) * AKAP(w,z,l) * AHWE(y,l) * AHWE(z,l) * PrB(v,w)
          enddo
        enddo
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrX))
enddo

LL = SUM(PrL)

end subroutine pairDHC

! #####################################################################

subroutine Clustering(PairID, PairType)
use Global
implicit none

integer, intent(IN) :: PairID(XP*nInd,2), PairType(XP*nInd)
integer :: k, x, n, m, ij(2), sx(2), topX, u, fcl, Par(2), topFS, topXi, chunk_printdot(10)
double precision :: LL(7,2), dLL, LLx(7, 2,2), dLLtmp(maxSibSize), dLLi
logical :: IsPair, FSM, DoLater 

chunk_printdot = mk_seq(nPairs,10)
if (quiet==-1 .and. nPairs==0)  call print_dot(10)

do x=1, nPairs
  if (MODULO(x,200)==0) call rchkusr()
  if (quiet==-1 .and. any(chunk_printdot==x)) call print_dot(count(chunk_printdot==x))
  LL = missing
  ij = PairID(x,:)

  do k=1,2
    if (k/=PairType(x) .and. PairType(x)/=3)  cycle                             
    if (any(Parent(ij,k)>0)) cycle
    if (Parent(ij(1),k) == Parent(ij(2),k) .and. Parent(ij(1),k) /= 0) cycle
    if (hermaphrodites==1 .and. ANY(parent(ij,3-k)>0) .and. ANY(parent(ij,k)==0))  cycle
    fcl = 0
    LLx = missing 
    call getFocal(ij(1), ij(2), 0, k, fcl)
    sx(1) = -Parent(ij(1),k)  
    sx(2) = -Parent(ij(2),k)
    
    if (k==1 .and. any(sx==0) .and. any(parent(ij,k)<0) .and. all(parent(ij,3-k)<0)) then
      DoLater = .FALSE.
      do n=1,2
        if (sx(n)==0)  cycle
        if (all(parent(SibID(1:ns(sx(n),k),sx(n),k), 3-k) < 0)) then
          call getFSpar(sx(n), k, .TRUE.,par(1))
          if (par(1) < 0) then
            if (nS(-par(1), 3-k) == nS(sx(n),k)) DoLater = .TRUE.  !else possibly add on wrong side i.o. merge
          endif 
        endif
      enddo
      if (DoLater)  cycle
    endif              
    
    if (sx(1)==0 .and. sx(2)==0) then
      if (AgeDiff(ij(1),ij(2))==missing) then
        call CheckRel(ij(1), k, ij(2), k, fcl, LLx(:,1,1), LLx(:,1,2))  
        call CheckRel(ij(2), k, ij(1), k, fcl, LLx(:,2,1), LLx(:,2,2))
        do u=1,7
          do n=1,2
            LL(u,n) = MaxLL(LLx(u,:,n))
          enddo
        enddo
      else if(AgeDiff(ij(1),ij(2))>=0) then
        call CheckRel(ij(1), k, ij(2), k, fcl, LL(:,1), LL(:,2))
      else
        call CheckRel(ij(2), k, ij(1), k, fcl, LL(:,1), LL(:,2))
      endif
      
      IsPair = .TRUE.
      do n=1,2  
        if (AgePhase==0 .and. n==2)  cycle
        if (AgePhase==2 .and. n==1)  cycle                                  
        if (LL(2,n)>0 .and. LL(3,n)>0) then
          IsPair = .FALSE.
          exit
        endif
        topX = 0
        dLL = 0D0
        call BestRel(LL(:,n), fcl, topX, dLL)        
        if (fcl==3 .and. (topX==2 .or. topX==3)) then
          IsPair = .TRUE.
        else if (fcl==2 .and. topX==2 .and. dLL > TA .and. &
         (LL(2,n) - MaxLL(LL((/1,4,5,6,7/),n)) > 2*TA .or. Complx==0 .or. &
           ALL(Parent(ij, 3-k)/=0))) then   
          IsPair = .TRUE.           
        else if (AgePhase==2 .and. n==1 .and. Complx>0) then   
          ! still do a basic check on non-age-dependend
          if (MaxLL(LL((/2,3/),n)) - MaxLL(LL((/1,6,7/),n))>TA .and. &
            MaxLL(LL((/4,5/),n)) - MaxLL(LL((/2,3/),n)) < TA) then
              IsPair = .TRUE.
          else
            IsPair = .FALSE.
            exit
          endif
        else
          IsPair = .FALSE.
          exit
        endif
      enddo
      if (.not. IsPair) cycle
      
      call NewSibship(ij(1), ij(2), k)   ! new sibship (pair)  
      if (fcl==2 .or. (topX==2 .and. dLL>2*TA)) then
        if (ALL(Parent(ij, 3-k)==0)) then  ! another new sibship
          call NewSibship(ij(1), ij(2), 3-k)
        else
          do m=1,2
            if (Parent(ij(m),3-k)/=0) then
              call setPar(ij(3-m), 3, Parent(ij(m),3-k), 3-k)
            endif
          enddo
        endif
      endif
      
    else if (sx(1)>0 .and. sx(2)>0 .and. sx(1) /= sx(2)) then
      IsPair = .TRUE.
      FSM = .FALSE.
      call CheckMerge(sx(1), sx(2), k,k, 1, LL(:,1), LL(:,2), FSM)
      do n=1,2   
        if (AgePhase==0 .and. n==2)  cycle
        if (AgePhase==2 .and. n==1)  cycle
        topX = 0
        dLL = 0D0        
        call BestRel(LL(:,n), 1, topX, dLL)
        if (topX /=1 .or. dLL < TA * dble(MIN(nS(sx(1),k), nS(sx(2),k))) &
         .or. (fcl==2 .and. (.not. FSM .or. &
          dLL < 2.0*TA * dble(MIN(nS(sx(1),k), nS(sx(2),k)))))) then
          IsPair = .FALSE.
          exit
        endif
      enddo
      if (.not. IsPair) cycle
      
      if (FSM .and. fcl==2) then
        call DoFSmerge(sx(1), sx(2), k)
      else 
        call DoMerge(sx(1), sx(2), k)
      endif
      
    else
      do m=1,2
        if (sx(m)>0 .and. sx(3-m)==0) then
          IsPair = .TRUE.
          FSM = .FALSE.
          call CheckRel(ij(3-m), 0, -sx(m), k, fcl, LL(:,1), LL(:,2))
          do n=1,2
            if (AgePhase==0 .and. n==2)  cycle
            if (AgePhase==2 .and. n==1)  cycle
            topX = 0
            dLL = 0D0 
            call BestRel(LL(:,n), fcl, topX, dLL)              
            if (.not. (topX==fcl .or. (fcl==3 .and. topX==2))) then
              IsPair = .FALSE.
              exit
            endif
          enddo
          if (.not. IsPair) cycle
          LLx = missing
          if (Complx>1 .and. ANY(SibID(1:ns(sx(m),k), sx(m), k) == Parent(ij(3-m),3-k))) then  ! inbreeding  
            call CheckPair(ij(3-m), Parent(ij(3-m),3-k), k, 3, LLx(:,1,1), LLx(:,1,2))
            call BestRel(LLX(:,1,2), 3, topXi, dLLi)  
            if (topXi /= 3) then
              IsPair = .FALSE. 
              cycle
            endif
          endif
          if (topX==2 .and. fcl/=2) then  
            call BestRel(LL(:,1), 2, topX, dLL)  
            if (topX==2 .and. dll > 2*TA) then 
              FSM = .TRUE.
            endif
          endif
          if (fcl==2 .or. FSM) then
            Par = 0
            topFS = 0
            dLLtmp = missing
            call getFSpar(sx(m), k, .TRUE., Par(m))           
             if (Par(m)/=0 .and. Parent(ij(3-m), 3-k) == Par(m)) then
              ! do nothing
            else if (Par(m)/=0 .and. all(Parent(SibID(1:ns(sx(m),k),sx(m),k),3-k)==Par(m))) then
              call setPar(ij(3-m), 3, Par(m), 3-k)
            else 
              call AddFS(ij(3-m), sx(m), k,0,k, LL(2,1), topFS, dLLtmp)
              if (Complx==0) then
                  call setPar(ij(3-m), 3, Parent(topFS, 3-k), 3-k)  
              else if (topFS>0) then
                if (parent(ij(3-m),3-k) == Parent(topFS, 3-k) .and. Parent(topFS, 3-k)/=0) then
                  ! do nothing
                else if (MAXVAL(dLLtmp, mask=dLLtmp<impossible)>2*TA) then
                  call CheckPair(ij(3-m), topFS, k, 2, LL(:,1), LL(:,2))
                  call BestRel(LL(:,2), 2, topX, dLL)
                  if (topX==2 .and. dll > 2*TA) then                
                    if (Parent(topFS, 3-k)/=0) then                
                      call setPar(ij(3-m), 3, Parent(topFS, 3-k), 3-k)
                    else if (Parent(ij(3-m), 3-k)/=0) then
                      call setPar(topFS, 3, Parent(ij(3-m), 3-k), 3-k)
                    else
                      call NewSibship(ij(3-m), topFS, 3-k)   ! new sibship (pair)      
                    endif
                  endif
                endif
              endif
            endif     
          endif  
          
          call setPar(ij(3-m), 3, -sx(m), k)
        endif
      enddo
    endif
  enddo 
enddo

end subroutine Clustering

! #####################################################################

subroutine Merging ! check if any of the existing clusters can be merged
use Global
implicit none

integer :: k, s, r, topX, xr, n, chunk_printdot(10), z
double precision :: LLm(7,2), dLL
logical :: FSM, OK

if (sum(nC)==0) then
  call print_dot(10)
  return
endif

chunk_printdot = mk_seq(sum(nC)-2, 10)      
z=0

do k=2,1,-1
!  if (Complx==0 .and. k==2) cycle
  do s=1,nC(k)-1
    if (modulo(s,20)==0)  call rchkusr()
    z=z+1
    if (quiet==-1 .and. any(chunk_printdot==z)) call print_dot(count(chunk_printdot==z))
    if (s >= nC(k)) exit
    r = s
    do xr=s+1, nC(k)
      r = r + 1
      if (r > nC(k)) exit   ! possible due to merged sibships      
      LLm = missing
      FSM = .FALSE.
      call CheckMerge(s, r, k, k, 1, LLm(:,1), LLm(:,2), FSM)
      if (LLM(1,2) > 0 .or. LLM(1,2) < LLM(7,2)) cycle
      if (.not. FSM .and. (Complx==0 .or. Hermaphrodites==1))  cycle   
      OK = .TRUE.
      topX = 0
      dLL = missing
      if (all(GpID(:,s,k) /= 0) .and. GpID(1,s,k)==GpID(1,r,k) .and. GpID(2,s,k)==GpID(2,r,k) .and. &
        ABS(LLM(7,1) - LLM(1,1)) < 0.1 .and. LLM(1,2) > LLM(7,2)) then
         OK = .TRUE.      ! identical LL s & r  FS vs. identical. occams razor --> merge
      else     
        do n=1,2    ! UseAge = (/.FALSE., .TRUE./)
          if (AgePhase==0 .and. n==2)  cycle
          if (AgePhase==2 .and. n==1)  cycle          
          call BestRel(LLm(:,n), 1, topX, dLL)
          if (topX /=1 .or. dLL < TA * dble(MIN(nS(s,k), nS(r,k)))) then
            OK = .FALSE.
            exit   
          endif
        enddo
      endif
      if (.not. OK)  cycle   
            
      if (FSM .and. (dLL > 2.0*TA * dble(MIN(nS(s,k), nS(r,k))) .or. &
        Complx==0 .or. Hermaphrodites==1)) then
        call DoFSmerge(s, r, k)
      else 
        call DoMerge(s, r, k)
      endif
      r = r-1  ! otherwise a cluster is skipped
      WHERE (chunk_printdot > (sum(nC)-2)) chunk_printdot = sum(nC)-2                       
    enddo
  enddo
enddo

! sometimes too few dots printed to screen
if (quiet==-1 .and. z < (sum(nC)-2)) then
  do r=z, sum(nC)-2
    if (any(chunk_printdot==r)) call print_dot(count(chunk_printdot==r))
  enddo
endif

end subroutine Merging

! #####################################################################

subroutine SibParent  
! for each sibship, check if a real indiv can replace the dummy parent
use Global
use OHfun
implicit none

integer :: k, s, xs, i, n, topX, CurNumC(2), Par, SClone, chunk_printdot(10), z, &
  j, nCandPar, CandPar(mxCP), h, SibTmp(maxSibSize), nSib, sib1, sxSib(maxSibSize)
double precision :: LL(7), dLL, LLtmp(7,2), ALR, LLO, LR, LLg(7)
logical :: NeedsOppMerge, Maybe, MaybeOpp, AncOK, FSM

if (sum(nC)==0) then
  call print_dot(10)
  return
endif

chunk_printdot = mk_seq(sum(nC), 10)      
z=0                                                    
CurNumC = nC
do k=1,2
  s = 0
  do xs=1, CurNumC(k)
    s = s+1
    if (modulo(s,20)==0)  call rchkusr()
    z=z+1 
    if (quiet==-1 .and. any(chunk_printdot==z)) call print_dot(count(chunk_printdot==z))
    if (s > nC(k)) exit   
    if (hermaphrodites==1) then
      call getFSpar(s, k, .TRUE., Par)
      if (Par < 0)  cycle
    endif   
    nCandPar = 0
    CandPar = 0
    NeedsOppMerge = .FALSE.
    
    do i=1,nInd
      Maybe = .TRUE.
      MaybeOpp = .FALSE.              
      if (nCandPar == mxCP) exit  !unlikely
      if (Sex(i)/=k .and. Sex(i)<3) cycle
      if (Parent(i,k)==-s) cycle
      if (ANY(GpID(:,s,k)==i)) cycle      
      if (DumClone(s,k)/=0 .and. sex(i)/=4)  cycle
      do n=1,nS(s,k)
        if (AgeDiff(SibID(n,s,k), i) <= 0) then
          maybe = .FALSE.
          exit
        endif
      enddo
      if (.not. Maybe) cycle         
      if (DoMtDif) then
        if (k==1 .and. mtDif(SibID(1,s,k), i))  cycle 
      endif        
      call CalcAgeLR(-s, k, i, k, k, -1, .TRUE., ALR)
      if (ALR==impossible)  cycle
      call ChkAncest(i,k, -s,k, AncOK)
      if (.not. AncOK)  cycle
      do n=1,nS(s,k)
        call CalcP2(SibID(n,s,k), 3, i, Parent(SibID(n,s,k),3-k), k, LR)
        if (LR == impossible .or. LR < 3*TF) then
          maybe = .FALSE.
          exit
        endif   
      enddo
      if (.not. Maybe) cycle
      call CalcP2(i, Sex(i), GpID(1,s,k), GpID(2,s,k), 1, LR)
      if (LR == impossible .or. LR < 3*TF)  cycle
      call chkPO(i, s, k, LR)
      if (LR < TF*nS(s,k))  cycle
      LLg = missing
      LL = missing
      topX = 0
      call CheckRel(-s, k, i, k, 1, LLg, LL)
      if (AgePhase <=1) then
        call BestRel(LLg, 1, topX, dLL)
      else
        call BestRel(LL, 1, topX, dLL)
      endif
      if (topX/=1) then
        Maybe = .FALSE.
        cycle
      endif
      if (Sex(i)>2 .and. DumClone(s,k)==0) then  ! check if parent of opposite sex instead
        MaybeOpp = .TRUE.
        Par = 0
        call getFSpar(s, k, .TRUE., Par)
        if (Par > 0) then
          MaybeOpp = .FALSE.
        else if (Par/=0 .and. Parent(i, 3-k) == Par) then
          MaybeOpp = .FALSE.  ! are HS
        else if (Par==0) then
          if (ANY(Parent(SibID(1:nS(s,k), s,k),3-k)>0)) then
            MaybeOpp = .FALSE.
          else  ! check if could all be FS
            do j=1, ns(s,k)-1
              do h=j+1, ns(s,k)
                call CalcAgeLR(sibID(j,s,k), 0, SibID(h,s,k), 0, 0, 2, .TRUE., ALR)
                if (ALR == impossible .or. ALR < 3.0*TF) then
                  MaybeOpp = .FALSE.
                  exit
                endif                 
              enddo
              if (.not. MaybeOpp) exit
            enddo
          endif
          if (MaybeOpp) then
            call OppMerge(s,k,LLO)
            if (LLO>NotImplemented .or. (CLL(s,k) - LLO) > ns(s,k)*TF) then
              MaybeOpp = .FALSE.
            endif
          endif                
        else if (Par < 0) then
          do n=1,nS(-Par,3-k)
            if (QLR_PO(i, SibID(n,-par,3-k)) < 5.0*TF) then   
              MaybeOpp = .FALSE.
              exit
            endif   
          enddo
          if (MaybeOpp) then
            do n=1,2
              if (GpID(n,-Par,3-k) <= 0) cycle
              if (QLR_PO(i, GpID(n,-Par,3-k)) < 5.0*TF) then
                MaybeOpp = .FALSE.
                exit
              endif   
            enddo
          endif                        
        endif
        if (MaybeOpp) then
          LLtmp = missing
          topX = 0
          if (Par < 0) then  ! may have more/fewer sibs
            call CheckRel(Par, 3-k, i, 3-k, 1, LLtmp(:,1), LLtmp(:,2))
            if (AgePhase <=1) then
              call BestRel(LLtmp(:,1), 1, topX, dLL)
            else
              call BestRel(LLtmp(:,2), 1, topX, dLL)
            endif
            if (topX/=1)  MaybeOpp = .FALSE.
          else if (Par == 0 .and. ns(s,k)>0) then
            sib1 = SibID(1,s,k)
            call PairPO(sib1, i, 3-k, 0, LLtmp(1,2))   
            call CalcU(sib1, k, i, 3-k, LLtmp(2,2))
            if (LLtmp(1,2)>0 .or. (LL(1)-LL(7)) - (LLtmp(1,2)-LLtmp(2,2)) > &
             TA*ns(s,k))  MaybeOpp = .FALSE.
          endif
        endif
        if (MaybeOpp) cycle
      endif
      
      if (Complx==0) then  ! ensure monogamous
        Par = 0       
        call getFSpar(s, k, .FALSE., Par)        
        if (Mate(i)/=Par .and. Mate(i)/=0 .and. Par/=0) then
          LLtmp = missing               
          if (Mate(i)<0 .and. Par<0) then
            call CheckMerge(-Par, -Mate(i), 3-k, 3-k, 8, LLtmp(:,1), LLtmp(:,2), FSM)
          else if (Mate(i)>0 .and. Par<0) then
            call CheckRel(Par, 3-k, Mate(i), 3-k, 1, LLtmp(:,1), LLtmp(:,2))
          endif
          if (LLtmp(1,2) < 0 .and. (LLtmp(1,2) - MaxLL(LLtmp(2:7,2)) > -TA)) then
            NeedsOppMerge = .TRUE.
          else
            cycle
          endif
        endif
      endif
      
      if (Maybe) then               
        nCandPar = nCandPar + 1
        CandPar(nCandPar) = i
      endif
    enddo  ! i
    
    if (nCandPar == 1) then
      Par = 0
      if (NeedsOppMerge) then
        call getFSpar(s, k, .FALSE., Par)
        if (Mate(CandPar(1)) < 0) then
          call DoMerge(-Mate(CandPar(1)), -Par, 3-k)  ! Par gets removed
        else
          call getOff(Par, 3-k, .TRUE., nSib, SibTmp, sxSib)
          do n=1,nSib
            call setPar(SibTmp(n), sxSib(n), Mate(CandPar(1)), 3-k)
          enddo
          call DoMerge(0, -Par, 3-k)
        endif
      else if (Complx==0 .and. Mate(CandPar(1))==0) then
        Mate(CandPar(1)) = Parent(SibID(1,s,k), 3-k)
      endif
      
      if (hermaphrodites/=0 .and. Sex(CandPar(1))==4) then
        SClone = DumClone(s,k)       
        if (SClone /= 0) then
          call getOff(-s, k, .FALSE., nSib, SibTmp, sxSib) 
          do n=1,nSib 
            call setParTmp(SibTmp(n), sxSib(n), CandPar(1), k)  ! else conflict w CheckSelfed()
          enddo 
          
          call getOff(-SClone, 3-k, .TRUE., nSib, SibTmp, sxSib)
          do n=1,nSib 
            call setPar(SibTmp(n), sxSib(n), CandPar(1), 3-k)
          enddo  
          call DoMerge(0, SClone, 3-k)  !removes Sclone cluster 
        endif
      endif

      call getOff(-s, k, .TRUE., nSib, SibTmp, sxSib)   ! includes dummy sibs
      do n=1,nSib 
        call setPar(SibTmp(n), sxSib(n), CandPar(1), k)
        if (SibTmp(n)>0)  ToCheck(SibTmp(n)) = .FALSE.  ! might not be most likely when considered separately 
      enddo  
      call DoMerge(0, s, k)  !removes cluster s 
      s = s-1  ! otherwise a cluster is skipped
      WHERE (chunk_printdot > (sum(nC))) chunk_printdot = sum(nC)
!    else   
!       TODO?
    endif
  enddo ! s
enddo ! k

end subroutine SibParent

! #####################################################################

subroutine MoreParent
! for each individual, check if a parent can be assigned now.
use Global
use OHfun
use CalcLik
use sort_module, ONLY: getRank_i                                
implicit none

integer :: x, i, j, k, s, curPar(2), nCP(2), CandPar(mxCP, 2), &
   TopTmp, BYrank(nInd)
double precision :: LLP(2), ALR(2), LL(7,2), LRQ, dLL, LRFS(mxCP,2,2), LLtmp(7,2)
logical :: DoNewPars, AncOK, DropS, KeepOld

if (hermaphrodites==1 .or. ALL(BY==BY(1))) then
  DoNewPars = .FALSE.
else
  DoNewPars = .TRUE.   ! do check for additional parent-offspring pairs. 
endif

call getRank_i(BYrank)

do x=1, nInd
  if (MODULO(x,20)==0)  call rchkusr()   
  if (quiet==-1 .and. any(chunk_printdot_i==x)) call print_dot()
  i = BYRank(x)                     
  if (ALL(Parent(i,:)/=0) .and. .not. ToCheck(i)) cycle
  call CalcLind(i)
  CurPar = Parent(i,:)
  nCP = 0
  CandPar = 0
  LRFS = missing
  
  do k=1,2
    if (curPar(k)/=0) then
      nCP(k) = nCP(k) +1
      CandPar(nCP(k), k) = curPar(k)
    endif
    call setParTmp(i, Sex(i), 0, k)  
    call SetEstBY(i, Sex(i))                        
    call SetEstBY(curPar(k), k)
  enddo
  call SetEstBY(i, Sex(i))
  
  do k=1,2
    if (Complx==0 .and. k==2 .and. all(parent <= 0)) exit
    if (nC(k)==0)  cycle
    do s=1, nC(k)
      if (nCP(k) == mxCP)  exit
      if (ANY(CandPar(:,k) == -s))  cycle   
      if (ANY(GpID(:,s,k)==i)) cycle
      call ChkAncest(-s, k, i, sex(i), AncOK)
      if (.not. AncOK)  cycle
      call CalcAgeLR(i,0,-s,k,0,1, .TRUE.,ALR(1))
      if (ALR(1)==impossible)  cycle
      if (DoMtDif) then
        if (k==1 .and. mtDif(i, SibID(1,s,k)))  cycle    
      endif  
      if (Complx==0) then
        call QFS(i,s,k,LRQ)
        if (LRQ < TF) cycle
      else
        call Qadd(i, s, k, LRQ)
        if (LRQ < MAX(ns(s,k),3)*TF) cycle   
      endif
      LL = missing
      LLtmp = missing
      if (DumClone(s,k)==0) then 
        call SibChk(i,s,k,3, 1, LL(:,1))
        if (LL(3,1)>0 .or. MaxLL(LL(2:3,1)) - LL(7,1) < TA .or. & 
         (MaxLL(LL(2:3,1)) - MaxLL(LL(4:6,1)) < TF .and. ANY(LL(4:6,1) < 0))) cycle 
      endif     
      
      if (Complx>1 .and. ANY(SibID(1:ns(s,k), s, k) == Parent(i,3-k))) then  ! inbreeding
        call CheckPair(i, Parent(i,3-k), k, 3, LLtmp(:,1), LLtmp(:,2))
        call BestRel(LLtmp(:,2), 3, topTmp, dLL)
        if (topTmp /= 3)  cycle
      endif
      
      nCP(k) = nCP(k) +1
      CandPar(nCP(k), k) = -s
      if (Complx==0 .and. nCP(3-k)<mxCP) then
        nCP(3-k) = nCP(3-k) +1
        if (DumMate(s,k) == 0) call Erstop("Mono error: no DumMate", .TRUE.)
        CandPar(nCP(3-k), 3-k) = DumMate(s,k)
      endif
    enddo
  enddo

  if (DoNewPars)  then
    do j=1, nInd  ! candidate parent.
      if (i==j) cycle
      if (ANY(CandPar == j) .and. Sex(j)<3)  cycle  ! already included 
      if (ANY(Parent(j,:)==i) .or. ANY(Parent(i,:)==j)) cycle            
      if (AgeDiff(i,j) <= 0)  cycle  ! note: unknown = missing > 0    
      if (DoMtDif) then      
        if (Sex(j)==1 .and. mtDif(i,j))  cycle   
      endif  
      call ChkAncest(j, sex(j), i, sex(i), AncOK)
      if (.not. AncOK)  cycle
      if (Sex(j) < 3) then
        if (Parent(i, Sex(j)) < 0) cycle   ! replacings done elsewhere
        if (nCP(Sex(j))==mxCP) cycle
      else 
        if (ANY(nCP == mxCP)) cycle
      endif
      ALR = missing
      LLP = missing
      if (QLR_PO(i,j) < 5.0*TF)  cycle
      call CalcAgeLR(i,sex(i), j, sex(j), 0, 1, .TRUE., ALR(1))
      if (ALR(1) == impossible)  cycle
      call CalcAgeLR(j, sex(j), i,sex(i), 0, 1, .TRUE., ALR(2))
      if (ALR(2) /= impossible .and. (ALR(1)-ALR(2)) < TF)  cycle 
       if (Sex(j) < 3) then
        call PairPO(i, j, sex(j), 0, LLP(1))
      else
        call PairPO(i, j, 1, 0, LLP(1))
      endif 
      if (LLP(1) > 0) cycle
      call CalcU(i,sex(i), j, sex(j), LLP(2))
      if ((LLP(1) - LLP(2)) < TF) cycle  
      do k=1,2
        if (Sex(j)<3 .and. Sex(j)/= k) cycle                             
        if (DoMtDif) then
          if (k==1 .and. mtDif(i,j))  cycle  ! when sex(j)>=3 
        endif 
        if (nCP(k) < mxCP .and. .not. any(candPar(:,k) == j))  then
          nCP(k) = nCP(k) +1
          CandPar(nCP(k), k) = j
        endif
        if (Complx==0 .and. Mate(j)/=0 .and. nCP(3-k)<mxCP) then
          nCP(3-k) = nCP(3-k) +1
          CandPar(nCP(3-k), 3-k) = Mate(j)
        endif
      enddo
    enddo  ! j
  endif
  
  KeepOld = .FALSE.
  if (ALL(nCP <=1) .and. ALL(candPar(1,:) == curPar)) then  ! no new/additional parents found
    KeepOld = .TRUE.
    do k=1,2
      if (curPar(k) < 0) then
        if (IsNewSibship(-curPar(k), k) .and. ns(-curpar(k),k) <=3)  KeepOld = .FALSE.
      endif
    enddo
    if (KeepOld) then  
      do k=1,2
        call setParTmp(i, Sex(i), curPar(k), k)  ! restore
        call SetEstBY(i, Sex(i))
        call SetEstBY(curPar(k), k)
      enddo
      call SetEstBY(i, Sex(i))
      ToCheck(i) = .FALSE.
      ! exception: GGpairs   added 2024-12-30
      do k=1,2
        if (curpar(k) < 0) then
          if (ns(-curpar(k),k)==1)  ToCheck(i) = .TRUE.  
        endif
      enddo
      cycle
    endif
  endif
  
  call SelectParent(i, Sex(i), nCP, CandPar, .FALSE.)

  if (ANY(Parent(i,:)/=CurPar)) then
    DropS = .FALSE.
    if (Complx == 0)  call UpdateMate(i, Sex(i), curPar, .FALSE.)   
    do k=1,2
      if (curPar(k)<0 .and. Parent(i,k)/=curPar(k)) then
        call CheckDropSibship(-curPar(k), k, DropS)
        if (hermaphrodites/=0 .and. .not. DropS) then
          call CheckSelfed(curPar(k),k)
        endif
      endif
    enddo
    call setEstBY(i,sex(i))
  else
    ToCheck(i) = .FALSE.
    do k=1,2
      if (curpar(k) < 0) then
        if (ns(-curpar(k),k)==1)  ToCheck(i) = .TRUE.  
      endif
    enddo                                    
  endif    
enddo

end subroutine MoreParent

! #####################################################################

subroutine UpdateMate(A, kA, OldPar, ParOnly)
use Global
implicit none

integer, intent(IN) :: A, kA, OldPar(2)
logical, intent(IN) :: ParOnly  ! T:only parents / F:also dummy parents
integer :: NewPar(2), m, nOff, sxOff(maxSibSize), Off(maxSibSize), x

if (Complx /= 0)  return

NewPar = getPar(A, kA)
if (all(NewPar == OldPar))  return

do m=1,2
  if (oldPar(m) > 0 .and. oldpar(m) /= newpar(m)) then
    call getOff(oldPar(m), m, .TRUE., nOff, Off, sxOff)
    if (nOff == 0) then
      if (any(Mate == OldPar(m))) then
        x = MINLOC(ABS(Mate - oldPar(m)), DIM=1)
        Mate(x) = 0
      else if (any(DumMate(:,3-m) == OldPar(m))) then
        x = MINLOC(ABS(DumMate(:,3-m) - oldPar(m)), DIM=1)
        DumMate(x, 3-m) = 0
      endif
    endif
  ! else if oldpar(m) < 0, CheckDropSibship will take care of it
  endif  
enddo

do m=1,2
  if (newpar(m) > 0) then
    if (Mate(newpar(m)) == 0) then
      Mate(newpar(m)) = newpar(3-m)
    else if (Mate(newpar(m)) /= newpar(3-m)) then
      call Erstop("Something going wrong with Mate", .TRUE.)
    endif
  else if (newpar(m) < 0) then
    if (DumMate(-newpar(m),m) == 0) then
      DumMate(-newpar(m),m) = newpar(3-m)
    else if (DumMate(-newpar(m),m) /= newpar(3-m)) then
      call Erstop("Something going wrong with DumMate", .TRUE.)
    endif
  endif 
enddo

if (A > 0 .and. any(newPar == 0) .and. any(NewPar/=0) .and. .not. ParOnly) then   ! create singleton sibship
  do m=1,2
    if (NewPar(m) == 0) then
      call NewSibship(A, 0, m)   ! takes care of (Dum)Mate
    endif
  enddo
endif

end subroutine UpdateMate

! #####################################################################

subroutine getfocal(A, B, s, k, focal)
use Global
implicit none

integer, intent(IN) :: A, B, s, k
integer, intent(OUT) :: focal
integer :: j, BB(MaxsibSize), nB, opParB, opParA

BB = 0
nB = 0
if (B/=0) then
  BB(1) = B
  nB = 1
else if (s/=0) then
  BB = SibID(:,s,k)
  nB = nS(s,k)
else
  call ErStop("getFocal: B=0 and s=0", .TRUE.)
endif

focal = 3
if (Complx==0) then
  focal = 2  ! FS
else if (nB == 0) then  ! empty sibship ; should only happen w hermaphrodites or Complx=0 
  focal = 3 
else if (Parent(A,k)==0 .and. Parent(A,3-k)/=0 .and. &
 any(Parent(BB(1:nB), 3-k) == Parent(A,3-k))) then
  focal = 2
else if (ALL(Parent(BB(1:nB),k)==0) .and. Parent(A,3-k)/=0 .and. &
 ALL(Parent(BB(1:nB), 3-k) == Parent(A,3-k))) then
  focal = 2
else if (ALL(Parent(A,:)==0) .and. ALL(Parent(BB(1),:)==0) .and. hermaphrodites/=2) then
  focal = 2
else if (all(Parent(A,:)<=0) .and. Parent(BB(1),k)<0 .and. ALL(Parent(BB(1:nB),3-k)<=0)) then
  call getFSpar(-Parent(BB(1),k), k, .TRUE., opParB)
  if (opParB < 0) then
    if (ALL(Parent(A,:)==0)) then
      focal = 2  ! FS add
    else if (Parent(A,k)<0) then
      call getFSpar(-Parent(A,k), k, .TRUE., opParA)
      if (opParA < 0) then
        focal = 2     ! FS merge
      endif 
    endif
  endif
endif

if (focal==2) then  ! exception: cannot be /unlikely FS based on age
  do j=1,nB 
    if (getAP(AgeDiff(A,BB(j)), 2, 0,0, Impossible) == Impossible) then
      focal = 3
      exit
    else if (AgePhase == 2) then
      if (getAP(AgeDiff(A,BB(j)), 3, 0, k, log10(zero)) - &
       MAX(getAP(AgeDiff(A,BB(j)), 2, 0,0, log10(zero)), &
           getAP(AgeDiff(A,BB(j)), 3, 0, 3-k, log10(zero))) > 2.0*ABS(TF)) then
        focal = 3
      endif
    endif
  enddo 
endif

end subroutine getfocal

! #####################################################################

subroutine SibGrandparents 
! for each sibship, find the most likely male and/or female grandparent
use Global
use sort_module, ONLY: sort_sibships                                    
implicit none

integer :: k, s, i, r, m, par, x, nCG(2), candGP(mxCP, 2), curGP(2), &
   ix, not4(5), n, chunk_printdot(10)
double precision :: LRG, ALRtmp(2), LLX(3,2), dx(maxSibSize), LLA(7)
logical :: skipCluster(maxval(nC),2), AncOK, DropS, Maybe
integer, allocatable :: s_sorted(:), k_sorted(:)                                                

if (sum(nC)==0) then
  call print_dot(10)
  return
endif

skipCluster = .FALSE.
do k=1,2
  if (nC(k)==0)  cycle
  do s=1, nC(k)
    if (ALL(Parent(SibID(1:nS(s,k), s, k), 3-k) < 0)) then
      call getFSpar(s, k, .TRUE.,par)
      if (par < 0) then
        if (nS(-par, 3-k) == nS(s,k) .and. DumClone(s,k)==0) then  !cannot tell if mat or pat
          skipCluster(s,k) = .TRUE.  ! done in FsibsGPs()  
        endif
      endif          
    endif
  enddo
enddo                    

call sort_sibships(s_sorted, k_sorted)

chunk_printdot = mk_seq(sum(nC), 10)

not4 = (/1,2,3,5,6/)
do x=1, SUM(nC)
  if (x > sum(nC))  exit   ! sibships dropped during run; do loop max not re-evaluated?
  if (MODULO(x,10)==0)  call rchkusr()
  s = s_sorted(x)
  k = k_sorted(x)
  
  if (quiet==-1 .and. any(chunk_printdot==x)) call print_dot(count(chunk_printdot==x))   
  if (ALL(GpID(:,s,k)/=0) .and. .not. IsNewSibship(s,k)) cycle
  if (ALL(GpID(:,s,k)>0) .and. ns(s,k)==1) cycle  
  if (skipCluster(s,k) .and. .not. (hermaphrodites==2 .and. k==1))  cycle
  nCG = 0 
  CandGP = 0
  CurGP = GpID(:,s,k)
  do m=1,2
    if (GpID(m,s,k)/=0) then
      nCG(m) = 1
      CandGP(1,m) = GpID(m,s,k)
    endif
    call setParTmp(-s, k, 0, m)
  enddo
  call setEstBY(-s,k)                   
  
  do i=1,nInd
    if (Parent(i,k)==-s) cycle
    if (ANY(CandGP==i)) cycle
    if (Sex(i)<3) then
      if (nCG(Sex(i))==mxCP)  cycle
    else
      if (ANY(nCG==mxCP))  cycle 
    endif
    if (DoMtDif) then
      if (k==1 .and. Sex(i)==1 .and. mtDif(SibID(1,s,k), i))  cycle   
    endif       
     Maybe = .TRUE.
    do n=1,ns(s,k)
      if (AgeDiff(SibID(n,s,k), i) <= 1) then
        Maybe = .FALSE.
        exit
      endif
    enddo
    if (.not. Maybe)  cycle
    ALRtmp = missing
    call CalcAgeLR(-s,k, i,Sex(i), 0,1, .TRUE., ALRtmp(1))
    if (ALRtmp(1) == impossible) cycle
    call CalcAgeLR(i,Sex(i), -s,k, 0,1, .TRUE., ALRtmp(2))
    if (ALRtmp(2)/=impossible .and. (ALRtmp(1)-ALRtmp(2)) < TF)  cycle
    call ChkAncest(i,0,-s,k, AncOK)
    if (.not. AncOK)  cycle
    call QGP(i, Sex(i), s, k, LRG)
    if (LRG < TF * MAX(dble(nS(s,k)),2D0))  cycle      ! 2TF for ns=1
    LLA = missing
    call GPfilter(i,s,k,LLA)
    if (LLA(4)>0 .or. LLA(4) - LLA(7) < -TA .or.  & 
     (LLA(4) - MaxLL(LLA(not4)) < 2*TF .and. ANY(LLA(not4) < 0))) cycle  
    do m=1,2
      if (Sex(i)<3 .and. Sex(i)/=m) cycle
      if (ncG(m) < mxCP) then  ! arbitrary threshold to limit comp. time
        ncG(m) = nCG(m) + 1
        CandGP(nCG(m), m) = i
      endif
    enddo
  enddo
  
  do m=1,2
    do r=1, nC(m) 
      if (ncG(m) == mxCP) exit
      if (m==k .and. s==r) cycle
      if (any(CandGP(:,m) == -r)) cycle  ! current GP
      call ChkAncest(-r,m, -s,k, AncOK)
      if (.not. AncOK)  cycle
      if (nS(r,m)==1 .and. ANY(SibID(1:nS(s,k),s,k) == SibID(1,r,m))) cycle
      if (DoMtDif) then
        if (k==1 .and. m==1 .and. mtDif(SibID(1,s,k), SibID(1,r,m)))  cycle 
      endif  
      if (m/=k .and. complx<2) then
        if (ALL(Parent(SibID(1:ns(s,k),s,k),m) == -r))  cycle
        if (ALL(Parent(SibID(1:ns(r,m),r,m),k) == -s))  cycle
      endif
      call CalcAgeLR(-s,k, -r,m, 0,1, .TRUE., ALRtmp(1))
      if (ALRtmp(1) == impossible .or. ALRtmp(1) < 3.0*TF) cycle
      call CalcAgeLR(-r,m, -s,k, 0,1, .TRUE., ALRtmp(2))
      if (ALRtmp(2)/=impossible .and. (ALRtmp(1)-ALRtmp(2)) < TF)  cycle 
      if (ALL(ABS(ALRtmp) < 0.001))  cycle  ! no age info                                                   
      if (hermaphrodites==0) then
        call QGP(-r, m, s, k,  LRG) 
        if (LRG < TF*dble(MIN(nS(s,k), nS(r,m)))) cycle  ! conservative.
      endif

      LLX = missing
      call PairUA(-s, -r, k, m, LLX(1,1))
      if (LLX(1,1)>0) cycle
      call CalcU(-s,k, -r,m, LLX(1,2)) 
      if ((LLX(1,1) - LLX(1,2)) < nS(s,k)*TF)  cycle
      call addFS(0, r, m, s, k, LLX(2,1), ix, dx) 
      if ((MaxLL(LLX(:,1)) - LLX(1,2)) < TA)  cycle
      if (ncG(m) < mxCP) then
        nCG(m) = nCG(m) + 1
        CandGP(nCG(m), m) = -r
      endif                
    enddo
  enddo
  
  if ((ALL(nCG <=1) .and. ALL(CandGP(1,:) == curGP) .and. .not. &
    (any(curGP==0) .and. any(curGP<0))) .or. &
    (ns(s,k)==1 .and. all(candGP>=0))) then    ! all options already considered by GGpair
    do m=1,2
      call setParTmp(-s, k, curGP(m), m)  
    enddo
    call setEstBY(-s,k)                   
    cycle
  endif 
  
  call SelectParent(-s, k, nCG, candGP, .FALSE.)
  
  if (any(GpID(:,s,k)/= curGP)) then
    if (Complx == 0) call UpdateMate(-s, k, curGP, .FALSE.)
    if (any(GpID(:,s,k) /= curGP .and. curGP/=0)) then
      ToCheck(SibID(1:ns(s,k), s, k)) = .TRUE.
    endif
  endif
  
  if (ALL(GpID(:,s,k)==0) .and. nS(s,k)==1) then  ! single sib left; remove sibship 
    i = SibID(1,s,k)                
    call CheckDropSibship(s, k, DropS)
    if (DropS) then
      call sort_sibships(s_sorted, k_sorted)
      do r=s, nC(k)
        skipCluster(r,k) = skipCluster(r+1,k)
      enddo
      WHERE (chunk_printdot > sum(nC)) chunk_printdot = chunk_printdot -1   
      if (Parent(i,3-k) < 0) then
        par = Parent(i,3-k)
        call CheckDropSibship(-Parent(i,3-k), 3-k, DropS)
        if (DropS) then
          call sort_sibships(s_sorted, k_sorted)
          do r=-par, nC(3-k)
            skipCluster(r,3-k) = skipCluster(r+1,k)
          enddo
          WHERE (chunk_printdot > sum(nC)) chunk_printdot = chunk_printdot -1  
        endif
      endif
    endif
    cycle
  endif
enddo  ! x

end subroutine SibGrandparents

! #####################################################################

subroutine GPfilter(A, SB, k, LLg)
use Global
implicit none

integer, intent(IN) :: A, SB, k
double precision, intent(OUT) :: LLg(7)
integer :: fsi, sibfcl, kA
double precision :: ALR, dx(maxSibSize), LLtmp

LLg = missing
call AddGP(A, SB, k, LLg(4))
if (LLg(4) > 0)  return
! U
call CalcU(A, k, -SB, k, LLg(7))  
if (LLg(4) - LLg(7) < TA)  return
if (Complx==0)  return   ! TODO: implement properly for monogamous?  
! GGP / 3rd degree rel   ! before or after FA?
call AddGGP(A, SB, k, LLg(6))  
if (LLg(4) - LLg(6) < TF .and. LLg(6)<0)  return
! FA
if (any(Parent(A,:)/=0)) then
  call CalcAgeLR(-SB,k, A,Sex(A), 0,2, .TRUE., ALR)
  if (ALR /= impossible)  call pairUA(-SB, A, k, 3, LLg(5))
  if (LLg(4) - LLg(5) < TF .and. LLg(5)<0) then
    if (Sex(A) < 3) then
      kA = Sex(A)
    else
      kA = 1
    endif    
    if (Parent(A,3-kA)/=0 .and. GpID(3-kA,SB,k)==0) then
      call setParTmp(-SB,k, Parent(A,3-kA),3-kA)
      call AddGP(A, SB, k, LLtmp)
      call setParTmp(-SB,k, 0,3-kA)
      if (LLtmp < 0 .and. LLtmp > LLg(4)) then
        LLg(4) = LLtmp
        if (LLg(4) - LLg(5) < TF)  return
      endif
    else
      return
    endif
  endif
endif

! FS/HS  
call getfocal(A, 0, SB, k, sibfcl)
if (sibfcl == 2) then
  call AddFS(A, SB, k,0,k, LLg(2), fsi, dx)
else if (any(Parent(A,:)/=0)) then
  call AddSib(A, SB, k, LLg(3))
endif

end subroutine GPfilter

! #####################################################################

subroutine SibChk(A, SB, k, focal, cat, LLg)  ! 1=filter, 2=confirm
use Global
implicit none

integer, intent(IN) :: A, SB, k, focal, cat  
double precision, intent(OUT) :: LLg(7)
integer :: fsi
double precision :: ALR, dx(maxSibSize), Threshold

if (cat==1) then
  Threshold = TF
else if (cat==2) then
  Threshold = TA
else
  Threshold = missing
  call Erstop("SibChk: cat must be 1 or 2", .TRUE.)
endif

LLg = missing
! FS/HS 
if (focal==2 .or. Complx==0) then
  call AddFS(A, SB, k,0,k, LLg(2), fsi, dx)  ! includes age check
else if (focal==3) then
  call AddSib(A, SB, k, LLg(3))
else
  call Erstop("SibChk: focal must be 2 or 3", .TRUE.)
endif

if (all(LLg > 0))  return
! U
if (cat==1) then
  call CalcU(A, k, -SB, k, LLg(7))  
  if (LLg(focal) - LLg(7) < TA) then
    if (focal==3 .and. Parent(A,3-k)==0 .and. LLg(focal)-LLg(7) > TF) then
      call AddFS(A, SB, k,0,k, LLg(2), fsi, dx)  
      if (LLg(2) - LLg(7) < TA)  return
    else
      return
    endif
  endif
endif

if (focal==3 .and. all(Parent(A,:)==0) .and. cat==1) then  ! 2nd rels indistinguishable
  LLg(4:5) = LLg(3)
  return
endif
! FA
call CalcAgeLR(A,Sex(A), -SB,k, 3,4, .TRUE., ALR) 
if (ALR /= impossible)  call addFA(A, SB, k, LLg(5))
if (LLg(focal) - LLg(5) < Threshold .and. LLg(5)<0)  return
! GP
call AddGP(A, SB, k, LLg(4))
if (LLg(focal) - LLg(4) < Threshold .and. LLg(4)<0)  return
! HA / 3rd degree rel
if (cat==2)  call pairUA(A, -SB, k, k, LLg(6))     !  .and. focal==3
!if (LLg(focal) - LLg(6) < Threshold .and. LLg(6)<0)  return

end subroutine SibChk

! #####################################################################

subroutine Calc4U(Par, B, kB,  A, kA, LLU, LLcor)  
use Global
implicit none

integer, intent(IN) :: Par(2), B, kB,  A, kA
double precision, intent(OUT) :: LLU(4), LLcor(3,2)
integer :: m, y, CY(4), kY(4), x, v, ParA(2)
double precision :: LLtmp(3), LLoverlap(4,4)
logical :: ConPar(4,4)

CY = (/ Par, B, A /)
kY = (/ 1, 2, kB, kA /)

ParA = getPar(A, kA)
do m=1,2
  if (ParA(m)==0)  cycle
  call setParTmp(A, kA, 0, m)   
enddo

LLU = 0D0  
! Individual likelihoods, excluding overlap with A if any
LLcor = 0D0
! correction factors: LL(all 4 indiv) = e.g. LL(A+B) from CheckRel + LLcor(3,m)
LLoverlap = 0D0
! likelihoods of overlaps/ change in likelihoods when considering vs ignoring connections
LLtmp = missing  

call CalcU(A, kA, 0, 0, LLU(4))
do y=1,3
  if (CY(y)==0) cycle
  call CalcU(CY(y),kY(y), A,kA, LLtmp(1))
  LLU(y) = LLtmp(1) - LLU(4)
enddo  
! pairs likelihoods, if no overlap                                  
do m=1,2
  LLcor(m,m) = LLU(3-m) + LLU(3) 
enddo
LLcor(3,:) = LLU(1) + LLU(2)  

ConPar = .FALSE.
if (ANY(CY(1:3)<0)) then
  do m=1,2
    if (Par(m)==0) cycle
    call Connected(Par(m), m, A, kA, ConPar(4,m))
    if (B/=0)  call Connected(Par(m), m, B, kB, ConPar(3,m))
  enddo
  call Connected(Par(1), 1, Par(2), 2, ConPar(2,1))
  if (B/=0)  call Connected(B, kB, A, kA, ConPar(4,3))
endif

if (ANY(ConPar)) then 
  ! likelihoods of overlaps/ change in likelihoods when considering vs ignoring connections
  do y=1,3
    do x=y+1, 4
      if (.not. ConPar(x,y))  cycle
      call CalcU(CY(x),kY(x), CY(y),kY(y), LLtmp(1))
      call CalcU(CY(x),kY(x), 0,0, LLtmp(2))
      call CalcU(CY(y),kY(y), 0,0, LLtmp(3))
      LLoverlap(x,y) = LLtmp(1) - LLtmp(2) - LLtmp(3)
      LLoverlap(y,x) = LLoverlap(x,y)
    enddo
  enddo
  
  do m=1,2        
    do y=1,3  ! focal
      if (y/=m .and. y/=3) cycle           
      if (y==1) then
        call CalcU(CY(2),kY(2), CY(3),kY(3), LLcor(y,m))
      else if (y==2) then
        call CalcU(CY(1),kY(1), CY(3),kY(3), LLcor(y,m))    
      else if (y==3) then
        call CalcU(CY(1),kY(1), CY(2),kY(2), LLcor(y,m))
      endif
      do x=1,3
        if (ConPar(4,x) .and. x/=y) then
          LLcor(y,m) = LLcor(y,m) + LLoverlap(x,4)
        endif
        do v=1,2
          if (ConPar(x,v) .and. (x==y .or. v==y)) then
            LLcor(y,m) = LLcor(y,m) + LLoverlap(x,v) 
          endif
        enddo
      enddo
    enddo
  enddo
endif

do m=1,2
  call setParTmp(A, kA, parA(m), m)
enddo

end subroutine Calc4U

! #####################################################################

subroutine GGpairs  ! find & assign grandparents of singletons
use Global
use sort_module, ONLY: getRank_i                                
implicit none

integer :: i, j, k, nCG(2,2), CandG(2,mxCP, 2), n, s, BYrank(nInd), x
double precision :: LRS, LRG, ALR, ALRx(2), LRx, LLx(7,2), LL(7,2), LLdGP
logical :: AncOK, MaybePair

call getRank_i(BYrank)

do x=1, nInd
  if (MODULO(x,200)==0)  call rchkusr()
  if (quiet==-1 .and. any(chunk_printdot_i==x)) call print_dot()
  
  i = BYRank(x)
  if (ALL(Parent(i,:)==0) .and. AgePhase <2 .and. hermaphrodites/=2)  cycle
  ! can't determine if mat or pat GP  ! TODO: more nuanced.
  if (ALL(Parent(i,:)/=0)) cycle
  nCG = 0  
  CandG = 0
  LL = missing
  
  do k=1,2
    if (Parent(i,k)/=0) cycle
    do j=1, nInd
      if (i==j)  cycle                
      if (ANY(nCG(k,:)>=mxCP)) cycle
      if (AgeDiff(i,j) <= 1)  cycle
      if (any(parent(j,:) == i))  cycle
      call ChkAncest(j, sex(j), i, sex(i), AncOK)
      if (.not. AncOK)  cycle
      call PairQHS(i, j, LRS)     
      if (LRS < TF)  cycle
      call CalcAgeLR(i,Sex(i), j,Sex(j), k,4,.TRUE.,ALR)
      if (ALR==impossible)  cycle  
      call LRGG(i,k,j,Sex(j),LRG)
      if (LRG < -TA)  cycle                                 
      if (DoMtDif) then
        if (k==1 .and. Sex(j)==1 .and. mtDif(i,j))  cycle   
      endif   
      if (Parent(i,3-k)<= 0  .and. hermaphrodites/=2) then      
        call CalcAgeLR(i,Sex(i), j,Sex(j), 3-k,4,.TRUE.,ALRx(1))
        if (ALRx(1)/=impossible .and. (ALR - ALRx(1)) < TA) then     
          if (Parent(i,3-k)==0) then   !  .and. .not. DoMtDif   !! it's complicated
            cycle   ! unclear if pat. or mat. GP     
          else if (Parent(i,3-k) < 0) then
            if (.not. any(GpID(:,-Parent(i,3-k),3-k) == j)) then  ! allow for double GP  
              call CalcAgeLR(Parent(i,3-k),3-k, j,Sex(j), 0,1, .TRUE., ALRx(2))
              if (ALRx(2)/=impossible .and. (ALR - ALRx(2)) < TA) then       
                call QGP(j, sex(j), -Parent(i,3-k),3-k, LRx)
                if (LRx > TF*ns(-Parent(i,3-k),3-k)) then
                  LLx = missing
                  call CheckAdd(j, -Parent(i,3-k),3-k, 4, LLx(:,1), LLx(:,2))
                   if ((LLx(4,1)- MaxLL(LLx((/1,2,6,7/),1))) > TF) then
                    call PairdGP(i,j,k, LLdGP)  ! double GP
                    if ((LLdGP - LLx(4,1)) < TA)  cycle   ! plausible that GP only of opposing sibship  
                  endif
                endif
              endif
            endif
          endif
        endif
      endif
      
      MaybePair = .TRUE.
      LL = missing
      call CheckPair(i,j,k, 4, LL(:,1), LL(:,2))
      do n=1,2
        if (AgePhase==0 .and. n==2)  cycle
        if (AgePhase==2 .and. n==1)  cycle
        if (LL(4,n)<0 .and. (LL(4,n)- MaxLL(LL((/1,2,6,7/),n))) > TF) then  
          MaybePair = .TRUE.
        else
          MaybePair = .FALSE.
          exit
        endif
      enddo
      if (.not. MaybePair) cycle
      do n=1,2
        if (Sex(j)/=n .and. sex(j)<3)  cycle
        if (nCG(k,n) == mxCP)  cycle
        nCG(k, n) = nCG(k, n) +1
        CandG(k, nCG(k,n), n) = j
      enddo
    enddo
  enddo
  
  if (ANY(nCG>0)) then
    do k=1,2
      if (ANY(nCG(k,:)>0)) then
        call NewSibship(i, 0, k)
        s = nC(k)
        call SelectParent(-s, k, nCG(k,:), CandG(k,:,:), .FALSE.)
        if (ALL(GpID(:,s,k)==0)) then 
          call RemoveSib(i, s, k)    
          call DoMerge(0, s, k)
        else if (Complx == 0) then
          call UpdateMate(-s, k, (/0,0/), .FALSE.)
        endif
      endif
    enddo
  endif  
  
  do k=1,2
    if (nc(k)==0)  cycle
    if (nS(nc(k), k) == 0 .or. SibID(1, nC(k), k) == 0) then
      call Erstop("grandparent pairs -- empty sibship!", .TRUE.)
    endif
  enddo 
enddo

contains
  subroutine PairdGP(A, B, k, LL)  ! double GP. Parent(A,3-k)<0; returns LL Parent(A,3-k) + B
  integer, intent(IN) :: A,B,k
  double precision, intent(OUT) :: LL
  integer :: l, x, y, curGP(2,2), m, z, v, i, Ei,j
  double precision :: PrL(nSnp), PrPA(3), PrG(3),PrXZ(3,3,3), PrB(3), PrGx(3), &
    PrE(3), PrV(3)

  ! assume all ancestor checks etc. have been done
  curGP = 0 
  do z=1,2
    curGP(:,z) = getPar(Parent(A,z), z)
  enddo
  
  if (Sex(B)<3) then
    m = Sex(B)
  else if (curGP(1,k) == 0) then
    m=1    ! doesn't really matter.
  else
    m = 2
  endif
  
  PrL = 0D0
  do l=1,nSnp
    call ParProb(l, curGP(3-m,k), 3-m, 0,0, PrG) 
    call ParProb(l, Parent(A,k), k, A, -4, PrPA) 
    call ParProb(l, B, 0, 0, 0, PrB)  
    call ParProb(l, curGP(3-k,3-m), 3-m, 0, 0, PrGx) 

    PrXZ = 0D0
    do x=1,3  ! PA(k)     
      do y=1,3  ! PA(3-k)
  !      if (Aselfed .and. x/=y)  cycle
        do z=1,3  !  PrG(3-m)
          do v=1,3
            PrV(v) = SUM(AKA2P(y,v,:) * PrGx) * AKA2P(x,v,z) * PrB(v)
          enddo
          PrXZ(x,y,z) = OKA2P(Genos(l,A),x,y) * PrPA(x) * PrG(z) * SUM(PrV)  
          do i=1, ns(-Parent(A,3-k),3-k)
            Ei = SibID(i,-Parent(A,3-k),3-k)
            if (nFS(Ei)==0)  cycle
            call ParProb(l, Parent(Ei,k), k, Ei, -1, PrE)  ! -1: exclude all FS of A
            do j=1,nFS(Ei)
              if (FSID(j,Ei)==A .or. FSID(j,Ei)==B)  cycle
              PrE = PrE * OKA2P(Genos(l,FSID(j,Ei)), :, y)
            enddo
            if (.not. all(PrE==1D0)) PrXZ(x,y,z) =  PrXZ(x,y,z) * SUM(PrE)             
          enddo    
        enddo
      enddo
    enddo
    PrL(l) = LOG10(SUM(PrXZ))
  enddo

  LL = SUM(PrL) + Lind(B)

  end subroutine PairdGP
  
end subroutine GGpairs

! ##############################################################################

subroutine FsibsGPs
! assign grandparents to full-sibling clusters
use Global
implicit none

integer :: x, j, fsx(maxSibSize), k,m,r, nCG(2,2), candGP(mxCP,2,2), SAB(2), ix, n
logical :: maybeGP(2,2), AncOK
double precision :: ALR, LRG, LLX(2,2), dx(maxSibSize), LRS

if (.not. (DoMtDif .or. any(AgePriorA(:,1,2) /= AgePriorA(:,1,3)) .or. &
  any(AgePriorA(:,2,2) /= AgePriorA(:,2,3)))) then   ! diff ageprior between MGM-PGM or MGP-PGP
  if (quiet==-1)  call rprint_status_tbl_no_dots()
  return   ! no way to distinguish between maternal and paternal grandparents
endif       

do x=1, nInd
  if (MODULO(x,200)==0)  call rchkusr()
  if (quiet==-1 .and. any(chunk_printdot_i==x)) call print_dot()
  if (nFS(x)==0)  cycle  ! not 'primary' sib of FS cluster
  if (.not. ALL(Parent(x,:) < 0) .or. ALL(parent(x,:)==0))  cycle
  if (all(Parent(x,:)==0) .and. BY(x) < 0)  cycle   ! high risk wrong way around   
  SAB = -parent(x,:)
  if (all(SAB/=0)) then
    if (ns(SAB(1),1) /= ns(SAB(2),2))  cycle  ! resolvable via SibGrandparent 
    if (any(GpID(:,SAB(1),1)/=0) .or. any(GpID(:,SAB(2),2)/=0))  cycle   ! resolvable via SibGrandparent 
!    if (ns(SAB(1),1)==1 .and. any(GpID(:,SAB(1),1)/=0) .and. any(GpID(:,SAB(2),2)/=0))  cycle
    if (all(BY(SibID(1:nS(SAB(1),1),SAB(1),1)) < 0))  cycle  ! high risk wrong way around + slow
  endif
 
  if (ALL(parent(x,:)==0)) then
    do k=1,2
      call NewSibship(x, 0, k)
    enddo
    SAB = -parent(x,:)
  endif
  
  fsx = 0
  fsx(1:nFS(x)) = FSID(1:nFS(x),x)
  nCG = 0  
  CandGP = 0
  
  do j=1, nInd
    maybeGP = .TRUE.
    if (ANY(fsx == j))  cycle
    do n=1,nFS(x)
      if (AgeDiff(fsx(n),j) <= 1) then
        maybeGP = .FALSE.
        exit
      endif
    enddo
    if (.not. any(MaybeGP))  cycle            
    do m=1,2
      if (m/=Sex(j) .and. Sex(j)<3)  maybeGP(m,:) = .FALSE.
    enddo
    if (Sex(j)/=2 .and. DoMtDif) then
      if (mtDif(x,j))  maybeGP(1,1) = .FALSE.   ! mat grandmother must have same mt haplo
    endif 
    if (.not. ANY(maybeGP))  cycle          
    do k=1,2
      call ChkAncest(j, sex(j), -SAB(k), k, AncOK)
      if (.not. AncOK)  maybeGP(:,k) = .FALSE.
    enddo    
    if (.not. ANY(maybeGP))  cycle 
    do k=1,2
      do m=1,2
        if (.not. maybeGP(m,k))  cycle
        call CalcAgeLR(-SAB(k), k, j,m, k,1,.TRUE.,ALR)
        if (ALR==impossible .or. ALR < 3.0*TF) then
          maybeGP(m,k) = .FALSE.
        endif
      enddo
    enddo
    if (.not. ANY(maybeGP))  cycle    
    if (ns(SAB(1),1)>1) then
      call QFSGP(j, Sex(j), SAB(1), 1, LRG)   ! TODO: sep subroutine for FS?
      if (LRG < TF)  cycle
    else
      call PairQHS(x, j, LRS)    
      if (LRS < TF)  cycle
      call LRGG(x,1,j,Sex(j),LRG)
      if (LRG < TA)  cycle
    endif
    
    do k=1,2
      do m=1,2
        if (nCG(m,k) == mxCP)  cycle
        if (maybeGP(m,k)) then
          nCG(m,k) = nCG(m,k) +1
          CandGP(nCG(m,k),m,k) = j
        endif
      enddo
    enddo
  enddo
  
  ! dummy GPs
  do m=1,2
    if (ns(SAB(m),m)==1)  exit                         
    do r=1, nC(m) 
      if (r==SAB(m)) cycle
      if (ANY(GpID(:,r,m) == SAB)) cycle
      maybeGP = .FALSE.   
      do k=1,2
        if (ncG(m,k) == mxCP) cycle
        if (nS(r,m)==1 .and. ANY(SibID(1:nS(SAB(k),k),SAB(k),k) == SibID(1,r,m))) cycle
        if (DoMtDif) then
          if (k==1 .and. m==1 .and. mtDif(SibID(1,SAB(k),k), SibID(1,r,m)))  cycle
        endif   
        if (m/=k .and. complx<2) then
          if (ALL(Parent(SibID(1:ns(SAB(k),k),SAB(k),k),m) == -r))  cycle
          if (ALL(Parent(SibID(1:ns(r,m),r,m),k) == -SAB(k)))  cycle
        endif
        call ChkAncest(-r,m, -SAB(k),k, AncOK)
        if (.not. AncOK)  cycle
        call CalcAgeLR(-SAB(k),k, -r,m, 0,1, .TRUE., ALR)
        if (ALR == impossible .or. ALR < 3.0*TF) cycle
        maybeGP(m,k) = .TRUE.
      enddo
      if (.not. ANY(maybeGP(m,:)))  cycle   
      if (hermaphrodites==0) then
        call QFSGP(-r, m, SAB(1), 1,  LRG) 
        if (LRG < TF*dble(MIN(nS(SAB(1),1), nS(r,m)))) cycle  ! conservative.
      endif      
      LLX = missing
      call PairUA(-SAB(1), -r, 1, m, LLX(1,1))
      if (LLX(1,1)>0) cycle
      call CalcU(-SAB(1),1, -r,m, LLX(1,2)) 
      if ((LLX(1,1) - LLX(1,2)) < nS(SAB(1),1)*TF) cycle
      call addFS(0, r, m, SAB(1), 1, LLX(2,1), ix, dx) 
      if ((MaxLL(LLX(:,1)) - LLX(1,2)) < TA)  cycle
      
      do k=1,2
        if (nCG(m,k) == mxCP)  cycle
        if (maybeGP(m,k)) then
          nCG(m,k) = nCG(m,k) +1
          CandGP(nCG(m,k),m,k) = -r
        endif
      enddo
    enddo
  enddo
  
  if (all(CandGP(:,1,1) == CandGP(:,1,2)) .and. all(CandGP(:,2,1) == CandGP(:,2,2))) then
    ! identical candidates for mat & pat side of full sibship, incl all 0
    do k=1,2
      if (ns(SAB(k),k)==1 .and. all(GpID(:,SAB(k),k)==0)) then
        call RemoveSib(x, SAB(k),k)
        call DoMerge(0, SAB(k),k)
      endif
    enddo     
    cycle   
  endif

  if (ANY(nCG>0)) then   
    do k=1,2  !2,1,-1
      if (all(nCG(:,k)==0))  cycle
      call SelectParent(-SAB(k), k, nCG(:,k), CandGP(:,:,k), .FALSE.)
    enddo
    call ChkGPs(SAB, CandGP)  ! drops GPs if could be GPs of 3-k
    
    do k=2,1,-1
      if (all(nCG(:,k)==0))  cycle
      if (any(GpID(:,SAB(k),k) /=0))  cycle
      call SelectParent(-SAB(k), k, nCG(:,k), CandGP(:,:,k), .FALSE.)
    enddo
    call ChkGPs(SAB, CandGP)
    
    do k=1,2
      if (ns(SAB(k),k)==1 .and. all(GpID(:,SAB(k),k)==0)) then
        call RemoveSib(x, SAB(k),k)
        call DoMerge(0, SAB(k),k)
      endif
    enddo
  endif  
enddo

end subroutine FsibsGPs

! ##############################################################################

subroutine ChkGPs(SAB, CandGP)
use Global
implicit none

integer, intent(IN) :: SAB(2), candGP(mxCP,2,2)
integer :: k,m
logical :: MaybeOpp(2,2)

! if both pairs could have been assigned as GP via parent 3-k, drop all
maybeOpp = .TRUE.
do k=1,2      
  do m=1,2
    if (GpID(m,SAB(k),k)/=0) then
      if (ANY(CandGP(:,m,3-k) == GpID(m,SAB(k),k))) then
        maybeOpp(k,m) = .TRUE.
      else
        maybeOpp(k,m) = .FALSE.
      endif
    endif
  enddo
enddo 
do k=1,2
  if (ALL(maybeOpp)) then
    do m=1,2
      call setPar(-SAB(k),k, 0,m)
    enddo
  endif 
enddo

end subroutine ChkGPs

! ##############################################################################

subroutine Qadd(A, SB, kB, LR)
use Global
implicit none

integer, intent(IN) :: A, SB, kB
double precision, intent(OUT) :: LR
integer :: l, x
double precision :: PrL(nSnp), PrX(3)

PrL = 0D0
do l=1,nSnp
  do x=1,3
    PrX(x) = OKAP(Genos(l,A), x, l) * DumP(x,l,SB,kB) / AHWE(x,l)
  enddo   ! simple LL identical for HS and GP
  PrL(l) = LOG10(SUM(PrX))
enddo
LR = SUM(PrL)

end subroutine Qadd

! #####################################################################

subroutine QGP(A, kA, SB, kB, LR)  ! A indiv or dummy, GP of SB
use Global
implicit none

integer, intent(IN) :: A, kA, SB, kB
double precision, intent(OUT) :: LR
integer :: l, x
double precision :: PrLR(nSnp), PrX(3,2), PrA(3)

if (ns(SB,kB)==1 .and. A>0) then
  call PairQHS(SibID(1,SB,kB), A, LR)
else
  PrLR = 0D0
  do l=1,nSnp
    call ParProb(l, A, kA, 0, 0, PrA)  ! no effect on time vs. LindX/DumP 1x
    do x=1,3               
      PrX(x,1) =XPr(1,x,l,SB,kB) * SUM(AKAP(x,:,l) * PrA)
      PrX(x,2) =XPr(1,x,l,SB,kB) * AHWE(x,l)
    enddo
    PrLR(l) = LOG10(SUM(PrX(:,1))) - LOG10(SUM(PrX(:,2))) 
  enddo
  LR = SUM(PrLR)
endif

end subroutine QGP

! #####################################################################

subroutine chkPO(A, SB, kB, LR)  ! A replaces dummy SB?
use Global
implicit none

integer, intent(IN) :: A, SB, kB
double precision, intent(OUT) :: LR
integer :: l, x, sib1
double precision :: PrLR(nSnp), PrX(3,2), LL(2), PrA(3)

if (ns(SB,kB)==1) then
  sib1 = SibID(1,SB,kB)
  call CalcU(sib1,kB,A,kB, LL(1))
  call PairPO(sib1, A, kB, 1, LL(2))
  LR = LL(2) - LL(1)
else
  PrLR = 0D0
  do l=1,nSnp
    call ParProb(l, A, kB, 0, 0, PrA)
    do x=1,3
      PrX(x,1) = XPr(1,x,l,SB,kB) * XPr(2,x,l,SB,kB)
      PrX(x,2) = XPr(1,x,l,SB,kB) * PrA(x)
    enddo
    PrLR(l) = LOG10(SUM(PrX(:,2))) - LOG10(SUM(PrX(:,1)))
  enddo
  LR = SUM(PrLR)
endif

end subroutine chkPO

! #####################################################################

subroutine QFS(A, SB, kB, LR)   ! only when Complx==0
use Global
implicit none

integer, intent(IN) :: A, SB, kB
double precision, intent(OUT) :: LR
integer :: l, x, y
double precision :: PrLR(nSnp), PrXY(3,3,2), PrY(3)

PrLR = 0D0
do l=1,nSnp
  call ParProb(l, Parent(SibID(1,SB,kB), 3-kB), 3-kB, -1, 0, PrY)
  do x=1,3
    do y=1,3
      PrXY(x,y,1) = OKA2P(Genos(l,A),x,y) * DumP(x,l,SB,kB) * PrY(y)
      PrXY(x,y,2) = OKA2P(Genos(l,A),x,y) * AHWE(x,l) * AHWE(y,l)
    enddo
  enddo  
  PrLR(l) = LOG10(SUM(PrXY(:,:,1))) - LOG10(SUM(PrXY(:,:,2)))
enddo
LR = SUM(PrLR)

end subroutine QFS

! #####################################################################

subroutine QFSGP(A, kA, SB, kB, LR)  ! A indiv or dummy, GP of SB, all B's are FS
use Global
use CalcLik          
implicit none

integer, intent(IN) :: A, kA, SB, kB
double precision, intent(OUT) :: LR
integer :: l, x,y, i
double precision :: PrLR(nSnp), PrXY(3,3,2), PrA(3), PrI(3,3)

i = FSID(maxSibSize+1, SibID(1,SB,kB))

if (ns(SB,kB)==1 .and. A>0) then
  call PairQHS(SibID(1,SB,kB), A, LR)
else
  PrLR = 0D0
  do l=1,nSnp
    call ParProb(l, A, kA, 0, 0, PrA)
    PrI = FSLik(l,i)                
!    call ParProb(l, Parent(i,3-kB),3-kB,-1,0, PrY)  ! GPs only
    do x=1,3
      do y=1,3
        PrXY(x,y,1) = PrI(x,y) * SUM(AKAP(x,:,l) * PrA) * AHWE(y,l)
        PrXY(x,y,2) = PrI(x,y) * AHWE(x,l) * AHWE(y,l)
      enddo
    enddo
    PrLR(l) = LOG10(SUM(PrXY(:,:,1))) - LOG10(SUM(PrXY(:,:,2))) 
  enddo
  LR = SUM(PrLR)
endif

end subroutine QFSGP

! #####################################################################

subroutine CheckRel(A, kA, B, kB, focalIN, LLg, LL)
use Global
implicit none

integer, intent(IN) :: A, kA, B, kB, focalIN
double precision, intent(OUT) :: LLg(7), LL(7)
logical:: FSJ  !do separately?
integer :: k, focal

focal = focalIN
FSJ = .FALSE.
LLg = missing
LL = missing
if (A==0 .or. B==0) then
  call Erstop("CheckRel A or B null ", .TRUE.)
else if (A==B .and. (A>0 .or. kA==kB)) then
  call Erstop("CheckRel A==B ", .TRUE.)
else if (A > 0 .and. B > 0) then
  if (kA == 0 .and. kB==0) then
    call Erstop("CheckRel kA == kB == 0!", .TRUE.)
  else if (kB /= 0) then 
    k = kB
  else if (kA /= 0) then
    k = kA
  endif
  call CheckPair(A, B, k, focal, LLg, LL)  
else if (A > 0 .and. B < 0) then
  if (kB<1 .or. kB>2)  call Erstop( "CheckRel A>0, B<0, invalid kB", .TRUE.)
  if (focal==0)  call Erstop("CheckRel focal == 0!", .TRUE.)
  if (focalIN==1)  focal =  3  ! -B parent of A -> B's HS of A
  call CheckAdd(A, -B, kB, focal, LLg, LL)
  if (focalIN==1 .or. focalIN==6) then    ! 7 default for getPairLL   #  .or. focalIN==7
    if (focalIN==6 .and. Parent(A,3-kB)==0 .and. Complx/=0) then  ! called by CalcCandParLL, want single vs parent-pair
      LLg(2) = 333D0
      LL(2) = 333D0
    endif
    LLg = ReOrderAdd(LLg)
    LL = ReOrderAdd(LL) 
  endif
else if (A < 0 .and. B > 0) then
  if (kA<1 .or. kA>2)  call Erstop("CheckRel A<0, B>0, invalid kA", .TRUE.)
  call CheckAdd(B, -A, kA, focal, LLg, LL)
else if (A < 0 .and. B < 0) then
  if (kA<1 .or. kA>2)  call Erstop("CheckRel A<0, B<0, invalid kA", .TRUE.)
  if (kB<1 .or. kB>2)  call Erstop( "CheckRel A<0, B<0, invalid kB", .TRUE.)
  ! note: focal=1: merge, focal=4: SB parent of SA. FSJ: full-sib merge
  call CheckMerge(-A, -B, kA, kB, focal, LLg, LL, FSJ)
endif

contains
  function ReOrderAdd(LL) result(LLtmp)
    ! reorder output from CheckAdd for compatibility with CheckPair (for POZ)
    double precision, intent(IN) :: LL(7)
    double precision :: LLtmp(7)

    LLtmp = missing
    LLtmp(1) = MaxLL(LL(2:3))
    LLtmp(2) = LL(5)
    if (Complx/=0)  LLtmp(3) = LL(6)  ! not: monogamous
    ! LLtmp(4) =   TODO !!!
    LLtmp(5) = LL(4)   ! ? not technically correct, but ... 
    if (complx==0) then
      LLtmp(6) = MaxLL(LL((/1,6/)))
    else
      LLtmp(6) = LL(1)
    endif
    LLtmp(7) = LL(7) 

  end function ReOrderAdd
end subroutine CheckRel

! #####################################################################

subroutine CheckAdd(A, SB, k, focal, LLg, LL)
use Global
use CalcLik          
implicit none

integer, intent(IN) :: A, SB, k, focal
double precision, intent(OUT) :: LLg(7), LL(7)
double precision :: LRHS, ALR(7), LLPH(2), ALRH(2), LLAU(2,3), ALRAU(2,3), LLUi, &
  LLC, LLz(7), ALRz(7), LLM(3), LLp(7), LLpg(7), LLFH(3), LLPX(2,2), dx(maxSibSize), &
  ALRq, LLHH(2,2), ALRtmp, LHH(3), LHH2, LLy(2,2), LLpo(ns(SB,k),2), ALRpo(ns(SB,k),2), &
  LLgp(ns(SB,k),3), ALRgp(ns(SB,k),3), LLfs(3,2), LLdGP(ns(SB,k)), LLOP, LLPA ! , LLHSPO(ns(SB,k),2)
integer :: x, y, FSPar, i, ParTmp(2), OpPar(maxSibSize), nop, fsi, ix, m, Bi, sib1, curpar(2)    
logical :: AncOK, fclsib, MaybeOpp, ParOK, ParAisBClone

LL = missing
LLg = missing
ALR = missing

! quick check
LRHS = missing
call Qadd(A, SB, k, LRHS)  ! 2nd degree relatives vs unrelated
if (LRHS < MIN(TF*2, TF*nS(SB,k)) .and. (focal/=4 .and. focal/=7 .and. focal/=6)) return
  
if (Sex(A)<3 .and. Sex(A)/=k) then
  LL(1) = impossible
  if (focal==1)  return
endif
  
if (focal==1) then              
  call CalcAgeLR(-SB,k, A,k, 0,-1, .TRUE., ALR(1))
  if (ALR(1)==impossible) then
    LL(1) = impossible
    return
  endif
else if (focal==2 .or. focal==3) then
  if (all(GpID(:,SB,k)==0) .and. ns(SB,k)>0) then
    call calcALR_addsib(A,SB,k,3,ALR(3))
  else
    call CalcAgeLR(A,Sex(A), -SB,k, 0,1, .TRUE., ALR(3))
  endif
  ALR(2) = ALR(3)   ! updated below if a FS found                                   
  if (ALR(3)==impossible) then
    LL(2:3) = impossible
    return
  endif
  do Bi=1, ns(SB,k)
    if (ns(SB,k)==0)  exit                      
    if ((Parent(A,3-k)==Parent(SibID(Bi,SB,k),3-k) .and. Parent(A,3-k)/=0) .or. &
      Complx==0) then
      call PairQFS(A, SibID(Bi,SB,k), LRHS)
    else
      call PairQHS(A, SibID(Bi,SB,k), LRHS)
    endif
    if (LRHS < 4*TF) then
      LL(2:3) = impossible
      return
    endif
  enddo
endif

call ChkAncest(A,k,-SB,k, AncOK)
if (.not. AncOK) then
  LL(1) = impossible
  LL(4) = impossible
  if (focal==1 .or. focal==4)  return
endif
if (Parent(A,k)/=0) then
  AncOK = .FALSE.
else
  call ChkAncest(-SB,k, A,k, AncOK)
endif
if (.not. AncOK) then
  LL(2:3) = impossible
  if (focal==2 .or. focal==3)  return
endif

! mt haplotype
if (DoMtDif) then
  if (k==1 .and. mtDif(SibID(1,SB,k), A)) then
    LL(1:3) = impossible
    if (Sex(A)==1)  LL(4) = impossible
    if (LL(focal) == impossible)  return 
  endif
endif 
        
call CalcU(A,k, -SB, k, LLg(7))   ! unrelated
LL(7) = LLg(7)

fsi=0
if (focal <4 .and. ns(SB,k) > 0) then
  if (focal==1)  call AddParent(A, SB, k, LLg(1))
  if (focal==2)  call AddFS(A, SB, k,0,k, LLg(2), fsi, dx)
  if (focal==3 .and. Complx/=0)  call AddSib(A, SB, k, LLg(3))
  do x=1,3
    if (focal==x) then
      if ((LLg(focal) > 0D0 .or. LLg(focal) - LL(7) < TA)) then
        LL(x) = addALR(LLg(x), ALR(x))
        if (focal ==3 .and. (Parent(A,3-k)==0 .or. Complx==0 .or. &
         ANY(Parent(SibID(1:ns(SB,k),SB,k), 3-k) == Parent(A,3-k)))) then
          call AddFS(A, SB, k,0,k, LLg(2), fsi, dx)
          if ((ALL(LLg(2:3) >0D0) .or. MaxLL(LLg(2:3)) - LL(7) < TA))  return
        else
          return
        endif
      endif
    endif
  enddo
endif 

fclsib = (focal==2 .or. focal==3 .or. focal==6)
call getFSpar(SB, k, .TRUE., FSpar) 

!=======

if (LL(1)/=impossible) then   
  if (ALR(1)==missing)  call CalcAgeLR(-SB,k, A,k, 0,-1, .TRUE., ALR(1))                         
  if (LLg(1)==missing .and. ALR(1)/=impossible)  call AddParent(A, SB, k, LLg(1))  ! A parent of SB 
endif

if (LL(3)/=impossible) then                           
  if (ALR(3)==missing)  call CalcAgeLR(A,Sex(A), -SB,k, 0,1, .TRUE., ALR(3))  ! SB parent of A
  if (LLg(3)==missing .and. (Complx>0 .or. ns(SB,k)==0) .and. ALR(3)/=impossible) then
      call AddSib(A, SB, k, LLg(3))
  endif                              
endif     

if (LLg(2)==missing .and. ns(SB,k)>0 .and. ALR(2)/=impossible)  call AddFS(A, SB, k,0,k, LLg(2), fsi, dx)
if (ns(SB,k)==1 .and. .not. any(Parent(A,:)==SibID(1,SB,k)))  fsi = SibID(1,SB,k)                                     

ALR(2) = ALR(3)
if (fsi/=0 .and. LL(2)/=impossible) then
  if (parent(A,3-k)==0 .or. Parent(A,3-k)==Parent(fsi,3-k)) then
    call CalcAgeLR(A,Sex(A), fsi,k, 0,2, .TRUE., ALR(2))
  endif
endif
  
if (ns(SB,k)==0 .and. LL(3)/=impossible) then
  LLg(2) = LLg(3)
  ALR(2) = ALR(3)
  if (Complx==0) then
    LLg(3) = missing
  endif
endif

if (nYears>2 .and. LL(4)/=impossible) then
  call CalcAgeLR(-SB,k, A,Sex(A), 0,1, .TRUE., ALR(4))  ! A parent of SB
  if (ALR(4)/=impossible) then
    call AddGP(A, SB, k, LLg(4))
  endif
endif

do x=1,4
  if (LL(x)==impossible) then
    LLg(x) = impossible
  else
    LL(x) = addALR(LLg(x), ALR(x))
  endif
enddo

! monogamous
if (Complx==0 .and. Mate(A)/=0 .and. any(GpID(:,SB,k)==Mate(A))) then
  return
endif 

!~~~~~~~~~~~~
!ParAisBClone = .FALSE.
!if (Hermaphrodites/=0 .and. Parent(A,3-k)<0) then
!  if (DumClone(SB,k) == -Parent(A,3-k))  ParAisBClone = .TRUE.  
!endif

if (Hermaphrodites/=0 .and. focal/=7) then     ! TODO: double check if this if-then is as intended
  if (ns(SB,k)==1) then
    do Bi=1, ns(SB,k)
      LLPH = missing
      ALRH = missing
      call CalcAgeLR(A, Sex(A), SibID(Bi,SB,k),k, 0,1, .TRUE., ALRH(1))
      call CalcAgeLR(A, Sex(A), SibID(Bi,SB,k),k, 3,4, .TRUE., ALRH(2))
      if (ALRH(1)/=impossible) then
        call PairPO(A, SibID(Bi,SB,k), k, 0, LLPH(1))
      endif
      if (ALRH(2)/=impossible) then
        call PairGP(A, SibID(Bi,SB,k), k, 0, LLPH(2))
      endif
      do x=1,2
        if (LLPH(x) < 0)  LLPH(x) = LLPH(x) - Lind(SibID(Bi,SB,k)) + CLL(SB,k)
      enddo
      if (focal==1) then
        LLg(6) = MaxLL((/LLg(6), LLPH(1)/))
        LL(6) = MaxLL((/LL(6), addALR(LLPH(1), ALRH(1))/))      
      else
        if ((LLPH(1) > LLg(1) .and. LLPH(1)<0) .or. LLg(1)>0) then
          LLg(1) = LLPH(1)
          LL(1) = addALR(LLPH(1), ALRH(1))
        endif
      endif
      if ((LLPH(2) > LLg(4) .and. LLPH(2)<0) .or. LLg(4)>0) then
        LLg(4) = LLPH(2)
        LL(4) = addALR(LLPH(2), ALRH(2))
      endif                        
    enddo
  endif
endif

!~~~~~~~~~~~~
LLAU = missing
ALRAU = missing
! FA 1: A FS of SB
! FA 2: SB GP of A, SB monogamous, SB's partner (thus) also GP of A
call CalcAgeLR(-SB,k, A,Sex(A), 0,2, .TRUE., ALRAU(1,3))
if (ALRAU(1,3)/=impossible .and. (Complx==0 .or. &
  (.not. (focal==4 .and. ALL(Parent(A,:)/=0)) &
   .and. .not. (focal==7 .and. GpID(3-k,SB,k)/=0))) .and. ns(SB,k)>0) then 
  call pairUA(-SB, A, k, 3, LLAU(1,3))
endif

if (Parent(A,k)/=0) then
  call CalcAgeLR(Parent(A,k),k, -SB,k, 0,1, .TRUE., ALRAU(2,3))
else
  call CalcAgeLR(A,Sex(A), -SB,k, k,4, .TRUE., ALRAU(2,3)) 
endif

if (ALRAU(2,3)/=impossible .and. .not. (focal==7 .and. fsi/=0 .and. Parent(A,3-k)==0) .and. &
  ns(SB,k)>0) then
  ! not considered during CalcParLLR: true other-parent is GP when m=1 (single parent LLR)
  if (ns(SB,k)==1) then
    call pairUA(A, SibID(1,SB,k), k, 3, LLAU(2,3))
  else if (Parent(A,k) < 0 .and. FSpar/=0 .and. all(parent(SibID(1:ns(SB,k),SB,k), 3-k)==FSpar)) then
    call pairUA(Parent(A,k), SibID(1,SB,k), k, 3, LLAU(2,3))
    if (LLAU(2,3) < 0) then
      call CalcU(A,k, SibID(1,SB,k),k, LLUi)
      LLAU(2,3) = LLAU(2,3) - LLUi + LLg(7)
    endif 
  else
    call addFA(A, SB, k, LLAU(2,3))
  endif
endif

LLg(5) = MaxLL(LLAU(:,3))
LL(5) = MaxLL((/ addALR(LLAU(1,3),ALRAU(1,3)), addALR(LLAU(2,3),ALRAU(2,3)) /))

LLC = missing 
if (complx==2 .and. ns(SB,k)>0 .and. (fclsib .or. focal==7) .and. LL(2)<0D0 .and. &
 Parent(A,3-k)==FSpar .and. (MaxLL(LLAU(:,3)) - MaxLL(LL(2:3)) > -TA)) then
  call FSHC(A, -SB, k, LLC)  ! Full sib & half-cousin
  if (LLC >LLg(2) .and. LLC<0) then
    LLg(2) = LLC
    LL(2) = addALR(LLg(2), ALR(2))  ! no cousin ageprior yet
  endif
endif

!~~~~~~~~~~~~
! LLg(6) HA (other 3rd degree rel: LLz further down)
if (Complx>0) then
  do x=1,2
    ! HA 1: A HS of SB:
    call CalcAgeLR(-SB,k, A,Sex(A), x,3, .TRUE., ALRAU(1,x))
    if (ALRAU(1,x)/=impossible .and. ns(SB,k)>0 .and. &
     .not. (focal==7 .and. x==3-k .and. Parent(A,3-k)==0) .and. &
         .not. (focal==4 .and. Parent(A,x)>0 .and. GpID(x,SB,k)==0)) then  ! else conflict with CalcCandParLL
      if (Parent(A,x)<0) then  ! shouldn't matter but does; TODO fix bug in pairUA
        call PairUA(-SB, Parent(A,x), k, x, LLAU(1,x))
        call CalcU(-SB,k, Parent(A,x),x, LLPA)
        if (LLAU(1,x)<0)  LLAU(1,x) = LLAU(1,x) - LLPA + LLg(7)
      else
        call pairUA(-SB, A, k, x, LLAU(1,x))
      endif
    endif   
    
    ! HA 2: SB GP of A 
    if (Parent(A,x)/=0) then
      call CalcAgeLR(Parent(A,x),x, -SB,k, 0,1, .TRUE., ALRAU(2,x))
    else
      call CalcAgeLR(A,Sex(A), -SB,k, x,4, .TRUE., ALRAU(2,x))
    endif
    if (ALRAU(2,x)/=impossible .and. .not. (focal==7 .and. x==3-k)) then    
      if (Parent(A,x)<0) then  ! shouldn't matter but does; TODO fix bug in pairUA
        call PairUA(Parent(A,x),-SB, x, k, LLAU(2,x))
        call CalcU(-SB,k, Parent(A,x),x, LLPA)
        if (LLAU(2,x)<0)  LLAU(2,x) = LLAU(2,x) - LLPA + LLg(7)
      else
        call pairUA(A, -SB, x, k, LLAU(2,x))
      endif   
    endif 
  enddo
  LLg(6) = MaxLL(RESHAPE(LLAU(:,1:2), (/2*2/) ))
  do x=1,2
    do y=1,2
      LLAU(y,x) = addALR(LLAU(y,x), ALRAU(y,x))
    enddo
  enddo
  LL(6) = MaxLL(RESHAPE(LLAU(:,1:2), (/2*2/) ))
endif     

LLz = missing
ALRz = missing
if ((LL(focal)<0D0 .and. LL(focal)>=LL(7)) .or. focal==4 .or. LL(6)>0D0 .or. LL(6)<LL(7)) then
  if (any(GpID(:,SB,k)/=0)) then
    do x=1,2
      if (GpID(x,SB,k)==0) then
        call CalcAgeLR(-SB,k, A,Sex(A), x,4, .TRUE., ALRz(1))
      endif
    enddo
  else
    call CalcAgeLR(-SB,k, A,Sex(A), 3,4, .TRUE., ALRz(1))
  endif      
  if (ALRz(1)/=impossible) then
    call AddGGP(A, SB, k, LLz(1))
  endif      
  if (nS(SB,k)>0) then
    do x=1,2
      if (focal==6 .and. parent(A,x)/=0)  cycle                                               
      call CalcAgeLR(A,k, -SB,k, x,5, .TRUE., ALRz(x+1))     
      if (ALRz(x+1)==impossible) then
        LLz(x+1) = impossible
      else
        call ParentHFS(A, 0,x, SB, k,3, LLz(x+1))
      endif
    enddo    
  endif
  if (Complx==2 .or. Complx==0) then
    do x=1,2   ! as checkmerge: full great-uncle  (2x 1/4)
      call CalcAgeLR(-SB,k, A,Sex(A), x,5, .TRUE., ALRz(3+x))
      if (ALRz(3+x) == impossible) then
        LLz(3+x) = impossible
      else 
        if (GpID(x,SB,k) <0 .and. .not. any(parent(SibID(1:ns(SB,k),SB,k),x) == GpID(x,SB,k))) then 
          call PairUA(GpID(x,SB,k), A, x, 3, LLz(3+x))
          if (LLz(3+x) < 0) then
            LLz(3+x) = LLz(3+x) - CLL(-GpID(x,SB,k), x) + CLL(SB,k)  
          endif
        else if (GpID(x,SB,k)==0) then   ! else cond. indep.
          call addGAU(A, SB, k, x, LLz(3+x))    
        endif
      endif
    enddo
  endif
  if (ns(SB,k)>0) then
    sib1 = SibID(1,SB,k)
    call PairCC(A, sib1, k, LLz(6))  ! full cousins
    if (LLz(6) < 0D0) then
      call CalcU(A,k, SibID(1,SB,k),k, LLUi)
      LLz(6) = LLz(6) - LLUi + LL(7)
    endif
    if (FSpar<0 .and. parent(A,k)==0 .and. all(GpID(:,SB,k)==0)) then
      if (ns(SB,k) == ns(-FSpar,3-k) .and. all(parent(SibID(1:ns(SB,k),SB,k), 3-k)==FSpar)) then
        do i=1,ns(SB,k)
          if (nFS(SibID(i,SB,k))==0)  cycle
          call pairDHC(A,k, SibID(i,SB,k), .TRUE., LLz(7))  ! double half cousins
        enddo
      endif
    endif
  endif
  ALRz(6:7) = 0D0  ! no ALR for cousins yet
  LLg(6) = MaxLL((/LLg(6), LLz/))
  do x=1,6
    LLz(x) = addALR(LLz(x), ALRz(x))
  enddo
  LL(6) = MaxLL((/LL(6), LLz/))
endif

! A parent of opp. sibship
LLOP = missing
if (FSpar<0 .and. (sex(A)>3 .or. Sex(A)/=k)) then
  call CalcAgeLR(FSpar,3-k, A,3-k, 0,-1, .TRUE., ALRq)
  if (ALRq /= impossible)  call AddParent(A, -FSpar, 3-k, LLOP)
  if (focal==1 .and. LLOP<0d0 .and. LLOP > LLg(6)) then
    LLg(6) = LLOP
    LL(6) = addALR(LLOP, ALRq)
  else if (LLOP<0d0 .and. (LLOP > LLg(1) .or. LLg(1) > 0d0)) then
    LLg(1) = LLOP
    LL(1) = addALR(LLOP, ALRq)
  endif  
endif    

LLM = missing    
LLp = missing
LLFH = missing   
LLPX = missing
MaybeOpp = .FALSE.
if (complx>0 .and. fclsib .and. hermaphrodites/=2 .and. &
 abs(MaxLL(LL(2:3)) - MaxLL(LL)) < 0.01 .and. Parent(A,3-k)==0 .and. ns(SB,k)>0) then 
  if (abs(MaxLL(LL)-LL(2)) < 0.01 .and. fsi/=0) then
    !fsi = ID of putative full sib of A within SB, returned by AddFS() 
    if (ns(SB,k)>1) then
!      call PairFullSib(A, fsi, LLM(1))  
      call CalcU(A, k, fsi, k, LLM(3)) 
    else
!      LLM(1) = LLg(2)
      LLM(3) = LLg(7)
    endif
    call PairHalfSib(A, fsi, 3-k, LLM(2))     
    if ((LLM(2) - LLM(3)) - (LLg(2) - LLg(7)) > TA) then
      LL(2) = MaybeOtherParent   ! more likely to be HS via 3-k
    endif
  endif
  
  MaybeOpp = .TRUE.
  if (FSpar > 0) then
    if (FSpar/=A) then
       call CheckPair(A, FSpar, 3-k, 1, LLpg, LLp)    !! DANGER !!!
      if (LLp(1)<0 .and. (LLp(1) - MaxLL(LLp)) > TF) then  
        LL(2:3) = MaybeOtherParent  ! FSpar plausible parent of A
      endif
    endif
  else if (FSpar==0 .and. ANY(Parent(SibID(1:ns(SB,k), SB, k), 3-k) < 0)) then 
    ! get unique opposite-sex dummy parents
    OpPar = 0
    nop = 0
    do i=1, nS(SB,k)
      if (ANY(OpPar == Parent(SibID(i,SB,k), 3-k)))  cycle
      nop = nop +1
      OpPar(nop) = Parent(SibID(i,SB,k), 3-k)
    enddo
    if (ANY(OpPar(1:nop) > 0) .or. nop > 2) then
      MaybeOpp = .FALSE.       
    endif
    if (MaybeOpp .and. nop==2) then
      call CalcU(OpPar(1), 3-k, OpPar(2), 3-k, LLM(1))
      call MergeSibs(-OpPar(1), -OpPar(2), 3-k, LLM(2))
      if ((LLM(2) - LLM(1)) < TF*nS(SB,k))  MaybeOpp = .FALSE.
    endif
    if (nop>0 .and. MaybeOpp) then
      do x=1, nop
        if (OpPar(x) > 0)  cycle
        call CalcU(A, 3-k, OpPar(x), 3-k, LLM(1))
        call AddSib(A, -OpPar(x), 3-k, LLM(2))
        if (LLM(2)<0D0 .and. (LLM(2) - LLM(1)) - (LLg(2) - LLg(7)) > TA*nS(SB,k)) then
          LL(2) = MaybeOtherParent   ! more likely to be added to opposing sibship only. 
        endif
        if (LLM(2)<0D0 .and. (LLM(2) - LLM(1)) - (LLg(3) - LLg(7)) > TA*nS(SB,k)) then
          LL(3) = MaybeOtherParent  
        endif
        if (LL(2)==MaybeOtherParent .and. LL(3)==MaybeOtherParent)  exit 
      enddo
    endif 
  else if (FSpar < 0) then
    call Qadd(A, -FSpar, 3-k, LLM(1))  ! 2nd degree relatives vs unrelated    
    if (LLM(1) < TF*nS(-FSpar,3-k))  MaybeOpp = .FALSE.
    call CalcAgeLR(A,Sex(A), FSpar,3-k, 0,1, .TRUE., ALRq)
    if (ALRq==impossible)  MaybeOpp = .FALSE.
  endif
  if (MaybeOpp .and. FSpar < 0) then
    LLM = missing
    call AddFS(A, -FSpar, 3-k,0,3-k, LLM(1), ix, dx)
    call AddSib(A, -FSpar, 3-k, LLM(2))
    call CalcU(A, 3-k, FSpar, 3-k, LLM(3))
    if (LLM(2) < 0D0) then
      if (complx>0) then
         if ((LLM(2) - LLM(3)) - (LLg(3) - LLg(7)) > TA*dble(MAX(nS(SB,k),nS(-FSpar,3-k)))) then
          LL(3) = MaybeOtherParent  
        endif
      endif
      if (LLM(1) < 0 .and. (LLM(1) - LLM(2)) > 2*TA) then
        if (Complx==2) then  ! HS + parents FS/PO?
          curPar = Parent(A,:)                    
          call setParTmp(A, Sex(A), -SB, k)
          call PairUA(A, FSpar, 3-k, 3-k, LLPX(1,1))  ! HS + HA
          call ParentHFS(A, 0,3-k,-FSpar, 3-k,3, LLPX(1,2))  ! HS + FC
          call setParTmp(A, Sex(A), curPar(k), k)
          call setParTmp(A, Sex(A), FSpar, 3-k)     ! check done by PairFullSib
          call PairUA(A, -SB, k, k, LLPX(2,1))  ! HA + HS
          call ParentHFS(A, 0,k,SB, k,3, LLPX(2,2))  ! FC + HS
          call setParTmp(A, Sex(A), curPar(3-k), 3-k)
          if ((LLg(2) - LLg(7)) - (MaxLL(LLPX(1,:)) - LLM(3)) < TA) then
            LL(2) = MaybeOtherParent
          endif
          if ((MaxLL(LLPX(1,:)) - LLM(3)) > (LLg(3) - LLg(7)) .and. &
           (MaxLL(LLPX(1,:)) - MaxLL(LLPX(2,:))) > TA .and. MaxLL(LLPX(1,:))<0D0) then
            LLg(3) = MaxLL(LLPX(1,:)) - LLM(3) + LLg(7)  
            LL(3) = addALR(LLg(3), ALR(3))
          else if (((MaxLL(LLPX(2,:)) - LLM(3)) - (LLg(3) - LLg(7))) > TA .and. &
            MaxLL(LLPX(2,:))<0D0) then  ! MAX(nS(SB,k),nS(-FSpar,3-k))
            LL(3) = MaybeOtherParent
          endif
        else
          LL(2) = LL(2)
        endif
      else if (LLM(3)<0 .and. (LLM(2) -LLM(3)) >2*TA .and. complx>0) then
        LL(2:3) = MaybeOtherParent  ! as likely to be added to opp. parent
      endif
    endif
  endif
  if ((FSpar <0 .or. (FSpar==0 .and. all(Parent(SibID(1:ns(SB,k),SB,k),3-k)==0))) .and. &
   LL(2)<0 .and. Complx==2 .and. ns(SB,k)>0) then  
    sib1 = SibID(1,SB,k)
    call calcU(A,k,sib1, k, LLFH(1))
    if (ALRAU(1,3)/=impossible) call pairFAHA(sib1, A, .TRUE., LLFH(2))
    if (ALRAU(2,3)/=impossible) call pairFAHA(A, sib1, .TRUE., LLFH(3)) 
    WHERE(LLFH(2:3)<0)  LLFH(2:3) = LLFH(2:3) - LLFH(1) + LLg(7) 
    if (ANY(LLFH(2:3)<0) .and. MaxLL(LLFH(2:3)) > LLg(5)) then
      LLg(5) = MaxLL(LLFH(2:3))
      LL(5) = MaxLL((/ addALR(LLFH(2),ALRAU(1,3)), addALR(LLFH(3),ALRAU(2,3)) /))
    endif
  endif
endif

LLHH = missing
if (ANY(LL(2:3)<0) .and. ANY(LLAU<0) .and. fclsib .and. complx==2 .and. &
  fsi/=0 .and. hermaphrodites/=2 .and. ns(SB,k)>0) then                                                  
    call CalcAgeLR(A,Sex(A), fsi,3-k, 3, 6, .TRUE., ALRtmp)
    if (ALRtmp /= impossible) then
      call PairHSHA(A, fsi, k, LLHH(1,k), (/.FALSE., .TRUE./))  !HS via k, & parent A is HS of B via 3-k
      call PairHSHA(fsi, A, k, LLHH(2,k), (/.TRUE., .FALSE./))    
    endif  
 if ((Parent(fsi,3-k)/=0 .or. Parent(A,3-k)/=0) .and. &  ! else symmetrical 
  .not. (focal==6 .and. Parent(A,3-k)==0 .and. Parent(fsi,3-k)/=0)) then  ! else mucks up CalcCandPar  
    call CalcAgeLR(A,Sex(A), fsi,k, 3, 6, .TRUE., ALRtmp)  
    if (ALRtmp /= impossible) then
      call PairHSHA(A, fsi, 3-k, LLHH(1,3-k), (/.FALSE., .TRUE./))   !HS via 3-k, & parent A is HS of B via k
      call PairHSHA(fsi, A, 3-k, LLHH(2,3-k), (/.TRUE., .FALSE./))
    endif
  endif 
  if (any(LLHH < 0d0)) then 
    if (FSpar==0 .or. any(parent(SibID(1:ns(SB,k),SB,k), 3-k)/=FSpar)) then  ! not all in SB are FS of fsi
      call CalcU(A, k, fsi, k, LLUi)
      WHERE(LLHH < 0d0)  LLHH = LLHH - LLUi + LLg(7)    
    endif
    if (any(LLHH(:,k) <0D0) .and. MaxLL(LLHH(:,k)) > LLg(3)) then  ! (MaxLL(LLHH(:,k)) - MaxLL(LLHH(:,3-k))) > TA .and. &
      LLg(3) = MaxLL(LLHH(:,k))  ! HS via k, with or without HA  (and not HS via 3-k; only run if fclsib)
      LL(3) = addALR( addALR(LLg(3),ALR(3)), Maxval(ALRAU(1:2,3-k)))   
    endif
    if (any(LLHH(:,3-k) <0D0) .and. MaxLL(LLHH(:,3-k)) > LLg(6)) then
      LLg(6) = MaxLL(LLHH(:,3-k))  ! HA via k, with or without HS via 3-k 
      call CalcAgeLR(A,Sex(A), fsi,sex(fsi), 3-k,3, .TRUE., ALRq)
      LL(6) = addALR( addALR(LLg(6),ALRq), Maxval(ALRAU(1:2,k)))  
    endif
  endif
endif

LHH = missing
LHH2 = missing              
if (complx==2 .and. nYears>1 .and. ns(SB,k)>0 .and. (fclsib .or. &  
  (focal>5 .and. GpID(3-k,SB,k)==0)) .and. MaxLL(LL(2:3))<0D0 .and. MaxLL(LL(2:3))>=LL(7)) then
  call AddSibInbr(A, SB, k, LHH)  
  ! 1: par(Parent(A,3-k),k)=SB, 2: Parent(A,3-k)=GpID(3-k,SB,k), 3: as 1, A FS of B's (PA == DB)
  if (ns(SB,k)==1 .and. all(GpID(:,SB,k)==0)) then
    call pairHSHAI(SibID(1,SB,k),A,k, LHH2)  ! B1 inbred  (needed for symmetry)
    if (LHH2 < 0D0 .and. (LHH2 > LHH(2) .or. LHH(2)>0))  LHH(2) = LHH2
  endif
  if (MaxLL(LHH(1:2)) > LLg(3) .and. MaxLL(LHH(1:2))<0D0) then
    LLg(3) = MaxLL(LHH(1:2))
    LL(3) = addALR(LLg(3), ALR(3))
  endif
  if (LHH(3) > LLg(2) .and. LHH(3)<0D0) then  ! MAX(LLg(3), LLg(2))
    LLg(2) = LHH(3)
    LL(2) = addALR(LLg(2), ALR(2))
  else if (LL(3)==impossible .and. LLg(2)<0 .and. MaxLL(LHH)<0D0 .and. &
   MaxLL(LHH) > LLg(2)) then  ! add inbred FS
    LLg(2) = MaxLL(LHH)
    LL(2) = addALR(MaxLL(LHH), ALR(2))
  endif
endif 

LLy = missing
if (Complx==2 .and. fclsib .and. ns(SB,k)>0 .and. MaxLL(LL(2:3))<0D0 .and. MaxLL(LLg(2:3))>MaxLL(LLg(5:7))) then
  do x=1,2
    if (ALRAU(1,x)/= impossible) then
      call HSmating(-SB, k, A, k, x, LLy(1,x))
    endif
    if (ALRAU(2,x)/= impossible) then
      call HSmating(A, k, -SB, k, x, LLy(2,x))
    endif
  enddo
  if (any(LLy < 0D0)) then
    LLg(6) = MaxLL((/LLg(6), LLy(1,:), LLy(2,:)/))
    do x=1,2
      do m=1,2
        LLy(m,x) = addALR(LLy(m,x), ALRAU(m,x))
      enddo
    enddo
    LL(6) = MaxLL((/LL(6), LLy(1,:), LLy(2,:)/))
  endif
endif 

LLpo = missing
ALRpo = missing
LLgp = missing
ALRgp = missing
LLdGP = Missing   
!LLHSPO = Missing           
if (focal/=7 .and. complx>0 .and. &   ! configurations not monogamous
 (FSpar<0 .or. FSpar==A .or. COUNT(nFS(SibID(1:ns(SB,k),SB,k))>0) <= 5)) then  ! no. full-sib groups  
  ! one of Bi parent of A?  (NOT when called by CalcCandParLL() !)
  call CalcAgeLR(A,Sex(A), -SB,k, 3-k,4, .TRUE., ALRq)
  if (ANY(Parent(A,:)<=0) .and. ALRq /=impossible .and. focal/=6) then 
    ParTmp = Parent(A,:)
    do m=1,2
      if (Parent(A,m)>0) cycle
      do i=1, ns(SB,k)
        Bi = SibID(i,SB,k)
        if (AgeDiff(A, Bi) < 0) cycle
        if (Sex(Bi)/=m .and. Sex(Bi)<3)  cycle
        call ChkValidPar(A, Sex(A), Bi,m, ParOK)
        if (.not. ParOK)  cycle
        call setParTmp(A, Sex(A), Bi, m)
        call CalcU(A, Sex(A), -SB, k, LLPO(i,1))
        ! if (fclsib .and. Complx==2 .and. LLPO(i,1)>LLg(3)) then
          ! call addsib(A,SB,k,LLHSPO(i,1))  ! HS via k, & PO via 3-k
        ! endif
        call setParTmp(A, Sex(A), ParTmp(m), m)
        call CalcCLL(SB, k)
        call CalcLind(A)
        call CalcAgeLR(A,Sex(A), Bi,m,0,1,.TRUE., ALRpo(i,1))
      enddo
    enddo
  endif

  ! A parent 3-k of one of Bi?
  if (focal/=1 .and. (Sex(A)==3-k .or. Sex(A)>2)) then  
    do i=1, ns(SB,k)
      Bi = SibID(i,SB,k)
      if (AgeDiff(Bi, A) < 0 .or. Parent(Bi,3-k)>0) cycle
      call ChkValidPar(Bi,Sex(Bi), A, 3-k, ParOK)
      if (.not. ParOK)  cycle
      ParTmp = Parent(Bi,:)
      call setParTmp(Bi, Sex(Bi), A, 3-k)
      call CalcU(A, Sex(A), -SB, k, LLPO(i,2))
      ! if (fclsib .and. Complx==2 .and. LLPO(i,2)>LLg(3)) then
        ! call addsib(A,SB,k,LLHSPO(i,2))
      ! endif
      call setParTmp(Bi, Sex(Bi), ParTmp(3-k), 3-k) 
      call CalcCLL(SB, k)
      call CalcAgeLR( Bi,Sex(Bi), A,3-k,0,1,.TRUE., ALRpo(i,2))
    enddo
  endif

  if (any(LLpo < 0D0)) then
    LLg(6) = MaxLL((/LLg(6), LLpo(:,1), LLpo(:,2)/))
    do x=1,2
      do i=1, ns(SB,k)
        LLpo(i,x) = addALR(LLpo(i,x), ALRpo(i,x))
      enddo
    enddo
    LL(6) = MaxLL((/LL(6), LLpo(:,1), LLpo(:,2)/))
  endif
  
  ! if (any(LLHSPO < 0d0)) then    !! Not sure if this makes things better or worse
    ! LLg(3) = MaxLL((/LLg(3), LLHSPO(:,1), LLHSPO(:,2)/))
    ! do x=1,2
      ! do i=1, ns(SB,k)
        ! LLHSPO(i,x) = addALR( addALR(LLHSPO(i,x), ALRpo(i,x)), ALR(3))
      ! enddo
    ! enddo
    ! LL(3) = MaxLL((/LL(3), LLHSPO(:,1), LLHSPO(:,2)/))  
  ! endif

  ! one of Bi grandparent of A?
  if (ANY(Parent(A,:)<=0) .and. focal/=6) then  ! messes up calccandpar 
    do i=1, ns(SB,k)
      Bi = SibID(i,SB,k)
      if (AgeDiff(A, Bi) < 1) cycle
      do x=1,2
        if (Parent(A,x) > 0)  cycle
        call CalcAgeLR(A,Sex(A), Bi,Sex(Bi),x,4,.TRUE., ALRgp(i,x))
        if (ALRgp(i,x) == impossible)  cycle
        call PairGP(A, Bi, x, focal, LLgp(i,x))
        if (LLgp(i,x) < 0D0) then
          call CalcU(A,k, SibID(i,SB,k),k, LLUi)
          LLgp(i,x) = LLgp(i,x) - LLUi + LLg(7)
        endif
      enddo
    enddo
  endif
  
  ! A GP of one of Bi via 3-k?
  do i=1, ns(SB,k)
    Bi = SibID(i,SB,k)
    if (nFS(Bi)==0 .or. Parent(Bi,3-k)>0 .or. nFS(Bi)==ns(SB,k))  cycle  ! last case: safety net elsewhere
    if (AgeDiff(Bi, A) < 1) cycle
    call CalcAgeLR(Bi,Sex(Bi), A,Sex(A), 3-k,4,.TRUE., ALRgp(i,3))
    if (ALRgp(i,3) == impossible)  cycle
    call PairGP(Bi, A, 3-k, focal, LLgp(i,3))
    if (LLgp(i,3) < 0D0) then
      call CalcU(A,k, SibID(i,SB,k),k, LLUi)
      LLgp(i,3) = LLgp(i,3) - LLUi + LLg(7)     
      if (focal==4 .and. LLg(4)<0 .and. LLg(4)-LLgp(i,3) < TA .and. Sex(A)<3 .and. ANY(GpID(:,SB,k)==0)) then  ! check if double GP
        call setParTmp(-SB,k, A,Sex(A))                               
        call PairGP(Bi, A, 3-k, focal, LLdGP(i))
        call setParTmp(-SB,k, 0,Sex(A))
        if (LLdGP(i)<0)  LLdGP(i) = LLdGP(i) - LLUi + LLg(7)
        if (LLdGP(i)<0 .and. LLdGP(i)>LLg(4)) then
          LLg(4) = LLdGP(i)
          LL(4) = addALR(LLg(4), ALR(4))
          LLgp(i,3) = MaybeOtherParent   ! do not consider for LLg(6)
        endif
      endif
    endif
  enddo
  
  if (any(LLgp < 0D0)) then
    LLg(6) = MaxLL((/LLg(6), LLgp(:,1), LLgp(:,2), LLgp(:,3)/))
    do x=1,3
      do i=1, ns(SB,k)
        LLgp(i,x) = addALR(LLgp(i,x), ALRgp(i,x))
      enddo
    enddo
    LL(6) = MaxLL((/LL(6), LLgp(:,1), LLgp(:,2), LLgp(:,3)/))
  endif
endif


LLfs = Missing
if (.not. fclsib .and. Parent(A,k)==0 .and. Parent(A,3-k)<0 .and. FSpar<0 .and. &
  LL(2)==impossible .and. LL(3)/=impossible .and. ns(SB,k)>0) then
  ! check if FS anyway, by merging parents 3-k
  ! not merged because unclear if pat or mat merge needed.  
  ParTmp = Parent(A,:)
  call setParTmp(A,0,0,3-k)
  call AddFS(A, SB, k,0,k, LLfs(1,1), fsi, dx)
  if (fsi/=0)  call CalcAgeLR(A,Sex(A), fsi,k, 0,2, .TRUE., ALR(2))
  if (LLfs(1,1) < 0 .and. ALR(2)/=impossible) then
    call AddSib(A, SB, k, LLfs(2,1))    ! w/o parent(A,3-k)
    call CalcU(A, k, -SB, k, LLfs(3,1)) 
    if (LLfs(1,1) - MAXVAL(LLfs(2:3,1)) > TA) then
      call MergeSibs(-FSpar, -ParTmp(3-k), 3-k, LLfs(1,2))
      call parenthfs(0,-FSpar,3-k, -ParTmp(3-k),3-k, 3, LLfs(2,2))
      call CalcU(FSpar, 3-k, ParTmp(3-k), 3-k, LLfs(3,2))
     if (LLfs(1,2) - MAXVAL(LLfs(2:3,2)) > TA) then
       LLg(2) = LLfs(1,1) - LLfs(3,1) + LLg(7)
       LL(2) = addALR(LLg(2), ALR(2))
     endif
    endif
  endif
  call setParTmp(A,0,ParTmp(3-k),3-k)
endif

do x=1,4
  if (LL(x) > 0) then
    LLg(x) = LL(x)
  endif
enddo

end subroutine CheckAdd 

! #####################################################################

subroutine Qmerge(SA, SB, k, LR)
use Global
implicit none

integer, intent(IN) :: SA, SB, k
double precision, intent(OUT) :: LR
integer :: l, x, y
double precision :: PrL(nSnp), PrX(3), PrXY(3,3)

PrL = 0D0
do l=1,nSnp
  do x=1,3
    PrX(x) = XPr(1,x,l,SA,k)*XPr(1,x,l,SB,k)* AHWE(x,l)
    do y=1,3
      PrXY(x,y) = XPr(1,x,l,SA,k)*XPr(1,y,l,SB,k)* AHWE(x,l) * AHWE(y,l)
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrX)) - LOG10(SUM(PrXY))   ! merge
enddo
LR = SUM(PrL)

end subroutine Qmerge

! #####################################################################

subroutine QFSmerge(SA, SB, k, LR)
use Global
use CalcLik          
implicit none

integer, intent(IN) :: SA, SB, k
double precision, intent(OUT) :: LR
integer :: l, x, z, par(2), i, j
double precision :: PrL(nSnp), PrXZ(3,3,2), PrI(3,3), PrJ(3,3)

par = 0
call getFSpar(SA,k, .TRUE., par(1))
call getFSpar(SB,k, .TRUE., par(2))
if (any(par==0))  return

i = FSID(maxSibSize+1, SibID(1,SA,k))
j = FSID(maxSibSize+1, SibID(1,SB,k))

PrL = 0D0
do l=1,nSnp
  PrI = FSLik(l,i)
  PrJ = FSLik(l,j)                                    
  do x=1,3
    do z=1,3
      PrXZ(x,z,:) = prI(x,z) * AHWE(x,l) * AHWE(z,l)
      PrXZ(x,z,1) = PrXZ(x,z,1) * prJ(x,z)
      PrXZ(x,z,2) = PrXZ(x,z,2) * SUM(PrJ(x,:) * AHWE(:,l))
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXZ(:,:,1))) - LOG10(SUM(PrXZ(:,:,2)))
enddo
LR = SUM(PrL)

end subroutine QFSmerge

! #####################################################################

subroutine CheckMerge(SA, SB, kA, kB, focal, LLg, LL, FSM) 
use Global
implicit none

integer, intent(IN) :: SA, SB, kA, kB, focal
double precision, intent(OUT) :: LLg(7), LL(7)
logical, intent(OUT) :: FSM
double precision ::  ALRtmp, LRHS, ALR(7), LLtmp(3), dx(maxSibSize), LLx(6), ALRx(6), &
  LLz(2,2), LLY(2,3), ALRy(2,3), LLHA(3), dLH(nS(SB,kB)), LLM(5), LLMo(5), LLHHA(2), &
  LLC, LLP, TAx
integer :: i, j, x, Par(2), ix, tmpGP, NSx(2,2), OpPars(maxSibSize, 2), DoQuickA, DoQuickB
logical :: ShareOpp, ShareSib, AncOK(2), ParOK

LLg = missing
LL = missing
FSM = .FALSE.  ! merge both k & 3-k
ShareOpp = .FALSE.
ShareSib = .FALSE.
if (kA /= kB) then
  LL(1) = impossible
  if (focal==1)  return
endif
if (focal==1 .and. kA==1 .and. DoMtDif) then
  if (mtDif(SibID(1,SA,kA), SibID(1,SB,kB))) then
    LL(1) = impossible
    return
  endif
endif
do i=1, nS(SA, kA)
  do x=1,2
    if (SibID(i, SA, kA)==GpID(x,SB,kB)) then
      LL(1) = impossible
      exit
    endif
  enddo
enddo
do j=1, nS(SB, kB)
  do x=1,2
    if (SibID(j, SB, kB)==GpID(x,SA,kA)) then
      LL(1) = impossible
      exit
    endif
  enddo
enddo
do i=1, nS(SA, kA)
  do j=1, nS(SB, kB)
    if (AgeDiff(SibID(i,SA,kA), SibID(j,SB,kB))==missing) cycle
    if (getAP( AgeDiff(SibID(i,SA,kA), SibID(j,SB,kB)), 3, 0, kA, Impossible) == Impossible) then
      LL(1) = impossible
      exit
    endif
    if (LL(1)==impossible) exit
    if (kA/=kB) then
      if (SibID(i, SA, kA)==SibID(j,SB,kB)) then
        ShareSib = .TRUE.
      endif
    endif
  enddo
enddo 
if (LL(1) == impossible .and. focal==1) return

if (focal==1) then
  call ChkAncest(-SA,kA, -SB,kB, AncOK(1))
  call ChkAncest(-SB,kB, -SA,kA, AncOK(2))
  if (any(.not. AncOK)) then
    LL(1) = impossible
    return
  endif

  do x=1,2
    if (GpID(x,SB,kB)/=0) then  !  .and. GpID(x,SB,kB)/=GpID(x,SA,kA)
      call CalcAgeLR(-SA,kA, GpID(x,SB,kB),x, 0,1, .TRUE., ALRtmp)
      if (ALRtmp == impossible) then
        LL(1) = impossible
        return
      endif
    endif
    if (GpID(x,SA,kA)/=0) then
      call CalcAgeLR(-SB,kB, GpID(x,SA,kA),x, 0,1, .TRUE., ALRtmp)
      if (ALRtmp == impossible) then
        LL(1) = impossible
        return
      endif
    endif
  enddo

  if (hermaphrodites==1 .and. ((DumClone(SA,kA)/=0 .or. DumClone(SB,kB)/=0) .and. &
   .not. (DumClone(SA,kA)/=0 .and. DumClone(SB,kB)/=0)))  then
    LL(1) = MaybeOtherParent      
    return
  endif
endif

OpPars = 0
if (kA==kB .and. Complx/=0) then
  do i=1, nS(SA, kA)
    OpPars(i,1) = Parent(SibID(i,SA,kA), 3-kA)
  enddo
  do j=1, nS(SB,kB)
    OpPars(j,2) = Parent(SibID(j,SB,kB), 3-kB)
    if (OpPars(j,2)/=0 .and. ANY(opPars(1:ns(SA,kA),1) == opPars(j,2))) then
      ShareOpp = .TRUE.
    endif
  enddo
endif

LRHS = missing
if (.not. ShareOpp .and. .not. ShareSib.and. Complx/=0) then
  call ChkDoQuick(SA, kA, DoQuickA)
  call ChkDoQuick(SB, kB, DoQuickB) 
  if (DoQuickA /= 2 .and. DoQuickB /= 2) then
    call Qmerge(SA, SB, kB,  LRHS)
    if (LRHS < 2.0*TF*dble(MAX(nS(SA,kA), nS(SB,kB)))) then
      LL(1) = impossible
    endif
    if (LL(1) == impossible .and. focal==1) return
  endif
endif

if (focal==1) then
  do i=1, ns(SA,kA)
    do j=1, nS(SB,kB)
      if ((OpPars(i,1)==OpPars(j,2) .and. OpPars(i,1)/=0) .or. Complx==0) then 
        call PairQFS(SibID(i,SA,kA), SibID(j,SB,kB), LRHS)
      else
        call PairQHS(SibID(i,SA,kA), SibID(j,SB,kB), LRHS)
      endif
      if (LRHS < 4*TF) then
        LL(1) = impossible
        return
      endif
    enddo
  enddo
endif                  

call CalcU(-SA,kA, -SB,kB, LLg(7))
LL(7) = LLg(7)

ALR = missing
if (LL(1)/=impossible)  call CalcALRmerge(SA, SB, kA, ALR(1))              
if (LL(1)/=impossible .and. ALR(1)/=impossible) then
  if (Complx/=0 .or. focal==8) then
    call MergeSibs(SA, SB, kA, LLg(1))   ! SB parent of A's
  else
    Par = 0
    LLM = missing
    call getFSpar(SA, kA, .TRUE., Par(1))
    call getFSpar(SB, kB, .TRUE., Par(2))
    if ((Par(1)==Par(2) .and. Par(1)/=0) .or. all(Parent(SibID(:,SA,kA),3-kA)==0) .or. &
     all(Parent(SibID(:,SB,kB),3-kB)==0)) then
      call MergeSibs(SA, SB, kA, LLg(1))  
    else
      call FSMerge(SA,SB,kA, LLM)
      LLg(1) = LLM(4)   ! merge both k & 3-k
    endif
  endif    
  LL(1) = addALR(LLg(1), ALR(1))
else
  LL(1) = impossible
  LLg(1) = LL(1)                
endif

if (focal==1 .and. (LLg(1) > 0D0 .or. LL(1)==impossible .or. &
  (LL(1) - LL(7) < TA .and. Complx/=0))) return

call CalcAgeLR(-SB,kB, -SA,kA, 0,1, .TRUE., ALR(2))
if (ALR(2) /= impossible) then
  call addFS(0, SA, kA, SB, kB, LLg(2), ix, dx)  ! SB FS with an A
  if(complx>0)  call PairUA(-SB, -SA, kB, kA, LLg(3))  ! SB HS with an A
  do x=2,3
    LL(x) = addALR(LLg(x), ALR(2))
  enddo
else
  LL(2:3) = impossible
endif

if (focal==1 .and. (LL(1) - MaxLL(LL(2:7)) < TA) .and. ANY(LL(2:7) < 0) .and. Complx/=0) return

LLtmp = missing
tmpGP = 0        
call CalcAgeLR(-SA,kA, -SB,kB, 0,1, .TRUE., ALR(4))
if (ALR(4)/=impossible) then
  if (focal/=4 .or. complx==0) then
    call addFS(0, SB, kB, SA, kA, LLtmp(1), ix, dx)  ! SB GP of A's
  endif
  if (focal==4) then  ! allow for replacement
    tmpGP = GpID(kB,SA,kA)
    call setParTmp(-SA, kA, 0, kB)
  endif
  if(complx>0)  call PairUA(-SA, -SB, kA, kB, LLtmp(2))  ! SB GP of A's
  if (hermaphrodites/=0)  call addHAselfed(SA,kA,SB,kB, LLtmp(3))  
  if (focal==4)  call setParTmp(-SA, kA, tmpGP, kB)
  LLg(4) = MaxLL(LLtmp)
  LL(4) = addALR(LLg(4), ALR(4))
else
  LL(4) = impossible
endif

if (.not. (focal==4 .and. any(GpID(:,SA,kA)==0 .and. GpID(:,SB,kB)/=0))) then  
  ! else GB assigned as GP of SA as side-effect, messes up CalcCandPar  
  call CalcAgeLR(-SA,kA, -SB,kB, 0,2, .TRUE., ALR(5))
  if (ALR(5) /= impossible) then 
    if(complx>0)  call ParentHFS(0, SA, kA, SB, kB,3, LLg(5))  ! SB & SA are FS
    LL(5) = addALR(LLg(5), ALR(5))
  else
    LL(5) = impossible
  endif
endif

LLx = missing
ALRx = 0D0
do x=1,4
  if (complx==0) cycle
  if (x==1 .or. x==2) then  
    if (focal==4 .and. GpID(x,SB,kB)>0)  cycle  ! else messes up CalcCandPar 
    if (GpID(x,SA,kA)==GpID(x,SB,kB) .and. GpID(x,SB,kB)/=0)  cycle  ! SB & SA would become FS, not HS  
    call CalcAgeLR(-SA,kA, -SB,kB, x,3, .TRUE., ALRx(x))
    if (ALRx(x) /= impossible) then
      call ParentHFS(0, SA, kA, SB, kB, x, LLx(x))
    endif
  else if (x==3) then
    call CalcAgeLR(-SA,kA, -SB,kB, 3,4, .TRUE., ALRx(x))
    if (ALRx(x) /= impossible) then
     ! if (focal==4 .and. LLX(1)<0d0) then  ! else link via an offspring of SB is skipped. 
     !   call dummyGP(SA, SB, kA, kB, 3, LLx(3))  
     ! else
        call dummyGP(SA, SB, kA, kB, focal, LLx(3))  ! SB GGP of A's
     ! endif
    endif      
  else if (x==4) then
    call CalcAgeLR(-SB,kB, -SA,kA, 3,4, .TRUE., ALRx(x))
    if (ALRx(x) /= impossible) then
      call dummyGP(SB, SA, kB, kA, 0, LLx(4))  ! SA GGP of B's
    endif 
  endif
enddo

LLz = missing
do x=1,2
  if (GpID(x, SA, kA) > 0) then   ! TODO: more general
    if (Parent(GpID(x, SA, kA), kB)==-SB) then
      LLz(x,2) = impossible
    else
      do i=1,2
        if (focal==4 .and. Parent(GpID(x, SA, kA), i)/=0) then
          LLz(x,2) = NotImplemented   ! becomes not about SA/SB
        else if (Parent(GpID(x, SA, kA), i)/=0 .and. GpID(i,SB,kB)/=0 .and. &
         Parent(GpID(x, SA, kA), i)/=GpID(i,SB,kB)) then
          LLz(x,2) = Impossible
        else 
          call ChkValidPar(GpID(x, SA, kA), x, GpID(i,SB,kB), i, ParOK)
          if (.not. ParOK)  LLz(x,2) = impossible
        endif
      enddo
    endif
    if (LLz(x,2)==missing) then
      call CalcU(-SB, kB, GpID(x,SA,kA), x, LLz(x,1))
      call PairUA(-SB, GpID(x,SA,kA), kB, 3, LLz(x,2))
      call CalcAgeLR(-SB,kB, GpID(x,SA,kA),x,0,2, .TRUE., ALRx(4+x))
    endif
    if (LLz(x,2) < 0D0 .and. ALRx(4+x)/=impossible) then
      LLx(4+x) = LL(7) + LLz(x,2) - LLz(x,1)
    endif
  endif
enddo

LLg(6) = MaxLL(LLx)  ! most likely 3rd degree relative
do x=1,6
  LLX(x) = addALR(LLX(x), ALRx(x))
enddo
LL(6) = MaxLL(LLx)

! avuncular relationships between dummies?
LLY = missing
ALRy = missing
do x=1,3
  if (x < 3)  call CalcAgeLR(-SA,kA, -SB, kB, x, 6, .TRUE., ALRy(1,x))
  if (x == 3)  call CalcAgeLR(-SA,kA, -SB, kB, 3, 5, .TRUE., ALRy(1,x))   
  if (ALRy(1,x) /= impossible) then
    call dummyHFA(SA,kA, SB,kB, x, LLy(1,x))
  endif
enddo
do x=1,3
  if (x < 3)  call CalcAgeLR(-SB, kB, -SA,kA,  x, 6, .TRUE., ALRy(2,x))
  if (x == 3)  call CalcAgeLR(-SB, kB, -SA,kA, 3, 5, .TRUE., ALRy(2,x))   
  if (ALRy(2,x) /= impossible) then
    call dummyHFA(SB,kB, SA,kA, x, LLy(2,x))
  endif
enddo

if (any(LLY < 0)) then
  LLg(6) = MaxLL((/LLg(6), LLY(1,:), LLY(2,:)/))
  do x=1,3
    do j=1,2
      LLY(j,x) = addALR(LLY(j,x), ALRy(j,x))
    enddo
  enddo
  LL(6) = MaxLL((/LL(6), LLY(1,:), LLY(2,:)/))
endif

LLHA = missing
dLH = missing           
if (complx>0 .and. LL(4)<0D0 .and. focal/=4 .and. LLtmp(1)<0D0 .and. &
  LLtmp(1)>=LLtmp(2) .and. LLtmp(1) > MaxLL((/LL(1:3), LL(5:7)/))) then   
  do j=1, nS(SB, kB)
    if (GpID(3-kB,SA,kB)==0) then
      shareOpp = .TRUE.
      par(1) = Parent(SibID(j, SB, kB), 3-kB)
    else if (Parent(SibID(j, SB, kB), 3-kB)==GpID(3-kB,SA,kB) .or. &
      Parent(SibID(j, SB, kB), 3-kB)==0) then
      shareOpp = .TRUE.
      par(1) = GpID(3-kB,SA,kB)
    else
      shareOpp = .FALSE.
    endif
    if (shareOpp) then
      call PairUA(-SA, SibID(j, SB, kB), kA, 3, LLHA(1))
      call PairUA(-SA, SibID(j, SB, kB), kA, 3-kB, LLHA(2))
      call CalcU(-SA, kA, SibID(j, SB, kB), kB, LLHA(3))
      if (LLHA(1)<0)  dLH(j) = LLHA(1) - MaxLL(LLHA(2:3))
    endif
  enddo
  if (MAXVAL(dLH, MASK=dLH<missing) < TA) then
    LLg(4) = LLtmp(2)    ! do not use LLtmp(1) from addFS
    LL(4) = addALR(LLg(4), ALR(4))
  endif
endif

LLM = missing
LLMo = missing
Par = 0
NSx = 0
LLHHA = missing
if (kA == kB .and. (Complx==0 .or. (focal==1 .and. (ABS(MaxLL(LL) - LL(1))<TA)))) then 
  call FSMerge(SA,SB,kA, LLM) ! 1:not, 2: via k, 3: via 3-k, 4:both, !!5: 3-k + par HS
  if (Complx/=2)  LLM(5) = 555D0   ! merge via 3-k + par HS
  LLM(1) = MaxLL((/LLM(1), LL(7)/))  ! do not merge  
  LLM(2) = MaxLL((/LLM(2), LLg(1), LLM(5)/))  ! merge via k
  call getFSpar(SA, kA, .FALSE., Par(1))
  call getFSpar(SB, kB, .FALSE., Par(2))                                       
  if (par(1)<0 .and. par(2)<0) then
    NSx(1,1) = nS(SA,kA)
    NSx(2,1) = nS(SB,kB)                      
    NSx(1,2) = ns(-par(1), 3-kA)
    NSx(2,2) = ns(-par(2), 3-kB)
    if (ANY(NSx(:,1) /= NSx(:,2))) then
      call FSMerge(-par(1),-par(2),3-kA, LLMo)
       if (Complx/=2)  LLMo(5) = 555D0   ! merge via 3-k + par HS
      call CalcU(Par(1), 3-kA, Par(2), 3-kA, LLMo(1))  ! more accurate
      call MergeSibs(-par(1), -par(2), 3-kA, LLMo(2))
       if (LLMo(2) <0) then
        LLM(1) = MaxLL((/LLM(1), LLMo(2) - LLMo(1) + LL(7) /))
      endif
    endif
  endif

   TAx = TA * dble(MIN(nS(SA,kA), nS(SB,kB)))
  if (Complx>0 .and. LLM(2)<0D0 .and. LLM(1) - LLM(2) > TAx) then
    LL(1) = MaybeOtherParent 
  endif

  if (Complx==2 .and. par(1)<0 .and. par(2)<0 .and. hermaphrodites/=2) then
    if (SUM(NSx(:,1)) == SUM(NSx(:,2))) then  
      if (LLM(4)<0D0) then
        call FSHC(-SA,-SB,kA,LLC)
        if (LLC > LLM(4) .and. LLC<0)  LLM(4) = LLC
      endif
      call clustHSHA(SA, SB, kA, LLHHA(1))
      if (LLHHA(1) - LLM(4) > TA)  LL(1) = MaybeOtherParent
      call clustHSHA(SB, SA, kA, LLHHA(2))
      if (LLHHA(2) - LLM(4) > TA)  LL(1) = MaybeOtherParent
    endif
  endif
  
  if (hermaphrodites > 0) then   
    if ((ns(SA,kA)==1 .and. SelfedIndiv(SibID(1,SA,kA)) .and. all(GpID(:,SA,kA)==0)) .or. &
      (ns(SB,kB)==1 .and. SelfedIndiv(SibID(1,SB,kB)) .and. all(GpID(:,SB,kB)==0))) then
      LLM(4) = LLM(2)  ! selfed singletons
      LLM(2:3) = impossible
    endif
  else if (Complx==0) then
    if (Par(1)>0 .and. Par(1)==Par(2)) then   
      LLM(4) = LLM(2)
    endif
  endif
  
  TAx = TA * dble(MIN(nS(SA,kA), nS(SB,kB)))
  if (MaxLL(LLM(1:4))==LLM(4) .and. (LLM(4)-LLM(2) >TAx .and. &
   (ALL(LLMo==missing) .or. LLMo(4)-LLMo(3) >TAx)) .or. &
   ((Complx==0 .or. hermaphrodites>0) .and. LLM(4)<0 .and. LLM(4) - LLM(1) > TAx)) then
    LLg(1) = LLM(4)  ! FS merge most likely - go ahead.
    LL(1) = addALR(LLg(1), ALR(1))  
    FSM = .TRUE.
  else if (Complx>0 .and. LLM(3)<0D0 .and. LLM(3)-LLM(1) > TA ) then   ! .and. hermaphrodites/=2
    if (LLM(3)-LLM(4) > TA) then
      LL(1) = MaybeOtherParent ! likely that opp. parent need to be merged, 
    else if (Par(1) < 0 .and. Par(2)<0) then
       if (NSx(1,1)==NSx(1,2) .and. NSx(2,1)==NSx(2,2)) then   ! 2 FS groups
        LL(1) = MaybeOtherParent
      endif
    endif
  endif
endif

if (focal==1 .and. (ABS(MaxLL(LL) - LL(1))<TA)) then  ! one of Bi is SA, or vv?
  do i=1, ns(SA,kA)
    if (AgeDiff(SibID(i,SA,kA), SibID(1,SB,kB)) < 0) cycle
    LLP = missing
    call AddParent(SibID(i,SA,kA), SB, kB, LLP)
    if (LLP < 0D0 .and. (LLP + CLL(SA,kA) - Lind(SibID(i,SA,kA))) > LL(1)) then
      LL(1) = MaybeOtherParent
      exit
    endif
  enddo
  if (LL(1)<0D0) then
    do j=1, ns(SB,kB)
      if (AgeDiff(SibID(j,SB,kB), SibID(1,SA,kA)) < 0) cycle
      call AddParent(SibID(j,SB,kB), SA, kA, LLP)
      if (LLP < 0D0 .and. (LLP + CLL(SB,kB) - Lind(SibID(j,SB,kB))) > LL(1)) then
        LL(1) = MaybeOtherParent
      endif
    enddo
  endif 
endif

do x=1,4
  if (LL(x) > 0) then
    LLg(x) = LL(x)
  endif
enddo

end subroutine CheckMerge 

! #####################################################################

subroutine getFSpar(SA, kA, strict, par)  
! all individuals in SA are FS to eachother
use Global
implicit none

integer, intent(IN) :: SA,  kA
logical, intent(IN) :: strict
integer, intent(OUT) :: Par
integer :: i, j, ParV(ns(SA,kA))

Par = 0
if (ns(SA,kA)==0)  return

ParV = 0
do i=1, nS(SA,kA)
  if (Parent(SibID(i,SA,kA), 3-kA)/=0) then
    Par = Parent(SibID(i,SA,kA), 3-kA)
    if (strict) then
      do j= i+1, nS(SA, kA)
        if (Parent(SibID(j,SA,kA), 3-kA) /= Par .and. &
         Parent(SibID(j,SA,kA), 3-kA)/=0) then
          Par = 0
          return
        endif
      enddo
    else 
      ParV(i) = Par
    endif
  endif
enddo

if (.not. strict) then ! > half by same opp. parent?
  Par = 0
  do i=1, nS(SA,kA)
   if (real(COUNT(ParV == ParV(i))) > nS(SA,kA)/2.0) then
      Par = ParV(i)
      return
    else if (real(COUNT(ParV == ParV(i))) == nS(SA,kA)/2.0 .and. ParV(i)<0) then
      Par = ParV(i)
      return
    endif
  enddo
endif

end subroutine getFSpar

! #####################################################################

subroutine OppMerge(SA, k, LL)  ! could opposing parents of SA all be the same dummy parent?
use Global
implicit none

integer, intent(IN) :: SA, k
double precision, intent(OUT) :: LL ! of SA
integer :: i, l, x, y, m,u, opPar(ns(SA,k)), GPY(2)
double precision :: PrL(nSnp), PrSA(3), PrXY(3, 3), PrGY(3,2)

opPar = 0
GPY = 0
do i=1, ns(SA,k)
  opPar(i) = Parent(SibID(i, SA, k), 3-k)
enddo
if (ANY(opPar > 0)) then
  LL = NotImplemented
  return
else
  do i=1, ns(SA,k)
    if (opPar(i) < 0) then
      if (ns(-opPar(i), 3-k) > COUNT(opPar == opPar(i))) then
        LL = missing  ! TODO?  
        return
      else
        do m=1,2
          if (GpID(m,-OpPar(i),3-k) /= GPY(m)) then
            if (GPY(m) == 0) then
              GPY(m) = GpID(m,-OpPar(i),3-k)
            else
              LL = impossible
              return
            endif
          endif
        enddo
      endif
    endif
  enddo
endif

PrL = 0D0
do l=1, nSnp
  call ParProb(l, -SA, k, -1, 0, PrSA)
  do m=1,2
    call ParProb(l, GPY(m), 3-k, 0, 0, PrGY(:,m))  
  enddo
  do x=1,3
    do y=1,3
      do u=1,3
        PrXY(x,y) = PrSA(x) * SUM(AKA2P(y, u, :) * PrGY(u,1) * PrGY(:,2))
      enddo
      do i=1, ns(SA,k)
        PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l, SibID(i,SA,k)), x, y)
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXY))
enddo
LL = SUM(PrL)

end subroutine OppMerge

! #####################################################################

subroutine FSmerge(SA,SB,k, LL)  
! calc LL if SA and SB merged via both pat & mat
use Global
implicit none

integer, intent(IN) :: SA, SB, k
double precision, intent(OUT) :: LL(5) ! 1:not, 2: via k, 3: via 3-k, 4:both, !!5: 3-k + par HS
integer :: l, x, y, i, u,v, G(2,2),z, m, Par(2), SX(2)
double precision :: ALR, PrL(nSnp,5), PrXY(3,3), PrUV(3,3), PrXV(3,3,3,3,5),&
  PrG(3,2,2), PrX(3,2), PrTmp(3), PrY(3,2)
logical :: DoParHS, MaybeOpp(2), AncOK(2)

! TODO: currently assumes no gps of sibship 3-k, no close inbreeding
LL = missing
! check if all FS
SX = (/SA, SB/)
MaybeOpp = .FALSE.
do i=1,2
  call getFSpar(SX(i), k, .TRUE., Par(i))
  if (Par(i)>0) cycle
  if (Par(i)==0 .and. ANY(Parent(SibID(1:nS(SX(i),k),SX(i),k),3-k)>0)) cycle
  MaybeOpp(i) = .TRUE.
enddo
if (ALL(MaybeOpp)) then   
  if (Par(1)==Par(2) .and. Par(1)/=0) then 
    if (ALL(Parent(SibID(1:nS(SA,k),SA,k),3-k)==Par(1)) .and. &
      ALL(Parent(SibID(1:nS(SB,k),SB,k),3-k)==Par(2))) then
      MaybeOpp = .FALSE.   ! already share same opp parent
    endif
  else if (Par(1)<0 .and. Par(2)<0) then
    ALR = missing
    call CalcAgeLR(Par(1),3-k, Par(2),3-k, 0,-1, .TRUE., ALR)
    if (ALR==impossible) MaybeOpp = .FALSE.
    if (DoMtDif) then
      if (k==2 .and. mtDif(SibID(1,-Par(1),3-k), SibID(1,-Par(2),3-k))) then
        MaybeOpp = .FALSE.
      endif
    endif
    if (nS(-Par(1),3-k) > ns(SA,k) .or. nS(-Par(2),3-k) > ns(SB,k)) then
      LL = NotImplemented   ! called separately from CheckMerge on 'other side'
    endif
  endif
endif
if (ANY(.not. MaybeOpp) .or. ALL(LL==NotImplemented)) return

call ChkAncest(Par(1),3-k, Par(2),3-k, AncOK(1))
call ChkAncest(Par(2),3-k, Par(1),3-k, AncOK(2))
if (any(.not. AncOK)) then
  LL = Impossible
  return
endif  

G = 0
do i=1,2
  if (GpID(i,SA,k)/=0) then
    if(GpID(i,SA,k)/=GpID(i,SB,k) .and. GpID(i,SB,k)/=0) then
      G(i,k) = 0  ! shouldn't happen
    else
      G(i,k) = GpID(i,SA,k)
    endif
  else
    G(i,k) = GpID(i,SB,k)
  endif
  if (Par(1)<0) then
    G(i,3-k) = GpID(i, -Par(1),3-k)
    if (Par(2) < 0) then
      if (GpID(i, -Par(2),3-k) /= G(i,3-k) .and. &
       GpID(i, -Par(2),3-k)/=0 .and. G(i,3-k)/=0) then
        LL = Impossible
        return
      else if (G(i,3-k)==0 .and. GpID(i, -Par(2),3-k)/=0) then
        G(i,3-k) = GpID(i, -Par(2),3-k)
      endif
    endif
  endif
enddo

if (ALL(GPID(:,SA,k)==0) .and. ALL(GPID(:,SB,k)==0)  .and. Complx==2) then
  DoParHS = .TRUE.
else
  DoParHS = .FALSE.
endif

PrL = 0D0
do l=1,nSnp 
  do m=1,2
    do i=1,2
      call ParProb(l, G(i,m), i, 0, 0, PrG(:,i,m))
    enddo
    do x=1,3
      do z=1,3
        PrTmp(z) = SUM(AKA2P(x,:,z) * PrG(:,1,m) * PrG(z,2,m))
      enddo
      PrX(x,m) = SUM(PrTmp)
    enddo
  enddo
  do i=1,2
    call ParProb(l, Par(i), 3-k, -1, 0, PrY(:,i))
  enddo
  do x=1,3  ! P1
    do y=1,3  ! P2
      PrXY(x,y) = 1D0  ! XPr(2,x,l, sA,k) * AHWE(y,l)
      do i=1,nS(SA,k)
        PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l, SibID(i,SA,k)), x, y)
      enddo
    enddo
  enddo
  do u=1,3
    do v=1,3
      PrUV(u,v) = 1D0  ! XPr(2,u,l, sB,k) * AHWE(v,l)
      do i=1,nS(SB,k)
        PrUV(u,v) = PrUV(u,v) * OKA2P(Genos(l, SibID(i,SB,k)), u, v)
      enddo
    enddo
  enddo

  PrXV = 0D0
  do x=1,3 
    do y=1,3
      do u=1,3
        do v=1,3
          PrXV(x,y,u,v,1) = PrXY(x,y) * XPr(2,x,l, sA,k) * PrY(y,1) * &
            PrUV(u,v) * XPr(2,u,l, sB,k) * PrY(v,2)
          if (Complx/=0) then
            PrXV(x,y,x,v,2) = PrXY(x,y) * PrX(x,k) * PrY(y,1) * &
              PrUV(x,v) * PrY(v,2)
          endif
        enddo
        if (Complx/=0) then
          PrXV(x,y,u,y,3) = PrXY(x,y) * XPr(2,x,l, sA,k) * PrX(y,3-k) * &
              PrUV(u,y) * XPr(2,u,l, sB,k)
        endif
        if (DoParHS) then
          do z=1,3
            PrTmp(z) = AKAP(x,z,l) * AKAP(u,z,l) * AHWE(z,l)
          enddo
          PrXV(x,y,u,y,5) = PrXY(x,y) * PrX(y,3-k) * PrUV(u,y) * &
            SUM(PrTMP)
        endif
      enddo
      PrXV(x,y,x,y,4) = PrXY(x,y) * PrX(x,k) * PrX(y,3-k) * PrUV(x,y)
    enddo
  enddo
  do x=1,5
    PrL(l,x) = LOG10(SUM(PrXV(:,:,:,:,x)))
  enddo
enddo
LL = SUM(PrL,DIM=1)
if (.not. DoParHS)  LL(5) = impossible
if (Complx==0)  LL(2:3) = impossible                                    

end subroutine FSmerge

! #####################################################################

subroutine NewSibship(A, B, k)  ! make new sibship
use Global
implicit none

integer, intent(IN) :: A, B, k
integer :: s

nC(k) = nC(k) + 1
s = nC(k)
DumBY(:,s,k,:) = LOG10(1.0D0/nYears)                                    
call SetPar(A, Sex(A), -s, k)
if (B/=0) then
  call SetPar(B, Sex(B), -s, k) 
  if (BY(A) < 0)  call setEstBY(A, Sex(A))
  call UpdateLL(Parent(A,3-k), 3-k)
endif
call CalcCLL(s, k)

IsNewSibship(s,k) = .TRUE.
if (hermaphrodites/=0)  call CheckSelfed(-s,k)

if (Complx==0) then
  if (Parent(A,3-k)/=0)  DumMate(s,k) = Parent(A,3-k)
  if (Parent(A,3-k) < 0) then
    DumMate(-Parent(A,3-k), 3-k) = -s
  else if (Parent(A,3-k) > 0) then
    Mate(Parent(A,3-k)) = -s
  endif
endif

end subroutine NewSibship

! #####################################################################

subroutine CheckDropSibship(s, k, Drop) 
use Global
implicit none

integer, intent(IN) :: s, k
logical, intent(OUT) :: Drop
integer :: i                  

if (s > nC(k))  return  ! already dropped. 

Drop = .FALSE.
if (ns(s, k) == 0) then
  Drop = .TRUE.
else if (ALL(GpID(:,s,k)==0) .and. ns(s,k)==1) then
  if (DumClone(s,k)/=0 .or. Complx==0) then
    Drop = .FALSE.
  else
    Drop = .TRUE.
  endif
endif                 
if (.not. Drop)  return

if (ns(s,k)==1) then
  i = SibID(1, s, k)
  call RemoveSib(i, s, k)  
endif
call DoMerge(0, s, k)      ! delete sibship      

end subroutine CheckDropSibship

! #####################################################################

recursive subroutine DoMerge(SA, SB, k)  ! if SA=0, delete SB
use Global
implicit none

integer, intent(IN) :: SA, SB, k
integer :: j, n, m, x, y, BB(ns(SB,k)), nB
logical :: valid(2)            

if (SA == SB) return

valid = .TRUE.
if (SA/=0) then
  call ChkAncest(-SA,k, -SB,k, valid(1))
  call ChkAncest(-SB,k, -SA,k, valid(2))
  if (.not. (all(valid))) then
    call Erstop("Pedigree loop created by merge", .TRUE.)
  endif
endif 

if (SA/=0) then
  if (nS(SA,k) + nS(SB,k) >= maxSibSize) then
    call Erstop("Reached Maximum Sibship Size (number of offspring per parent), please increase '--maxsibsize'", .FALSE.)
  endif 
  nB = ns(SB,k)
  BB = SibID(1:ns(SB,k), SB, k)
  do j=1,nB
    if (nB==0)  exit
    call setPar(BB(j), 3, -SA, k)
  enddo
  ! dummy offspring fixed below

  do m=1,2
    if (GpID(m, SA, k)==0 .and. GpID(m, SB, k)/=0) then  ! checked for mismatches earlier
      call setPar(-SA,k, GpID(m,SB,k), m)   ! takes care of updating CLL, SClone, etc.
    endif  ! else keep GpID(i,SA,k)
  enddo
endif

do x=SB, nC(k)  !remove cluster SB, shift all subsequent ones
  SibID(:, x, k) = SibID(:, x+1, k)
  nS(x, k) = nS(x+1, k)
  GpID(:, x,k) = GpID(:, x+1,k)
  do n=1, nS(x,k)
    Parent(SibID(n,x,k),k) = -x
  enddo
  CLL(x,k) = CLL(x+1, k)
  XPr(:,:,:,x,k) = XPr(:,:,:,x+1,k)
  DumP(:,:,x,k) = DumP(:,:,x+1,k)
  DumBY(:,x,k,:) = DumBY(:,x+1,k,:)
  IsNewSibship(x,k) = IsNewSibship(x+1, k)
  DumMate(x,k) = DumMate(x+1, k)
  DumClone(x,k) = DumClone(x+1, k)
enddo
SibID(:,nC(k),k) = 0
GpID(:,nC(k),k) = 0
nS(nC(k), k) = 0
DumMate(nC(k), k) = 0  
DumClone(nC(k), k) = 0  

do x=SB, nC(k)
  if (Complx == 0) then
    if (any(Mate == -x .and. Sex==3-k)) then
      y = MINLOC(ABS(Mate + x), DIM=1, MASK = Sex==3-k)
      if (x==-SB) then
        Mate(y) = 0
      else
        Mate(y) = -x+1   ! shift towards zero.
      endif
    else if (any(DumMate(:,3-k) == -x)) then
      y = MINLOC(ABS(DumMate(:,3-k) + x), DIM=1)
      if (x==-SB) then
        DumMate(y, 3-k) = 0
      else
        DumMate(y, 3-k) = -x+1
      endif
    endif
  endif
  if (hermaphrodites/=0) then
    if (any(DumClone(:,3-k) == x)) then
      y = MINLOC(ABS(DumClone(:,3-k) - x), DIM=1)
      if (x==SB) then
        DumClone(y, 3-k) = 0  
      else
        DumClone(y, 3-k) = x-1
      endif
    endif
  endif
enddo                     

do m=1,2  !fix GPs
  do n=1, nC(m)
    if (GpID(k, n, m) == -SB) then
      GpID(k, n, m) = -SA
      if (all(GpID(:,n,m)==0) .and. ns(n,m)==1) then
        j = SibID(1,n,m)
        call RemoveSib(j,n,m)      
        call DoMerge(0, n, m)   ! Recursive
      endif
    endif
    do x=SB+1, nC(k)  
      if (GpID(k, n, m) == -x)  GpID(k, n, m) = -x+1 
    enddo
  enddo
enddo
nC(k) = nC(k) -1

if (SA/=0) then
  IsNewSibship(SA,k) = .TRUE.
  ToCheck(SibID(1:ns(SA,k),SA,k)) = .TRUE. 
  call setEstBY(-SA, k)                       
  if (hermaphrodites/=0) then
    call CheckSelfed(-SA,k) 
    if (all(SelfedIndiv(SibID(1:ns(SA,k),SA,k)))) then
      ToCheck(SibID(1:ns(SA,k),SA,k)) = .FALSE. 
    endif
  endif
endif

end subroutine DoMerge

! #####################################################################

subroutine DoFSMerge(SA, SB, k)   ! merge via k .and. k-3
use Global
implicit none

integer, intent(IN) :: SA, SB, k
integer :: i, ParA, ParB

! assume all checks have been done beforehand
! not implemented yet: Par(1) > 0, Par(2) <= 0 or vv
 call getFSpar(SA, k, .TRUE., ParA)
 call getFSpar(SB, k, .TRUE., ParB)
 
if (ParA==0 .and. any(parent(SibID(1:ns(SA,k),SA,k), 3-k)/=0)) then
  ParA = 9999
else if (ParA/=0 .and. any(parent(SibID(1:ns(SA,k),SA,k), 3-k)==0)) then
  do i=1, nS(SA,k)
    if (parent(SibID(i,SA,k), 3-k)==0) then
      call SetPar(SibID(i, SA, k), 3, parA, 3-k)
    endif
  enddo
endif
if (ParB==0 .and. any(parent(SibID(1:ns(SB,k),SB,k), 3-k)/=0)) then
  ParB = 9999
else if (ParB/=0 .and. any(parent(SibID(1:ns(SB,k),SB,k), 3-k)==0)) then
  do i=1, nS(SB,k)
    if (parent(SibID(i,SB,k), 3-k)==0) then
      call SetPar(SibID(i, SB, k), 3, parB, 3-k)
    endif
  enddo
endif
 
if (ParA < 0 .and. ParB < 0) then
    call DoMerge(-ParA, -ParB, 3-k)
else if (ParA==0 .and. ParB==0) then
    call NewSibship(SibID(1,SA,k), 0, 3-k)
    do i=2, ns(SA,k)
      call setPar(SibID(i,SA,k), 3, -nC(3-k), 3-k)
    enddo
    do i=1, nS(SB, k)
      call setPar(SibID(i, SB, k), 3, -nC(3-k), 3-k)
    enddo 
    call CalcCLL(nC(3-k), 3-k)
    call CalcCLL(SA, k)
    call CalcCLL(SB, k)
    call CalcCLL(nC(3-k), 3-k)
else if (ParA < 0 .and. ParB == 0) then
    do i=1, nS(SB, k)
    call setPar(SibID(i,SB,k), 3, ParA, 3-k)
  enddo
else if (ParB < 0 .and. ParA == 0) then
  do i=1, nS(SA, k)
    call setPar(SibID(i,SA,k), 3, ParB, 3-k)
  enddo
! else not implemented yet
endif

 call DoMerge(SA, SB, k)  ! takes care of MakeFS

end subroutine DoFSMerge

! #####################################################################

subroutine getOff(P, kP, dums, nOff, Off, sxOff)  ! list all offspring for parent P
use Global
implicit none

integer, intent(IN) :: P, kP
logical, intent(IN) :: dums  ! include dummy offspring
integer, intent(OUT) :: nOff, Off(maxSibSize), sxOff(maxSibSize)
integer :: i, k, m, s

nOff = 0
Off = 0
sxOff = 3
if (P==0) return                                
do k=1,2
  if (P>0 .and. kP/=1 .and. kP/=2) then
    if (Sex(P)<3 .and. Sex(P)/=k) cycle
  else if (k/=kP) then 
    cycle
  endif
  do i=1, nInd
    if (Parent(i,k) == P) then
      nOff = nOff + 1
      Off(nOff) = i
      sxOff(nOff) = Sex(i)
    endif
    if (nOff == maxSibSize) then
      call Erstop("Reached Maximum Sibship Size (number of offspring per parent), please increase '--maxsibsize'", .FALSE.)
    endif
  enddo
  if (dums) then
    do m=1,2
      do s=1,nC(m)
        if (GpID(k,s,m) == P) then
          nOff = nOff + 1
          Off(nOff) = -s
          sxOff(nOff) = m 
        endif
        if (nOff == maxSibSize) then
          call Erstop("Reached Maximum Sibship Size (number of offspring per parent), please increase '--maxsibsize'", .FALSE.)
        endif
      enddo
    enddo
  endif
enddo

end subroutine getOff

! #####################################################################

subroutine CalcU(A, kAIN, B, kBIN, LL)  ! A, SB, k, SA, LL
use Global
use CalcLik          
implicit none

integer, intent(IN) :: A, kAIN, B, kBIN
double precision, intent(OUT) :: LL
integer :: kA, kB, Ai, Bj, SA, SB, cat, m, n, par(2), i, tmpGP
logical :: con, swap, OpG, conP
double precision :: AgeD

LL = missing
con = .FALSE.
kA = 0
kB = 0
if (A>0) then
  call CalcLind(A)
else if (A<0) then
  kA = kAIN
  call CalcCLL(-A, kA)
endif
if (B>0) then
  call CalcLind(B)
else if (B<0) then
  kB = kBIN
  call CalcCLL(-B, kB)
endif
!==================================

if (A==0) then
  if (B==0) then
    LL = 0D0
  else if (B>0) then
    LL = Lind(B)
  else if (B<0) then
    LL = CLL(-B, kB)
  endif
  return
else if (B==0) then
  if (A>0) then
    LL = Lind(A)
  else if (A<0) then
    LL = CLL(-A,kA)
  endif
  return
else if (A>0 .and. B<0) then
  if (Parent(A,kB)==B) then
    LL = CLL(-B,kB)
    return
  else if (ANY(GpID(:,-B,kB) == A)) then
    LL = CLL(-B,kB) + Lind(A)  ! CLL already conditional on A
    return
  else if (ALL(Parent(A,:)>=0)) then  
    LL = Lind(A) + CLL(-B, kB)
    return
  else
    call Connected(A,1,B,kB, con)
    if (.not. con) then
      LL = Lind(A) + CLL(-B, kB)
      return
    endif
  endif
else if (B>0 .and. A<0) then
  if (Parent(B,kA)==A) then
    LL = CLL(-A, kA)
    return
  else if (ANY(GpID(:,-A,kA) == B)) then
    LL = CLL(-A,kA) + Lind(B)  
    return
  else if (ALL(Parent(B,:)>=0)) then 
    LL = CLL(-A, kA) + Lind(B) 
    return
  else
    call Connected(B,1,A,kA, con)
    if (.not. con) then
      LL = CLL(-A, kA) + Lind(B) 
      return
    endif
  endif
endif

!==================================
! determine relationship between focal individuals/clusters

Ai = 0
Bj = 0
SA = 0
SB = 0
cat = 0
swap = .FALSE.              

if (A>0 .and. B>0) then  ! == pairs ==
  do m=1,2
    if (Parent(A,m)/=0 .and. Parent(A,m)==Parent(B,m)) then
       par(m) = Parent(A,m)
    else
      par(m) = 0  ! unknown or unequal
    endif
  enddo

  if (par(1)/=0 .and. par(2)/=0) then
    cat = 2  ! FS
  else if (par(1)/=0 .or. par(2)/=0) then
    cat = 3  ! HS
  else 
    do m=1,2
      if (parent(A,m) < 0) then
        do n=1,2
          if (GpID(n, -parent(A,m), m) == B) then
            cat = 0 !4  ! already conditioned on.
          else
            if (GpID(n, -parent(A,m), m)==Parent(B, n) .and. &
             Parent(B, n)<0) then
              cat = 5
            endif
          endif
        enddo
      else if (parent(B,m) < 0) then
        if (ANY(GpID(:, -parent(B,m), m) == A)) then
          cat = 0 !4
          swap = .TRUE.
        else
          do n=1,2
            if (GpID(n, -parent(B,m), m) == Parent(A, n) .and. &
             Parent(A, n)<0) then
              cat = 5
              swap = .TRUE.
            endif
          enddo
        endif
      endif
    enddo
  endif

  if (cat==0 .or. cat==5) then  ! TODO? cat=5
    LL = Lind(A) + Lind(B)
    return
  else if (cat==2 .and. par(1)<0 .and. par(2)<0) then
    Ai = A
    Bj = B
    SA = -par(1)
    kA = 1
    SB = -par(2)
    kB = 2
    cat = 0
  else
    call Upair(A, B, cat, LL)
    return
  endif

else if (A>0 .and. B<0) then
  SB = -B
  Ai = A
  if (ALL(Parent(A,:) < 0)) then
    SA = -Parent(A,3-kB)
    kA = 3-kB
    do m=1,2
      do i=1,ns(-Parent(A,m),m)
        if (Parent(SibID(i,-Parent(A,m),m), 3-m)/=Parent(A,3-m)) then
          call Connected(SibID(i,-Parent(A,m),m),m,B,kB, conP)
          if (conP) then
            SA = -parent(A,m)
            kA = m
            exit
          endif
        endif
      enddo
    enddo
  else if (Parent(A,3-kB) < 0) then
    SA = -Parent(A,3-kB)
    kA = 3-kB    
  else if (Parent(A,kB) < 0) then
    SA = -Parent(A,kB)
    kA = kB
  endif ! else: Lind + CLL (earlier) 
else if (B>0 .and. A<0) then
  SA = -A
  Bj = B
  if (ALL(Parent(B,:) < 0)) then
    SB = -Parent(B,3-kA)
    kB = 3-kA
    do m=1,2
      do i=1,ns(-Parent(B,m),m)
        if (Parent(SibID(i,-Parent(B,m),m), 3-m)/=Parent(B,3-m)) then
          call Connected(SibID(i,-Parent(B,m),m),m,A,kA, conP)
          if (conP) then
            SB = -parent(B,m)
            kB = m
            exit
          endif
        endif
      enddo
    enddo
  else if (Parent(B,3-kA) < 0) then
    SB = -Parent(B,3-kA)
    kB = 3-kA
  else if (Parent(B,kA) < 0) then
    SB = -Parent(B,kA)
    kB = kA 
  endif
else if (A<0 .and. B<0) then
  SA = -A  
  SB = -B
endif


if (Hermaphrodites/=0) then   !!! NEW 2025-12-26 !!!
  if (kA/=kB .and. DumClone(SA,kA) == SB) then
    if (Bj==0 .and. Ai/=0 .and. ns(SA,kA)==(ns(SB,kB)+1)) then
      LL = CLL(SA,kA)
      return
    else if (Bj/=0 .and. Ai==0 .and. (ns(SA,kA)+1)==ns(SB,kB)) then
      LL = CLL(SB,kB)
      return
    endif
    ! ELSE: TODO: fix bug in UCLUST
  endif
endif

cat = 0
if (GpID(kB, SA, kA) == -SB) then
  cat = 1  ! PO
else if (GpID(kA, SB, kB) == -SA) then
  cat = 1
  swap = .TRUE.
else 
  do m=1,2
    if (GpID(m, SA, kA)==GpID(m, SB, kB) .and. GpID(m, SA, kA)/=0) then
     if (GpID(3-m,SA,kA)==GpID(3-m,SB,kB) .and. GpID(3-m,SA,kA)/=0) then  
        cat = 2  ! FS
      else
        cat = 3  ! HS
      endif
    else 
      if (GpID(m, SA, kA)<0) then
        if (GpID(kB, -GpID(m, SA, kA), m) == -SB) then
          cat = 4  ! GP
        endif
      endif
      if (GpID(m, SB, kB)<0) then
        if (GpID(kA, -GpID(m, SB, kB),m) == -SA) then
          cat = 4
          swap = .TRUE.
        endif
      endif
    endif
  enddo  ! FA between SA, SB not currently considered.
endif

OpG = .FALSE.
if (con .and. cat==0) then
  if (A<0 .and. B>0) then
    do i=1, ns(-A,kA)
      if (Parent(SibID(i,-A,kA), 3-kA) < 0) then
        if (ANY(GpID(:,-Parent(SibID(i,-A,kA), 3-kA), 3-kA)==B)) then
          SB = -Parent(SibID(i,-A,kA), 3-kA)
          kB = 3-kA
          Bj = 0
          OpG = .TRUE.
        endif
      endif
    enddo
else if (A>0 .and. B<0) then
  do i=1, ns(-B,kB)
    if (Parent(SibID(i,-B,kB), 3-kB) < 0) then
        if (ANY(GpID(:,-Parent(SibID(i,-B,kB), 3-kB), 3-kB)==A)) then
          SA = -Parent(SibID(i,-B,kB), 3-kB)
          kA = 3-kB
          Ai = 0
          OpG = .TRUE.
          swap = .TRUE.
        endif
      endif
    enddo
  endif
endif

if (cat==0) then ! swap if BY(A) < BY(B)   (A older)   
  call EstAgeDif(A, kA, B, kB, AgeD)  ! = AgeDiff(A,B) if BY(A) & BY(B) both exactly known
  if (AgeD < 0D0) then
    swap = .TRUE.
  endif
endif

if (con .and. A<0 .and. B<0 .and. kA==kB) then
  do i=1,ns(-B,kB)
    if (Parent(SibID(i,-B,kB), 3-kB) < 0) then
      if (GpID(3-kA, -Parent(SibID(i,-B,kB), 3-kB), 3-kB) < 0) then
        tmpGP = GpID(3-kA, -Parent(SibID(i,-B,kB), 3-kB), 3-kB)
        if (ANY(Parent(SibID(1:ns(-A,kA),-A,kA), 3-kA) == tmpGP) .and. &
         .not. ANY(Parent(SibID(1:ns(-B,kB),-B,kB), 3-kB) == tmpGP)) then
          swap = .TRUE.
        endif
      endif
    endif
  enddo
endif

if (.not. swap) then
  call UClust(-SA, -SB, kA, kB, cat, Ai, Bj, LL)
else
  call UClust(-SB, -SA, kB, kA, cat, Bj, Ai, LL)
endif

if (opG) then
  if (B>0) then
    LL = LL - CLL(SB,kB) + Lind(B)
  else if (A>0) then
    LL = LL - CLL(SA,kA) + Lind(A)
  endif
endif

end subroutine CalcU

! #####################################################################

subroutine Upair(A, B, cat, LL)
use Global
implicit none

integer, intent(IN) :: A, B, cat
double precision, intent(OUT) :: LL
integer :: m, l, n, x, y, par(2)
double precision :: PrL(nSnp), PrP(3,2), PrPA(3), PrPB(3), PrXY(3,3)

LL = missing
do m=1,2
  if (Parent(A,m)/=0 .and. Parent(A,m)==Parent(B,m)) then
     par(m) = Parent(A,m)
  else
    par(m) = 0  ! unknown or unequal
  endif
enddo
      
PrL = 0D0
do l=1, nSnp  
  if (cat==2) then
    do m=1,2
      call ParProb(l, Par(m), m, A, B, PrP(:,m))
    enddo
    do x=1,3
      do y=1,3
        PrXY(x,y) = OKA2P(Genos(l,A),x,y) * OKA2P(Genos(l,B),x,y) * &
          PrP(x,1) * PrP(y,2)
      enddo
    enddo
  else if (cat==3) then  ! HS
    do m=1,2
      if (Par(m)==0) cycle
      call ParProb(l, Par(m), m, A, B, PrP(:,m))
      call ParProb(l, Parent(A, 3-m), 3-m, A, 0, PrPA)
      call ParProb(l, Parent(B, 3-m), 3-m, B, 0, PrPB)
      do x=1,3  ! shared parent
        do y=1,3  ! parent A
          PrXY(x,y) = OKA2P(Genos(l,A),x,y) * PrP(x,m) * PrPA(y) * &
             SUM(OKA2P(Genos(l,B),x,:) * PrPB)
        enddo
      enddo
    enddo
  else if (cat==4) then
    do m=1,2
      if (Parent(A,m)<0) then
        do n=1,2
          if (GpID(n, -parent(A,m), m) == B) then
            call ParProb(l, parent(A,m), m, A, -4, PrP(:,m))  
            call ParProb(l, parent(A,3-m), 3-m, A, 0, PrPA)
            call ParProb(l, GpID(3-n, -parent(A,m), m), 3-n, 0, 0, PrPB)
            call ParProb(l, B, n, 0, 0, PrP(:,3-m))
            do x=1,3  ! in-between parent
              do y=1,3  ! other parent of A
                PrXY(x,y) = SUM(OKA2P(Genos(l,A),x,:) *PrPA) *PrP(x,m)*&
                   SUM(AKA2P(x, y,:) * PrP(y,3-m) * PrPB)
              enddo
            enddo
          endif
        enddo
      endif
    enddo     
  endif
  PrL(l) = LOG10(SUM(PrXY))
enddo
LL = SUM(PrL)

end subroutine Upair

! #####################################################################

subroutine UClust(A, B, kA, kB, cat, Ai, Bj, LL)
use Global
use CalcLik           
implicit none

integer, intent(IN) :: A, B, kA, kB, cat, Ai, Bj
double precision, intent(OUT) :: LL
integer :: nA, AA(ns(-A,kA)), GA(2), nB, BB(ns(-B,kB)), GB(2), &
  AB(2*maxSibSize), UseEE(ns(-A,kA)+ns(-B,kB)), TypeEE(ns(-A,kA)+ns(-B,kB)), &
  MateABpar(ns(-A,kA)+ns(-B,kB)), catA(maxSibSize), catB(maxSibSize), &
  l,x,y, v, i, j, z,m, f, e, u, DoneA(maxSibSize), g, Ei
double precision :: PrL(nSnp,2), PrGA(3,2), PrGB(3,2), PrGGP(3), PrFS(3,3), &
  PrUZ(3,3, 3,3,3,3,2), PrE(3), PrH(3), PrW(3), PrEE(3,(ns(-A,kA)+ns(-B,kB)))
logical :: ParAisClone(maxSibSize), ParBisClone(maxSibSize), AisBclone, &
  SIMPL, DoRSibs(maxSibSize, 2)

LL = missing

nA = nS(-A, kA)
AA = SibID(1:nA, -A, kA)
GA = GpID(:, -A, kA)

nB = nS(-B, kB)
BB = SibID(1:nB, -B, kB)
GB = GpID(:, -B, kB)

!============================================

AB = 0
UseEE = 0
TypeEE = 0
MateABpar = 0
if (kA==kB) then
  AB(1:nB) = BB
  AB((nB+1):(nB+nA)) = AA
  call FindEE(AB(1:(nB+nA)), nB, nA, kB, UseEE, MateABpar) 
  BB = AB(1:nB)
  AA = AB((nB+1):(nB+nA))
  TypeEE = 3-kB
  do i=1, nB  ! safety net
    if (UseEE(i) > nB) then
      UseEE(i) = 0  ! else use before store
    endif
  enddo
else if (kA/=kB) then
  call FindEE(BB, nB, 0, kB, UseEE(1:nB), MateABpar(1:nB))  ! may reorder BB
  call FindEE(AA, nA, 0, kA, UseEE((nB+1):(nB+nA)), MateABpar((nB+1):(nB+nA)))
  do i=1, nA
    if (UseEE(nB+i)/=0) then
      UseEE(nB+i) = nB + UseEE(nB+i)
    endif
  enddo
  TypeEE(1:nB) = 3-kB
  TypeEE((nB+1):(nB+nA)) = 3-kA
endif

!============================================
catA = 0
catB = 0
do i = 1, nA
  if (kA /= kB .and. Parent(AA(i), kB) == B) then
    catA(i) = 1
    UseEE(nB+i) = 0
  else if (GA(3-kA) == Parent(AA(i), 3-kA) .and. GA(3-kA) /= 0) then
    catA(i) = 2
  else if (GB(3-kA) == Parent(AA(i), 3-kA) .and. GB(3-kA) /= 0) then
    catA(i) = 3
  else if (Parent(AA(i), 3-kA) < 0) then
    if (GpID(kA, -Parent(AA(i), 3-kA), 3-kA)==A) then
      catA(i) = 6
      UseEE(nB+i) = 0
    else if (GpID(kB, -Parent(AA(i), 3-kA), 3-kA)==B) then
      catA(i) = 8
      UseEE(nB+i) = 0
    endif
  endif
enddo

do j=1, nB
  if (kA /= kB .and. Parent(BB(j), kA) == A) then
    catB(j) = 1
    UseEE(j) = 0
  else if (GA(3-kB) == Parent(BB(j), 3-kB) .and. GA(3-kB) /= 0) then
    catB(j) = 2
  else if (GB(3-kB) == Parent(BB(j), 3-kB) .and. GB(3-kB) /= 0) then
    catB(j) = 3  
  else if (Parent(BB(j), 3-kB) < 0) then
    if (GpID(kA, -Parent(BB(j),3-kB), 3-kB)==A) then
      catB(j) = 6
      UseEE(j) = 0
    else if (GpID(kB, -Parent(BB(j),3-kB), 3-kB)==B) then
      catB(j) = 8
      UseEE(j) = 0
    endif 
  endif
enddo

do i = 1, nA
  do j = 1, nB
    if (kA == kB) then
      if (Parent(AA(i), 3-kA) == Parent(BB(j), 3-kB) .and. Parent(BB(j), 3-kB)<0) then  
        catA(i) = 7
        if (catB(j)==0)  catB(j) = 7
      endif
    endif
  enddo
enddo

ParAisClone = .FALSE.
ParBisClone = .FALSE.
AisBclone = .FALSE.
if (hermaphrodites /= 0) then
  if (kA/=kB .and. DumClone(-A,kA) == -B) then
    AisBclone = .TRUE.
  else
    if (DumClone(-A,kA)/=0) then
      do i = 1, nA
        if (Parent(AA(i), 3-kA) == -DumClone(-A,kA))   ParAisClone(i) = .TRUE.
      enddo 
    endif
    if (DumClone(-B,kB)/=0) then
      do j=1,nB
        if (Parent(BB(j), 3-kB) == -DumClone(-B,kB))   ParBisClone(j) = .TRUE.
      enddo
    endif
  endif
endif


!==================================
if (cat==0 .and. ALL(catA==0) .and. ALL(CatB==0) .and. ALL(UseEE==0) .and. &
  .not. AisBclone) then
  LL = CLL(-A,kA) + CLL(-B,kB)
  return
endif
!==================================

SIMPL = ALL(catA==0) .and. ALL(catB==0) .and. Ai==0 .and. Bj==0 .and. &
    ALL(UseEE==0) .and. .not. AisBclone .and. .not. any(ParAisClone) .and. &
    .not. any(ParBisClone)

DoRsibs = .TRUE. 
if (.not. SIMPL) then
  call ChkTooManySibs(AA, nA, kA, DoRsibs(:,1))
  call ChkTooManySibs(BB, nB, kB, DoRsibs(:,2))
endif

PrL = 0D0
do l=1, nSnp
  PrUZ = 0D0
  
  do m=1, 2
    if ((ANY(catA==2) .and. m/=kA) .or. (ANY(catB==2) .and. m/=kB)) then
      call ParProb(l, GA(m), m, -1, 0, PrGA(:, m))
    else
      call ParProb(l, GA(m), m, 0, 0, PrGA(:, m))
    endif
    if ((ANY(catA==3) .and. m/=kA) .or. (ANY(catB==3) .and. m/=kB)) then
      call ParProb(l, GB(m), m, -1, 0, PrGB(:, m))
    else
      call ParProb(l, GB(m), m, 0, 0, PrGB(:, m))
    endif
  enddo
  if (cat==4) then
    do m=1,2
      if (GA(m)<0) then
        if (GpID(kB, -GA(m), m) == B) then
          call ParProb(l, GpID(3-kB, -GA(m), m), 3-kB, 0, 0, PrGGP) 
        endif
      endif
    enddo
  endif
  
  ! == grandparents ==
  do x=1,3 
    do y=1,3
      if (AisBclone .and. y/=x)  cycle                                
      do u=1,3  ! GP A, kB
        do z=1,3  ! GP A, 3-kB
          do v=1,3  ! GP B, kB
            if (cat == 1) then
             if (kA==kB .and. GA(3-kB)==GB(3-kB) .and. GA(3-kB)/=0) then 
                PrUZ(x,y,y,z,v,z,1) = AKA2P(x,y,z) * AKA2P(y,v,z) *&
                 PrGA(z,3-kB) * PrGB(v,kB)
              else
                PrUZ(x,y,y,z,v,:,1) = AKA2P(x,y,z) * AKA2P(y,v,:) * &
                  PrGA(z,3-kB) * PrGB(v,kB) * PrGB(:,3-kB)
              endif
            else if (cat==2) then
              PrUZ(x,y,u,z,u,z,1) = AKA2P(x,u,z) * AKA2P(y,u,z) *&
               PrGA(u,kB) * PrGA(z,3-kB)
            else if (cat==3) then
              do m=1,2
                if (GA(m)/=0 .and. GA(m) == GB(m)) then
                  if (m==kB) then
                    PrUZ(x,y,u,z,u,:,1) = AKA2P(x,u,z) * AKA2P(y,u,:) *&
                       PrGA(u,m) * PrGA(z,3-m) * PrGB(:,3-m)   
                  else
                    PrUZ(x,y,u,z,v,z,1) = AKA2P(x,u,z) * AKA2P(y,v,z) *&
                       PrGA(u,3-m) * PrGA(z,m) * PrGB(v,3-m)
                  endif
                endif
              enddo
            else if (cat==4) then 
              do m=1,2
                if (GA(m)<0) then
                  if (GpID(kB, -GA(m), m) == B) then
                    if (m==kB) then
                      PrUZ(x,y,u,z,v,:,1) =AKA2P(x,u,z) *&
                       SUM(AKA2P(u,y,:) *PrGGP) *PrGA(z,3-kB) *&
                       AKA2P(y,v,:) * PrGB(v,kB) * PrGB(:, 3-kB) 
                    else
                      PrUZ(x,y,u,z,v,:,1) = AKA2P(x,u,z) *&
                     SUM(AKA2P(z,y,:) *PrGGP) *PrGA(u,kB) *AKA2P(y,v,:)&
                     * PrGB(v,kB) * PrGB(:, 3-kB)
                    endif
                  endif
                endif
              enddo
            else
              PrUZ(x,y,u,z,v,:,1) = AKA2P(x,u,z) * AKA2P(y,v,:) * &
                PrGA(u,kB) * PrGA(z,3-kB) * PrGB(v,kB) * PrGB(:, 3-kB)
            endif
            PrUZ(x,y,u,z,v,:,2) = PrUZ(x,y,u,z,v,:,1)
          enddo
        enddo
      enddo
    enddo
  enddo
   
  ! == siblings ==   
  if (SIMPL) then
    do x=1,3 
      do y=1,3
        PrUZ(x,y,:,:,:,:,2) = PrUZ(x,y,:,:,:,:,2) * XPr(1,x,l,-A,kA) *&
         XPr(1,y,l,-B,kB)
      enddo  ! TODO: needs special for cat<4 ?
    enddo
  
  else

  do x=1,3  ! SA
    doneA = 0
    do y=1,3  ! SB
      if (AisBclone .and. y/=x)  cycle                                
      PrEE = 0D0
      do j=1, nB
        if (nFS(BB(j))==0) cycle
        if (catB(j)==1 .or. catB(j)==2 .or. catB(j)==3 .or. ParBisClone(j)) then
          PrE = 1D0
        else if (catB(j)==6) then
          call ParProb(l, GpID(3-kA,-Parent(BB(j),3-kB),3-kB),3-kA, 0, 0, PrH)
          do e=1,3
            PrE(e) = SUM(AKA2P(e,x,:) * PrH) 
          enddo
        else if (catB(j)==8) then  
          call ParProb(l, GpID(3-kB,-Parent(BB(j),3-kB),3-kB),3-kB, 0, 0, PrH)
          do e=1,3
            PrE(e) = SUM(AKA2P(e,y,:) * PrH) 
          enddo
        else if (UseEE(j)/=0) then
          call ParProb(l, MateABpar(j), 3-TypeEE(j), 0,0,PrH)
          do e=1,3
            do u=1, 3
              PrW(u) = SUM(AKA2P(e,u,:) * PrEE(u,UseEE(j)) * PrH)
            enddo
            PrE(e) = SUM(PrW)
          enddo
          PrE = PrE/SUM(PrE)
        else if (catB(j)==1 .or. catB(j)==2 .or. catB(j)==3) then
          PrE = 1D0
        else
          if (DoRsibs(j,2)) then
            call ParProb(l, Parent(BB(j),3-kB), 3-kB, -1, 0, PrE)
          else
            call ParProb(l, Parent(BB(j),3-kB), 3-kB, BB(j), -1, PrE)
          endif
        endif
        
        if (Parent(BB(j),3-kB) < 0 .and. catB(j)/=1 .and. DoRsibs(j,2)) then           
          do e=1,3
            do g=1, nS(-Parent(BB(j),3-kB), 3-kB)
              Ei = SibID(g, -Parent(BB(j),3-kB), 3-kB)
              if (nFS(Ei) == 0) cycle
              if (Parent(Ei, kB) == B) cycle  
              if (Parent(Ei, kA) == A) cycle
              PrFS = FSLik(l,Ei)
              call ParProb(l, Parent(Ei, kB), kB, Ei, -1, PrH) 
              PrH = PrH * PrFS(:,e)
              if (.not. all(PrH==1D0))  PrE(e) = PrE(e) * SUM(PrH)
            enddo
          enddo
        endif
        
        if (Bj/=0 .or. catB(j)==7 .or. (catB(j)==1 .and. Ai/=0)) then 
          do f=1, nFS(BB(j))
            if (Bj==0 .or. FSID(f, BB(j))==Bj) cycle
            if (Parent(BB(j),kA)==A .and. (Ai==0 .or. FSID(f, BB(j))==Ai)) cycle
            PrE = PrE * OKA2P(Genos(l,FSID(f,BB(j))), y, :)
          enddo
        endif
        
        if (catB(j)==7 .and. Ai/=0) then 
          do i=1,nA
            if (Parent(AA(i), 3-kB) /= Parent(BB(j), 3-kB)) cycle
            if (Parent(AA(i), kB) == B) cycle
            if (AA(i)==Ai)  cycle
            PrE = PrE * OKA2P(Genos(l,AA(i)), x, :)
          enddo
        endif
        
        if (catB(j)==1) then  ! Parent(BB(j), 3-kB)==PA
          PrUZ(x,y,:,:,:,:,1) = PrUZ(x,y,:,:,:,:,1) * PrE(x)
        else if (CatB(j)==2) then
          do z=1,3
            PrUZ(x,y,:,z,:,:,1) = PrUZ(x,y,:,z,:,:,1) * PrE(z)   
          enddo
        else if (CatB(j)==3) then
          do z=1,3
            PrUZ(x,y,:,:,:,z,1) = PrUZ(x,y,:,:,:,z,1) * PrE(z)   
          enddo                       
        else if (ParBisClone(j)) then     
          PrUZ(x,y,:,:,:,:,1) = PrUZ(x,y,:,:,:,:,1) * PrE(y)
        else if (.not. all(PrE==1D0)) then
          PrUZ(x,y,:,:,:,:,1) = PrUZ(x,y,:,:,:,:,1) * SUM(PrE)
        endif
        
        do f=1, nFS(BB(j)) ! includes some AA if cat=1 
          if (Bj==0 .or. FSID(f, BB(j))==Bj .or. &
            (Parent(BB(j),kA)==A .and. (Ai==0 .or. FSID(f, BB(j))==Ai))) then
            PrE = PrE * OKA2P(Genos(l,FSID(f, BB(j))), y, :)
          endif
        enddo

        if (catB(j)==7) then 
          do i=1,nA
            if (Parent(AA(i), 3-kB) /= Parent(BB(j), 3-kB)) cycle
            if (Parent(AA(i), kB) == B) cycle
            if (Ai/=0 .and. AA(i)/=Ai)  cycle
            PrE = PrE * OKA2P(Genos(l,AA(i)), x, :)
            DoneA(i) = 1
          enddo
        endif
        
        if (catB(j)==1) then  ! Parent(BB(j), 3-kB)==PA
          PrUZ(x,y,:,:,:,:,2) = PrUZ(x,y,:,:,:,:,2) * PrE(x)
        else if (CatB(j)==2) then
          do z=1,3
            PrUZ(x,y,:,z,:,:,2) = PrUZ(x,y,:,z,:,:,2) * PrE(z)   
          enddo
        else if (CatB(j)==3) then
          do z=1,3
            PrUZ(x,y,:,:,:,z,2) = PrUZ(x,y,:,:,:,z,2) * PrE(z)   
          enddo                       
        else if (ParBisClone(j)) then
          PrUZ(x,y,:,:,:,:,2) = PrUZ(x,y,:,:,:,:,2) * PrE(y)
        else if (.not. all(PrE==1D0)) then
          PrUZ(x,y,:,:,:,:,2) = PrUZ(x,y,:,:,:,:,2) * SUM(PrE)
        endif
        PrEE(:,j) = PrE
      enddo  ! B_j
    
      do i=1, nA
        if (DoneA(i)==1) cycle
        if (nFS(AA(i))==0) cycle
        if (Parent(AA(i),kB)==B) cycle
        if ((catA(i)>1 .and. catA(i)<4) .or. ParAisClone(i)) then  ! catA==1 already done
          PrE = 1D0
        else if (catA(i)==6) then ! Parent(AA(i), 3-kA) <0
          call ParProb(l, GpID(3-kA,-Parent(AA(i),3-kA),3-kA),3-kA,0,0,PrH)
          do e=1,3
            PrE(e) = SUM(AKA2P(e,x,:) * PrH) 
          enddo
        else if (catA(i)==8) then ! Parent(AA(i), 3-kA) <0
          call ParProb(l, GpID(3-kB,-Parent(AA(i),3-kA),3-kA),3-kB,0,0,PrH)
          do e=1,3
            PrE(e) = SUM(AKA2P(e,y,:) * PrH) 
          enddo
        else if (UseEE(nB+i)/=0) then
          call ParProb(l, MateABpar(nB+i), 3-TypeEE(nB+i), 0,0,PrH)
          do e=1,3
            do u=1, 3
              PrW(u) = SUM(AKA2P(e,u,:) * PrEE(u,UseEE(nB+i)) * PrH)
            enddo
            PrE(e) = SUM(PrW)
          enddo
          PrE = PrE/SUM(PrE)
        else
          if (DoRsibs(i,1)) then
            call ParProb(l, Parent(AA(i), 3-kA), 3-kA, -1, 0, PrE)  
          else
            call ParProb(l, Parent(AA(i), 3-kA), 3-kA, AA(i), -1, PrE)  
          endif
        endif
        
        if (Parent(AA(i), 3-kA) < 0 .and. DoRsibs(i,1)) then  
          do e=1,3
            do g=1, nS(-Parent(AA(i), 3-kA), 3-kA)
              Ei = SibID(g, -Parent(AA(i), 3-kA), 3-kA)
              if (nFS(Ei) == 0) cycle
              if (Parent(Ei, kA) == A) cycle 
              if (Parent(Ei, kB) == B) cycle
              PrFS = FSLik(l,Ei)
              call ParProb(l, Parent(Ei, kA), kA, Ei, -1, PrH)            
              PrH = PrH * PrFS(:,e)
              if (.not. all(PrH==1D0))  PrE(e) = PrE(e) * SUM(PrH)
            enddo
          enddo
        endif
        
        if (Ai/=0) then
          do f=1, nFS(AA(i))
            if (FSID(f, AA(i))==Ai) cycle
            PrE = PrE * OKA2P(Genos(l,FSID(f,AA(i))), x, :)
          enddo
        endif
        
        if (catA(i)==2) then
          do z=1,3
            if (kA==kB) then
              PrUZ(x,y,:,z,:,:,1) = PrUZ(x,y,:,z,:,:,1) * PrE(z)
            else
              PrUZ(x,y,z,:,:,:,1) = PrUZ(x,y,z,:,:,:,1) * PrE(z)         
            endif
          enddo
        else if (catA(i)==3) then  
          do z=1,3
            if (kA==kB) then
              PrUZ(x,y,:,:,:,z,1) = PrUZ(x,y,:,:,:,z,1) * PrE(z)
            else
              PrUZ(x,y,:,:,z,:,1) = PrUZ(x,y,:,:,z,:,1) * PrE(z)
            endif
          enddo
        else if (ParAisClone(i)) then    ! TODO: .or. AisBclone?
          PrUZ(x,y,:,:,:,:,1) = PrUZ(x,y,:,:,:,:,1) * PrE(x)
        else if (.not. all(PrE==1D0)) then
          PrUZ(x,y,:,:,:,:,1) = PrUZ(x,y,:,:,:,:,1) * SUM(PrE)    
        endif
        
        do f=1, nFS(AA(i)) 
          if (Ai/=0 .and. FSID(f, AA(i))/=Ai) cycle
!          DoneA(i)=2    ! for debuging only
          PrE = PrE * OKA2P(Genos(l,FSID(f,AA(i))), x, :)
        enddo
        
        if (catA(i)==2) then
          do z=1,3
            if (kA==kB) then
              PrUZ(x,y,:,z,:,:,2) = PrUZ(x,y,:,z,:,:,2) * PrE(z)
            else
              PrUZ(x,y,z,:,:,:,2) = PrUZ(x,y,z,:,:,:,2) * PrE(z)         
            endif
          enddo
        else if (catA(i)==3) then  
          do z=1,3
            if (kA==kB) then
              PrUZ(x,y,:,:,:,z,2) = PrUZ(x,y,:,:,:,z,2) * PrE(z)
            else
              PrUZ(x,y,:,:,z,:,2) = PrUZ(x,y,:,:,z,:,2) * PrE(z)
            endif
          enddo
        else if (ParAisClone(i)) then
          PrUZ(x,y,:,:,:,:,2) = PrUZ(x,y,:,:,:,:,2) * PrE(x)
        else if (.not. all(PrE==1D0)) then
          PrUZ(x,y,:,:,:,:,2) = PrUZ(x,y,:,:,:,:,2) * SUM(PrE)    
        endif
        PrEE(:,nB+i) = PrE
      enddo  ! i
    enddo  ! x
  enddo  ! y
  endif
  do f=1,2
    PrL(l,f) = LOG10(SUM(PrUZ(:,:,:,:,:,:,f)))
  enddo
enddo
LL = SUM(PrL(:,2)) - SUM(PrL(:,1))

end subroutine UClust

! #####################################################################

subroutine FindEE(AB, nA, nB, k, UseEE, MatePar)  ! find PO pairs among mates
use Global
use sort_module
implicit none

integer, intent(IN) :: nA, nB, k
integer, intent(INOUT) :: AB(nA+nB)
integer, intent(OUT) :: UseEE(nA+nB), MatePar(nA+nB)
integer :: i, j,x, nAB(2), ABM(MAX(nA,nB),2), MateE(MAX(nA,nB), 2), GGK(2), &
  UseM(MAX(nA,nB),2), Order(2*maxSibSize), MateI
logical :: reorder, OrderAgain
double precision :: EEtmp(2*maxSibSize)

UseEE = 0
MatePar = 0

nAB = (/nA, nB/)
ABM = 0
ABM(1:nA, 1) = AB(1:nA)
ABM(1:nB, 2) = AB((nA+1):(nA+nB))
MateE = 0
do x=1,2
  do i=1, nAB(x)
    if (nFS(ABM(i,x))==0 .and. nAB(x)>1)  cycle
    MateE(i,x) = Parent(ABM(i,x), 3-k)
  enddo
enddo
if ((nAB(1)==1 .and. nAB(2)<2) .or. COUNT(MateE < 0) < 2) return

GGK = 0
do x=1,2
  if (ABM(1,x)==0)  cycle
  if (Parent(ABM(1,x),k) < 0) then  ! else not called?
    GGK(x) = GpID(3-k, -Parent(ABM(1,x),k), k)   
  endif
enddo

! re-order AA and BB, so that PrE calculated before used
UseM = 0
reorder = .FALSE.
do x=1,2
  if (nAB(x)<=1) cycle
  do i=1, nAB(x)
    if (MateE(i,x) < 0) then 
      if (GpID(3-k, -MateE(i,x), 3-k) < 0 .and. &
       .not. ANY(GGK == GpID(3-k, -MateE(i,x), 3-k))) then
        do j=1, nAB(x)
          if (MateE(j,x) == GpID(3-k, -MateE(i,x), 3-k)) then
            UseM(i,x) = j
            if (j > i) reorder = .TRUE.
            exit
          endif
        enddo
      endif
    endif
  enddo
  
  if (reorder) then
    EEtmp(1:nAB(x)) = dble(UseM(1:nAB(x),x))
    Order = (/ (i, i=1, nAB(x), 1) /)
    call QsortC(EEtmp(1:nAB(x)), Order(1:nAB(x)))
    ABM(1:nAB(x),x) = ABM(Order(1:nAB(x)), x)
    UseM(1:nAB(x),x) = UseM(Order(1:nAB(x)),x)
    OrderAgain = .FALSE.
    do i=1, nAB(x)
      if (UseM(i,x) /= 0) then
        do j=1, nAB(x)
          if (UseM(i,x) == Order(j)) then
            UseM(i,x) = j
            if (j>i)  OrderAgain = .TRUE.
            exit
          endif
        enddo
      endif
    enddo
    if (OrderAgain) then
      EEtmp(1:nAB(x)) = dble(UseM(1:nAB(x),x))
      Order = (/ (i, i=1, nAB(x), 1) /)
      call QsortC(EEtmp(1:nAB(x)), Order(1:nAB(x)))
      ABM(1:nAB(x),x) = ABM(Order(1:nAB(x)), x)
    endif
  endif
enddo

AB = 0
AB(1:nA) = ABM(1:nA, 1)
AB((nA+1):(nA+nB)) = ABM(1:nB, 2)

do i=1, nA+nB
  if (nFS(AB(i))==0 .and. ((i<=nA .and. nA>1) .or. (i>nA .and. nB>1)))  cycle
  MateI = Parent(AB(i), 3-k)
  if (MateI < 0) then 
    if (GpID(3-k, -MateI, 3-k) < 0 .and. &
     .not. ANY(GGK == GpID(3-k, -MateI, 3-k))) then
      do j=1, i
      if (nFS(AB(j))==0 .and. ((j<=nA .and. nA>1) .or. (j>nA .and. nB>1)))  cycle
        if (Parent(AB(j), 3-k) == GpID(3-k, -MateI, 3-k)) then
          UseEE(i) = j
          MatePar(i) = GpID(k, -MateI, 3-k)
          exit
        endif
      enddo
    endif
  endif
enddo

end subroutine FindEE

! #####################################################################

subroutine AddSib(A, SB, k, LL)  
use Global
use CalcLik           
implicit none

integer, intent(IN) :: A, SB, k
double precision, intent(OUT) :: LL
integer :: l, x, Bj, f, j, y, Ei, z, v, DoQuick
double precision :: PrL(nSnp), PrX(3), PrY(3), LRQ, LLtmp(2), LLU, PrZ(3), &
  PrYZ(3,3), PrXb(3,2), PrE(3), PrB(3,3), PrFS(3,3)
logical :: Inbr, ParOK 

LL = missing
if (Parent(A,k)==-SB) then
  LL = AlreadyAss
else if (Parent(A,k)/=0) then
  LL = impossible
endif
if (LL/=missing) return

if (ns(SB,k)==0 .and. any(GpID(:,SB,k)==0) .and. any(GpID(:,SB,k)>0)) then
  do x=1,2
    if (GpID(x,SB,k)>0) then
      call pairGP(A, GpID(x,SB,k), k, 4, LL)  ! requires 'focal'
      if (LL < 0d0)  LL = LL - Lind(GpID(x,SB,k))
      return
    endif
  enddo
endif

do f=1, nS(SB,k)
  if (ns(SB,k)==0)  exit                      
  Bj = SibID(f, SB, k)
  if (Parent(A, 3-k) /= 0) then
    if (Parent(Bj, 3-k) == Parent(A, 3-k)) then
      LL = impossible  ! use addFS() instead
    endif
  endif
  if (getAP (AgeDiff(A, Bj), 3, 0, k, Impossible) == Impossible) then  
    LL=impossible
  endif 
enddo
if (LL/=missing) return

 call ChkValidPar(A,Sex(A), -SB,k, ParOK)
if (.not. ParOK) then
  LL = impossible
  return
endif

call Qadd(A, SB, k, LRQ)
if (LRQ < -HUGE(0D0)) then
  LL = impossible
  return
endif

Inbr = .FALSE.
if (Parent(A,3-k) < 0) then
  if (Parent(A,3-k) == GpID(3-k, SB, k)) then
    Inbr = .TRUE.  ! inbreeding loop created
  else if (GpID(k,-Parent(A,3-k),3-k) == -SB) then
    Inbr = .TRUE.
  endif
endif
do f=1, nS(SB,k)
  if (ns(SB,k)==0)  exit                      
  Bj = SibID(f, SB, k)
  if (Parent(A,3-k) == Bj)  Inbr = .TRUE.
  if (Parent(Bj,3-k) == A)  Inbr = .TRUE.
enddo

call ChkDoQuick(SB,k,DoQuick)

if ((Parent(A,3-k)<0 .and. DoQuick/=-2) .or. Inbr .or. DoQuick>1 .or. DoQuick==-3) then 

  if (Parent(A,3-k) < 0) then
    call CalcU(-SB, k, A, 3-k, LLU)
    call CalcU(-SB,k, Parent(A,3-k),3-k, LLtmp(1))
  else if (GpID(3-k,SB,k) < 0) then
    call CalcU(-SB, k, 0, 0, LLU)
    call CalcU(-SB,k, GpID(3-k,SB,k),3-k, LLtmp(1))
  endif
  call setParTmp(A,0,-SB,k)
  if (Parent(A,3-k) < 0) then
    call CalcU(-SB,k, Parent(A,3-k),3-k, LLtmp(2)) 
    LL = LLU + (LLtmp(2) - LLtmp(1))
  else if (GpID(3-k,SB,k) < 0) then
    call CalcU(-SB,k, GpID(3-k,SB,k),3-k, LLtmp(2))
    LL = LLU + (LLtmp(2) - LLtmp(1))
  else
    LL = CLL(SB,k)
  endif
  call setParTmp(A,0,0,k)
  if (GpID(3-k,SB,k) < 0)  call CalcCLL(-GpID(3-k,SB,k), 3-k)
  
else 
  PrL = 0D0
  
  if (DoQuick == -2) then    ! all FS
    do l=1,nSnp
      call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrY)
      do j=1, nS(SB,k)
        Bj = SibID(j,SB,k)
        if (nFS(Bj)==0) cycle
        call ParProb(l, Parent(Bj,3-k), 3-k, -1, 0, PrZ)
        PrB = FSLik(l,Bj)                 
        do x=1,3
          do y=1,3
            do z=1,3
              PrYZ(y,z) = PrB(x,z) * PrZ(z) * OKA2P(Genos(l,A), x, y) * PrY(y) 
            enddo
          enddo
          PrX(x) = XPr(2,x,l, SB,k) * SUM(PrYZ)
        enddo
      enddo
      PrL(l) = LOG10(SUM(PrX))   
    enddo
  
  else if (abs(DoQuick) == 1) then
    do l=1,nSnp
      call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrY)
      do x=1,3
        PrX(x) = XPr(3,x,l, SB,k) * SUM(OKA2P(Genos(l,A), x, :) * PrY)   
      enddo
      PrL(l) = LOG10(SUM(PrX))   
    enddo

  else 
    do l=1,nSnp 
      do x=1,3
        PrXb(x,:) = XPr(2,x,l,SB,k)  ! GPs
        do j=1, nS(SB,k)
          Bj = SibID(j,SB,k)
          if (nFS(Bj)==0) cycle
          call ParProb(l, Parent(Bj,3-k), 3-k, -1, 0, PrZ)
          do z=1,3
            if (Parent(Bj,3-k)<0) then
              do v=1, nS(-Parent(Bj, 3-k), 3-k)
                Ei = SibID(v, -Parent(Bj, 3-k), 3-k)  
                if (NFS(Ei) == 0) cycle
                if (Parent(Ei, k) == -SB) cycle
                PrFS = FSLik(l,Ei)
                call ParProb(l, Parent(Ei, k), k, Ei,-1, PrE)                
                PrE = PrE * PrFS(:,z)
                if (.not. ALL(PrE==1D0))  PrZ(z) = PrZ(z) * SUM(PrE)  
              enddo  
            endif
          enddo
          if (.not. ALL(PrZ==1D0))  PrXb(x,1) = PrXb(x,1) * SUM(PrZ)
          PrB = FSLik(l,Bj)
          PrZ = PrZ * PrB(:,x)
          PrXb(x,2) = PrXb(x,2) * SUM(PrZ)
        enddo
      enddo

      call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrY)
      do x=1,3
        PrXb(x,2) = PrXb(x,2) * SUM(OKA2P(Genos(l,A), x, :) * PrY)
      enddo
      PrL(l) = LOG10(SUM(PrXb(:,2))) - LOG10(SUM(PrXb(:,1)))  
    enddo
  endif
  
  LL = SUM(PrL)
endif

end subroutine AddSib

! #####################################################################

subroutine AddSibInbr(A,SB,k,LL)
use Global
use CalcLik             
implicit none

integer, intent(IN) :: A, SB, k
double precision, intent(OUT) :: LL(3)
integer :: l, x, y, GA(2), GG, Par, i, u, Bj, j, GX, Ei
double precision :: PrL(nSnp,3), PrXY(3,3), PrZ(3), PrPA(3), LLtmp(3), &
  ALR(3), LLU(4), PrXYU(3,3,3), PrLU(nSnp,3), PrE(3), PrFS(3,3), PrB(3)
logical :: maybe(3)

! 1: Par(Parent(A,3-k),k)=SB 
! 2: Parent(A,3-k)=GpID(3-k,SB,k)
! 3: as 1, A FS of B's (PA == DB)

LL = missing
maybe = .TRUE.
GA = getPar(Parent(A,3-k), 3-k)
if (GA(k) == -SB) then
  maybe(2) = .FALSE.
else if (GA(k) /= -SB .and. GA(k)/=0) then
  maybe(1) = .FALSE.
endif

if (hermaphrodites/=0) then
  LL = NotImplemented
  return
endif

GG = GpID(k,SB,k)
if (GpID(3-k,SB,k)/=Parent(A,3-k) .and. GpID(3-k,SB,k)/=0 .and. Parent(A,3-k)/=0) then
  maybe(2) = .FALSE.
else if (ANY(SibID(1:ns(SB,k),SB,k)==Parent(A,3-k))) then
  maybe(2) = .FALSE.
endif
if (.not. ANY(maybe(1:2))) then
  LL = impossible
  return
endif

Par = 0       
if (maybe(1)) then
  call getFSpar(SB, k, .TRUE., Par)
  if (Par==0) then
    maybe(3) = .FALSE.
  else if (Parent(A,3-k)/=Par .and. Parent(A,3-k)/=0) then
    maybe(3) = .FALSE.
  else if (Par>0) then
    if (Parent(Par, k)/=SB .and. Parent(Par, k)/=0) then
      maybe(3) = .FALSE.
    endif
  else if (Par<0) then
    if (GpID(k, -Par, 3-k)/=SB .and. GpID(k, -Par, 3-k)/=0) then
      maybe(3) = .FALSE.
    endif
  endif
  if (maybe(3) .and. GA(3-k)==0) then
    GA = getPar(Par, 3-k)
  endif
else
  maybe(3) = .FALSE.
endif

call CalcAgeLR(Parent(A,3-k),3-k, -SB,k, 0,1, .TRUE., ALR(1))
call CalcAgeLR(-SB,k, Parent(A,3-k),3-k,  0,1, .TRUE., ALR(2))
call CalcAgeLR(Par,3-k, -SB,k, 0,1, .TRUE., ALR(3))
do x=1,3
  if (ALR(x) == impossible .or. ALR(x) < 3.0*TF)  maybe(x) = .FALSE.
enddo

if (.not. ANY(maybe)) then
  LL = impossible
  return
endif

GX = 0
if (maybe(2)) then
  if (Parent(A,3-k)/=0) then
    GX = Parent(A,3-k)
  else
    GX = GpID(3-k, SB, k)
  endif
endif

PrL = 0D0
PrLU = 0D0          
do l=1,nSnp
  if (maybe(1)) then
    if (Parent(A,3-k)>0) then
      call ParProb(l, GA(3-k), 3-k, Parent(A,3-k), 0, PrZ) 
      call ParProb(l, Parent(A,3-k),3-k,0,0,PrPA)
    else
      call ParProb(l, GA(3-k), 3-k, 0, 0, PrZ) 
      call ParProb(l, Parent(A,3-k),3-k,A,-4,PrPA)
    endif
    if (Parent(A,3-k)==0)   PrPA = 1D0
    PrB = XPr(3,:,l, SB,k)
    if (GA(k)==-SB .and. Parent(A,3-k)>0) then  ! avoid double counting
      do x=1,3
        PrE = OKA2P(Genos(l,Parent(A,3-k)),x,:) * PrZ
        if (SUM(PrE) > 0D0)  PrB(x) = PrB(x) / SUM(PrE)
      enddo  
    endif
    do x=1,3
      do y=1,3
        PrXY(x,y) = OKA2P(Genos(l,A), x, y) * PrB(x) * PrPA(y) * &
          SUM(AKA2P(y,x,:) * PrZ)
        do u=1,3
          PrXYU(x,y,u) = OKA2P(Genos(l,A), u, y) * PrB(x) * PrPA(y) * &
             SUM(AKA2P(y,u,:) * PrZ) * AHWE(u,l) 
        enddo
      enddo
    enddo
    PrL(l,1) = LOG10(SUM(PrXY))   ! Parent(A,3-k) offspring of SB 
    if (GA(k)/=-SB)  PrLU(l,1) = LOG10(SUM(PrXYU))  
  endif
  
 !===
  if(maybe(3)) then
    do i=1, ns(SB, k)
      if (nFS(SibID(i,SB,k))==0 .or. Parent(SibID(i,SB,k),3-k)/=Par)  cycle
      call ParProb(l, Par,3-k,SibID(i,SB,k),-5,PrPA)  ! exclude both GPs & Bi & FS of Bi
    enddo
    if (Par < 0) then
      if (Parent(A,3-k)==Par) then
        PrPA = PrPA/OKAP(Genos(l,A),:,l)
        PrPA = PrPA/SUM(PrPA)
      endif
      call ParProb(l, GA(3-k), 3-k, 0, 0, PrZ) 
    else
      call ParProb(l, GA(3-k), 3-k, Par, 0, PrZ) 
    endif
    do x=1,3
      do y=1,3
        PrXY(x,y) = PrPA(y) * SUM(AKA2P(y,x,:) * PrZ) * XPr(2,x,l, SB,k)  !Xpr(2,) = GP's
        do i=1, ns(SB, k)
          PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,SibID(i,SB,k)), x, y)
        enddo
        PrXYU(x,y,:) = PrXYU(x,y,:) * OKAP(Genos(l,A), :, y)
        PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,A), x, y)
      enddo
    enddo
    PrL(l,3) = LOG10(SUM(PrXY)) 
    PrLU(l,3) = LOG10(SUM(PrXYU))    
  endif
  
  !===
  if (maybe(2)) then
    call ParProb(l, GG, k, 0, 0, PrZ)
    call ParProb(l, GX,3-k, -1,0, PrPA)   
    PrXY = 1D0
    do x=1,3  ! SB
      do y=1,3  ! GX
        if (GX < 0) then
          do i=1, ns(-GX,3-k)
            Ei = SibID(i,-GX,3-k)
            if (Ei==A .or. Parent(Ei,k)==-SB .or. nFS(Ei)==0)  cycle
            call ParProb(l, Parent(Ei,k), k, Ei, -1, PrE)
            PrFS = FSLik(l,Ei)
            PrE = PrE * PrFS(:,y)
            if (.not. all(PrE==1D0))  PrXY(x,y) = PrXY(x,y) * SUM(PrE)
          enddo
        endif
        do j=1, nS(SB,k)
          Bj = SibID(j,SB,k)
          if (nFS(Bj)==0) cycle
          PrFS = FSLik(l,Bj)
          if (Parent(Bj,3-k)==GX .and. GX/=0) then
            PrXY(x,y) = PrXY(x,y) * PrFS(x,y)
          else
            call ParProb(l, Parent(Bj,3-k), 3-k, Bj, -1, PrE)
            PrE = PrE * PrFS(x,:)
            if (.not. all(PrE==1D0))  PrXY(x,y) = PrXY(x,y) * SUM(PrE)
          endif
        enddo

        if (GpID(3-k,SB,k)==0) then
          do u=1,3
            PrXYU(x,y,u) = PrXY(x,y) * SUM(AKAP(x,:,l) * PrZ) * PrPA(y) *  &
             OKA2P(Genos(l,A), u, y) * SUM(AKA2P(u,y,:) * AHWE(:,l)) 
          enddo
          PrXY(x,y) = PrXY(x,y) * SUM(AKA2P(x,y,:) * PrPA(y) * PrZ)
        else
          PrXY(x,y) = PrXY(x,y) * SUM(AKA2P(x,y,:) * PrPA(y) * PrZ)
          do u=1,3
            if (Parent(A,3-k)==GX .and. GX/=0) then
              PrXYU(x,y,u) = PrXY(x,y) * OKA2P(Genos(l,A), u, y) * SUM(AKA2P(u,y,:) * AHWE(:,l)) 
            else !if (Parent(A,3-k)==0) then
              PrXYU(x,y,u) = PrXY(x,y) * SUM(AKAP(u,:,l) * OKA2P(Genos(l,A), u, :) *  AHWE(:,l))
            endif
          enddo
        endif
        PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,A), x, y)       
      enddo
    enddo
    PrL(l,2) = LOG10(SUM(PrXY))   ! SB offspring of Parent(A,3-k)
    PrLU(l,2) = LOG10(SUM(PrXYU)) ! A inbred via another parent
  endif
enddo

LLtmp = SUM(PrL, dim=1)
if (maybe(1) .and. Parent(A,3-k)>0 .and. GA(k)/=-SB) then
  LLtmp(1) = LLtmp(1) - Lind(Parent(A,3-k))
endif

LL = impossible
LLU(1:3) = SUM(PrLU, dim=1)
call CalcU(A,k, -SB, k, LLU(4))   ! unrelated & A non-inbred
do x=1,3
  if (.not. maybe(x))  cycle
  if (LLU(x) > LLU(4) .and. LLU(x)/=0d0) then 
    LL(x) = LLtmp(x) - LLU(x) + LLU(4)
  else
    LL(x) = LLtmp(x)
  endif
enddo

end subroutine AddSibInbr

! #####################################################################

! subroutine AddGrandSib(A,SB,k,LL)   ! SB is GP of A, via unobserved parent
! ! NOTE: THIS IS A SPECIAL CASE OF PAIRUA
! use Global
! implicit none

! integer, intent(IN) :: A, SB, k
! double precision, intent(OUT) :: LL

! integer :: l, x,y, m
! double precision :: PrL(nSnp), PrXY(3,3), PrPA(3)

! if (all(Parent(A,:)/=0)) then
  ! LL = NotImplemented
  ! return
! else if (Parent(A,1)==0) then
  ! m=1
! else
  ! m=2
! endif

! LL = missing
! PrL = 0D0
! do l=1,nSnp
  ! call ParProb(l, Parent(A,3-m), 3-m, A, 0, PrPA)
  ! do x=1,3  ! SB
    ! do y=1,3  ! parent of A
      ! PrXY(x,y) = SUM(OKA2P(Genos(l,A), y, :) * PrPA) * AKAP(y,x,l) * XPr(3,x,l, SB,k)     
    ! enddo
  ! enddo
  ! PrL(l) = LOG10(SUM(PrXY))   
! enddo
! LL = SUM(PrL)

! end subroutine AddGrandSib

! #####################################################################

subroutine MergeSibs(SA, SB, k, LL)  
use Global
use CalcLik             
implicit none

integer, intent(IN) :: SA, SB, k
double precision, intent(OUT) :: LL
integer :: G(2), nAB(2), AB(2,maxSibsize), catG, catGG(2), GGP(2), ParPar(2), &
  catA(ns(SA,k)), catB(ns(SB,k)), l,x,y, r,v, Bj, Ai, i,  m, Ei, e,f, j, z
double precision :: PrL(nSnp, 2), PrG(3,2), PrXYZ(3,3,3,2), PrX(3,2), PrE(3), PrH(3), PrFS(3,3)
logical :: AncOK(2), ParIsClone(2,maxSibsize)

LL = missing
G = 0  
do m=1,2  
  if (GpID(m,SA,k) /= 0) then
    if (GpID(m,SB,k) /= 0 .and. GpID(m,SA,k)/=GpID(m,SB,k)) then
      LL = impossible  ! incompatible grandparents
    else
      G(m) = GpID(m,SA,k)  ! including if GP(B) is dummy
    endif
  else if (GpID(m,SA,k) == 0) then
    G(m) = GpID(m,SB,k)
  endif
enddo
if (GpID(k, SA,k)==-SB .or. GpID(k, SB, k)==-SA) then
  LL = impossible
endif
if (LL==impossible) return

AncOK = .TRUE.
call ChkAncest(-SA,k, -SB,k, AncOK(1))
call ChkAncest(-SB,k, -SA,k, AncOK(2))
if (any(.not. AncOK)) then
  LL = impossible
  return
endif

nAB(1) = ns(SA, k)
nAB(2) = ns(SB, k)
AB(1, 1:ns(SA,k)) = SibID(1:ns(SA,k), SA, k)
AB(2, 1:ns(SB,k)) = SibID(1:ns(SB,k), SB, k)

catG = 0
catGG = 0
GGP = 0
if (ANY(G/=0)) then
  do j=1,2
    do i=1,nAB(j)
      if (nFS(AB(j,i))==0) cycle
      if (Parent(AB(j,i), 3-k)==0) cycle
      if (Parent(AB(j,i), 3-k) == G(3-k)) then
         if (catG==0) then
          catG = AB(j,i)
          exit  
        endif   
      endif
      ParPar = getPar(Parent(AB(j,i), 3-k), 3-k)
      do m=1,2
        if (ParPar(m) == G(m) .and. G(m)/=0) then
          catGG(m) = AB(j,i)
          GGP(3-m) = ParPar(3-m)
        endif
      enddo
    enddo
    if (catG /= 0) exit
  enddo
endif

catA = 0
catB = 0
do r = 1, nS(SB, k)
  Bj = SibID(r, SB, k) 
  if (nFS(Bj)==0 .or. Parent(Bj,3-k)==0) cycle
  do v=1,nS(SA,k)
    Ai = SibID(v, SA, k)   
    if (nFS(Ai)==0) cycle
    if (Parent(Ai,3-k) == Parent(Bj,3-k) .and. (Parent(Bj,3-k) < 0 .or. catG==Ai)) then  
      catA(v) = r
      catB(r) = 1 
    endif
  enddo
enddo

ParIsClone = .FALSE.
if (hermaphrodites/=0) then
  if (DumClone(SA,k)/=0) then
    do i=1, nS(SA,k)
      if (Parent(SibID(i,SA,k), 3-k) == -DumClone(SA,k)) then
        ParisClone(1,i) = .TRUE.
      endif
    enddo
  endif
  if (DumClone(SB,k)/=0) then
    do i=1, nS(SB,k)
      if (Parent(SibID(i,SB,k), 3-k) == -DumClone(SB,k)) then
        ParisClone(2,i) = .TRUE.
      endif
    enddo
  endif
endif

PrL = 0D0
do l=1,nSnp
  do m=1,2
    if (m/=k .and. catG/=0) then
      call ParProb(l, G(m), m, -1, 0, PrG(:,m))
    else if (catGG(m)/=0) then  
      if (Parent(catGG(m),3-k) > 0) then
        call ParProb(l, G(m), m, Parent(catGG(m),3-k), 0, PrG(:,m))
      else
        call ParProb(l, G(m), m, 0, 0, PrG(:,m))
      endif
    else
      call ParProb(l, G(m), m, 0, 0, PrG(:,m))
    endif
  enddo
  do x=1,3
    do y=1,3
      do z=1,3
        PrXYZ(x,y,z,:) = AKA2P(x, y, z) * PrG(y,3-k) * PrG(z,k)
      enddo     
    enddo
    PrX(x,1) = SUM(PrXYZ(x,:,:,1))
    PrX(x,2) = PrX(x,1)
  enddo
  
  do z=1,3   ! GP k
    do y=1,3  ! GP 3-k
      if ((y>1 .or. z>1) .and. ALL(catGG==0) .and. catG==0) cycle
  do x=1,3
    do j=1,2
      do r=1, nAB(j)
       if (j==2) then
         if (catB(r)==1) cycle ! done as FS of an A
       endif
        Ai = AB(j,r)
        if (NFS(Ai) == 0) cycle
        if (catG==Ai .or. ParIsClone(j,r)) then
          PrE = 1D0
        else if (ANY(catGG == Ai)) then
          if (ALL(catGG == Ai)) then  ! parent(Ai, 3-k) FS with SB
            PrE = AKA2P(:,z,y)
          else      
            do m=1,2
              if (catGG(m)==Ai) then
                if (Parent(Ai,3-k)>0) then
                  call ParProb(l, GGP(3-m), 3-m, Parent(Ai,3-k),0,PrH)
                else
                  call ParProb(l, GGP(3-m), 3-m, 0,0,PrH)
                endif
                do e=1,3
                  if (m==k) then
                    PrE(e) = SUM(AKA2P(e, z, :) * PrH)
                  else
                    PrE(e) = SUM(AKA2P(e, y, :) * PrH)
                  endif
                enddo
              endif
            enddo
          endif
          if (Parent(Ai,3-k)>0) then
            PrE = PrE * OcA(:,Genos(l, Parent(Ai,3-k)))
          endif
!          PrE = PrE/SUM(PrE)
        else
          call ParProb(l, Parent(Ai,3-k), 3-k, -1, 0, PrE)
        endif

        if (Parent(Ai,3-k) < 0) then 
          do e=1,3
            if (catG==Ai .and. y/=e)  cycle                               
            do f=1, nS(-Parent(Ai,3-k), 3-k)
              Ei = SibID(f, -Parent(Ai,3-k), 3-k)
              if (nFS(Ei) == 0) cycle
              if (Parent(Ei, k)==-SB .or. Parent(Ei,k)==-SA) cycle  
              if (catGG(3-k)>0) then
                if (Parent(catGG(3-k),3-k) == Ei)  cycle
              endif  
              call ParProb(l, Parent(Ei, k), k, Ei, -1, PrH)
              PrFS = FSLik(l,Ei)
              PrH = PrH * PrFS(:,e)
              if (.not. ALL(PrH==1D0))  PrE(e) = PrE(e) * SUM(PrH)
            enddo
          enddo
        endif

        if (.not. ALL(PrE==1D0)) then
          if (catG==Ai) then
            PrXYZ(x,y,z,1) = PrXYZ(x,y,z,1) * PrE(y)
          else if (ANY(catGG/=0) .or. catG/=0) then
            PrXYZ(x,y,z,1) = PrXYZ(x,y,z,1) * SUM(PrE)
          else if (ParIsClone(j,r)) then
            PrX(x,1) = PrX(x,1) * PrE(x)
          else 
            PrX(x,1) = PrX(x,1) * SUM(PrE)
          endif
        endif
        
        PrFS = FSLik(l,Ai)
        PrE = PrE * PrFS(:,x)

        if (j==1) then
          if (catA(r)/=0) then
            Bj = AB(2, catA(r))
            PrFS = FSLik(l,Bj)
            PrE = PrE * PrFS(:,x)
          endif
        endif

        if (.not. ALL(PrE==1D0)) then
          if (catG==Ai) then
            PrXYZ(x,y,z,2) = PrXYZ(x,y,z,2) * PrE(y)
          else if (ANY(catGG/=0) .or. catG/=0) then
            PrXYZ(x,y,z,2) = PrXYZ(x,y,z,2) * SUM(PrE)
          else if (ParIsClone(j,r)) then
            PrX(x,2) = PrX(x,2) * PrE(x)
          else 
            PrX(x,2) = PrX(x,2) * SUM(PrE)
          endif 
        endif
      enddo  ! r
    enddo  ! j
  enddo  ! x
  enddo  ! y (catGG>0 only)
  enddo  ! z (catGG>0 only)
  do m=1,2
    if (ANY(catGG/=0) .or. CatG/=0) then
      PrL(l,m) = LOG10(SUM(PrXYZ(:,:,:,m)))
    else
      PrL(l,m) = LOG10(SUM(PrX(:,m)))
    endif
  enddo
enddo
LL = SUM(PrL(:,2)) - SUM(PrL(:,1))

end subroutine MergeSibs

! #####################################################################

subroutine AddFS(A, SB, kB, SA, kAx, LL, TopSib, dLL)  ! A/SA FS with any B?
use Global
use CalcLik            
implicit none

integer, intent(IN) :: A, SB, kB, SA, kAx
integer, intent(OUT) :: TopSib    ! most likely FS of A within SB
double precision, intent(OUT) :: LL, dLL(maxSibSize)  ! dLL
integer :: PA, kA, AncA(2,mxA), GA(2), InbrX, Par(nS(SB,kB)), MaybeFS(nS(SB,kB)), &
  f, i, Bj, Inbr(nS(SB,kB)),  GB(ns(SB,kB), 2), DoQuick, l,x,y,g,h,Ei,z,j,Rj, ParAtmp 
double precision :: ALR, LRQ, LLtmp(2), LLUX, PrL(nSnp, nS(SB,kB),2), PrXY(3,3,2), &
  PrG(3), PrX(3,2), PrY(3,2), PrW(3), PrZ(3), PrV(3), PrB(3,3), PrF(3,3) 
logical :: ParOK, ParBisClone(ns(SB,kB))

LL = missing
TopSib = 0
dLL = missing

if (nS(SB,kB)==0) then
  LL = impossible
  return   ! nobody to be FS with
endif 

PA = 0
if (A /= 0) then
  if (kAx==1 .or. kAx==2) then
    kA = kAx
  else
    kA = 1
  endif
  PA = Parent(A, 3-kB)
  call GetAncest(A, kA, AncA)
else !if (SA /= 0) then   ! TODO: does it matter if kA=kB?
  kA = kAx
  PA = GpID(3-kB, SA, kA)
  call GetAncest(-SA, kA, AncA)
endif
GA = getPar(PA, 3-kB)

if (A/=0) then
  if (Parent(A,kB)/=0 .and. Parent(A,kB)/=-SB) then
    LL = impossible
  else
    call ChkValidPar(A,Sex(A), -SB,kB, ParOK)
    if (.not. ParOK)  LL = impossible
  endif
else !if (SA/=0) then
  if (GpID(kB, SA, kA)/=0 .and. GpID(kB, SA, kA)/=-SB) then
    LL = impossible
  else
    call ChkValidPar(-SA,kA, -SB,kB, ParOK)
    if (.not. ParOK)  LL = impossible
  endif
endif
if (LL /= missing) return

InbrX = 0
if (ANY(AncA(kB, 3:mxA) == -SB) .and. A < 0) then  
  LL = NotImplemented   ! TODO: check  
else if (A > 0 .and. PA/=0) then
  if (AncA(kB,5-kB)==-SB .and. PA<0) then
    InbrX = -1  ! P-O mating
  else if (ANY(Parent(SibID(1:nS(SB,kB),SB,kB),3-kB)==GA(3-kB)) .and. GA(3-kB)/=0) then
    InbrX = -2
  else if (GpID(3-kB, SB, kB)==PA) then
    InbrX = -3
  else if (ANY(AncA(kB, 3:mxA) == -SB)) then
    if (.not. any(Parent(SibID(1:ns(SB,kB),SB,kB),3-kB)==PA)) then
      LL = NotImplemented
    endif
  endif
endif
if (LL /= missing) return

if (SA/=0 .and. GpID(3-kA, SB, kB)/=0) then
  do i=1, ns(SA,kA)
    if (Parent(SibID(i,SA,kA), 3-kA) == GpID(3-kA, SB, kB)) then
      LL = NotImplemented
      return
    endif
  enddo
endif

Par = 0  ! shared parent 3-kB  (cand. parent(kB) == SB)
MaybeFS = 1
if (PA/=0 .and. ANY(Parent(SibID(1:nS(SB,kB),SB,kB), 3-kB) == PA)) then
  MaybeFS = 0
  do f=1, nS(SB,kB)
    if (NFS(SibID(f, SB, kB))==0) then
      MaybeFS(f) = -1
    else if (Parent(SibID(f, SB, kB), 3-kB) == PA) then
      MaybeFS(f) = 1
      Par(f) = PA
    else
      Par(f) = Parent(SibID(f, SB, kB), 3-kB)
    endif
  enddo

else
 do f=1, nS(SB,kB)
  if (NFS(SibID(f, SB, kB))==0) then
    MaybeFS(f) = -1
    cycle
  endif   
  do i=1,nFS(SibID(f, SB, kB))
    Bj = FSID(i, SibID(f, SB, kB))
    if (A == Bj) then
      LL = AlreadyAss
    else if (A >0) then
      if (Parent(A,3-kB) == Bj) then  
         MaybeFS(f) = 0     ! can't be FS with own parent    
      else if (Parent(Bj, 3-kB) == A) then
        MaybeFS(f) = 0
      else 
        call ChkValidPar(A, Sex(A), Parent(Bj,3-kB), 3-kB, ParOK)
        if (.not. ParOK) then
          MaybeFS(f) = 0
        else
          call CalcAgeLR(A, Sex(A), Bj, Sex(Bj), kB, 2, .TRUE., ALR)
          if (ALR==impossible)  MaybeFS(f) = 0
        endif
      endif
      
    else if (SA/=0) then
      if (kA/=kB .and. Parent(Bj, 3-kB) == -SA) then
        MaybeFS(f) = 0  ! cannot be FS with own parent
        LL = NotImplemented   ! TODO: implement. 
        cycle
      else
        call ChkValidPar(-SA,kA, Parent(Bj,3-kB), 3-kB, ParOK)
        if (.not. ParOK) then
          MaybeFS(f) = 0
        else
          call CalcAgeLR(-SA,kA, Bj, Sex(Bj), kB, 2, .TRUE., ALR)
          if (ALR==impossible)  MaybeFS(f) = 0
        endif  
      endif
    endif
    if (Bj == PA .or. (A/=0 .and. A == Parent(Bj, 3-kB))) then
      MaybeFS(f) = 0
      cycle
    endif
    if (PA>0) then
      if (any(Parent(PA,:)==Bj)) then
        MaybeFS(f) = 0
        cycle
      endif
    endif
    
    Par(f) = Parent(Bj, 3-kB)
    if (PA/=0 .and. PA/=Par(f) .and. Par(f)/=0) then
      MaybeFS(f) = 0
    else if (Par(f)==0) then
      call CalcP2(Bj, Sex(Bj), -SB, PA, kB, LRQ)  
      if (LRQ == impossible) then  
        MaybeFS(f) = 0
        cycle
      endif 
      Par(f) = PA
    endif
  enddo
 enddo
endif
if (LL /= missing) return
if (ALL(MaybeFS==0 .or. MaybeFS==-1)) then
  LL = impossible
  return
endif

Inbr = 0
GB = 0
do f=1, nS(SB,kB)    
  GB(f,:) = getPar(Par(f), 3-kB)
  if (GB(f,kB) == -SB)  Inbr(f) = 1 
  if (Par(f) == GpID(3-kB, SB, kB) .and. Par(f)/=0) then   ! DoQuick = -1
    Inbr(f) = 2
  endif
enddo

do f=1, nS(SB,kB)
  if (nFS(SibID(f, SB, kB))==0) cycle
  if (MaybeFS(f)<1 .or. Par(f)==0 .or. Par(f)==PA) cycle
  if (A/=0)  call ChkValidPar(A, Sex(A), Par(f), 3-kB, ParOK)
  if (SA/=0) call ChkValidPar(-SA, kA, Par(f), 3-kB, ParOK)
  if (.not. ParOK)  MaybeFS(f) = 0  
enddo
if (ALL(MaybeFS==0 .or. MaybeFS==-1)) then
  LL = impossible
  return
endif

call ChkDoQuick(SB,kB,DoQuick)

ParBisClone = .FALSE.
if (hermaphrodites/=0 .and. DoQuick==-3) then   ! -3: SB has a dummyclone
  do f=1, nS(SB,kB)
    if (DumClone(SB,kB) == -Parent(SibID(f,SB,kB), 3-kB)) then
      ParBisClone(f) = .TRUE.
    endif
  enddo
endif

if (A/=0 .and. nYears>1 .and. ns(SB,kB)>1) then    ! check if A is more likely to be parent of sib instead  
  do f=1, nS(SB,kB)
    Bj = SibID(f, SB, kB)
    if (MaybeFS(f)<1 .or. Par(f)/=0 .or. Parent(Bj, 3-kB)/=0) cycle  
    if (AgeDiff(Bj, A) <= 0) cycle
    if (Sex(A)<3 .and. Sex(A)/=3-kB)  cycle  ! TODO check both sexes?
    call ChkValidPar(Bj, 3, A, 3-kB, ParOK)  
    if (.not. ParOK)  cycle
    call CalcU(-SB, kB, A, kB, LLtmp(1))
    call setParTmp(Bj, 3, A, 3-kB)
    call CalcU(-SB, kB, A, kB, LLtmp(2))
    call setParTmp(Bj, 3, 0, 3-kB)
    call CalcCLL(SB,kB)
    if (LLtmp(1) - LLtmp(2) < TA) then
      MaybeFS(f) = 0
    endif
  enddo
endif
if (ALL(MaybeFS==0 .or. MaybeFS==-1)) then
  LL = impossible
  return
endif

LLtmp = missing               
if (A>0 .and. PA/=0 .and. any(Parent(SibID(1:ns(SB,kB),SB,kB),3-kB)==PA) .and. &  ! A already HS via 3-k
   (InbrX==-2 .or. InbrX==-3 .or. DoQuick>1))  then  
  call CalcU(-SB, kB, A, Sex(A), LLUX)
  if (Parent(A,3-kB) < 0) then
    call CalcU(-SB, kB, Parent(A,3-kB), 3-kB, LLtmp(1))
  else
    LLtmp(1) = LLUX
  endif
  ParAtmp = Parent(A,kB)  ! 0 or -SB                           
  call setParTmp(A, sex(A), -SB, kB)
  if (Parent(A,3-kB) < 0) then
    call CalcU(-SB, kB, Parent(A,3-kB), 3-kB, LLtmp(2))
  else
    LLtmp(2) = CLL(SB,kB)
  endif
  call setParTmp(A, sex(A), ParAtmp, kB) 
  call CalcCLL(SB,kB)
  LL = LLUX + LLtmp(2) - LLtmp(1)  

  do f=1, ns(SB,kB)
    if (nFS(SibID(f, SB, kB))==0) cycle
    if (Parent(SibID(f, SB, kB), 3-kB) == Parent(A, 3-kB)) then
      TopSib = SibID(f, SB, kB)
      dLL(f) = LLtmp(2) - LLtmp(1)
    endif
  enddo
  
else 
  PrL = 0D0
  
  if (A>0 .and. DoQuick==-2) then  ! SB are all FS
    do f=1, nS(SB,kB)
      if (MaybeFS(f) < 1) cycle
      Bj = SibID(f, SB, kB)
      do l=1, nSnp
        PrB = FSLik(l,Bj)
        do x=1,3
          do y=1,3
            PrXY(x,y,:) = XPr(2,x,l,SB,kB) * XPr(2,y,l,-par(f),3-kB) * PrB(x,y)
            PrXY(x,y,1) = PrXY(x,y,1) * OKA2P(Genos(l,A),x,y)
            if (PA/=0) then
              PrXY(x,y,2) = PrXY(x,y,2) * OKAP(Genos(l,A), y, l)
            else
              PrXY(x,y,2) = PrXY(x,y,2) * OHWE(Genos(l,A), l)
            endif
          enddo
        enddo
        do i=1,2
          PrL(l,f,i) = LOG10(SUM(PrXY(:,:,i)))
        enddo
      enddo
    enddo
  
  else
    do l=1,nSnp
      if (any(inbr == 2)) then
        call ParProb(l, GpID(kB,SB,kB), kB, -1,0, PrG)
      endif 
      do f=1, nS(SB,kB)
        if (MaybeFS(f) < 1) cycle
        if (any(inbr == 2)) then
          PrX = 1D0
        else
          PrX(:,1) = XPr(2,:,l, SB, kB)
          PrX(:,2) = PrX(:,1)
        endif
        do x=1,3   ! SB 
          do g=1,nS(SB,kB)
            Bj = SibID(g, SB, kB)
            if (NFS(Bj) == 0) cycle
            if (Inbr(g)==1) then
              if (g==f) then
                PrY(:,1) = 1D0
              else
                call ParProb(l, Parent(Bj, 3-kB), 3-kB, Bj,-5, PrY(:,1))  
                ! no GPs & no Bj & no FS of Bj 
              endif
              call ParProb(l, GB(g,3-kB), 3-kB, 0,0,PrW)
              do y=1,3
                PrY(y,1) = PrY(y,1) * SUM(AKA2P(y,:,x) * PrW)
              enddo
            else if (ParBisClone(g)) then
              PrY(:,1) = 1D0  
            else
              call ParProb(l, Par(g), 3-kB, -1,0, PrY(:,1))
              if (Inbr(g)==2) then
                do y=1,3
                  PrY(y,1) = PrY(y,1) * SUM(AKA2P(x,y,:) * PrG)
                enddo
              endif
            endif
            
            PrY(:,2) = PrY(:,1)  ! 1: FS, 2: HS via 3-k or U 
            do y=1,3
              if (ParBisClone(g) .and. y/=x)  cycle
              if (Par(g) < 0) then 
                do h = 1, nS(-Par(g), 3-kB)
                  Ei = SibID(h, -Par(g), 3-kB)
                  if (Parent(Ei, kB) == -SB .or. Ei==A) cycle  
                  if (NFS(Ei) == 0) cycle 
                  if (g==f .and. Parent(Ei,kB)<0) then 
                    call ParProb(l, Parent(Ei,kB), kB, -1,0, PrZ)
                    do z=1,3
                      do j=1, ns(-parent(Ei,kB),kB)
                        Rj = SibID(j,-parent(Ei,kB),kB)
                        if (Parent(Rj,3-kB) == Par(g) .and. Par(g)/=0)  cycle
                        if (nFS(Rj)==0)  cycle
                        call ParProb(l, Parent(Rj,3-kB),3-kB, Rj,-1, PrV)
                        PrF = FSLik(l,Rj)
                        PrV = PrV * PrF(:,z)
                        if (.not. all(PrV==1D0))  PrZ(z) = PrZ(z) * SUM(PrV)
                      enddo
                    enddo
                  else
                    call ParProb(l, Parent(Ei,kB), kB, Ei,-1, PrZ)
                  endif
                  PrF = FSLik(l,Ei)
                  PrZ = PrZ * PrF(:,y)
                  if (.not. all(PrZ==1D0))  PrY(y,:) = PrY(y,:) * SUM(PrZ)
                enddo
              endif
              
              PrB = FSLik(l,Bj)
              PrY(y,:) = PrY(y,:) * PrB(x,y)
              
              if (g==f) then  
                if (A/=0) then
                  PrY(y,1) = PrY(y,1) * OKA2P(Genos(l,A), x, y)
                  PrY(y,2) = PrY(y,2) * OHWE(Genos(l,A), l)   ! vs unrelated !
                else if (SA/=0) then
                  PrY(y,1) = PrY(y,1) * SUM(XPr(1,:,l, SA,kA) *AKA2P(:,x,y))
                  if (PA/=0) then
                    PrY(y,2) =PrY(y,2) *SUM(XPr(1,:,l, SA,kA) *AKAP(:,y, l))
                  else
                    PrY(y,2) = PrY(y,2) * SUM(XPr(1,:,l, SA,kA) * AHWE(:,l))
                  endif
                endif 
              endif          
            enddo  ! y
            do i=1,2
              if (ParBisClone(g)) then
                PrX(x,i) = PrX(x,i) * PrY(x,i)
              else if (.not. all(PrY(:,i)==1D0)) then
                PrX(x,i) = PrX(x,i) * SUM(PrY(:,i))
              endif
            enddo
          enddo  ! g
        enddo  ! x
        PrL(l,f,:) = LOG10(SUM(PrX, DIM=1))
      enddo  ! f   
    enddo  ! l
  endif       

  dLL = impossible
  do f = 1, nS(SB, kB)
    if (MaybeFS(f)<1) cycle
    Bj = SibID(f, SB, kB)
    if (NFS(Bj) == 0) then
      cycle
    else if (nFS(Bj) == 1) then
      dLL(f) = SUM(PrL(:, f,1)) - SUM(PrL(:,f,2))
    else
      do g=1, nS(SB,kB)
        if (Parent(SibID(g, SB, kB), 3-kB) == Parent(Bj, 3-kB)) then
          dLL(g) = SUM(PrL(:, f,1)) - SUM(PrL(:,f,2))
        endif
      enddo
    endif
  enddo

  if (A/=0) then
    ! calc LLUX with parA=0, else LL(FS)|PA/=0 /= LL(FS)|PA==0
    if(PA/=0)  call setParTmp(A,kA,0,3-kB)
    call CalcU(A,kA, -SB, kB, LLUX)
    if(PA/=0)  call setParTmp(A,kA,PA,3-kB)         
    LL = MAXVAL(dLL, MASK=dLL/=impossible) + LLUX
    TopSib = MAXLOC(dLL, MASK=dLL/=impossible, DIM=1)
    if(TopSib>0)  TopSib = SibID(TopSib, SB, kB)
  else if (SA/=0) then
    call CalcU(-SA, kA, -SB, kB, LLUX)
    do f = 1, nS(SB, kB)
      if (dLL(f)==impossible) cycle
      if (Par(f)==0 .and. nS(SA,kA)>1) then
        dLL(f) = dLL(f) + LLUX
      else  ! consider changes in SA (e.g. inbreeding loops) 
        call setParTmp(-SA,kA, Par(f), 3-kB)
        call PairUA(-SA, -SB, kA, kB, dLL(f))  
      endif
    enddo
    call setParTmp(-SA,kA, PA, 3-kB)
    call CalcCLL(SA,kA)
    LL = MaxLL(dLL) 
  endif  
endif

end subroutine AddFS

! #####################################################################

subroutine AddParent(A, SB, k, LL)  ! is A parent of sibship SB?  (replace dummy)
use Global
use CalcLik             
implicit none

integer, intent(IN) :: A, SB, k
double precision, intent(OUT) :: LL
integer :: GG(2), GB(2), m, n, Inbr, DoQuick, l,x,y,Bi, Ei, i,j,z
double precision :: LRQ, PrL(nSnp), PrGG(3,2), PrXY(3,3,2),  PrZ(3), PrE(3), PrB(3,3), PrPA(3)
logical :: AncOK, ParBisClone(ns(SB,k)), Aselfed

LL = missing
do i=1,ns(SB,k)
  if (AgeDiff(SibID(i,SB,k), A) <= 0) then
    LL = impossible
    return
  endif
enddo

call ChkAncest(A,0, -SB,k, AncOK) ! e.g. if age A unknown, or all age B's unknown
if (.not. AncOK) then
  LL = impossible
  return
endif

if (ns(SB,k)==0 .and. any(GpID(:,SB,k)==0) .and. any(GpID(:,SB,k)>0)) then
  do m=1,2
    if (GpID(m,SB,k)>0) then
      call pairPO(A, GpID(m,SB,k), k, 1, LL)  ! requires 'focal'
      if (LL < 0d0)  LL = LL - Lind(GpID(m,SB,k))
      return
    endif
  enddo
endif

GB = GPID(:,SB,k)   ! current grandparents
GG = 0   ! grandparents if A replaces dummy
do m=1,2
  if (Parent(A,m)/= 0) then
    if (GB(m)/= 0 .and. GB(m)/=Parent(A,m)) then
      LL = impossible
      return
    else
      GG(m) = Parent(A,m)
    endif
  else if(GB(m)/=0) then
    GG(m) = GB(m)
  endif
enddo

LRQ = missing
do n=1, ns(SB,k)
  call CalcP2(SibID(n,SB,k), 3, A, Parent(SibID(n,SB,k), 3-k), k, LRQ)
  if (LRQ == impossible) then
    LL = impossible
    return
  endif
enddo
call CalcP2(A, k, GB(1), GB(2), 1, LRQ)
if (LRQ == impossible) then
  LL = impossible
  return
endif

Inbr = 0
if (GG(3-k)/=0) then
  do n=1, nS(SB, k)
    if (nFS(SibID(n,SB,k))==0) cycle
    if (Parent(SibID(n,SB,k), 3-k)==GG(3-k)) then
      Inbr = n
    endif
  enddo 
endif

call ChkDoQuick(SB,k,DoQuick)

ParBisClone = .FALSE.
Aselfed = .FALSE.
if (hermaphrodites/=0) then
  if (DoQuick==-3) then
    do i=1, nS(SB,k)
      if (DumClone(SB,k) == -Parent(SibID(i,SB,k), 3-k)) then
        ParBisClone(i) = .TRUE.
      endif
    enddo
  endif
  if (all(GG > 0) .and. GG(1) == GG(2)) then
    Aselfed = .TRUE.
  else if (all(GG < 0)) then
    if (DumClone(-GG(1),1) == -GG(2))  ASelfed = .TRUE.
  endif
endif

PrL = 0D0
do l=1,nSnp
  if (Inbr==0 .and. DoQuick>0) then
    do m=1,2
      call ParProb(l, GG(m), m, A,0, PrGG(:,m))   
    enddo
    do x=1,3
      do y=1,3
        PrXY(x,y,2) = XPr(1,x,l, SB,k) * OcA(x,Genos(l,A)) * &
          SUM(AKA2P(x, y, :) *  PrGG(y,3-k) * PrGG(:, k))
      enddo
    enddo
    PrL(l) = LOG10(SUM(PrXY(:,:,2)))
  else
    call ParProb(l, GG(k), k, A,0, PrGG(:,k))     
    call ParProb(l, GG(3-k), 3-k, -1,0, PrGG(:,3-k))
    call ParProb(l, Parent(A,k), k, A,0, PrPA)     
    do x=1,3  ! A
      do y=1,3  ! G 3-k
        if (Aselfed) then
          PrXY(x,y,:) = OcA(x,Genos(l,A)) * AKA2P(x, y, y) * PrGG(y, k)  ! TODO?
        else
          if (Parent(A,3-k)==GG(3-k) .and. GG(3-k)/=0) then
            PrXY(x,y,1) = XPr(2,x,l, SB,k) * SUM(OKA2P(Genos(l,A),y,:) * PrPA)
          else
            PrXY(x,y,1) = XPr(2,x,l, SB,k) * SUM(LindX(:,l,i))
          endif
          PrXY(x,y,2) = OcA(x,Genos(l,A)) * SUM(AKA2P(x, y, :) *PrGG(y,3-k) *PrGG(:, k)) 
        endif
        do n=1, nS(SB,k)
          Bi = SibID(n, SB, k)
          if (nFS(Bi)==0) cycle
          if (Inbr == n .or. ParBisClone(n)) then
            PrZ = 1D0
          else
            call ParProb(l, Parent(Bi,3-k), 3-k, -1,0, PrZ)
          endif
          
          if (Parent(Bi,3-k) < 0) then
            do z=1,3
              if (Inbr == n .and. z/=y) cycle
              if (ParBisClone(n) .and. z/=x)  cycle 
              do i=1, ns(-Parent(Bi,3-k),3-k)
                Ei = SibID(i,-Parent(Bi,3-k),3-k)
                if (Parent(Ei,k)==-SB .or. nFS(Ei)==0)  cycle
                call ParProb(l, Parent(Ei,k), k, Ei, -1, PrE)
                do j=1, nFS(Ei)
                  if (FSID(j,Ei)==A)  cycle
                  PrE = PrE * OKA2P(Genos(l,FSID(j,Ei)),z,:)
                enddo
                if (.not. ALL(PrE==1D0))  PrZ(z) = PrZ(z) * SUM(PrE)  
              enddo
            enddo
          endif
          
          PrB = FSLik(l,Bi)
          if (Inbr == n) then
            PrXY(x,y,:) = PrXY(x,y,:) * PrZ(y)
            PrXY(x,y,2) = PrXY(x,y,2) * PrB(x,y)
          else if (ParBisClone(n)) then
            PrXY(x,y,:) = PrXY(x,y,:) * PrZ(x)
            PrXY(x,y,2) = PrXY(x,y,2) * PrB(x,x)
          else
            if (.not. all(PrZ==1D0)) PrXY(x,y,1) = PrXY(x,y,1) * SUM(PrZ)
            PrZ = PrZ * PrB(:,x)
            PrXY(x,y,2) = PrXY(x,y,2) * SUM(PrZ)
          endif
        enddo
      enddo
    enddo
    PrL(l) = LOG10(SUM(PrXY(:,:,2))) - LOG10(SUM(PrXY(:,:,1)))
  endif
enddo

LL = SUM(PrL)
if (LL < -HUGE(0D0))  LL = impossible
!if (Inbr==0 .and. DoQuick>0)  LL = LL + Lind(A)

end subroutine AddParent

! #####################################################################

subroutine AddGP(A, SB, k, LL)  ! add A as a grandparent to sibship SB
use Global
use CalcLik
implicit none 

integer, intent(IN) :: A, SB, k
double precision, intent(OUT) :: LL
integer :: cat, catG, DoQuick, curGP(2), l, x,y, m, i,  z, g, Bi,  Ei,w,v,j, pari
double precision :: PrL(nSnp), LLtmp(2), LLU, PrG(3), PrPA(3,2), &
  PrXYZ(3,3,3,2), PrP(3), PrW(3), PrE(3), PrF(3,3), PrA(3), PrXY(3,3)
logical :: ParOK, ParBisClone(ns(SB,k))

LL = missing
if (Sex(A)<3) then
  m = Sex(A)
else if (GpID(1,SB,k)==0) then
  m = 1
else if (GpID(2,SB,k)==0) then
  m = 2
else
  LL = impossible
  return
endif

call ChkValidPar(-SB,k, A, m, ParOK)
if (.not. ParOK)  then
  LL = impossible
  return
endif

cat = 0
catG = 0
curGP = GPID(:, SB,k)
pari = 0                 
if (Parent(A, 3-k)==GpID(3-k,SB,k) .and. Parent(A, 3-k) /= 0) then
  cat = 1
else
  do i=1,nS(SB,k)
    if (nFS(SibID(i,SB,k))==0) cycle
    if (Parent(SibID(i, SB, k), 3-k) == 0) cycle
    if (Parent(SibID(i, SB, k), 3-k) == A) then
      cat = 4
    else if (Parent(SibID(i, SB, k), 3-k) == Parent(A, 3-k)) then
      catG = i
    else if (Parent(SibID(i, SB, k), 3-k) == GpID(3-k,SB,k)) then
      cat = 2
      exit
    else if (Parent(SibID(i,SB,k),3-k) < 0) then
      pari = Parent(SibID(i,SB,k),3-k) 
      if (GpID(k, -pari, 3-k) == -SB) then
        cat = 3
        exit
      else if (any(GpID(:, -pari, 3-k) == A)) then  ! double GP
        cat = 5
        exit
      else if (any(GpID(:, -pari, 3-k) == curGP .and. curGP/=0)) then
        cat = 6       
        exit
      endif
    endif     
  enddo
endif 

call ChkDoQuick(SB,k,DoQuick)    

if (DoQuick==2) then  ! inbreeding: Parent(Bj,3-k) = Offspr(i); offspr(i) < 0
  Bi = 0
  do i=1, nS(SB,k)
    if (Parent(SibID(i,SB,k), 3-k) < 0) then
      Bi = Parent(SibID(i,SB,k), 3-k)
      if (GpID(k, -Bi, 3-k) == -SB)  exit
    endif
  enddo
endif 

ParBisClone = .FALSE.
if (hermaphrodites/=0 .and. DoQuick==-3) then
  do i=1, nS(SB,k)
    if (DumClone(SB,k) == -Parent(SibID(i,SB,k), 3-k)) then
      ParBisClone(i) = .TRUE.
    endif
  enddo
endif

if (Complx==0 .and. Mate(A)/=0) then  
  if ((Mate(A) /= curGP(3-m) .and. curGP(3-m)/=0) .or. &
  (k==3-m .and. Mate(A)==-SB)) then
    LL = impossible
    return
  else if (curGP(3-m) == 0) then
    curGP(3-m) = Mate(A)
  endif
endif

PrL = 0D0
LLU = missing
LLtmp = missing

if (DoQuick==1 .and. cat==0 .and. catG==0 .and. all(parent(A,:)==0) .and. all(curGP==0)) then
  do l=1,nSnp
    PrXY = 0D0
    do x=1,3  ! SB
      do y=1,3
        PrXY(x,y) = XPr(1,x,l, SB,k) * AKAP(x,y,l) * AHWE(y,l) * OcA(y,Genos(l,A))
      enddo
    enddo
    PrL(l) = LOG10(SUM(PrXY))
  enddo
  LL = SUM(PrL)   ! +Lind(A) ?

else if (cat/=0 .or. DoQuick>1 .or. (Parent(A,3-k)<0 .and. catG/=0)) then  
! inbreeding loop present / will be created:   
  call CalcU(-SB, k, A, 3-k, LLU)
  if (curGP(3-m) < 0) then
    call CalcU(-SB, k, curGP(3-m), 3-m, LLtmp(1))
  else if (Parent(A,3-k)<0) then
    call CalcU(-SB, k, Parent(A,3-k), 3-k, LLtmp(1))
  else if (cat==5 .or. cat==6) then  ! .or. cat==3 ?
    call CalcU(-SB, k, pari, 3-k, LLtmp(1))
  else if (DoQuick==2) then
    call CalcU(-SB, k, Bi, 3-k, LLtmp(1))  
  endif
  
  call setParTmp(-SB, k, A, m)
  if (curGP(3-m) < 0) then
    call CalcU(-SB, k, curGP(3-m), 3-m, LLtmp(2))
    LL = LLU + (LLtmp(2) - LLtmp(1))
  else if (Parent(A,3-k)<0) then
    call CalcU(-SB, k, Parent(A,3-k), 3-k, LLtmp(2))
    LL = LLU + (LLtmp(2) - LLtmp(1))
  else if (cat==5 .or. cat==6) then
    call CalcU(-SB, k, pari, 3-k, LLtmp(2))
    LL = LLU + (LLtmp(2) - LLtmp(1))
   else if (DoQuick==2) then
    call CalcU(-SB, k, Bi, 3-k, LLtmp(2)) 
    LL = LLU + (LLtmp(2) - LLtmp(1))                             
  else
    LL = CLL(SB,k) + Lind(A)
  endif
  call setParTmp(-SB, k, curGP(m), m)
  if (curGP(3-m) < 0)  call CalcCLL(-curGP(3-m), 3-m)
  if (Parent(A,3-k) < 0)  call CalcCLL(-Parent(A,3-k), 3-k)
  
else
  do l=1,nSnp
    call ParProb(l, curGP(3-m), 3-m, 0, 0, PrG)
     if (catG/=0) then  
       call ParProb(l, Parent(A,3-k), 3-k, -1, 0, PrPA(:,3-k))
       call ParProb(l, Parent(A,k), k, -1, 0, PrPA(:,k))
     else
      call ParProb(l, A, m, 0, 0, PrA)   
     endif
    
    PrXYZ = 0D0
    do x=1,3  ! SB
      do y=1,3
        if (catG == 0) then
          PrXYZ(x,y,1,:) = SUM(AKA2P(x, :, y) * PrG * PrA(y))
        else
          do z=1,3  ! Parent(A,3-k) 
            do g=1,3
              PrP(g) = AKA2P(x,g,y) * PrG(g) * OcA(y,Genos(l,A)) * &
               SUM(AKA2P(y,z,:) * PrPA(:,k) * PrPA(z,3-k))  ! approx
            enddo
            PrXYZ(x,y,z,:) = SUM(PrP)
          enddo
        endif
      enddo
      if (DoQuick > 0 .and. catG==0) then
        PrXYZ(x,:,1,2) = PrXYZ(x,:,1,2) * XPr(1,x,l, SB,k)
      else
        do y=1,3  ! A
          do z=1,3  ! Parent(A,3-k)   
            if (z>1 .and. catG==0)  cycle
            do i=1, nS(SB,k)
              Bi = SibID(i, SB, k)
              if (nFS(Bi)==0) cycle
              if (catG == i .or. ParBisClone(i)) then
                PrW = 1D0
              else
                call ParProb(l, Parent(Bi,3-k), 3-k, -1, 0, PrW)
              endif
              if (Parent(Bi,3-k)<0) then
                do w=1,3
                  if (catG==i .and. w/=z)  cycle
                  if (ParBisClone(i) .and. w/=x)  cycle                                     
                  do v=1, nS(-Parent(Bi, 3-k), 3-k)
                    Ei = SibID(v, -Parent(Bi, 3-k), 3-k)  
                    if (NFS(Ei) == 0) cycle
                    if (Parent(Ei, k) == -SB) cycle
                    call ParProb(l, Parent(Ei, k), k, Ei,-1, PrE)
                    if (catG==i) then 
                      do j=1, nFS(Ei)
                        if (FSID(j,Ei)==A)  cycle
                        PrE = PrE * OKA2P(Genos(l,FSID(j,Ei)),w,:)
                      enddo
                    else
                      PrF = FSLik(l,Ei)
                      PrE = PrE * PrF(:,w)
                    endif
                    if (.not. ALL(PrE==1D0))  PrW(w) = PrW(w) * SUM(PrE)  
                  enddo  
                enddo
              endif

              if (catG == i) then
                PrXYZ(x,y,z,1) = PrXYZ(x,y,z,1) * PrW(z)
               else if (catG /= 0) then
                if (.not. all(PrW==1D0))  PrXYZ(x,y,z,1) = PrXYZ(x,y,z,1) * SUM(PrW)
              else if (ParBisClone(i)) then
                PrXYZ(x,y,:,1) = PrXYZ(x,y,:,1) * PrW(x)
              else if (.not. all(PrW==1D0)) then
                PrXYZ(x,y,:,1) = PrXYZ(x,y,:,1) * SUM(PrW)
              endif 
              
              PrF = FSLik(l,Bi)
              PrW = PrW * PrF(:,x)              
              
              if (catG == i) then
                PrXYZ(x,y,z,2) = PrXYZ(x,y,z,2) * PrW(z)
              else if (catG /= 0) then
                if (.not. all(PrW==1D0))  PrXYZ(x,y,z,2) = PrXYZ(x,y,z,2) * SUM(PrW)
              else if (ParBisClone(i)) then
                PrXYZ(x,y,:,2) = PrXYZ(x,y,:,2) * PrW(x)
              else if (.not. all(PrW==1D0)) then
                PrXYZ(x,y,:,2) = PrXYZ(x,y,:,2) * SUM(PrW)
              endif
            enddo 
          enddo
        enddo
      endif
    enddo
    PrL(l) = LOG10(SUM(PrXYZ(:,:,:,2))) - LOG10(SUM(PrXYZ(:,:,:,1)))
  enddo
  
  LL = SUM(PrL) + Lind(A)
endif  
  
end subroutine AddGP

! #####################################################################

subroutine AddGGP(A, SB, k, LL)
use Global
use CalcLik             
implicit none
! A a GGP of sibship SB? (only calculating over non-gp-assigned)

integer, intent(IN) :: A, SB, k
double precision, intent(OUT) :: LL
integer :: m, n, catG, GG, AncG(2,mxA), DoQuick, l, x, y,z,i, v, Bj, w,Ei
double precision :: PrL(nSnp), PrXYZ(3,3,3,2), PrZ(3),PrA(3),PrP(3),PrV(3), &
  PrW(3), PrE(3), PrG(3), PrF(3,3)

LL = missing
if (GpID(1, SB,k)/=0) then
  if (GpID(2, SB,k)/=0) then  ! should be assigned as parent-of-gp
    LL = impossible   !(or AlreadyAss)
    return
  else
    m = 2
  endif
else
  m = 1  ! doesn't really matter (?); GpID(m, SB, k) == 0.
endif
GG = GpID(3-m, SB, k)

if (Sex(A)<3) then
  n = Sex(A)
else
  n = 1
endif

catG =0
AncG = 0
if (GG/=0) then
  if (GG==Parent(A,3-m)) then
    catG = 1
  endif
  call GetAncest(GG, 3-m, AncG)   
  if (ANY(AncG == A)) then
    if ((GG>0 .and. AncG(n,2)==A) .or. (GG<0 .and. &
     AncG(n,3-m+2)==A)) then  
      catG = 2  ! already GGP via 3-m; check if double ggp
    else
      LL = NotImplemented  ! possible; not yet implemented
    endif
  else
    do v=1,2
      if (Parent(A,v)==0) cycle
      if ((GG<0 .and. ANY(AncG(v, 2:4)==Parent(A,v))) .or. &
       (GG>0 .and. ANY(AncG(v, 1:2)==Parent(A,v)))) then
        LL = NotImplemented   ! TODO: stricter implementation?
      endif
    enddo
  endif
endif
if (GpID(3-k,SB,k) < 0) then
  do i=1,nS(SB,k)
    if (Parent(SibID(i, SB, k), 3-k) == GpID(3-k,SB,k)) then 
      LL = NotImplemented
      return
    endif
  enddo
endif
if (Parent(A,3-k)<0) then
  do i=1,nS(SB,k)
    if (Parent(SibID(i, SB, k), 3-k) == Parent(A,3-k)) then 
      LL = NotImplemented
      return
    endif
  enddo
endif
do i=1,nS(SB,k)
  if (Parent(SibID(i, SB, k), 3-k) == A) then   ! TODO: implement?
    LL = NotImplemented
    return
  else if (Parent(SibID(i,SB,k), 3-k) < 0) then
    if (any(getPar(Parent(SibID(i,SB,k), 3-k), 3-k) == A)) then
      LL = NotImplemented
      return 
    endif
  endif
enddo    

if (catG/=0) then  ! age check
  if (GG>0) then
    if (AgeDiff(GG, A) >= 0 .and. catG==1 .and. AgeDiff(GG, A)/=missing) &
     LL = impossible  ! A older than GG
    if (AgeDiff(GG, A) <= 0 .and. catG==2)  LL = impossible  ! GG older than A
  else if (GG<0) then
    do v=1, nS(-GG, 3-m)  ! TODO? Age, ancestors
      if (AgeDiff(SibID(v,-GG,3-m), A) >= 0 .and. catG==1 .and. &
        AgeDiff(SibID(v,-GG,3-m), A)/=missing)  LL = impossible
      if (AgeDiff(SibID(v,-GG,3-m), A) <= 0 .and. catG==2)  LL = impossible
    enddo
  endif
endif
if (LL == impossible) return

if (Complx==0 .and. Mate(A)/=0) then  
  if (Mate(A) == GPID(3-n, SB,k)) then
    LL = impossible
    return
  else
    catG = 3
  endif
endif


call ChkDoQuick(SB,k,DoQuick)
if (DoQuick == -1)  then
  LL = NotImplemented  ! inbreeding loops, approx. below invalid
  return
endif

PrL = 0D0
do l=1,nSnp
  call ParProb(l, A, 0, 0, 0, PrA)
  if (catG==1) then
    call ParProb(l, GG, 3-m, A, 0, PrZ)
    call ParProb(l, Parent(A,m), m, -1, 0, PrP)
    PrA = OcA(:, Genos(l,A))                            
  else if (catG==2) then
    call ParProb(l, GG, 3-m, -4, 0, PrZ)  ! offspring contribution only
    if (GG > 0) then
      call ParProb(l, Parent(GG,3-n), 3-n, GG, 0, PrP) 
    else if (GG < 0) then
      call ParProb(l, GpID(3-n, -GG,3-m), 3-n, 0, 0, PrP)
    else
      PrP = AHWE(:,l)
    endif
  else
    call ParProb(l, GG, 3-m, 0, 0, PrZ)
  endif
  if (catG==3) then
    call ParProb(l, Mate(A), 3-n, 0, 0, PrG)
  endif                                          
  do x=1,3  ! SB
    do y=1,3  ! in between SB and A
      do z=1,3  ! other GP 
        do v=1,3  ! A
          if (catG==1) then
            PrV(v) = AKAP(y, v, l) * SUM(AKA2P(v,z,:) * PrP)
          else if (catG==2) then
            PrV(v) = AKAP(y, v, l) * SUM(AKA2P(z,v,:) * PrP)
          else if (catG==3) then
            PrV(v) = SUM(AKA2P(y,v,:) * PrG)
          else
            PrV(v) = AKAP(y, v, l)
          endif
        enddo
        PrXYZ(x,y,z,:) = AKA2P(x, y, z) * PrZ(z) * SUM(PrV * PrA)
      enddo
    enddo
    if (DoQuick > 0) then   
      PrXYZ(x,:,:,2) = PrXYZ(x,:,:,2) * XPr(1,x,l, SB,k)
    else
      do i=1, nS(SB,k)
        Bj = SibID(i, SB, k)
        if (nFS(Bj)==0)  cycle
        call ParProb(l, Parent(Bj,3-k), 3-k, -1, 0, PrW)
        if (Parent(Bj,3-k)<0) then
          do w=1,3
            do v=1, nS(-Parent(Bj, 3-k), 3-k)
              Ei = SibID(v, -Parent(Bj, 3-k), 3-k)  
              if (NFS(Ei) == 0) cycle
              if (Parent(Ei, k) == -SB) cycle
              call ParProb(l, Parent(Ei, k), k, Ei,-1, PrE)       
              PrF = FSLik(l,Ei)
              PrE = PrE * PrF(:,w)
              if (.not. ALL(PrE==1D0))  PrW(w) = PrW(w) * SUM(PrE)  
            enddo  
          enddo
          if (.not. all(PrW==1D0))  PrXYZ(x,:,:,1) = PrXYZ(x,:,:,1) * SUM(PrW)
        endif          
        PrF = FSLik(l,Bj)
        PrW = PrW * PrF(:,x)                                      
        if (.not. all(PrW==1D0))  PrXYZ(x,:,:,2) = PrXYZ(x,:,:,2) * SUM(PrW)
      enddo 
    endif
  enddo
  PrL(l) = LOG10(SUM(PrXYZ(:,:,:,2))) - LOG10(SUM(PrXYZ(:,:,:,1)))          
enddo
LL = SUM(PrL) + Lind(A)

end subroutine AddGGP

! #####################################################################

subroutine addGAU(A, SB, k, m, LL)  ! A great-full-avuncular of B's
use Global
use CalcLik             
implicit none

integer, intent(IN) :: A, SB, k, m
double precision, intent(OUT) :: LL
integer :: nOff, Offspr(maxSibSize), sxOff(maxSibSize), DoQuick, &
  l, x, y, z, v,g, i, Bi, Ei, w,j, AncG(2,mxA), Ax
double precision :: PrL(nSnp), PrGGG(3,2), PrXV(3,3,3,3,2), PrG(3), PrW(3), PrE(3), PrF(3,3)

LL = missing
if (GpID(m,SB,k)/=0) then
  LL = NotImplemented
else if (Parent(A,3-m) == GpID(3-m,SB,k) .and. Parent(A,3-m)/=0) then
  LL = AlreadyAss  
else 
  call getOff(-SB,k, .TRUE., nOff, Offspr, sxOff)
  do g=1,2
    if (ANY(Offspr(1:nOff)==Parent(A,g))) then
      LL = impossible  
    else if (g/=k .and. Parent(A,g)/=0 .and. & 
     ANY(Parent(SibID(1:ns(SB,k),SB,k),g)==Parent(A,g))) then
      LL = NotImplemented
    endif
  enddo
endif
if (LL /= missing) return

AncG = 0
call getAncest(GpID(3-m,SB,k),3-m, AncG)
if (ANY(AncG == A)) then
  LL = NotImplemented
  return
endif

call ChkDoQuick(SB,k,DoQuick)
if (DoQuick == -1)  then
  LL = NotImplemented  ! inbreeding loops, approx. below invalid
  return
endif

Ax = 0
!if (any(parent(A,:)==0)) Ax = A

PrL = 0D0
do l=1,nSnp
  do g=1,2
    call ParProb(l, Parent(A,g), g, Ax, 0, PrGGG(:,g)) 
  enddo
  call ParProb(l, GpID(3-m,SB,k),3-m,0,0,PrG)
  do x=1,3  ! sibship parent
    do y=1,3
      do z=1,3
        do v=1,3
          PrXV(x,y,z,v,:) = SUM(AKA2P(x,y,:) *PrG) *AKA2P(y,z,v) *PrGGG(z,1) *PrGGG(v,2) 
          !if (any(parent(A,:)==0))  PrXV(x,y,z,v,2) = PrXV(x,y,z,v,2) * OKA2P(Genos(l,A), z, v)
          if (DoQuick > 0) then
            PrXV(x,y,z,v,2) = PrXV(x,y,z,v,2) * XPr(1,x,l, SB,k)
          else
            do i=1, nS(SB,k)
              Bi = SibID(i, SB, k)
              if (nFS(Bi)==0)  cycle
              call ParProb(l, Parent(Bi,3-k), 3-k,-1,0, PrW)
              if (Parent(Bi,3-k)<0) then
                do w=1,3
                  do j=1, nS(-Parent(Bi, 3-k), 3-k)
                    Ei = SibID(j, -Parent(Bi, 3-k), 3-k)  
                    if (NFS(Ei) == 0) cycle
                    if (Parent(Ei, k) == -SB) cycle
                    call ParProb(l, Parent(Ei, k), k, Ei,-1, PrE)       
                    PrF = FSLik(l,Ei)
                    PrE = PrE * PrF(:,w)
                    if (.not. ALL(PrE==1D0))  PrW(w) = PrW(w) * SUM(PrE)  
                  enddo  
                enddo
                PrXV(x,y,z,v,1) = PrXV(x,y,z,v,1) * SUM(PrW)
              endif
              PrF = FSLik(l,Bi)
              PrW = PrW * PrF(:,x)             
              PrXV(x,y,z,v,2) = PrXV(x,y,z,v,2) * SUM(PrW)
            enddo
          endif
        enddo
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXV(:,:,:,:,2))) -  LOG10(SUM(PrXV(:,:,:,:,1)))          
enddo
LL = SUM(PrL)

LL = LL + Lind(A)    ! if (.not. any(parent(A,:)==0))  


end subroutine addGAU

! #####################################################################

subroutine ParentHFS(A, SA, kA, SB, kB, hf, LL)  
! parents of SA and SB HS/FS?
use Global
use CalcLik             
implicit none

integer, intent(IN) :: A, SA, kA, SB, kB, hf
double precision, intent(OUT) :: LL
integer :: PA, nA, AncA(2,mxA), AA(maxSibSize), G(2), AncB(2,mxA), GA(2), GB(2), &
  catA(maxSibSize), catB(nS(SB,kB)), catG, GGP(2), DoQuickA, DoQuickB, &
  m, l, x, y, u,v, i, j,z, r, Ei, e, DoneA(MaxSibSize), Bj
double precision :: LLm(2), ALR, PrL(nSnp), PrG(3,2), PrXV(3,3,3,3,3,2), PrPA(3, 2),&
 PrGA(3), PrGB(3), PrE(3), PrH(3), PrGG(3), PrF(3,3)

PA = 0
nA = 0
LLm = missing  
if (A/=0) then
  PA = Parent(A,kA)
  call GetAncest(A, kA, AncA)
else
  PA = -SA
  call GetAncest(-SA, kA, AncA)
endif
if (PA/=0 .and. GpID(kA,SB,kB)==PA) then
  LL = impossible
  return
endif

do m=1,2
  if (m/=hf .and. hf/=3) cycle
  if (PA < 0) then
    if (GpID(m,-PA,kA)/=GpID(m,SB,kB) .and. GpID(m,-PA,kA)/=0 .and. GpID(m,SB,kB)/=0) then
      LLm(m) = impossible
    endif
    nA = nS(-PA,kA)
    AA(1:nA) = SibID(1:nA, -PA, kA)
  else 
    nA = 1
    AA(1) = A
    if (PA > 0) then
      if (Parent(PA,m)/=GpID(m,SB,kB) .and. parent(PA,m)/=0 .and. GpID(m,SB,kB)/=0) then
        LLm(m) = impossible
      endif
    endif
  endif
enddo
if (ALL(LLm == impossible)) then
  LL = impossible
  return
endif

G = 0
call GetAncest(-SB, kB, AncB)  
GA = AncA(:,kA+2)
GB = AncB(:,kB+2)
do m=1,2
  if (m/=hf .and. hf/=3) cycle
  if (GA(m)/=0) then
    if (GA(m) == -SB .and. m==kB) then
      LLm(m) = impossible
    else if (GB(m)/=0 .and. GA(m)/=GB(m)) then
      LLm(m) = impossible
    else if (GB(m)==0) then
      G(m) = GA(m)
    else if (GB(m)==GA(m)) then
      G(m) = GA(m)
      LLm(m) = AlreadyAss  ! already are sibs
    else
      LLm(m) = impossible
    endif
  else 
    G(m) = GB(m)
  endif
enddo
if (ALL(LLm == impossible .or. LLm==AlreadyAss)) then
  LL = impossible
  return
endif     

if (hf==3) then  ! FS
  if (ANY(AncA(kB, 3:mxA) == -SB)) then
    LLm = impossible
  else if (A>0) then
    if (ANY(AncB == A)) then
      LLm = impossible
    endif
  else if (SA/=0) then
    if (ANY(AncB(kA,3:mxA) == -SA)) then
      LLm = impossible
    endif
  endif
endif 
if (ALL(LLm == impossible)) then
  LL = impossible
  return
endif

ALR = missing
if (hf==3) then
  if (A>0)   call CalcAgeLR(  A, kA, -SB, kB,kA, 5, .TRUE., ALR)
  if (SA/=0) call CalcAgeLR(-SA, kA, -SB, kB, 0, 2, .TRUE., ALR) 
  if (ALR == impossible)  LLm = impossible    
else
  do m=1,2
    if (A>0)   call CalcAgeLR(  A, kA, -SB, kB, m, 6, .TRUE., ALR)
    if (SA/=0) call CalcAgeLR(-SA, kA, -SB, kB, m, 3, .TRUE., ALR) 
    if (ALR == impossible)  LLm(m) = impossible    
  enddo
endif
if (ALL(LLm == impossible)) then
  LL = impossible
  return
endif

if (hf==3) then
  if (LLm(1)==impossible .or. LLm(2)==impossible) then
    LL = impossible
  else if (LLm(1)==AlreadyAss .and. LLm(2)==AlreadyAss) then
    LL = AlreadyAss  ! already are FS
  endif
else 
  if (LLm(hf)==impossible) then 
    LL = impossible
  else if (LLm(hf)==AlreadyAss) then
    LL = AlreadyAss
  else if (GB(3-hf)==GA(3-hf)) then
    LL = impossible   ! already HS, would become FS
  endif
endif

if (ANY(AncA(kB, 5:mxA)==-SB)) then
  LL = NotImplemented  ! highly unlikely (but not strictly impossible: TODO)
else if (AncA(kA,2)/=0 .and. ANY(AncB(kA, 5:mxA) == AncA(kA,2))) then 
  LL = impossible
endif
if (LL /=missing) return  
   
catA = 0  
catB = 0
do i=1, nA
  if (hermaphrodites/=0 .and. PA<0) then
    if (Parent(AA(i), 3-kA) == -DumClone(-PA,kA) .and. DumClone(-PA,kA)/=0) then
      catA(i) = 12
      cycle
    endif
  endif 
  if (kA/=kB) then
    if (Parent(AA(i), kB) == AncB(kB, 2) .and. AncB(kB, 2)<0) then
      catA(i) = 1
    endif
  else if (kA == kB .and. Parent(AA(i), 3-kA) /= 0) then  
    do j=1, nS(SB, kB)
      if (Parent(AA(i), 3-kA) == Parent(SibID(j,SB,kB), 3-kB)) then
        catA(i) = 2
        catB(j) = 2
      endif
    enddo
  endif
  if (Parent(AA(i), 3-kA) /= 0) then
    if (G(3-kA) == Parent(AA(i), 3-kA)) then  ! incl. hf==3
      if (kA==kB) then
        catA(i) = 3  ! (u) 3-kA = 3-kB == hf 
      else if (kA/=kB) then
        catA(i) = 4  ! (z)
      endif 
    else if (kA==hf .and. GA(3-kA) == Parent(AA(i), 3-kA)) then
      catA(i) = 4  ! (z)
    else if (kA==hf .and. GB(3-kA) == Parent(AA(i), 3-kA)) then
      catA(i) = 5  ! (v)
    endif
  endif
enddo    

if (Complx<2 .and. (any(catA/=0) .or. any(catB/=0))) then   ! TODO DOUBLE CHECK IF SOME VALID
  LL = NotImplemented
  return
endif

do i=1, nS(SB, kB)
  if (hermaphrodites/=0 .and. DumClone(SB,kB)/=0) then
    if (Parent(SibID(i,SB,kB), 3-kB) == -DumClone(SB,kB)) then
      catB(i) = 12
      cycle
    endif
  endif
  if (kA/=kB) then
    if (Parent(SibID(i,SB,kB), kA) ==AncA(kA,2) .and. AncA(kA,2)<0) then
      catB(i) = 1
    endif
  endif
  if (Parent(SibID(i,SB,kB), 3-kB) /= 0) then
    if (G(3-kB) == Parent(SibID(i,SB,kB), 3-kB)) then
      catB(i) = 3  ! (u)  (for hf<3 .and. hf==3)
    else if (kB==hf .and. GA(3-kB) == Parent(SibID(i,SB,kB), 3-kB)) then
      catB(i) = 4  ! (z) (GA of type 3-kB if hf==kB) 
    else if (kB==hf .and. GB(3-kB) == Parent(SibID(i,SB,kB), 3-kB)) then
      catB(i) = 5  ! (v)
    endif
  endif
enddo 

catG = 0
GGP = 0
if (hf<3) then
  GGP = getPar(G(hf), hf)
  if (GGP(3-hf) == GA(3-hf) .and. GA(3-hf)/=0) then
    catG = 1
  else if (GGP(3-hf) == GB(3-hf) .and. GB(3-hf)/=0) then
    catG = 2
  endif
  if (catG == 0)  GGP = 0
endif

DoQuickA = 1
DoQuickB = 1
if (SA/=0)  call ChkDoQuick(SA,kA,DoQuickA)
call ChkDoQuick(SB,kB,DoQuickB)

PrL = 0D0
do l=1,nSnp
  do m=1,2
    if (m/=hf .and. hf/=3) cycle
    if (ANY(CatA==3) .or. ANY(CatB==3)) then
      call ParProb(l, G(m), m, -1,0, PrG(:,m)) 
    else if (catG/=0) then
      call ParProb(l, G(m), m, -4, 0, PrG(:,m))
      if (G(m) > 0) then
        call ParProb(l, GGP(hf), hf, G(m), 0, PrGG) 
      else
        call ParProb(l, GGP(hf), hf, 0, 0, PrGG)
      endif
    else
      call ParProb(l, G(m), m, 0,0, PrG(:,m)) 
    endif
  enddo
  if (hf < 3) then
    if (ANY(CatA==4) .or. ANY(CatB==4)) then
      call ParProb(l, GA(3-hf), 3-hf, -1,0, PrGA)
    else
      call ParProb(l, GA(3-hf), 3-hf, 0,0, PrGA)
    endif
    if (ANY(CatA==5) .or. ANY(CatB==5)) then
      call ParProb(l, GB(3-hf), 3-hf, -1,0, PrGB)
    else
      call ParProb(l, GB(3-hf), 3-hf, 0,0, PrGB)
    endif
  endif
  if (A>0) then
    call ParProb(l, Parent(A,kA), kA, A,-4, PrPA(:,kA))
    call ParProb(l, Parent(A,3-kA), 3-kA, A,0, PrPA(:,3-kA))
  endif
  
  PrXV = 0D0
  do x=1,3  ! SA/PA
    do y=1,3  ! SB
      do u=1,3  ! G_hf / G_3-kB (hf==3)
        do z=1,3  ! G_A (hf/=3) / G_kB (hf==3)
          do v=1,3 ! G_B (hf/=3)
            if (hf==3) then  ! 0 for z/=v
              PrXV(x,y,u,z,z,:) = AKA2P(x,u,z) * AKA2P(y,u,z) *&
               PrG(u,3-kB) * PrG(z, kB)
            else
              if (GA(3-hf) < 0 .and. GA(3-hf) == -SB .and. hf==3-kB) then
                PrXV(x,y,u,y,v,:) = AKA2P(x,u,y) * AKA2P(y,u,v) *&
                 PrG(u,hf) * PrGB(v)
              else if (GB(3-hf) < 0 .and. GB(3-hf) == -SA .and. hf==3-kA) then
                PrXV(x,y,u,z,x,:) = AKA2P(x,u,z) * AKA2P(y,u,x) *&
                 PrG(u,hf) * PrGA(z)
              else if (catG == 1) then
                PrXV(x,y,u,z,v,:) =AKA2P(x,u,z)*AKA2P(y,u,v)*PrG(u,hf)*&
                  SUM(AKA2P(u,z,:) * PrGG) * PrGA(z) * PrGB(v)
              else if (catG == 2) then
                PrXV(x,y,u,z,v,:) =AKA2P(x,u,z)*AKA2P(y,u,v)*PrG(u,hf)*&
                  SUM(AKA2P(u,v,:) * PrGG) * PrGA(z) * PrGB(v)
              else
                PrXV(x,y,u,z,v,:) = AKA2P(x,u,z)*AKA2P(y,u,v)*PrG(u,hf)&
                * PrGA(z) * PrGB(v)
              endif
            endif
            if (A /=0) then
              if (Parent(A, kA)/=0) then
                PrXV(x,y,u,z,v,:) = PrXV(x,y,u,z,v,:) * PrPA(x, kA)
              endif
            endif
          enddo
        enddo
      enddo
      
      DoneA = 0            
      if ((DoQuickB > 0 .and. ALL(catA==0) .and. ALL(catB==0)) .or. &
        DoQuickB>1 .or. DoQuickA>1) then   
        if (SA/=0) then
          PrXV(x,y,:,:,:,2) = PrXV(x,y,:,:,:,2) * XPr(1,x,l, SA,kA) *&
           XPr(1,y,l, SB,kB)
        else if (A>0) then
         PrXV(x,y,:,:,:,2) =PrXV(x,y,:,:,:,2)*SUM(OKA2P(Genos(l,A),x,:)&
          * PrPA(:,3-kA)) * XPr(1,y,l, SB,kB)
        endif            

      else
        do r=1, nS(SB,kB)
          Bj = SibID(r, SB, kB) 
          if (NFS(Bj) == 0) cycle 
          if (catB(r)==0 .or. catB(r)==2) then
            call ParProb(l, Parent(Bj,3-kB), 3-kB, -1, 0, PrE)
          else
            PrE = 1D0
          endif                                           

          if (Parent(Bj,3-kB) <0 .and. CatB(r)/=1) then
            do e=1,3
              if (catB(r)==12 .and. e/=y)  cycle                                  
              do v=1, nS(-Parent(Bj,3-kB), 3-kB)
                Ei = SibID(v, -Parent(Bj,3-kB), 3-kB)
                if (nFS(Ei) == 0) cycle
                if (Parent(Ei, kB) == -SB) cycle
                if (Parent(Ei, kA) == Parent(AA(1),kA) .and. Parent(AA(1),kA)/=0) cycle  ! FS of A if A>0
                call ParProb(l, Parent(Ei, kB), kB, Ei, -1, PrH)
                do i=1,nFS(Ei)
                  if (A>0 .and. (FSID(i,Ei)==A .or. FSID(i,Ei)==Parent(AA(1),kA))) cycle
                  PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)),:,e)
                enddo
                if (.not. all(PrH==1D0))  PrE(e) = PrE(e) * SUM(PrH)
              enddo
            enddo
          endif

          if ((catB(r)==0 .or. catB(r)==2) .and. .not. all(PrE==1D0)) then
            PrXV(x,y,:,:,:,1) = PrXV(x,y,:,:,:,1) * SUM(PrE)
          else if (catB(r)==1) then  ! Parent(Bj,3-kB) = PA, kA/=kB
            PrXV(x,y,:,:,:,1) = PrXV(x,y,:,:,:,1) * PrE(x)
          else if (catB(r)==3) then  ! hf==kB, Parent(Bj,3-kB) = GA
            do u=1,3
              PrXV(x,y,u,:,:,1) = PrXV(x,y,u,:,:,1) * PrE(u)
            enddo
          else if (catB(r)==4) then  ! Parent(Bj,3-kB) = G
            do e=1,3
              PrXV(x,y,:,e,:,1) = PrXV(x,y,:,e,:,1) * PrE(e)
            enddo
          else if (catB(r)==5) then  ! hf==kB, Parent(Bj,3-kB) = GB
            do v=1,3
              PrXV(x,y,:,:,v,1) = PrXV(x,y,:,:,v,1) * PrE(v)
            enddo
          else if (catB(r)==12) then ! selfing
            PrXV(x,y,:,:,:,1) = PrXV(x,y,:,:,:,1) * PrE(y)
          endif
     
          PrF = FSLik(l,Bj)
          PrE =  PrE * PrF(:,y)

          if (any(catA == 2) .and. Parent(Bj,3-kB)/=0) then  ! kA==kB, share parent 3-kB
            do v = 1, nA
              if (A/=0 .and. AA(v)/=A) cycle
              if (Parent(AA(v), 3-kA)/=Parent(Bj,3-kB)) cycle
              PrE =  PrE * OKA2P(Genos(l,AA(v)), x, :)
              doneA(v) = 1
            enddo
          endif
          
          if ((catB(r)==0 .or. catB(r)==2) .and. .not. all(PrE==1D0)) then
            PrXV(x,y,:,:,:,2) = PrXV(x,y,:,:,:,2) * SUM(PrE)
          else if (catB(r)==1) then  ! Parent(Bj,3-kB) = PA, kA/=kB
            PrXV(x,y,:,:,:,2) = PrXV(x,y,:,:,:,2) * PrE(x)
          else if (catB(r)==3) then  ! hf==kB, Parent(Bj,3-kB) = GA
            do u=1,3
              PrXV(x,y,u,:,:,2) = PrXV(x,y,u,:,:,2) * PrE(u)
            enddo
          else if (catB(r)==4) then  ! Parent(Bj,3-kB) = G
            do e=1,3
              PrXV(x,y,:,e,:,2) = PrXV(x,y,:,e,:,2) * PrE(e)
            enddo
          else if (catB(r)==5) then  ! hf==kB, Parent(Bj,3-kB) = GB
            do v=1,3
              PrXV(x,y,:,:,v,2) = PrXV(x,y,:,:,v,2) * PrE(v)
            enddo
          else if (catB(r)==12) then ! selfing
            PrXV(x,y,:,:,:,2) = PrXV(x,y,:,:,:,2) * PrE(y)
          endif
        enddo
        
        do r = 1, nA
          if (DoneA(r)==1) cycle
          if (SA/=0 .and. NFS(AA(r)) == 0) cycle         
          if (kA/=kB .and. Parent(AA(r),3-kA)==-SB) cycle  ! done
          if (catA(r)==0) then
            call ParProb(l, Parent(AA(r),3-kA), 3-kA, -1, 0, PrE)
          else
            PrE = 1D0
          endif

          if (Parent(AA(r), 3-kA) < 0 .and. (SA/=0 .or. &
           ANY(FSID(1:nFS(AA(r)), AA(r))==A))) then
            do e=1,3
              if (catA(r)==12 .and. e/=x)  cycle                                  
              do i=1, nS(-Parent(AA(r), 3-kA), 3-kA)
                Ei = SibID(i, -Parent(AA(r), 3-kA), 3-kA)
                if (nFS(Ei) == 0) cycle
                if (Parent(Ei, kB) == -SB) cycle
                if (A>0 .and. Ei==A) cycle
                if (SA/=0 .and. Parent(Ei, kA) == -SA) cycle
                call ParProb(l, Parent(Ei, kA), kA, Ei, -1, PrH)  
                do j=1, nFS(Ei)
                  if (A/=0 .and. FSID(j, Ei)==A) cycle
                  PrH = PrH * OKA2P(Genos(l,FSID(j,Ei)), :, e)
                enddo
                if (.not. all(PrH==1D0))  PrE(e) = PrE(e) * SUM(PrH)
              enddo
            enddo
          endif
          
          if (catA(r)<3 .and. .not. all(PrE==1D0)) then
            PrXV(x,y,:,:,:,1) = PrXV(x,y,:,:,:,1) * SUM(PrE)
          else if (catA(r)==3) then 
            do u=1,3
              PrXV(x,y,u,:,:,1) = PrXV(x,y,u,:,:,1) * PrE(u)
            enddo
          else if (catA(r)==4) then 
            do z=1,3
              PrXV(x,y,:,z,:,1) = PrXV(x,y,:,z,:,1) * PrE(z)
            enddo
          else if (catA(r)==5) then 
            do v=1,3
              PrXV(x,y,:,:,v,1) = PrXV(x,y,:,:,v,1) * PrE(v)
            enddo
          else if (catA(r)==12) then   
            PrXV(x,y,:,:,:,1) = PrXV(x,y,:,:,:,1) * PrE(x)
          endif
          
          do i=1, MAX(1, nFS(AA(r)))
            if (SA/=0 .or. (FSID(i, AA(r))==A .and. catA(r)/=2)) then 
              PrE =  PrE * OKA2P(Genos(l,FSID(i,AA(r))), x, :)
              doneA(r) = 2
            else if (catA(r)/=2) then
              PrE =  PrE * OKA2P(Genos(l,FSID(i,AA(r))), x, :)
            endif
          enddo
          
          if (catA(r)<3 .and. .not. all(PrE==1D0)) then
            PrXV(x,y,:,:,:,2) = PrXV(x,y,:,:,:,2) * SUM(PrE)
           else if (catA(r)==3) then 
            do u=1,3
              PrXV(x,y,u,:,:,2) = PrXV(x,y,u,:,:,2) * PrE(u)
            enddo
          else if (catA(r)==4) then 
            do z=1,3
              PrXV(x,y,:,z,:,2) = PrXV(x,y,:,z,:,2) * PrE(z)
            enddo
          else if (catA(r)==5) then 
            do v=1,3
              PrXV(x,y,:,:,v,2) = PrXV(x,y,:,:,v,2) * PrE(v)
            enddo
          else if (catA(r)==12) then   
            PrXV(x,y,:,:,:,2) = PrXV(x,y,:,:,:,2) * PrE(x)
          endif
        enddo
      endif
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXV(:,:,:,:,:,2))) - LOG10(SUM(PrXV(:,:,:,:,:,1)))
enddo

LL = SUM(PrL)
  
end subroutine ParentHFS

! #####################################################################

subroutine DummyGP(SA, SB, kA, kB, focal, LL)  
! SB GP of SA? via observed or unobserved
use Global
use CalcLik             
implicit none

integer, intent(IN) :: SA, SB, kA, kB, focal
double precision, intent(OUT) :: LL
integer :: i, r, AncB(2,mxA), Bi, catA, catB, GA(2), catG, &
  DoQuickA, DoQuickB, GGP(2), m, l, x, y, z, v, Ai, ix, e,f, Ei, parGA(2)
double precision :: LLGX(2,2), LLtmp(maxSibSize, 2), PrL(nSnp), PrZ(3), PrV(3), &
  PrG(3), PrPG(3), PrXYZ(3,3,3,3,2), PrW(3), dx(maxSibSize), PrH(3), PrF(3,3), LLZ(2,2)
logical :: MaybeGP(ns(SB,kB),2), ParAisClone(ns(SA,kA)), ParBisClone(ns(SB,kB)), &
  GA_SB_connected(2), ParOK

LL = missing
do i=1, nS(SB,kB)
  Bi = SibID(i, sB, kB)
  if (any(GpID(:,SA,kA) == Bi)) then
    LL = AlreadyAss
  else if (kA /= kB) then
    if (Parent(Bi, kA) == -SA) then
      LL = NotImplemented
    endif
  else if (kA == kB) then
    do r= 1, nS(SA, kA)
      if (Parent(Bi, 3-kB)==Parent(SibID(r,SA,kA), 3-kA) &
       .and. Parent(Bi, 3-kB)/=0) then
        LL = NotImplemented   ! TODO
      endif
    enddo
  endif
enddo
if (LL == NotImplemented)  return 

call GetAncest(-SB, kB, AncB)
if (ANY(AncB(kA,3:mxA) == -SA)) then
  LL = impossible 
  return
else if (ANY(AncB(3-kA, 3:mxA) < 0)) then
  do r= 1, nS(SA, kA)
    if (ANY(AncB(3-kA, 3:mxA)/=0 .and. AncB(3-kA, 3:mxA) == &
     Parent(SibID(r,SA,kA), 3-kA))) then
      LL = NotImplemented 
      return
    endif
  enddo
endif

catA = 0
catB = 0
GA = GpID(:, SA, kA)
do r = 1, nS(sB,kB)   ! check if overlap
  Bi = SibID(r, sB, kB)
  if (NFS(Bi) == 0) cycle
  if (Parent(Bi, 3-kB) == GA(3-kB) .and. GA(3-kB)<0) then
    catB = Bi
  endif
enddo
do r = 1, nS(sA,kA)   ! check if inbreeding loop
  Ai = SibID(r,SA,kA)
  if (NFS(Ai)==0) cycle
  if (Parent(Ai, 3-kA) == GA(3-kA) .and. GA(3-kA)<0) then 
    catA = Ai
  endif
enddo    

catG = 0
GA_SB_connected = .FALSE.
do m=1,2
  if (GA(m) == GpID(m, SB, kB) .and. GA(m)/=0) then
    catG = m
  endif
  call connected(GA(m),m, -SB,kB, GA_SB_connected(m))
  if (GA(m) < 0) then
    if (any(GpID(:,-GA(m),m) == GpID(:,SB,kB) .and. GpID(:,SB,kB)/=0)) then
      GA_SB_connected(m) = .TRUE.
    endif
  endif  
enddo

call ChkDoQuick(SA,kA,DoQuickA)
call ChkDoQuick(SB,kB,DoQuickB) 
 
ParAisClone = .FALSE.                             
ParBisClone = .FALSE.
if (hermaphrodites/=0) then
  if (DumClone(SA,kA)/=0) then   ! DoQuickA==-3
    do i=1, nS(SA,kA)
      if (Parent(SibID(i,SA,kA), 3-kA) == -DumClone(SA,kA)) then
        ParAisClone(i) = .TRUE.
      endif
    enddo
  endif
  if (DumClone(SB,kB)/=0) then   ! DoQuickB==-3
    do i=1, nS(SB,kB)
      if (Parent(SibID(i,SB,kB), 3-kB) == -DumClone(SB,kB)) then
        ParBisClone(i) = .TRUE.
      endif
    enddo
  endif
endif

MaybeGP = .FALSE.   ! Bi potential GP of SA?
do m=1,2
  if (GA(m)/=0) cycle
  do r=1, ns(SB,kB)
    Bi = SibID(r, sB, kB)
    if (sex(Bi)/=m .and. sex(Bi)<3)  cycle
    call chkvalidpar(-SA,kA, Bi,m, ParOK)
    if (ParOK)  MaybeGP(r,m) = .TRUE.
  enddo
enddo

LLGX = missing  ! D2: 1: a B is GP; 2: an unsampled offspring of SB is GP
LLtmp = missing
LLZ = missing
GGP = 0
parGA = 0
PrPG = 0d0          
do m=1,2
  if (focal==4 .and. all(GA==0) .and. complx/=0) then
    LLGX(m,1) = NotImplemented  
    ! else: SB is GP of SA --> Bi is parent of SA --> mate-of-SB is also GP of SA
  else if (m==kB .and. GA(kB) == -SB) then
    LLGX(m,:) = impossible
    cycle
  else if (GA(m) == 0) then
    do r=1, nS(sB,kB)  ! TODO: also loop over dummy offspring?
      if (.not. MaybeGP(r,m))  cycle
      Bi = SibID(r, sB, kB) 
      call AddGP(Bi, SA, kA, LLtmp(r,m))
      if (LLtmp(r,m) < 0)  LLtmp(r,m) = LLtmp(r,m) - Lind(Bi) + CLL(SB,kB)
    enddo
    LLGX(m,1) = MaxLL(LLtmp(:,m))
  else 
    parGA = getPar(GA(m),m)
    if (parGA(kB) /= 0) then
      LLGX(m,:) = impossible
      cycle
    endif
    call ChkValidPar(GA(m),m, -SB,kB, ParOK)
    if (.not. ParOK) then
      LLGX(m,:) = impossible
      cycle
    endif  
    if (GA(m) > 0) then
      call AddSib(GA(m), SB, kB, LLZ(1,m))
      call AddFS(GA(m), SB, kB,0,m, LLZ(2,m), ix, dx)
      LLGX(m,1) = MaxLL(LLZ(:,m))
      if (LLGX(m,1) < 0D0)  LLGX(m,1) = LLGX(m,1) - Lind(GA(m)) + CLL(SA, kA)
    else
      call PairUA(GA(m),-SB,m,kB,LLGX(m,1))
      if (LLGX(m,1) < 0d0) then
        call CalcU(GA(m),m, -SB,kB,LLZ(1,m))
        call CalcU(-SA,kA,  -SB,kB,LLZ(2,m))
        LLGX(m,1) = LLGX(m,1) - LLZ(1,m) + LLZ(2,m)
      endif
    endif
    LLGX(m,2) = impossible
    cycle
  endif
    
  if (any(GA_SB_connected)) then
    LLGX(m,2) = NotImplemented  
    cycle 
  endif
  
  if (Complx==0 .and. DumMate(SB,kB)/=0) then
    GGP(m) = DumMate(SB,kB)
  else
    GGP(m) = parGA(3-kB)
  endif

  PrL = 0D0
  do l=1,nSnp
    if ((catB /= 0 .and. m==kB) .or. catA/=0) then
      call ParProb(l, GA(3-m), 3-m, -1, 0, PrZ)
    else
      call ParProb(l, GA(3-m), 3-m, 0, 0, PrZ)
    endif
    if (catB/=0 .and. m==3-kB .and. GA(m)<0) then
      PrG = 1D0
    else
      call ParProb(l, GA(m), m, -4, 0, PrG)  ! GA(m)'s offspring if<0, else 1D0
    endif
    if (GGP(m) /= DumMate(SB,kB) .or. GGP(m)==0) then
      call ParProb(l, GGP(m), 3-kB, 0, 0, PrPG)
    endif

    PrXYZ = 0D0
    do x=1,3  ! SA
      do y=1,3  ! parent of SA, offspr of SB
        do z=1,3  ! other parent of SA
          PrXYZ(x,y,z,:,:) = AKA2P(x, y, z) * PrG(y) * PrZ(z)
        enddo

        if ((catA==0 .and. DoQuickA>0) .or. DoQuickA>1) then
          PrXYZ(x,y,:,:,2) = PrXYZ(x,y,:,:,2) * XPr(1,x,l, SA,kA)   
        else
          do z=1,3  ! mate of SA
            do r=1, nS(sA,kA)
              Ai = SibID(r, sA, kA)  
              if (NFS(Ai) == 0) cycle
              if (Ai == GA(m)) cycle
              if (catA==Ai  .or. ParAisClone(r)) then
                PrW = 1D0
              else
                call ParProb(l, Parent(Ai, 3-kA), 3-kA, -1, 0, PrW)
              endif
              
              if (Parent(Ai,3-kA) < 0) then
                do e=1,3
                  do f=1, nS(-Parent(Ai, 3-kA), 3-kA)
                    Ei = SibID(f, -Parent(Ai, 3-kA),3-kA)
                    if (nFS(Ei) == 0) cycle
                    if (Parent(Ei,kA) == -SA)  cycle
                    call ParProb(l,Parent(Ei,kA),kA,Ei,-1,PrH) 
                    PrF = FSLik(l,Ei)
                    PrH = PrH * PrF(:,x)          
                    if (.not. all(PrH==1D0))  PrW(e) = PrW(e) * SUM(PrH)
                  enddo
                enddo
              endif
              
              if (catA==Ai) then
                if (m==kA) then
                  PrXYZ(x,y,z,:,1) = PrXYZ(x,y,z,:,1) * PrW(z)
                else if (m/=kA) then
                  PrXYZ(x,y,z,:,1) = PrXYZ(x,y,z,:,1) * PrW(y)
                endif
              else if (ParAisClone(r)) then
                PrXYZ(x,y,z,:,1) = PrXYZ(x,y,z,:,1) * PrW(x)
              else if (.not. all(PrW==1D0)) then
                PrXYZ(x,y,z,:,1) = PrXYZ(x,y,z,:,1) * SUM(PrW)
              endif
              
              PrF = FSLik(l,Ai)
              PrW = PrW * PrF(:,x)
              
              if (catA==Ai) then
                if (m==kA) then
                  PrXYZ(x,y,z,:,2) = PrXYZ(x,y,z,:,2) * PrW(z)
                else if (m/=kA) then
                  PrXYZ(x,y,z,:,2) = PrXYZ(x,y,z,:,2) * PrW(y)
                endif
              else if (ParAisClone(r)) then
                PrXYZ(x,y,z,:,2) = PrXYZ(x,y,z,:,2) * PrW(x)
              else if (.not. all(PrW==1D0)) then
                PrXYZ(x,y,z,:,2) = PrXYZ(x,y,z,:,2) * SUM(PrW)
              endif
            enddo  ! r
          enddo  ! z
        endif
      
        do v=1,3  ! SB
          if (DoQuickB==-2 .or. Complx==0) then  ! all FS
            do r=1, nS(sB,kB)
              Bi = SibID(r, sB, kB)  
              if (NFS(Bi) == 0) cycle                            
              call ParProb(l, Parent(Bi, 3-kB), 3-kB, -1, 0, PrW)
              call ParProb(l, -SB, kB, -1, 0, PrV)  ! GPs of SB
              PrXYZ(x,y,:,v,:) = PrXYZ(x,y,:,v,:) * PrV(v)
              PrF = FSLik(l,Bi)
              
              if (Parent(Bi,3-kB)== GGP(m) .and. GGP(m)/=0) then  ! Complx==0 .or. 
                PrXYZ(x,y,:,v,1) = PrXYZ(x,y,:,v,1) * SUM(AKA2P(y,v,:) * PrW)
                PrXYZ(x,y,:,v,2) = PrXYZ(x,y,:,v,2) * SUM(AKA2P(y,v,:) * PrW * PrF(:,v))
              else 
                PrXYZ(x,y,:,v,:) = PrXYZ(x,y,:,v,:) * SUM(AKA2P(y,v,:) * PrPG)
                PrXYZ(x,y,:,v,2) = PrXYZ(x,y,:,v,2) * SUM(PrW * PrF(:,v))
              endif
            enddo

         else if (catG==0  .and. DoQuickB/=-3) then
           PrXYZ(x,y,:,v,:) = PrXYZ(x,y,:,v,:) * SUM(AKA2P(y, v, :) * PrPG)
           PrXYZ(x,y,:,v,1) = PrXYZ(x,y,:,v,1) * DumP(v,l, SB,kB) 
           PrXYZ(x,y,:,v,2) = PrXYZ(x,y,:,v,2) * XPr(3,v,l, SB,kB)  ! GPs & sibs
          
          else  ! SA inbred or selfing
          
            PrXYZ(x,y,:,v,:) = PrXYZ(x,y,:,v,:) * SUM(AKA2P(y, v, :) * PrPG)   
            call ParProb(l, GpID(m,SB,kB), m, 0, 0, PrH)
            do z=1,3
              PrXYZ(x,y,z,v,:) = PrXYZ(x,y,z,v,:) * SUM(AKA2P(v,z,:) *PrH)
            enddo
              
            do r=1, nS(sB,kB)
              Bi = SibID(r, sB, kB)  
              if (NFS(Bi) == 0) cycle
              if (catB==Bi .or. ParBisClone(r)) then
                PrW = 1D0
              else
                call ParProb(l, Parent(Bi, 3-kB), 3-kB, -1, 0, PrW)
              endif
              
              if (Parent(Bi,3-kB) < 0) then
                do e=1,3
                  do f=1, nS(-Parent(Bi, 3-kB), 3-kB)
                    Ei = SibID(f, -Parent(Bi, 3-kB),3-kB)
                    if (nFS(Ei) == 0) cycle
                    if (Parent(Ei,kB) == -SB)  cycle
                    call ParProb(l,Parent(Ei,kB),kB,Ei,-1,PrH) 
                    PrF = FSLik(l,Ei)
                    PrH = PrH * PrF(:,v)                                   
                    if (.not. all(PrH==1D0))  PrW(e) = PrW(e) * SUM(PrH)
                  enddo
                enddo
              endif
              
              if (catB==Bi .and. m==kB) then
                do z=1,3
                  PrXYZ(x,y,z,v,1) = PrXYZ(x,y,z,v,1) * PrW(z)
                enddo
              else if (catB==Bi .and. m/=kB) then
                PrXYZ(x,y,:,v,1) = PrXYZ(x,y,:,v,1) * PrW(y)
              else if (ParBisClone(r)) then
                PrXYZ(x,y,:,v,1) = PrXYZ(x,y,:,v,1) * PrW(v)
              else if (.not. all(PrW==1D0)) then
                PrXYZ(x,y,:,v,1) = PrXYZ(x,y,:,v,1) * SUM(PrW)
              endif
              
              PrF = FSLik(l,Bi)
              PrW = PrW * PrF(:,v) 
              
              if (catB==Bi .and. m==kB) then
                do z=1,3
                  PrXYZ(x,y,z,v,2) = PrXYZ(x,y,z,v,2) * PrW(z)
                enddo
              else if (catB==Bi .and. m/=kB) then
                PrXYZ(x,y,:,v,2) = PrXYZ(x,y,:,v,2) * PrW(y)
              else if (ParBisClone(r)) then
                PrXYZ(x,y,:,v,2) = PrXYZ(x,y,:,v,2) * PrW(v)
              else if (.not. all(PrW==1D0)) then
                PrXYZ(x,y,:,v,2) = PrXYZ(x,y,:,v,2) * SUM(PrW)
              endif

            enddo  ! r
          endif     
        enddo  ! v
      enddo  ! y          
    enddo  ! x
    PrL(l) = LOG10(SUM(PrXYZ(:,:,:,:,2))) - LOG10(SUM(PrXYZ(:,:,:,:,1)))      
  enddo
  LLGX(m,2) = SUM(PrL)
enddo
LL = MaxLL((/LLGX(:,1), LLGX(:,2)/))

end subroutine DummyGP

! ######################################################################

subroutine dummyHFA(SA,kA,SB,kB, hf, LL)   ! SB (not Bi's) is HA/FA of SA (not Ai's)
use Global
use CalcLik             
implicit none

integer, intent(IN) :: SA, kA, SB, kB, hf
double precision, intent(OUT) :: LL
integer :: i, r, AncB(2,mxA), m, DoQuickA, DoQuickB, l, x, y, z, u,v
double precision :: PrL(nSnp), PrGG(3,2), PrGA(3), PrXY(3,3,3,3,3)

if (nS(SA,kA)==0 .or. ns(SB,kB)==0) then
  LL = NotImplemented
  return
endif

LL = missing
do i=1, nS(SB,kB)
  if (Parent(SibID(i,SB,kB), 3-kB) == GpID(3-kB,SA,kA) .and. GpID(3-kB,SA,kA)/=0) then
    LL = NotImplemented
    return     
  endif
  if (kA /= kB) then
    if (Parent(SibID(i,SB,kB), kA) == -SA) then
      LL = NotImplemented
      return     
    endif
  else if (kA == kB .and. Parent(SibID(i,SB,kB), 3-kB)/=0) then   
    do r= 1, nS(SA, kA)
      if (Parent(SibID(i,SB,kB), 3-kB)==Parent(SibID(r,SA,kA), 3-kA)) then
        LL = NotImplemented
        return      
      endif
    enddo
  endif
enddo
do r= 1, nS(SA, kA)
  if (Parent(SibID(r,SA,kA), 3-kA) == GpID(3-kA,SB,kB) .and. GpID(3-kA,SB,kB)/=0) then
    LL = NotImplemented
    return     
  endif
enddo
 

call GetAncest(-SB, kB, AncB)
if (ANY(AncB(kA,3:mxA) == -SA)) then
  LL = impossible 
  return
endif

if (all(GpID(:,SA,kA)/=0)) then
  LL = NotImplemented
  return
else 
  do m=1,2
    if (GpID(m,SA,kA)/=0 .and. GpID(m,SA,kA)==GpID(m,SB,kB)) then
      LL = NotImplemented
      return
    endif
  enddo
endif

if (GpID(1,SA,kA)==0) then
  m = 1
else
  m = 2
endif

call ChkDoQuick(SA,kA,DoQuickA)
call ChkDoQuick(SB,kB,DoQuickB)
if (DoQuickA==-1 .or. DoQuickA==-3 .or. DoQuickB==-1 .or. DoQuickB==-3) then
  LL = NotImplemented
  return
endif

PrL = 0D0
do l=1, nSnp
  do x=1,2
    call ParProb(l, GpID(x,SB,kB), x, 0, 0, PrGG(:,x))
  enddo
  call ParProb(l, GpID(3-m,SA,kA), 3-m, 0,0, PrGA)
  
  do x=1,3  ! SA
    do y=1,3  ! GpID(m,SA,kA)
      do z=1,3  ! SB
        do u=1,3  ! GpID(1,SB,kB)
          do v=1,3  ! GpID(2,SB,kB)
            PrXY(x,y,z,u,v) = XPr(1,x,l,SA,kA)  * SUM(AKA2P(x,y,:) * PrGA) * &
              XPr(1,z,l,SB,kB) * AKA2P(z,u,v) * PrGG(u,1) * PrGG(v,2)
            if (hf==3) then
              PrXY(x,y,z,u,v) = PrXY(x,y,z,u,v) * AKA2P(y,u,v)
            else if (hf==1) then
              PrXY(x,y,z,u,v) = PrXY(x,y,z,u,v) * AKAP(y,u,l)
            else if (hf==2) then
              PrXY(x,y,z,u,v) = PrXY(x,y,z,u,v) * AKAP(y,v,l)
            endif
          enddo
        enddo
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXY))
enddo

LL = SUM(PrL)
  
end subroutine dummyHFA

! ######################################################################
!       Age priors  
! ######################################################################

subroutine getEstBY(A, kA, lvl, BYLR)
use Global
implicit none
! lvl: 1=own; 2=own + est. from relatives exact BY + parent yearlast; 
! 3= 2 + est. BY from par + GP, 4= 2 + est. BY from offspr, 5 = all (NOT est BY from sibs)

integer, intent(IN) :: A, kA, lvl
double precision, intent(OUT) :: BYLR(nYears)

BYLR = LOG10(zero)                  
if (A > 0) then
  if (BY(A) > 0) then
    BYLR(BY(A)) = zero
  else
    BYLR = IndBY(:, A, lvl)
  endif
else if (A < 0) then
  BYLR = DumBY(:, -A, kA, lvl)
endif

end subroutine getEstBY

! ######################################################################

subroutine setEstBY(A, k)
use Global
implicit none
   
! birth year probability distribution, based on offspring, parent & GP BY. (& own min/max) 
! updates IndBY(:,A,:) or DumBY(:,-A, k,:) as side effect, with D3 resp D4:
! 1: own exact BY or BYrange
! 2: own + contributions from Par, GP, Off, sibs with exact BY + YearLast from Par + GP
! 3: [2] + contr from est. BY from par + GP
! 4: [2] + contr from est. BY from offspr
! 5: [3] + contr from est. BY from offspr  (= all)

integer, intent(IN) :: A, k
integer :: w
double precision :: LPBY(nYears, 5), BYup(nYears,2), BYdown(nYears,2), BYsibs(nYears)

if (A == 0) then
  return
else if (A > 0) then
  if (BY(A) > 0)  return  ! birth year known exactly
endif

LPBY = zero
if (A > 0) then
  LPBY(:,1) = IndBY(:, A, 1)   ! own exact BY / BYrange only
else !if (A < 0) then
  LPBY(:,1) = LOG10(1.0D0/(nYears-1))   ! dummy --> by definition no own BY info
endif

BYup = 0D0
BYdown = 0D0
BYsibs = 0D0           
call CalcBYup(A, k, BYup)    ! info from parents + grandparents
call CalcBYdown(A, k, BYdown)  ! info from offspring
call CalcBYsibs(A, k, BYsibs) 

LPBY(:,2) = LPBY(:,1) + BYup(:,1) + BYdown(:,1) + BYsibs  ! relatives with exact BY
LPBY(:,3) = LPBY(:,1) + BYup(:,2) + BYdown(:,1) + BYsibs
LPBY(:,4) = LPBY(:,1) + BYup(:,1) + BYdown(:,2) + BYsibs
LPBY(:,5) = LPBY(:,1) + BYup(:,2) + BYdown(:,2) + BYsibs  ! everything

! scale to sum to unity 
LPBY = 10**LPBY
do w=2,5  
  if (SUM(LPBY(:,w)) > 0D0) then
    LPBY(:,w) = LPBY(:,w) / SUM(LPBY(:,w)) 
  endif
enddo
LPBY = LOG10(LPBY)

do w=2,5
  if (A > 0) then
    IndBY(:, A, w) = LPBY(:,w)  
  else !if (A < 0) then
    DumBY(:, -A, k, w) = LPBY(:,w)
  endif
enddo

end subroutine setEstBY

! ######################################################################

subroutine CalcBYup(A, kA, BYA)
! BY probabilities of indiv A based on its parents + grandparents
! NOTE: not scaled!
use Global
implicit none

integer, intent(IN) :: A, kA
double precision, intent(OUT) :: BYA(nYears, 2)   ! D2: exact BY only; all BY
integer :: Par(2), GP(2,2), y, m, g, x, Yfirst, Ylast, tmpBY
double precision :: BYP(nYears, 2), BYG(nYears, 2,2), tmpX(nYears)
logical :: parent_info(2), gp_info(2,2)

BYA = zero
Par = getPar(A, kA)
if (ALL(par == 0))   return

GP = 0
GP(:,1) = getPar(Par(1), 1)
GP(:,2) = getPar(Par(2), 2)

! get current value of parent's & grandparents BY estimates
BYP = LOG10(zero)
BYG = LOG10(zero)
parent_info = .FALSE.
GP_info = .FALSE.
do m=1,2
  if (Par(m)==0)  cycle
  call getEstBY(Par(m), m, 3, BYP(:,m))  ! self + exact + parents + GP
  if (any(ABS(BYP(:,m) - LOG10(1d0/nYears)) > quite_tiny))  parent_info(m) = .TRUE.
  do g=1,2
    if (GP(g,m)==0)  cycle                       
    call getEstBY(GP(g,m), g, 3, BYG(:,g,m))
    if (any(ABS(BYG(:,g,m) - LOG10(1d0/nYears)) > quite_tiny))  GP_info(g,m) = .TRUE.
  enddo
enddo

if (any(GP_info)) then
  Yfirst = 3
else if (any(parent_info)) then
  Yfirst = 2  ! has parents with age info --> cannot be born in year 1 
else
  Yfirst = 0
endif
do m=1,2
  if (Par(m)>0) then
    if (BY(Par(m)) > 0) then
      Yfirst = MAX(Yfirst, BY(Par(m)) +1)   ! born at the earliest 1 year after parent was born
    endif
  endif
enddo

! last possible BY based on parents' & GP's last year of reproduction
Ylast = 999
do m=1,2
  if (Par(m) == 0)  cycle
  if (Par(m) > 0)  Ylast = MIN(Ylast, YearLast( Par(m) ) )
  do g=1,2
    if (GP(g,m) > 0)  Ylast = MIN(Ylast, YearLast( GP(g,m) ) +MaxAgePO )  
  enddo
enddo

if (.not. any(parent_info) .and. .not. any(GP_info) .and. Ylast==999)  return


do m=1,2
  if (.not. parent_info(m) .and. .not. any(GP_info(:,m)))  cycle
  do y=1,nYears    
    if (y < Yfirst) then
      BYA(y,:) = LOG10(zero)
      cycle
    else if (y > Ylast) then
      BYA(y,:) = LOG10(zero)  ! year y after parents last possible year of reprod
      exit
    else if (ANY(BYP(:,m)>=0D0)) then  ! parent has exact BY
      tmpBY = MAXLOC(BYP(:,m), DIM=1)   ! Par(m) may be dummy --> cannot use BY(Par(m)) 
      BYA(y,:) = BYA(y,:) + getAP(y - tmpBY, 1, m, 0, LOG10(zero))
    else 
      ! weighed sum over all possible parent birth years      
      tmpX = 0D0
      do x=1, y-1 
        tmpX(x) = 10**(BYP(x,m) + getAP(y-x, 1, m, 0, LOG10(zero)))  ! parent born in year x
      enddo
      BYA(y,2) = BYA(y,2) + LOG10(SUM(tmpX))
    endif
    
    ! grandparents 
    do g=1,2
      if (.not. GP_info(g,m))  cycle
      if (ANY(BYG(:,g,m)>=0D0)) then   ! GP has exact BY
        tmpBY = MAXLOC(BYG(:,g,m), DIM=1)
        BYA(y,:) = BYA(y,:) + getAP(y - tmpBY, 4, m, g, LOG10(zero))
      else ! Yfirst ensures that y>2
        tmpX = 0D0
        do x=1, y-2 
          tmpX(x) = 10**(BYG(x,g,m) + getAP(y-x, 4, m, g, LOG10(zero)))  ! grandparent born in year x
        enddo 
        BYA(y,2) = BYA(y,2) + LOG10(SUM(tmpX))
      endif
    enddo

  enddo
enddo

end subroutine CalcBYup

! ######################################################################

subroutine CalcBYdown(A, kA, BYA)
! BY probabilities of indiv A based on its offspring
! NOTE: not scaled!
use Global
implicit none

integer, intent(IN) :: A, kA
double precision, intent(OUT) :: BYA(nYears, 2)   ! D2: exact BY only; all BY
integer :: nOff, Offspr(maxSibSize), sxOff(maxSibSize), y, x, i, tmpBY
double precision :: tmpX(nYears)
double precision, allocatable :: BYO(:,:)
logical :: offspring_info(maxSibSize)

allocate(BYO(nYears, maxSibSize))
BYA = zero

call getOff(A,kA, .TRUE., nOff, Offspr, sxOff)
if (nOff == 0)   return

BYO = LOG10(zero)  ! number of offspring born in year y
offspring_info = .FALSE. 
if (nOff > 0) then
  do i=1, nOff
    call getEstBY(Offspr(i), sxOff(i),4, BYO(:,i))   ! self + exact + offspring
    if (any(ABS(BYO(:,i) - LOG10(1d0/nYears)) > quite_tiny))  offspring_info(i) = .TRUE.
  enddo
endif

if (any(offspring_info))  BYA(nYears,:) = LOG10(zero)  ! has offspring --> cannot be born in last year    
do y=1, nYears-1   ! A's BY
  do i=1, nOff
    if (.not. offspring_info(i))  cycle
    if (ANY(BYO(:,i)>=0D0)) then  ! offspring i has exact BY
      tmpBY = MAXLOC(BYO(:,i), DIM=1)   ! Off(i) may be dummy --> cannot use BY(Off(i))
      BYA(y,:) = BYA(y,:) + getAP(tmpBY - y, 1, kA, 0, LOG10(zero))
    else
      ! weighed sum over all possible offspring birth years x
      tmpX = 0D0
      do x=y+1, nYears  ! offspring BY
        tmpX(x) = 10**(BYO(x,i) + getAP(x-y, 1, kA, 0, LOG10(zero)))
      enddo 
      BYA(y,2) = BYA(y,2) + LOG10(SUM(tmpX))
    endif
    if (BYA(y,2) < -HUGE(0.0D0)) exit  ! e.g. i born in/prior to year y - no need to look at other offspr
  enddo
enddo

deallocate(BYO)

end subroutine CalcBYdown

! ######################################################################

subroutine CalcBYsibs(A, kA, BYA)
! BY probabilities of indiv A based on its siblings (mat + pat + full)
! NOTE: not scaled!  
! NOTE2: sibs with exact BY only, to avoid double contributions of shared parents/shared offspring
use Global
implicit none

integer, intent(IN) :: A, kA
double precision, intent(OUT) :: BYA(nYears)   
integer :: Par(2), nSibs(3), Sibs(maxSibSize, 3), sxSibs(maxSibSize), m, y, i, AgeD
double precision :: BYtmp(nYears, 3)                                    

BYA = zero           
Par = getPar(A, kA)
if (ALL(par == 0)) then  ! no parents --> no siblings.  
  return
endif

nSibs = 0
sibs = 0        
do m=1,2
  call getOff(Par(m),m, .FALSE., nSibs(m), Sibs(:,m), sxSibs)  ! non-dummy sibs only
enddo

if (ALL(nSibs(1:2) <= 1)) then   ! no siblings (only offspring of parents = focal indiv, if not dummy)
  return
endif

! FS: intersect between mat + pat & exclude self
do m=1,2
  do i=1, nSibs(m)
    if (Sibs(i,m) == A .or. BY(sibs(i,m)) <0) then  ! do not use if unknown/uncertain BY
      Sibs(i,m) = 0
    else if (Parent(Sibs(i,m),1) == Par(1) .and. Parent(Sibs(i,m),2) == Par(2) .and. &
      Par(1)/=0 .and. Par(2)/=0) then
      if (m==1) then
        nSibs(3) = nSibs(3) +1
        Sibs(nSibs(3), 3) = Sibs(i,m) 
      endif
      Sibs(i,m) = 0   ! only half siblings for pat/mat to avoid double counting
    endif
  enddo
enddo

BYtmp = zero               
do m=1,3
  if (nSibs(m) == 0) then
    BYtmp(:,m) = LOG10(1.0D0/nYears)
  else
    do y=2, nYears
      do i=1, nSibs(m)
         if (Sibs(i,m) == 0)  cycle
        if (BY(Sibs(i,m)) > 0) then  ! sibling i has known BY
          AgeD = ABS(BY(Sibs(i,m)) - y)  ! absolute age difference
          if (m < 3) then  ! half sibs
            if (Par(m) < 0) then
              BYtmp(y,m) = BYtmp(y,m) + getAP(AgeD, 3, 0, m, LOG10(zero))
            else ! genotyped parent: sib genotypes not in likelihood --> unbalanced if multiplying
              BYtmp(y,m) = BYtmp(y,m) + 10**getAP(AgeD, 3, 0, m, LOG10(zero))
            endif
          else  ! full sibs
            if (any(Par < 0)) then
              BYtmp(y,m) = BYtmp(y,m) + getAP(AgeD, 2, 0, 0, LOG10(zero))
            else
              BYtmp(y,m) = BYtmp(y,m) + 10**getAP(AgeD, 3, 0, m, LOG10(zero))
            endif
          endif
        endif
      enddo
    enddo
  endif
enddo

! scale 
do m=1,3
  if (m < 3) then
    if (Par(m) <= 0 .or. COUNT(sibs(:,m)/=0)==0)  BYtmp(:,m) = 10**BYtmp(:,m)
  else
    if (any(Par <= 0) .or. COUNT(sibs(:,m)/=0)==0)  BYtmp(:,m) = 10**BYtmp(:,m)
  endif
  if (SUM(BYtmp(:,m)) > 0D0) then
   BYtmp(:,m) = BYtmp(:,m) / SUM(BYtmp(:,m)) 
  endif
enddo
BYtmp = log10(BYtmp)

BYA = BYtmp(:,1) + BYtmp(:,2) + BYtmp(:,3)

BYA = 10**BYA
if (SUM(BYA) > 0D0) then
  BYA = BYA / SUM(BYA) 
endif
BYA = LOG10(BYA)

end subroutine CalcBYsibs

! #####################################################################

subroutine EstAgeDif(A, kA, B, kB, AgeD)   ! estimate age difference, incl. from BYrange. Only called by CalcU()
use Global
implicit none

integer, intent(IN) :: A, kA, B, kB
double precision, intent(OUT) :: AgeD
integer :: y, x
double precision :: pBY(nYears, 2)
double precision, allocatable :: ADtmp(:,:)

allocate(ADtmp(nYears, nYears))

if (A>0 .and. B>0) then
  if (AgeDiff(A,B) < 999) then
    AgeD = REAL(AgeDiff(A,B), 8)
  endif
endif

pBY = LOG10(zero)
call getEstBY(A, kA, 5, pBY(:, 1))  ! all contributions from all relatives
call getEstBY(B, kB, 5, pBY(:, 2))
pBY = 10**pBY  ! log -> regular scale

ADtmp = 0D0
do x=1, nYears  ! A 
  if (pBY(x,1) < TINY(0.0D0)) cycle  ! no chance that A is born in year x
  do y=1,nYears  ! B
    if (pBY(y,2) < TINY(0.0D0)) cycle
    ADtmp(x,y) = pBY(x,1) * pBY(y,2) * (x-y)   ! if B older, than AgeD > 0
  enddo
enddo

AgeD = SUM(ADtmp)

deallocate(ADtmp)

end subroutine EstAgeDif

! #####################################################################

subroutine CalcAgeLR(A, kA, B, kB, m, focal, AllDumRel, ALR) ! m: mat/pat relatives
use Global
implicit none

integer, intent(IN) :: A, kA, B, kB, m, focal
logical, intent(IN) :: AllDumRel
double precision, intent(OUT) :: ALR
integer :: AB(2), kAB(2), x, y, i, n, fcl, YearLast_B
double precision :: BYLR(nYears, 2), ALRm(2), ALRtmp(nYears, nYears)
!double precision, allocatable :: ALRtmp(:,:)    ! TODO: separate module with age stuff; allocate 1x
logical :: has_BY_info(2)

if (.not. ANY((/-1,1,2,3,4,5,6/) == focal))  call Erstop('CalcAgeLR: illegal focal', .TRUE.)    

AB = (/ A, B /)
kAB = (/ kA, kB /)  
ALR = zero
if (A==0 .or. B==0) then
  return  
else if (A>0 .and. B>0) then
  if (focal==1 .and. BY(A)>0) then
    if (BY(A) > YearLast(B)) then  ! YearLast: unknown = +999
      ALR = impossible
      return
    endif
  endif
  
  if ((focal==2 .or. focal==3) .and. (any(Parent(A,:)/=0) .or. any(Parent(B,:)/=0))) then
    do i=1,2
      do n=1,2
        if (n/=m .and. focal/=2)  cycle
        if (Parent(AB(i),n) > 0) then
          ALR = getAP(AgeDiff(AB(3-i),Parent(AB(i),n)), 1, n, 0, Impossible)
          if (ALR == Impossible)  return   ! NOT else ALR = ALR * ALRp  - keep it simple. 
        endif
      enddo
    enddo
  endif

  if (m < 3) then  ! incl. m=0
    ALR = getAP(AgeDiff(A,B), focal, kB, m, LOG10(zero))
  else
    ALRm = LOG10(zero)
    do n=1,2
      ALRm(n) = getAP(AgeDiff(A,B), focal, kB, n, LOG10(zero))
    enddo
    ALR = MAXVAL(ALRm)
  endif 
  if (ALR < -HUGE(0.0D0)) then
    ALR = impossible
    return
  endif
  
  if (AgeDiff(A,B) /= 999)  return
endif  

ALR = Missing
fcl = focal
if (focal==1) then   ! short-cut instead of via dummy BYLR, faster
  if (A>0 .and. B<0) then
    if (BY(A)>0 .and. ns(-B,kB)==0) then
      do n=1,2
        if (GpID(n,-B,kB)>0 .and. GpID(3-n,-B,kB)==0) then
          if (BY(GpID(n,-B,kB))>0) then
            ALR = getAP(AgeDiff(A, GpID(n,-B,kB)), 4, n, kB, Impossible)
          endif
        endif
      enddo
    else if (ns(-B,kB)==1 .and. all(GpID(:,-B,kB)==0)) then
      if (Parent(A,3-kB)/=0 .and. Parent(SibID(1,-B,kB),3-kB)==Parent(A,3-kB)) then
        fcl=2
      else
        fcl=3
      endif
      if (BY(A)>0 .and. BY(SibID(1,-B,kB))>0) then
        ALR = getAP(AgeDiff(A,SibID(1,-B,kB)), fcl, 0, kB, Impossible)
      else
        AB(2) = SibID(1,-B,kB)     
      endif
    else if (BY(A)>0) then
      fcl = 3
      do i=1,ns(-B,kB)
        if (BY(SibID(i,-B,kB))< 0) cycle
        if (Parent(A,3-kB)/=0 .and. Parent(SibID(i,-B,kB),3-kB)==Parent(A,3-kB))  fcl=2
        ALR = getAP(AgeDiff(A,SibID(i,-B,kB)), fcl, 0, kB, Impossible)
        if (ALR == Impossible)  return
      enddo
      ALR = Missing  ! reset
      fcl = 1
    endif
  else if (A<0 .and. B>0) then
    if (ns(-A,kA)==1 .and. all(GpID(:,-A,kA)==0)) then
      if (BY(B)>0 .and. BY(SibID(1,-A,kA))>0) then
        ALR = getAP(AgeDiff(SibID(1,-A,kA),B),4,kB,kA, Impossible)
      else
        AB(1) = SibID(1,-A,kA)
        fcl = 4 
      endif   
    else if (BY(B)>0) then
      do i=1,ns(-A,kA)
        if (BY(SibID(i,-A,kA))< 0) cycle
        ALR = getAP(AgeDiff(SibID(i,-A,kA),B), 4,kB,kA, Impossible)
        if (ALR == Impossible)  return
      enddo
      ALR = Missing             
    endif
  endif  
endif
if (ALR /= Missing)  return

! temporary add-on; need to re-implement whole subroutine... TODO
if (A>0 .and. B<0 .and. (focal==5 .or. focal==6)) then  ! check: GP(SB) will become GP of A
  ALRm = 0d0
  do n=1,2
    if (focal==6 .and. n/=m)  cycle
    if (GpID(n,-B,kB) <=0) cycle
    if (BY(GpID(n,-B,kB)) <0)  cycle 
    ALRm(n) = getAP(AgeDiff(A,GpID(n,-B,kB)), 4,kA,n, Impossible)    
    if (ALRm(n) == Impossible) then
      ALR = impossible
      return
    endif
  enddo
endif

BYLR = LOG10(zero)  ! likelihood ratio to be born in year X
if (AllDumRel) then
  do i=1,2
    call getEstBY(AB(i), kAB(i), 5, BYLR(:, i))  ! all contributions from all relatives
  enddo
else 
  call getEstBY(AB(1), kA, 4, BYLR(:, 1))  ! excl contributions from est. BY from par + GP
  call getEstBY(AB(2), kB, 3, BYLR(:, 2))  ! excl contributions from est. BY from offspring
endif

has_BY_info = .FALSE.
do i=1,2
  if (any(ABS(BYLR(:, i) - LOG10(1d0/nYears)) > quite_tiny))  has_BY_info(i) = .TRUE.
enddo

if (any(.not. has_BY_info)) then
  ALR = 0D0  ! no age info available for either or both individuals
   return
endif    

if (fcl==1 .or. fcl==4) then  ! quick check
  do y=2, nYears  ! B 
    if (BYLR(y,2) < -HUGE(0.0D0)) cycle
    ! at oldest possible BY of B:
    if (ALL(BYLR((y-1):nYears, 1) < -HUGE(0.0D0))) then
      ALR = impossible
      return
    else
      exit
    endif
  enddo
endif

YearLast_B = 9999    
if (AB(2)>0)  YearLast_B = YearLast(AB(2))
 
ALRm = LOG10(zero)
do n=1,2
  if (n/=m .and. (m==1 .or. m==2 .or. focal>5))  cycle
  if (m==0 .and. n==2)  cycle   ! PO or FS (does not use n)
  ALRm(n) = calc_ALRm(kA,kB, n,fcl)
enddo
ALR = MAXVAL(ALRm)
if (ALR < -HUGE(0.0D0) .or. ALR/=ALR .or. ALR < 2*ALR_tiny)   ALR = impossible


contains
  function calc_ALRm(kA,kB, m,focal)  
    integer, intent(IN) :: kA,kB, m,focal
    double precision :: calc_ALRm
    
    ALRtmp = LOG10(zero)  ! -Inf
    do y=1,nYears  ! B
      if (BYLR(y,2) < -HUGE(0.0D0)) cycle
      do x=1, nYears  ! A 
        if (BYLR(x,1) < -HUGE(0.0D0)) cycle 
        if (focal==1 .and. x > YearLast_B)  cycle  ! maybe-BY for A after B's last year of reproduction
        if (focal==4 .and. x >(YearLast_B + MaxAgePO))  cycle
        if ((x-y) < -MaxAgePO .or. (x-y) > nYears)  cycle
        if (focal==-1) then  ! A==B
          if (x==y)  ALRtmp(x,y) = BYLR(x,1) + BYLR(y,2)
        else if (focal<=5) then
          ALRtmp(x,y) = BYLR(x,1) + BYLR(y,2) + getAP(x-y, focal, kB, m, LOG10(zero))  
        else
          ALRtmp(x,y) = BYLR(x,1) + BYLR(y,2) + getAP(x-y, focal, m, kA, LOG10(zero))
        endif                    
      enddo
    enddo
    calc_ALRm = LOG10(SUM(10**ALRtmp))  ! sum across age differences
    
  end function calc_ALRm
  
end subroutine CalcAgeLR

! ######################################################################

subroutine CalcALRmerge(SA, SB, k, ALR)  ! change in ALR when merging
use Global
implicit none

integer, intent(IN) :: SA, SB, k
double precision, intent(OUT) :: ALR
double precision :: ALRj                        
integer :: i, j

ALR = 0D0
do i = 1, nS(SA,k)
  if (BY(SibID(i,SA,k))<0) cycle
  do j=1, nS(SB, k)
    ALRj = getAP(AgeDiff( SibID(i,SA,k), SibID(j,SB,k)), 3, 0, k, Impossible)
    if (ALRj == Impossible) then
      ALR = impossible
      return
    else
      ALR = ALR + ALRj
    endif
  enddo
enddo

ALR = ALR / (ns(SA,k) * ns(SB,k))   ! else not comparable across sibship sizes

end subroutine CalcALRmerge

! ######################################################################
  
subroutine calcALR_addsib(A,SB,k,focal, dALR)   ! change in ALR when adding A to SB
use Global
implicit none
  integer, intent(IN) :: A, SB, k, focal
  double precision, intent(OUT) :: dALR
  integer :: j
  double precision :: ALRj
  
  dALR = 0D0
  do j=1,ns(SB,k)
    call CalcAgeLR(A,3, SibID(j,SB,k),3, k, focal, .TRUE., ALRj)
    if (ALRj == Impossible) then
      dALR = impossible
      return
    else
      dALR = dALR + ALRj
    endif
  enddo
  
  dALR = dALR / ns(SB,k)

end subroutine calcALR_addsib

! #####################################################################

! #####################################################################

subroutine BestRel(LLIN, focal, X, dLL)
use Global
implicit none
! return which relationship is most likely, by threshold TA
! assuming order PO,FS,HS,GG,FAU,HAU,U in LL vector

double precision, intent(IN) :: LLIN(7)
integer, intent(IN) :: focal
integer, intent(OUT) :: X
double precision, intent(OUT) :: dLL   ! diff best vs next best
integer :: i,j, Y(7)
double precision :: LL(7)
logical:: Maybe(6)

X = 0
dLL = 0D0
LL = LLIN

if (ALL(LL(1:6) > 0)) then
  X = 8
  return
endif

if (focal==3 .and. LL(3)<0d0 .and. LL(2)<0d0) then   ! want sib vs non-sib
  if (focal==3 .and. LL(2) - LL(3) < TA) then   ! FS less likely, or sliiiightly more likely
    LL(2) = 333D0  
  else if (focal /= 2 .and. LL(2)>=LL(3)) then
    LL(3) = 333D0
  endif
endif

if ((LL(7) - MAXVAL(LL(1:6), MASK=LL(1:6)<0d0)) > TA) then  
  X = 7  ! unrelated
else
  maybe = .TRUE.
  do i=1,6
    if (LL(i)>0) then
      maybe(i) = .FALSE.
    else
      do j=1,7
        if (i==j) cycle
        if (LL(j)>0d0) cycle
        if ((LL(i) - LL(j)) < TA) then
          maybe(i) = .FALSE.   ! i has no longer highest LL
        endif
      enddo
    endif
  enddo
  if (COUNT(maybe)==0) then
    X = 8  ! unclear
  else if (COUNT(maybe)==1) then
    X = MAXLOC(LL(1:6), MASK=maybe, DIM=1)
  endif       
endif

Y = (/(i, i=1,7, 1)/)
if (X<8 .and. X>0) then
  dLL = LL(X) - MAXVAL(LL, MASK=(LL<0d0 .and. Y/=X))
endif  
      
end subroutine BestRel

! #####################################################################

subroutine BestRel2(LLIN, X, dLL)
use Global
implicit none
! as BestRel, but no threshold, and consider all 1st & 2nd degree rel

double precision, intent(IN) :: LLIN(7)
integer, intent(OUT) :: X
double precision, intent(OUT) :: dLL !(2)   ! diff best vs next best
integer :: i,j, Y(7)
double precision :: LL(7), small
logical:: Maybe(6)

X = 0
dLL = 0D0
LL = LLIN
small = 0.01d0

if (MAXVAL(LL(1:6), MASK=LL(1:6)<0d0) - LL(7) < MAX(TA, small) .or. &
  ALL(LL(1:6) > 0d0)) then  
  X = 7  ! unrelated 
else
  maybe = .TRUE.
  do i=1,6
    if (LL(i)>0d0) then
      maybe(i) = .FALSE.
    else
      do j=1,7
        if (i==j) cycle
        if (LL(j)>0) cycle
        if ((LL(i) - LL(j)) < small) then
          maybe(i) = .FALSE.   ! i has no longer highest LL
        endif
      enddo
    endif
  enddo
  if (COUNT(maybe)==1) then
    X = MAXLOC(LL(1:6), MASK=Maybe, DIM=1)
  else if (ABS(MaxLL(LL(3:5))-MaxLL(LL))<small .and. COUNT(LL(3:5)<0d0) >1) then
  !  MAXVAL(LL(3:5), MASK=LL(3:5)<0d0,DIM=1) - MINVAL(LL(3:5), MASK=LL(3:5)<0d0,DIM=1) < small) then
    X = 9  ! any 2nd degree relative
  else
    X = 8  ! unclear
  endif   
endif

Y = (/(i, i=1,7, 1)/)
if (count(LL < 0d0) < 2) then
  dLL = -777D0   ! NOT positive !!
else if (X<8) then
  dLL = LL(X) - MAXVAL(LL, MASK=(LL<0d0 .and. Y/=X))
!  dLL(2) = LL(X) - MaxLL(LL(6:7))
else if (X==9) then
  dLL = MaxLL(LL(3:5)) - MaxLL(LL((/1,2,6,7/)))
!  dLL(2) = MaxLL(LL(3:5)) - MaxLL(LL(6:7))
endif  

end subroutine BestRel2

! #####################################################################

subroutine UpdateAllProbs
use Global
use CalcLik           
use sort_module               
implicit none

integer :: i, k, s, x, y, r, n=30, BYrankI(nInd)
double precision :: Lind_IN(nInd)
double precision, allocatable :: XPr_IN(:,:,:,:,:)
integer, allocatable :: s_sorted(:), k_sorted(:)           

allocate(XPr_IN(3,3,nSnp,nInd/2,2))

do k=1,2
  do s=1,nC(k)
    if (nC(k)==0)  cycle                    
    if (ALL(GpID(:,s,k)==0) .and. ALL(SibID(:,s,k)==0)) then
      call Erstop("Empty sibship!", .TRUE.)
    endif
  enddo
enddo

! individuals ~~~~                  
call getRank_i(BYrankI)

do x=1, nInd
  i = BYRankI(x)
  call CalcLind(i)   
enddo

do r=1,2
  do x=1, nInd
    i = BYRankI(x)
    call setEstBY(i, Sex(i))     
  enddo

  do x=nInd, 1, -1
    i = BYRankI(x)
    call setEstBY(i, Sex(i))     
  enddo
enddo

! sibships ~~~     
  call sort_sibships(s_sorted, k_sorted)

  do y=1,n  
    do r=1,n
      XPr_IN = XPr
      do x=1, SUM(nC)
        s = s_sorted(x)
        k = k_sorted(x)                       
        call CalcCLL(s,k)  
      enddo        
      if (all(abs(XPr_IN - XPr) < 0.1))  exit
    enddo
    Lind_IN = Lind
    do x=1, nInd
      i = BYRankI(x)
      call CalcLind(i)
    enddo
    if (all(abs(Lind_IN - Lind) < 0.1))   exit
  enddo

  do x=1, SUM(nC)
    s = s_sorted(x)
    k = k_sorted(x)                    
    call setEstBY(-s, k)
  enddo

deallocate(XPr_IN, s_sorted, k_sorted)

end subroutine UpdateAllProbs

! #####################################################################

subroutine CalcCLL(s, k) 
use Global
use CalcLik             
implicit none
! returns XPr: likelihood;  DumP: probability, scaled  (no age prior.),
! split into 1: sibs only 2: gp effect only, 3: all

integer, intent(IN) :: s, k ! S: sibship number, k: mat(1),pat(2),unk(3)
integer :: Sibs(ns(s,k)), UseEE(ns(s,k)), MatePar(ns(s,k)),  cat, catG, &
  IsInbr(ns(s,k)), AncR(2,mxA), FSX, &
  l, x, i, Ei, r, y, z, g, Ri, v, e
double precision :: PrL(nSnp), PrY(3), PrYp(3,ns(s,k)), PrGG(3,2),&
 PrZ(3),PrXZ(3,3,2), PrE(3), PrEE(3, ns(s,k)), LPrX(3,2), PrRF(3,3,ns(s,k)), PrF(3,3)
logical :: ParIsClone(ns(s,k)), DoRsibs(maxSibSize)
integer, allocatable :: HasInbr(:,:)   ! TODO: allocate 1x & reuse

allocate(HasInbr(ns(s,k), ns(s,k)))

if (ALL(GpID(:,s,k)==0) .and. ALL(SibID(:,s,k)==0)) then
  CLL(s,k) = 0D0
  XPr(1,:,:,s,k) = 1D0 
  do l=1,nSnp
    XPr(2,:,l,s,k) = AHWE(:,l)
    XPr(3,:,l,s,k) = AHWE(:,l)
    DumP(:,l,s,k) = AHWE(:,l)
  enddo
  return 
endif

Sibs = SibID(1:ns(s,k), s, k)
UseEE = 0
MatePar = 0
call FindEE(Sibs, ns(s,k), 0, k, UseEE, MatePar)   ! may shuffle sibs

!================= 
cat = 0
catG = 0
IsInbr = 0
HasInbr = 0
AncR = 0
ParIsClone = .FALSE.                    
do r=1,nS(s,k)
  Ri = Sibs(r)  
  if (nFS(Ri) > ns(s,k)) then
    call Rprint("nFS > nS! ", (/k, s, SibID(1:5, s, k)/), (/0D0/), "INT")
    call ErStop("something wrong with sibship cluster", .TRUE.)
  endif
  if (Parent(Ri, 3-k)==0) cycle
  if (Parent(Ri, 3-k)==GpID(3-k,s,k) .and. nFS(Ri)/=0) then  
    cat = Ri
    UseEE(r) = 0
  endif
  do v=1, nS(s,k)
    if (r==v) cycle
    if (nFS(Sibs(v))==0) cycle
    do i=1, nFS(Sibs(v))
      if (Parent(Ri, 3-k) == FSID(i, Sibs(v))) then
        IsInbr(r) = FSID(i, Sibs(v))
        HasInbr(v,i) = r !-1
      endif
    enddo
  enddo
  if (IsInbr(r)/=0) cycle
  call GetAncest(Ri,k,AncR)
  if (AncR(k, 5-k) == -s) then 
    IsInbr(r) = Parent(Ri, 3-k)   ! via dummy
  endif
  if (hermaphrodites/=0 .and. DumClone(s,k)/=0) then
    if (DumClone(s,k) == -Parent(Ri,3-k))  ParIsClone(r) = .TRUE.
  endif
enddo  

FSX = 0  
if (ALL(GpID(:,s,k)<0) .and. cat==0) then  ! check if sibship par inbred
  if (GPID(1,s,k) == GPID(1, -GPID(2,s,k),2)) then
    catG = 2
  else if (GPID(2,s,k) == GPID(2, -GPID(1,s,k),1)) then
    catG = 1
  else 
    do i=1, ns(-GpID(1,s,k), 1)
      if (Parent(SibID(i, -GpID(1,s,k), 1), 2) == GpID(2,s,k)) then  ! FS of dummy par
        catG = 3
        if (nFS(SibID(i, -GpID(1,s,k), 1))/=0) then
          FSX = SibID(i, -GpID(1,s,k), 1)
        endif
      endif
    enddo 
  endif
endif

call ChkTooManySibs(Sibs, ns(s,k), k, DoRsibs)  ! prevent numerical issues w huge sibships

PrL = 0D0       
do l=1,nSnp
  do g=1,2   !grandparents
    if (g/=k .and. cat>0) then
      call ParProb(l, GpID(g,s,k), g, -1,0, PrGG(:,g))
    else if (catG==g) then
      PrGG(:,g) = XPr(3,:,l, -GpID(g,s,k),g)
    else if (catG==3) then
      call ParProb(l, GpID(g,s,k), g, FSX,-1, PrGG(:,g))
    else
      call ParProb(l, GpID(g,s,k), g, 0,0, PrGG(:,g))
    endif
  enddo
  
  do x=1,3  ! genotype dummy parent
    do z=1,3
      if (catG==k) then
        PrXZ(x,z,:) = SUM(AKA2P(x, z, :) * PrGG(:,k))
      else if (catG==3-k) then
        PrXZ(x,z,:) = SUM(AKA2P(x, z, :) * PrGG(z,3-k))
      else if (catG==3) then
        PrE = PrGG(:,k)
        do i=1, nFS(FSX)
          PrE = PrE * OKA2P(Genos(l,FSID(i,FSX)), :, z)
        enddo          
        PrXZ(x,z,:) = SUM(AKA2P(x, z, :) * PrGG(z,3-k) * PrE)
      else  ! catG==0
        PrXZ(x,z,:) = SUM(AKA2P(x, z, :) * PrGG(z,3-k) * PrGG(:,k))  ! GPs
      endif
    enddo
  enddo
  if (catG>0) then
    do v=1,2
      PrXZ(:,:,v) = PrXZ(:,:,v)/SUM(PrXZ(:,:,v)) 
    enddo
  endif
  do x=1,3
    XPr(2,x,l, s,k) = SUM(PrXZ(x,:,2))  ! GP 
  enddo
  LPrX = log10(SUM(PrXZ, DIM=2))  ! sum over z
  
  PrYp = 0D0
  PrRF = 1D0               
  do r=1, nS(s,k)
    Ri = Sibs(r)
    if (NFS(Ri) == 0) cycle
    if (IsInbr(r) /= 0 .or. cat==Ri .or. ParIsClone(r)) then
      PrYp(:,r) = 1D0
    else if (DoRsibs(r)) then
      call ParProb(l, Parent(Ri, 3-k), 3-k, -1,0, PrYp(:,r))
    else
      call ParProb(l, Parent(Ri, 3-k), 3-k, Ri, -1, PrYp(:,r)) 
    endif
    if (all(HasInbr(r,:)==0)) PrRF(:,:,r) = FSLik(l,Ri)   
  enddo

  do z=1,3
    if (z>1 .and. cat==0) cycle
  do x=1,3
    PrEE = 0D0          
    do r=1, nS(s,k)
      Ri = Sibs(r)  ! array with IDs
      if (NFS(Ri) == 0) cycle  ! moved to its FS
      if (IsInbr(r) > 0) then
        cycle
      else if (IsInbr(r) < 0) then
        call ParProb(l, GpID(3-k,-Parent(Ri, 3-k),3-k), 3-k, 0,0, PrZ) 
        do y=1,3
          PrYp(y,r) = SUM(AKA2P(y,x,:) * PrZ)
        enddo
      else if (UseEE(r) /= 0) then
        call ParProb(l, MatePar(r), k, 0,0, PrZ)  !  GpID(k,-Parent(Ri, 3-k),3-k)
        do y=1,3
          do e=1,3
            PrE(e) = SUM(AKA2P(y,e,:) * PrEE(e,UseEE(r)) * PrZ)
          enddo
          PrYp(y,r) = SUM(PrE)
        enddo
        PrYp(:,r) = PrYp(:,r) / SUM(PrYp(:,r))
      endif
      PrY = PrYp(:,r)
        
      if (Parent(Ri, 3-k)<0 .and. DoRsibs(r))  then
        do y=1,3   ! parent 3-k 
          do v=1, nS(-Parent(Ri, 3-k), 3-k)
            Ei = SibID(v, -Parent(Ri, 3-k), 3-k)  
            if (NFS(Ei) == 0) cycle
            if (Parent(Ei, k) == -s) cycle
            call ParProb(l, Parent(Ei, k), k, Ei,-1, PrE)
            PrF = FSLik(l,Ei)
            PrE = PrE * PrF(:,y)
            if (.not. ALL(PrE==1D0))  PrY(y) = PrY(y) * SUM(PrE)
          enddo
        enddo
      endif
      
      if (.not. ALL(PrY==1D0)) then     
        if (cat==Ri) then
          PrXZ(x,z,1) = PrXZ(x,z,1) * PrY(z)
        else if (ParIsClone(r)) then
          PrXZ(x,z,1) = PrXZ(x,z,1) * PrY(x)          
        else if (cat/=0) then
          PrXZ(x,z,1) = PrXZ(x,z,1) * SUM(PrY)
        else
          PrXZ(x,:,1) = PrXZ(x,:,1) * SUM(PrY)
        endif
      endif
      if (cat==0 .and. .not. all(PrY==1D0)) then
        if (ParIsClone(r)) then
          LPrX(x,1) = LPrX(x,1) + log10(PrY(x))
        else
          LPrX(x,1) = LPrX(x,1) + log10(SUM(PrY))
        endif
      endif
     
      if (all(HasInbr(r,:)==0)) then
        PrY = PrY * PrRF(:,x,r)
      else
        do i=1, nFS(Ri)  ! default: nFS = 1    
          if (HasInbr(r,i)==0) then
            PrY = PrY * OKA2P(Genos(l,FSID(i,Ri)), x, :)
          else
            do y=1,3
              do e=1,3
                PrE(e) = AKA2P(e, x, y) * OcA(e,Genos(l,FSID(i,Ri)))
                do v=1, nS(s,k)
                  if (IsInbr(v)==FSID(i,Ri)) then  
                    PrE(e) = PrE(e) * OKA2P(Genos(l,Sibs(v)), x, e)
                  endif
                enddo
              enddo
              if (.not. ALL(PrE==1D0))  PrY(y) = PrY(y) * SUM(PrE)
            enddo
          endif
        enddo
      endif
      
      if (.not. ALL(PrY==1D0)) then        
        if (cat==Ri) then
          PrXZ(x,z,2) = PrXZ(x,z,2) * PrY(z)
        else if (ParIsClone(r)) then
          PrXZ(x,z,2) = PrXZ(x,z,2) * PrY(x)          
        else if (cat/=0) then
          PrXZ(x,z,2) = PrXZ(x,z,2) * SUM(PrY)
        else
          PrXZ(x,:,2) = PrXZ(x,:,2) * SUM(PrY)
        endif
        if (cat==0 .and. .not. all(PrY==1D0)) then
          if (ParIsClone(r)) then
            LPrX(x,2) = LPrX(x,2) + log10(PrY(x))
          else
            LPrX(x,2) = LPrX(x,2) + log10(SUM(PrY))
          endif
        endif
      endif
      PrEE(:,r) = PrY    
    enddo ! r 
  enddo ! x
  enddo ! z (cat/=0 only)
  if (cat==0) then
    XPr(3,:,l, s,k) = 10**LPrX(:,2) / SUM(10**LPrX(:,1))
  else
    do x=1,3  ! account for GP, dumm offspr & connected sibships
      XPr(3,x,l, s,k) = SUM(PrXZ(x,:,2))/ SUM(PrXZ(:,:,1))    ! offspring + grandparents
    enddo
  endif
  WHERE (XPr(3,:,l,s,k) /= XPr(3,:,l,s,k)) XPr(3,:,l,s,k) = 0D0  ! when Err=0
  do x=1,3
    DumP(x,l, s,k) = XPr(3,x,l, s,k)/ SUM(XPr(3,:,l, s,k))
    XPr(1,x,l, s,k) = XPr(3,x,l, s,k) / XPr(2,x,l, s,k)   ! offspring only
  enddo 
  PrL(l) = LOG10(SUM(XPr(3,:,l, s,k)))
enddo
CLL(s,k) = SUM(PrL) 
 
WHERE (XPr(1,:,:,s,k) /= XPr(1,:,:,s,k)) XPr(1,:,:,s,k) = 0D0  ! 0/0 when MAF=0
WHERE (DumP(:,:,s,k) /= DumP(:,:,s,k)) DumP(:,:,s,k) = 0D0
 
 
if (CLL(s,k)> .001 .or. CLL(s,k)/=CLL(s,k)) then   !.or. CLL(s,k)< -HUGE(1D0) 
  call Rprint("Problem: ", (/k, s, ns(s,k)/), (/0.0D0/), "INT")
  call Erstop("Invalid sibship LL - try increasing Err", .FALSE.)
endif

deallocate(HasInbr)

end subroutine CalcCLL

! #####################################################################

subroutine ChkTooManySibs(Sibs, n, k, DoRsibs)
use Global
implicit none

integer, intent(IN) :: n
integer, intent(IN) :: sibs(n), k
logical, intent(OUT) :: DoRsibs(maxSibSize)
integer :: r, i

! prevent numerical issues when sibships are very large
DoRsibs = .FALSE.
do r=1,n
  i = Sibs(r)
  if (nFS(i)==0)  cycle
  if (Parent(i,3-k) >=0) cycle
  if (ns(-parent(i,3-k),3-k) >50 .and. nFS(i) < ns(-parent(i,3-k),3-k)/5) then    ! What thresholds??
    DoRsibs(r) = .FALSE.
  else
    DoRSibs(r) = .TRUE.
  endif
enddo

end subroutine ChkTooManySibs

! #####################################################################

subroutine ParProb(l, i, k, A, B, prob)  
use Global
use CalcLik           
implicit none

integer, intent(IN) :: l, i, k, A,B
double precision, intent(OUT) :: prob(3)
integer :: x,j, AB(2), A1, parA
double precision :: PrP(3, 2), PrY(3), PrAF(3,3)
logical :: AllIN

prob = AHWE(:, l)    
A1 = 0             

if (i == 0) then  ! no parent
  if (A==-4 .or. B==-4 .or. B==-5) then
    prob = 1D0
  else
    prob = AHWE(:, l)
  endif

else if (i > 0) then  ! real parent
  if (A==-4 .or. B==-4 .or. B==-5) then
    prob = OcA(:,Genos(l,i))
  else
    prob = LindX(:, l, i)  ! =AHWE if Lind(i) not yet calculated  ! unscaled; scaled below.
  endif

else if (i < 0) then  ! dummy parent
  if (A==0) then   ! probability
    prob = DumP(:,l, -i,k)    
  else if (A == -1) then  ! grandparent contribution only
    prob = XPr(2,:,l, -i, k)
  else if (A==-4) then  ! offspring contribution only
    prob = XPr(1,:,l, -i, k)
  else if (A<0) then  ! shouldn't happen
    call Erstop("Invalid call to ParProb!", .TRUE.)

  else if (A>0) then   ! exclude indiv A from calc & standardise
    if (Parent(A,k)/=i .and. B<=0) then
      prob = DumP(:,l, -i,k)
    else if (ns(-i,k)<=1 .and. B> -4) then
      prob = XPr(2,:,l, -i, k)  ! grandparent contribution only
    else
     
      AB = (/ A, B /)
      A1 = FSID(maxSibSize+1, A)
      do j=1,2
        if (j==1 .and. (Parent(A,k)/=i)) cycle
        if (j==2) then
          if(B<=0) cycle
          if (Parent(B,k)/=i)  cycle
        endif
        if (Parent(AB(j), 3-k)==0) then  
          PrP(:,j) = AHWE(:,l)
        else if (Parent(AB(j), 3-k)>0) then
          PrP(:,j) = LindX(:, l, Parent(AB(j), 3-k))
          PrP(:,j) = PrP(:,j)/SUM(PrP(:,j))                                   
        else if (Parent(AB(j), 3-k)<0) then  
          PrP(:,j) = DumP(:,l, -Parent(AB(j),3-k), 3-k)
        endif
      enddo
      
      AllIN = .FALSE.
      if (B==-1) then               
        parA = Parent(A, 3-k)
        if (ns(-i,k) <= 1) then
          AllIN = .TRUE.
        else if (ParA/=0 .and. all(Parent(SibID(1:ns(-i,k),-i,k),3-k) == ParA)) then 
          AllIN = .TRUE.
        endif
      endif
   
      prob = 0D0
      PrAF = 1D0
      if ((B==-1 .and. .not. AllIN) .or. (B==-5 .and. nFS(A1) /= ns(-i,k))) then
        PrAF = FSLik(l,A1)
      endif
      do x = 1, 3
        if (B>=0) then    ! .or. (B==-1 .and. nFS(A1)<=1)
          prob(x) = XPr(3,x,l, -i, k)
          if (Parent(A,k)==i) then
            PrY = OKA2P(Genos(l,A),x,:) * PrP(:,1)
            if (SUM(PrY) > 0D0)  prob(x) = prob(x) / SUM(PrY)
          endif
          if (B>0) then
            if (Parent(B,k)==i) then
              PrY = OKA2P(Genos(l,B), x, :) * PrP(:,2)
              if (SUM(PrY) > 0D0)  prob(x) = prob(x) / SUM(PrY)  
            endif
          endif
        else if (B==-1) then  ! exclude all FS of A
          if (ALL(Xpr(3,:,l,-i,k)==1D0)) then   ! at initiate 
            prob(x) = AHWE(x, l)
          else if (AllIN) then
            prob(x) = XPr(2,x,l, -i, k)   ! GP only
          else
            PrY = PrP(:,1) * PrAF(:,x)
            if (SUM(PrY) > 0D0)  prob(x) = XPr(3,x,l, -i, k) / SUM(PrY)
          endif     
        else if (B==-4) then ! exclude both GPs & A
          if (ns(-i,k)==1 ) then  !  
            prob(x) = 1D0   ! AHWE(:,l)  !
          else
            PrY = OKA2P(Genos(l,A),x,:)* PrP(:,1)
            if (SUM(PrY) > 0D0)  prob(x) = XPr(1,x,l,-i,k) /SUM(PrY)
          endif     
        else if (B==-5) then ! exclude both GPs & A & FS of A
          if (nFS(A1) == ns(-i,k)) then  ! ns(-i,k)==1 
            prob(x) = 1D0   ! AHWE(:,l)  !
          else
            PrY = PrP(:,1) * PrAF(:,x)
            if (SUM(PrY) > 0D0)  prob(x) = XPr(1,x,l,-i,k) / SUM(PrY)   ! *AHWE(x,l)
          endif
        else
          call Erstop("ParProb: invalid B!", .TRUE.)
        endif
      enddo
    endif
  endif
endif
if (SUM(prob)>0D0 .and. .not. ALL(prob==1D0)) then
  prob = prob/SUM(prob)
endif
 
if (ANY(prob< 0D0) .or. ANY(prob/=prob) .or. ANY(prob>1.01D0)) then
  call Rprint( "Indiv, k,A,B: ", (/i, k, A, B/), (/0.0D0/), "INT") 
  if (A1/=0)  call Rprint("A1, nFS, ns, par: ", (/A1, nFS(A1), ns(-i,k), Parent(A,:)/), (/0.0D0/), "INT") 
  call Rprint("prob: ", (/0/), prob, "DBL")  
  call Erstop("Invalid ParProb!", .TRUE.)
endif

end subroutine ParProb

! #####################################################################

subroutine OffProb(l,i,k, prob)
use Global
implicit none

integer, intent(IN) :: l, i, k
double precision, intent(OUT) :: prob(3)

if (i > 0) then
  prob = OcA(:,Genos(l,i))
else if (i < 0) then
  prob = XPr(1,:,l,-i,k)
else
  prob = 1D0
endif

end subroutine OffProb

! #####################################################################

subroutine Connected(A, kA, B, kB, Con)
use Global
implicit none

integer, intent(IN) :: A, kA, B, kB
logical, intent(OUT) :: Con
integer :: i, j, m, nA, nB, AA(maxSibsize), BB(maxSibsize), n

Con = .FALSE.
if (A==0 .or. B==0)  return

AA = 0
BB = 0
if (A>0) then
  nA = 1
  AA(1) = A
else
  nA = nS(-A,kA)
  AA(1:nA) = SibID(1:nA, -A, kA)
endif

if (B>0) then
  nB = 1
  BB(1) = B
else
  nB = nS(-B,kB)
  BB(1:nB) = SibID(1:nB, -B, kB)
endif

do j=1, nB
  do i=1, nA
    do m=1,2  
      if (Parent(AA(i), m) < 0) then
        if (Parent(AA(i),m) == Parent(BB(j),m)) then
          Con = .TRUE.
          return
        else if(ANY(GpID(:,-Parent(AA(i), m),m) == BB(j))) then
          Con = .FALSE.  ! TODO: update Uclust(). already conditioned on?
!          return
        else if (A<0 .and. m==kA) then
          if(ANY(GpID(:,-Parent(AA(i), m),m) < 0)) then
            do n=1,2
              if (GpID(n,-Parent(AA(i), m),m) == Parent(BB(j),n) .and. &
                Parent(BB(j),n)<0) then
                Con = .TRUE.
                return
              endif 
            enddo
          endif
        endif
      endif
      if (Parent(BB(j),m)<0) then
        if (ANY(GpID(:,-Parent(BB(j),m),m) == AA(i))) then
          Con = .FALSE.  ! TODO
!          return
        else if (B<0 .and. m==kB) then
          if(ANY(GpID(:,-Parent(BB(j),m),m) < 0)) then
            do n=1,2
              if (GpID(n,-Parent(BB(j), m),m) == Parent(AA(i),n) .and. &
                Parent(AA(i),n)<0) then
                Con = .TRUE.
                return
              endif 
            enddo
          endif
        endif
      endif
    enddo
  enddo
enddo

end subroutine Connected

! #####################################################################

subroutine GetAncest(A, kIN, Anc)
use Global
implicit none

integer, intent(IN) :: A, kIN
integer, intent(OUT) :: Anc(2, mxA)  ! 32 = 5 generations
integer :: m, j, i, k, Par(2)

! Anc(1,:)  female ancestors
! Anc(2,:)  male ancestors
! Anc(:,2)  parents
! Anc(:,3)  mat. gp
! Anc(:,4)  pat. gp

Anc = 0
if (A==0) return

k = kIN
if (A > 0) then  ! real indiv
  if (kIN < 1 .or. kIN > 2)  k = 1
  Anc(k,1) = A
else !if (A < 0) then  ! dummy indiv
  if (kIN < 1 .or. kIN > 2) then
    call Erstop("getAncest: k must be 1 or 2 if A<0", .TRUE.)
  else
    Anc(k,2) = A  !! 
  endif
endif

Par = getPar(A,k)
if (ALL(Par == 0))  return
if (A > 0)  Anc(:, 2) = Par

do j = 2, mxA/2  
  do m = 1, 2
    i = 2 * (j-1) + m
    Anc(:,i) = getPar(Anc(m,j), m)
  enddo
  if (j==2 .and. ALL(Anc(:, 3:4) == 0))  return
  if (j==4 .and. ALL(Anc(:, 5:8) == 0))  return
  if (j==8 .and. ALL(Anc(:, 9:16) == 0))  return
  if (j==16 .and. ALL(Anc(:, 17:32) == 0))  return                                              
enddo

if ((A>0 .and. ANY(Anc(:, 2:mxA)==A)) .or. (A<0 .and. ANY(Anc(k,3:mxA)==A))) then
  call Rprint( "Female ancestors: ", Anc(1,1:8), (/0.0D0/), "INT")
  call Rprint( "Male ancestors: ", Anc(2,1:8), (/0.0D0/), "INT")
  call Erstop("An individual is its own ancestor! Need more birth years or better SNP data", .FALSE.)
endif

end subroutine GetAncest

! #####################################################################

subroutine ChkAncest(A, kA, B, kB, OK)  ! check that B is not an ancestor of A
use Global
implicit none

integer, intent(IN) :: A, kA, B, kB
logical, intent(OUT) :: OK
integer :: AncA(2, mxA), j

OK = .TRUE.
if (A==0 .or. B==0) then
  return 
else if (A==B .and. kA==kB) then
  OK = .FALSE.
  return
endif
          
call GetAncest(A, kA, AncA)

if (B > 0) then
  if (ANY(AncA == B))  OK = .FALSE.  
else if (kB == 1 .or. kB==2) then
  if (ANY(AncA(kB,:) == B))  OK = .FALSE.
  if (hermaphrodites/=0 .and. DumClone(-B,kB)/=0) then
    if (ANY(AncA(3-kB,:) == -DumClone(-B,kB)))  OK = .FALSE.
  endif
else
  call ErStop("ChkAncest: kB must be 1 or 2 if B<0", .TRUE.)
endif

if (OK .and. B < 0 .and. A < 0) then   ! check 1 extra generation
  if (ns(-B,kB)==0)  return
  do j=1, ns(-B,kB)
    if (ANY(AncA == SibID(j,-B,kB))) then
      OK = .FALSE.
      exit
    endif
  enddo
endif

end subroutine ChkAncest

! #####################################################################

subroutine CalcParentLLR(LLR_parent, LLR_GP)
! Calc parental LLR (vs next most likely relationship)
use Global
implicit none

double precision, intent(OUT) :: LLR_Parent(nInd,3), LLR_GP(3, nInd/2, 2)
integer :: i, CurPar(2), k, m, nonG(6), CurGP(2), s, g, chunk_printdot(40), z
double precision :: LLtmp(2,2,2), LLg(7), LLa(7)
logical :: AllSibsSelfed(2), NoGP, FSM  

call UpdateAllProbs()
LLR_parent = missing
LLR_GP = missing

if (quiet<1)  call rprint_progbar_header()
chunk_printdot = mk_seq(nInd+nC(1)+nc(2),40)

do i=1, nInd
  if (MODULO(i,100)==0) call rchkusr()
  if (quiet<1 .and. any(chunk_printdot==i))  call rprint_progbar_dot()  
  if (Parent(i,1)==0 .and. Parent(i,2)==0) cycle
  
  CurPar = Parent(i,:)
  do k=1,2  ! remove i from sibgroup
    call setParTmp(i, Sex(i), 0, k)
  enddo

  AllSibsSelfed = .FALSE.
  NoGP = .FALSE.
  if (hermaphrodites/=0 .and. all(curPar <0)) then
    if (DumClone(-curPar(1),1) == -curpar(2)) then
      do k=1,2
        if (ns(-curPar(k),k)==0) then
          AllSibsSelfed(k) = .TRUE.
        else if (all(SelfedIndiv(SibID(1:ns(-curPar(k),k),-curPar(k),k)))) then
          AllSibsSelfed(k) = .TRUE.
        endif
      enddo
      if (all(GpID(:,-curpar(1),1)==0))  NoGP = .TRUE.
    endif
  endif

  LLtmp = missing
  do m=1,2  ! m=1: no opp. sex parent;  m=2: with opp. sex parent
    if (m==2 .and. (CurPar(1)==0 .or. CurPar(2)==0)) cycle
    do k=1,2  ! mother, father
      if (CurPar(k)==0) cycle
      if (k==2 .and. Complx==0 .and. .not. any(curPar==0))  cycle
      if (AllSibsSelfed(k) .and. .not. (k==1 .and. AllSibsSelfed(2) .and. .not. NoGP))  cycle
     if (m==2) then  ! temp. assign parent 3-k    .and. .not. AllSibsSelfed(3-k)
        call setParTmp(i, Sex(i), curPar(3-k), 3-k)
      endif
      
      if (CurPar(k) > 0) then
        call CheckPair(i, CurPar(k), k, 7, LLg, LLa)  
        LLtmp(1,k,m) = LLg(1)
        LLtmp(2,k,m) = MaxLL(LLg(2:7))
      else if (CurPar(k) < 0) then
        call CheckAdd(i, -CurPar(k), k, 7, LLg, LLa)
        if (m==1) LLg(2) = 333D0   ! FS does not count here
        LLtmp(1,k,m) =  MaxLL(LLg(2:3))
        LLtmp(2,k,m) =  MaxLL((/LLg(1), LLg(4:7)/))
      endif
      
      if (m==2 .and. k==1) then
        call setParTmp(i, Sex(i), 0, 3-k)
      endif
    enddo
  enddo 
  
  do k = 1,2   ! restore
    if (Parent(i,k)/=curPar(k))  call setParTmp(i, Sex(i), CurPar(k), k) 
  enddo
  ! curPar(1) mostly restored in round m=2, k=2

  if (all(AllSibsSelfed) .and. NoGP) then
    call IsSelfed(i, .FALSE., LLR_parent(i,3))
  else 
    call ParLLtoLR(LLtmp, LLR_parent(i,:), AllSibsSelfed)
  endif  
enddo

nonG = (/1,2,3,5,6,7/)
z = nInd     
do k = 1,2  
  if (nC(k)==0)  cycle
  do s=1, nC(k)
    if (MODULO(s,10)==0) call rchkusr()
    z = z+1                                               
    if (quiet<1 .and. any(chunk_printdot==z)) call rprint_progbar_dot()
    CurGP = GpID(:, s, k)
    do g=1,2
      call setParTmp(-s, k, 0, g)
    enddo      
    LLtmp = missing
    do m=1,2
      if (m==2 .and. (CurGP(1)==0 .or. CurGP(2)==0)) cycle
      do g=1,2
        if (CurGP(g) == 0) cycle
        if (g==2 .and. Complx==0 .and. .not. any(curGP==0))  cycle
        if (m==2) then  ! temp. assign GP 3-g
          call setParTmp(-s, k, CurGP(3-g), 3-g)
        endif
        
        if (curGP(g) > 0) then
          call checkAdd(CurGP(g),s,k, 7, LLg, LLa)  ! B=GP + CurGP(m)_7
        else if (curGP(g) < 0) then
          call checkMerge(s, -CurGP(g), k, g, 4, LLg, LLa, FSM)   !TODO: use FSM?
          if (m==1) then
            call PairUA(-s, CurGP(g), k, g, LLg(4))  
          endif
        endif       
        LLtmp(1,g,m) = LLg(4)
        LLtmp(2,g,m) = MaxLL(LLg(nonG))
        if (m==2 .and. g==1) then  ! reset to 0
          call setParTmp(-s, k, 0, 3-g)
        endif 
      enddo
    enddo

    do g=1,2
      call setParTmp(-s, k, CurGP(g), g)  ! restore
    enddo
    
    call ParLLtoLR(LLtmp, LLR_GP(:,s,k), (/.FALSE.,.FALSE./))
  enddo
enddo

end subroutine CalcParentLLR

! ######################################################################

subroutine ParLLtoLR(LLtmp, LLR, AllSibsSelfed)
use Global
implicit none

double precision, intent(IN) :: LLtmp(2,2,2)
double precision, intent(OUT) :: LLR(3)
logical, intent(IN) :: AllSibsSelfed(2)
integer :: k
double precision :: LLX(2)

! LLtmp dims: focal/not-focal ; dam/sire ; without/with other parent 
LLR = missing
if (Complx > 0) then
  do k=1,2  ! max with - max w/o 
    if (LLtmp(1,k,1) < 0) then
      LLR(k) = LLtmp(1,k,1) - LLtmp(2,k,1)
    else
      LLR(k) = LLtmp(1,k,1)  ! something wrong  / AllSibsSelfed  / monogamous
    endif
  enddo
endif

if (hermaphrodites/=0) then
  do k=1,2
    if (AllSibsSelfed(k)) then
      LLR(k) = LLR(3-k)
    endif
  enddo
endif

do k=1,2
  if (LLtmp(1,k,2) < 0) then
    LLX(k) = LLtmp(1,k,2) -MaxLL((/LLtmp(2,k,2), LLtmp(:,k,1)/))
  else if (Complx==0 .and. LLtmp(1,k,1) < 0) then   ! single grandparent
    LLX(k) = LLtmp(1,k,1) - LLtmp(2,k,1)
  else
    LLX(k) = LLtmp(1,k,2)
  endif
enddo
LLR(3) = MINVAL(LLX)

end subroutine ParLLtoLR

! ######################################################################

subroutine setPar(A, kA, P, kP)    ! Assigns parent P to A, incl. sex & age update
use Global
implicit none

integer, intent(IN) :: A, kA, P, kP
integer :: curPar(2)

if (A==0)  return

curPar = getPar(A, kA)
if (curPar(kP) /= P)  then
  call setParTmp(A, kA, P, kP)
  call SetEstBY(curPar(kP), kP)
endif

call UpdateLL(P, kP)
call UpdateLL(curPar(3-kP), 3-kP)
call UpdateLL(A,kA)

call SetEstBY(A, kA)
call SetEstBY(P, kP)

if (P > 0) then
  if (Sex(P) == 3)   Sex(P) = kP
endif

if (A>0 .and. P/=0)  ToCheck(A) = .TRUE.
if (A<0 .and. P/=0) then
  if (ns(-A,kA) <= 3) then 
    ToCheck(SibID(1:ns(-A,kA),-A,kA)) = .TRUE.
  endif  
endif

if (hermaphrodites/=0) then
  if (A>0) then
    call CheckSelfed(A, sex(A))   ! sets Selfed(P,kP) if P<0
  else if (DumClone(-A,kA)/=0) then
    call setParTmp(-DumClone(-A, kA), 3-kA, P, kP)
    call SetEstBY(-DumClone(-A, kA), 3-kA) 
  endif
endif

end subroutine setPar

! ######################################################################

subroutine setParTmp(A, kA, P, kP)   ! Temporary assigns parent P to A
use Global
use CalcLik           
implicit none

integer, intent(IN) :: A, kA, P, kP
integer :: curPar(2), nOffP, OffP(maxSibSize), sxOffP(maxSibSize), i

curPar = getPar(A, kA)

if (A==0)  call Erstop("SetParTmp: A=0", .TRUE.)
if (kP/=1 .and. kP/=2)  call Erstop("SetParTmp: kP must be 1 or 2", .TRUE.)
if (A<0 .and. kA/=1 .and. kA/=2)  call Erstop("SetParTmp: kA must be 1 or 2 if A<0", .TRUE.)
if (P==0 .and. curPar(kP)==0)  return

if (P < 0) then
  if (-P > nC(kP)) then
    call Erstop("setParTmp: Sibship number out of bounds", .TRUE.)
  endif
endif

! remove old par
if (A > 0) then
  if (curPar(kP) > 0) then
    call RemoveFS(A)
    Parent(A,kP) = 0
    call CalcLind(A)

  else if (curPar(kP) < 0) then
    call RemoveSib(A, -curPar(kP), kP)   ! NOTE: doesn't drop sibship if singleton w/o GP
  endif
  
  if (P > 0 .and. curPar(3-kP)/=0) then  ! check for FS
    nOffP = 0
    if (ANY(Parent(:,kP) == P)) then   ! P already has some offspring
      call getOff(P, kP, .FALSE., nOffP, OffP, sxOffP)   ! currently only >0 FS considered
    endif
    Parent(A, kP) = P 
    call CalcLind(A)                
    if (nOffP > 0) then
      do i=1, nOffP
        if (Parent(A,3-kP) == Parent(OffP(i), 3-kP)) then
          call MakeFS(A, OffP(i))
          call CalcLind(OffP(i))
        endif
      enddo
    endif

  else if (P > 0) then
    Parent(A, kP) = P
    call CalcLind(A)                

  else if (P < 0) then
    call DoAdd(A, -P, kP)     
  endif

  if (curPar(3-kP) < 0) then
    call CalcCLL(-curPar(3-kP), 3-kP)
    if (ns(-curpar(3-kP),3-kP) > 0) then
      do i=1, ns(-curpar(3-kP),3-kP)
        call CalcLind(SibID(i, -curpar(3-kP),3-kP)) 
      enddo
    endif
    if (P < 0) then
      call CalcCLL(-P, kP)
      do i=1, ns(-P,kP)
        call CalcLind(SibID(i,-P,kP))
      enddo
    endif
  endif

else  ! A < 0
  GpID(kP, -A, kA) = P
  call CalcCLL(-A,kA)
  if (ns(-A,kA)>0) then
    do i=1, ns(-A,kA)
      call CalcLind(SibID(i,-A,kA))
    enddo
  endif
endif

end subroutine setParTmp

! #####################################################################

subroutine DoAdd(A, s, k)
use Global
use CalcLik           
implicit none

integer, intent(IN) :: A, s, k
integer :: i, u

if (nS(s,k) +1 >= maxSibSize) then
  call Erstop("Reached Maximum Sibship Size (number of offspring per parent), please increase 'MaxSibshipSize'", .FALSE.)
endif

Parent(A, k) = -s
if (.not. ANY(SibID(1:nS(s,k),s,k)==A)) then
  SibID(nS(s,k)+1, s, k) = A  ! add A to sibship
  nS(s,k) = nS(s,k) + 1
endif

do u=1, nS(s,k)  ! check for FS   
  i = SibID(u,s,k)
  if (i==A .or. nFS(i)==0) cycle 
  if (Parent(A, 3-k)/=0 .and. Parent(A, 3-k)==Parent(i, 3-k)) then
    call MakeFS(A, i)
  endif
enddo

call calcCLL(s,k)
if (Parent(A,3-k) < 0) call CalcCLL(-Parent(A,3-k), 3-k)
do u=1,ns(s,k)
  call CalcLind(SibID(u,s,k))
enddo

end subroutine DoAdd

! #####################################################################

subroutine RemoveSib(A, s, k)  ! removes individual A from sibship s. Does NOT drop sibship if ns=0
use Global
use CalcLik           
implicit none

integer, intent(IN) :: A, s, k
integer :: u, h, reps

call RemoveFS(A)

do u=1,ns(s,k)
  if (SibID(u,s,k)==A) then
    if (u<ns(s,k)) then  ! shift sibs
      do h=u, nS(s, k)-1  ! drop HS
        SibID(h, s, k) = SibID(h+1, s, k)
      enddo
    endif
    SibID(nS(s,k), s, k) = 0
    nS(s,k) = nS(s,k) -1
    exit
  endif
enddo

Parent(A, k) = 0

reps = 1
if (any(parent(SibID(1:ns(s,k),s,k),3-k)==A))  reps=2  ! removing inbreeding loop
do h=1,reps
  call CalcCLL(s,k)
  if (Parent(A,3-k) < 0) call CalcCLL(-Parent(A,3-k), 3-k)
  if (ns(s,k) >0) then
    do u=1,ns(s,k)
      call CalcLind(SibID(u,s,k))
    enddo
  endif
  call CalcLind(A)
enddo

end subroutine RemoveSib

! ######################################################################

subroutine UpdateLL(A, k)   ! update LL of A & if A<0 of its mates
use Global
use CalcLik           
implicit none

integer, intent(IN) :: A, k
integer :: Mates(maxSibSize), x, i, j, nOff, sxOff(maxSibSize), Off(maxSibSize)
double precision :: XPR_x(3, nSnp), tol
logical :: OK

if (A==0)  return

if (A >0) then
  call CalcLind(A)
  return
endif

if (ns(-A,k)==0)  return

tol = 0.01  ! What tolerance?                             
Mates = 0
Mates(1:ns(-A,k)) = Parent(SibID(1:ns(-A,k), -A,k), 3-k)

do x=1,10
  OK = .TRUE.
  do i=1, ns(-A,k)
    if (nFS(SibID(i,-A,k)) == 0)  cycle
    if (Mates(i) < 0) then
      if (ANY(Mates == GpID(3-k, -Mates(i),3-k)) .and. GpID(3-k, -Mates(i),3-k) < 0) then
        call CalcCLL(-GpID(3-k, -Mates(i),3-k), 3-k)   ! Used by UseEE
      endif
!      if (ns(-Mates(i),3-k) <= 50 .or. nFS(SibID(i,-A,k)) >= ns(-Mates(i),3-k)/5) then
        XPR_x = XPr(3,:,:,-Mates(i), 3-k)
        call CalcCLL(-Mates(i), 3-k)
        if (any(abs(XPR_x - XPr(3,:,:,-Mates(i), 3-k)) > tol))  OK = .FALSE.
!      endif
    endif
  enddo
  XPR_x = XPr(3,:,:,-A,k)
  call CalcCLL(-A,k)
  if (OK .and. all(abs(XPR_x - XPr(3,:,:,-A,k)) < tol))  exit    ! What tolerance?
enddo

call getOff(A, k, .TRUE., nOff, Off, sxOff)  ! includes dummy offspring
do i=1, nOff
  if (nOff==0)  exit                  
  if (Off(i) > 0) then
    if (nFS(Off(i)) > 0 .and. Mates(i) < 0) then
!      if (ns(-Mates(i),3-k) <= 50 .or. nFS(Off(i)) >= ns(-Mates(i),3-k)/5) then
        if (ns(-Mates(i), 3-k) > 0) then
          do j=1, ns(-Mates(i), 3-k)
            call CalcLind(SibID(j, -Mates(i), 3-k))
          enddo
        endif
!      endif
    endif
    call CalcLind(Off(i))
  else
    call CalcCLL(-Off(i), sxOff(i))
  endif
enddo

end subroutine UpdateLL

! #####################################################################

subroutine MakeFS(A, B)
use Global
implicit none

integer, intent(IN) :: A,B
integer :: x, i, j, Ai, Bj

if (nFS(A)>0) then
  Ai = A
else
  Ai = FSID(maxSibSize+1, A)
endif
if (nFS(B)>0) then
  Bj = B
else
  Bj = FSID(maxSibSize+1, B)
endif

if (ANY(FSID(1:nFS(Ai),Ai)==B) .or. ANY(FSID(1:nFS(Bj),Bj)==A)) then
  return ! already are FS.
endif

i = MIN(Ai,Bj)
j = MAX(Ai,Bj)
do x=1, nFS(j)   
  FSID(nFS(i)+x, i) = FSID(x, j)
  FSID(maxSibSize+1, FSID(x,j)) = i
enddo
nFS(i) = nFS(i) + nFS(j)
FSID(maxSibSize+1,i) = i    ! 'primary' sib
FSID(:,j) = 0
FSID(1,j) = j
FSID(maxSibSize+1,j) = i
nFS(j) = 0

end subroutine MakeFS

! ######################################################################

subroutine RemoveFS(A)
use Global
implicit none

integer, intent(IN) :: A
integer :: op, np, i, j

if (nFS(A) == 1) then
  return
else if (nFS(A) > 1) then
  op = A
  np = MINVAL(FSID(1:nFS(A), A), MASK=(FSID(1:nFS(A),A)/=A)) 
else !if (nFS(A) == 0) then
  op = FSID(maxSibSize+1, A)  ! 'primary' sib
  np = op
endif

i = 2  ! 1st one stays op
do j=1, nFS(op)
  if (FSID(j,op)==A) then
    FSID(j,op) = 0  ! if nFS(A)=0 .and. nFS(op)=2
    cycle
  endif
  if (FSID(j,op)==np) cycle
  FSID(i, np) = FSID(j, op)
  if (op /= np) then
    FSID(maxSibSize+1, FSID(j, op)) = np
  endif
  i = i+1
enddo

nFS(np) = nFS(op)-1
FSID(maxSibSize+1, np) = np
nFS(A) = 1
FSID(:,A) = 0
FSID(1,A) = A
FSID(maxSibSize+1, A) = A

end subroutine RemoveFS

! ######################################################################

subroutine CheckSelfed(A, kA)
use Global
implicit none

integer, intent(IN) :: A, kA
integer :: AS(2), m, j, Aj, x
double precision :: LRself

if (hermaphrodites==0)  return

if (A > 0) then
  call IsSelfed(A, .FALSE., LRself)
  if (all(Parent(A,:) == 0)) then
    if (LRself > TA) then   ! threshold? be consistent with selectparent()
      SelfedIndiv(A) = .TRUE.  
    ! else leave as is?
    endif
    return
  else if (any(Parent(A,:) > 0)) then  
    if (Parent(A,1)==Parent(A,2)) then
      if (LRself < 5*TF) then    
        call Rprint("Indiv causing error:", (/A, Parent(A,:)/), (/0D0/), "INT")
        call Rprint ("LRself:", (/0/), (/LRself/), "DBL")
        call Erstop("CheckSelfed: dam = sire, but LRself < 5*TF", .TRUE.)
      else
        SelfedIndiv(A) = .TRUE.
      endif
    else if (all(Parent(A,:)/=0) .and. Parent(A,1)/=Parent(A,2)) then
      if (LRself > TA) then   
        call Rprint("Indiv causing error:", (/A, Parent(A,:)/), (/0D0/), "INT")
        call Rprint ("LRself:", (/0/), (/LRself/), "DBL")
        call Erstop("CheckSelfed: dam /= sire, but LRself > TA", .TRUE.)
      else
        SelfedIndiv(A) = .FALSE.  
      endif
    ! else assignment in progress, leave as is?
    endif
  endif
  if (all(Parent(A,:) >= 0))  return
endif

AS = 0
if (A > 0)  AS = -Parent(A,:)
if (A < 0)  AS(kA) = -A

do m=1,2
  if (AS(m) <= 0)  cycle
  do j=1, ns(AS(m),m)
    Aj = SibID(j,AS(m),m)
    if (nFS(Aj)==0)  cycle
    call IsSelfed(Aj, .TRUE., LRself)
    if (LRself > TA) then
      if (Parent(Aj,3-m) < 0) then
        do x = 1, nFS(Aj)
          SelfedIndiv(FSID(x,Aj)) = .TRUE.
        enddo
        DumClone(AS(m), m) = -Parent(Aj,3-m)
        DumClone(-Parent(Aj,3-m), 3-m) = AS(m)
      else if (Parent(Aj,3-m) > 0) then
      call Erstop("SetPar: parents incompatible with LRself > TA", .TRUE.)
    else   ! transitory during assignment only (?)
       do x = 1, nFS(Aj)
         SelfedIndiv(FSID(x,Aj)) = .TRUE.
       enddo
      endif
    endif
  enddo
enddo
  
end subroutine CheckSelfed

! #####################################################################

! @@@@   INPUT & PRECALC PROB.   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! #####################################################################

subroutine PrepAgeData(nrows_AP, AP_IN, BYrange)
use Global
implicit none

integer, intent(IN) :: nrows_AP
double precision, intent(IN) :: AP_IN(nrows_AP,5)
integer, intent(INOUT) :: BYrange(nInd, 2)
integer :: i, BYLast, x, y

!===  determine first & last possible birth year  ==============   
BYzero = MINVAL(BY, MASK=BY>=0) -1   ! defaults to HUGE(ARRAY) 
BYlast = MAXVAL(BY, MASK=BY>=0)      ! defaults to -HUGE(ARRAY)
if (ANY(BYrange >= 0)) then
  BYzero = MIN(BYzero, MINVAL(BYrange(:,1), MASK=BYrange(:,1)>=0) -1)
  BYlast = MAX(BYlast, MAXVAL(BYrange(:,2), MASK=BYrange(:,2)>=0))
endif

!===  determine MaxAgePO  ==============   
maxAgePO = 1  ! maximum PO age difference (needed for dummy parents)
do y = 2, nrows_AP
  if (ANY(AP_IN(y, 1:2)>TINY(0D0))) then
    maxAgePO = y - 1  ! first row is agediff of 0
  endif
enddo

if (BYzero < 99999) then
  BYzero = BYzero - MaxAgePO !+1  
  if ((BYlast - BYzero) < 2*MaxAgePO .or. BYlast==0) then
    BYzero = BYzero - MaxAgePO    ! dummy parents + real grandparents w unknown BY
  endif  
  nYears = BYlast - BYzero   ! defines nYears!
else   ! all birth years unknown
  BYzero = 0
  nYears = 2*MaxAgePO
endif


!===  shift BY  ==============
WHERE (BY >=0) BY = BY - BYzero  
do x=1,2
  WHERE (BYRange(:,x) >=0) BYRange(:,x) = BYRange(:,x) - BYzero
enddo
WHERE(YearLast <= 0)  YearLast = 999
WHERE (YearLast /= 999)  YearLast = YearLast - BYzero       

! min/max BY based on YearLast
do i=1, nInd
  if (YearLast(i)==999)  cycle
  if (BYrange(i,1) < 0)  BYrange(i,1) = MAX(YearLast(i) - MaxAgePO, 1)
  if (BYrange(i,2) < 0)  BYrange(i,2) = MIN(YearLast(i), nYears)
enddo

WHERE (BYRange(:,1) <0) BYrange(:,1) = 1
WHERE (BYRange(:,2) <0) BYrange(:,2) = nYears  

!===  Initiate indiv BY prob distr  ==============
allocate(IndBY(1:nYears, nInd, 5))  ! year - indiv - own/wo/w dummy off+par
IndBY = LOG10(1.0D0/nYears)
do i=1, nInd   
  if (BY(i) >0) then 
    IndBY(:, i, :) = LOG10(zero)
    IndBY(BY(i), i, :) = zero  
  else if (BYrange(i,1)/=1 .or. BYrange(i,2)/=nYears) then
    if (BYrange(i,1) < 1)  call ErStop('BY.min < 1', .TRUE.)
    if (BYrange(i,2) > nYears)  call ErStop('BY.max > nYears', .TRUE.)
    IndBY(:, i, :) = LOG10(zero)
    IndBY(BYrange(i,1) : BYrange(i,2), i, :) = LOG10(1.0D0/(BYrange(i,2) - BYrange(i,1) +1))
  endif
enddo

allocate(DumBY(1:nYears, nInd/2, 2,5)) 
DumBY = 0D0


!===  initiate AgePrior array  ==============
call mk_APextra()  
! incl.  allocate(AgePriorA(-MaxAgePO : nYears, 5, 3)); incl.  calc for GP & AU
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! AgePriorA D2 + D3:
!  1    2     3
!  M   MGM   PGM
!  P   MGP   PGP
! FS   MFA   PFA
! MS  MMHA  PMHA
! PS  MPHA  PPHA
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~


contains

  ! ===================
  subroutine mk_APextra()
  integer, allocatable :: AgeDifCnt(:)
  double precision, allocatable :: AgeDifProp(:), BYP1(:), BYP2(:), APtmp(:,:), &
   APz(:), AgeDifCnt_dbl(:)
  integer :: a,i,j,r, r_ik, r_kj
  logical :: discrete_generations
  
  discrete_generations = .FALSE.
  if (COUNT(AP_IN(:,1)>0d0) == 1 .and. COUNT(AP_IN(:,2)>0d0) == 1) then
    if (MAXLOC(AP_IN(:,1), DIM=1) == MAXLOC(AP_IN(:,2), DIM=1))   discrete_generations = .TRUE.
    ! NOTE: No check if AP for sibling age difference is consistent; for that use R pkg
  endif
    
  ! ageprior array:
  ! D1: age difference
  ! D2: mum/dad/FS/MS/PS
  ! D3: of self / of mum / of dad
  allocate(AgePriorA(-MaxAgePO : MaxAgePO*2, 5, 3)) 
  AgePriorA = 0.0D0
  
  ! 1st slice of D3 = input ageprior (read from file)
  forall (r=1:5)  AgePriorA(0:MaxAgePO, r, 1) = AP_IN(1:(MaxAgePO+1), r) 
  AgePriorA(0, 1:2, 1) = 0.0D0   ! PO CANNOT have age diff of 0. (overrule input)
  do r = 3,5  ! mirror FS,MS,PS
    forall (a=1:MaxAgePO)  AgePriorA(-a, r, 1) = AgePriorA(a, r, 1)
  enddo
  
  if (ALL(AP_IN(:,1:2) < 0.001)) then  ! AP_IN with single row for sibs only
    AgePriorA(1:MaxAgePO, 1:2, 1) = 1.0D0
  endif
  
  if (discrete_generations) then
    AgePriorA(2*MaxAgePO, 1:2, 2:3) = 1.0D0  ! GPs
    AgePriorA(MaxAgePO,   3:5, 2:3) = 1.0D0  ! avuncular 

  else
    ! distribution of age differences among all sample pairs (reference for scaling)
    allocate(AgeDifCnt(0:MaxAgePO*2))
    AgeDifCnt = 0
    ! known birth years
    do i=1, nInd
      do j= i+1, nInd
        if (BY(j) < 0)  cycle
        a = ABS(BY(i) - BY(j))
        if (a > 2*MaxAgePO)  cycle
        AgeDifCnt(a) = AgeDifCnt(a) +1
      enddo
    enddo
    ! add pairs with one/both unknown birth year
    allocate(AgeDifCnt_dbl(0:MaxAgePO*2))
    AgeDifCnt_dbl = dble(AgeDifCnt)
    do i=1, nInd
      do j= i+1, nInd
        if (BY(i) >0 .and. BY(j)>0)  cycle  ! already done
        AgeDifCnt_dbl = AgeDifCnt_dbl + get_agedif(i,j)
      enddo
    enddo     
        
    ! counts --> proportions
    allocate(AgeDifProp(-MaxAgePO:MaxAgePO*2))
    AgeDifProp = 0.0D0
    AgeDifProp(0:MaxAgePO*2) = AgeDifCnt_dbl / sum(AgeDifCnt_dbl)
    ! flatten agedifprop, to deal with poorly sampled age differences that are possible for GP-GO
    WHERE(AgeDifProp < TINY(0d0))  AgeDifProp = TINY(0d0)  ! avoid division by zero errors on some compilers
    do a=1, MaxAgePO*2  ! only up to max GP-GO age difference
      if (AgeDifProp(a) < 0.001)  AgeDifProp(a) = 0.001d0  ! pretend at least 1k pairs for each age dif
    enddo
    AgeDifProp = AgeDifProp / sum(AgeDifProp)      ! rescale to sum to unity again
    ! mirror to negative age differences (for avuncular pairs)
    do a=1, MaxAgePO
      AgeDifProp(-a) = AgeDifProp(a)
    enddo
      
    ! calc GP & AU (for individual i born in year 0)
    allocate(APtmp(1:MaxAgePO, -MaxAgePO : MaxAgePO*2))
    allocate(BYP1(-MaxAgePO : MaxAgePO*2))
    allocate(BYP2(-MaxAgePO : MaxAgePO))
    allocate(APz(-MaxAgePO : MaxAgePO*2))
    do r_ik = 1,2  ! mum/dad k of i
      BYP1 = AgePriorA(:, r_ik, 1) * AgeDifProp  ! scaled ratio --> unscaled proportions (approx)
      do r_kj = 1,5
        APtmp = 0.0D0
        BYP2 = AgePriorA(-MaxAgePO : MaxAgePO, r_kj, 1) * AgeDifProp(-MaxAgePO : MaxAgePO)
        do a = 1, MaxAgePO   ! age difference i-k
          APtmp(a, (a-MaxAgePO) : (a+MaxAgePO) ) = BYP1(a) * BYP2
        enddo
        APz = SUM(APtmp, dim=1)
        APz = APz / SUM(APz)  ! scale to sum to unity
        AgePriorA(:, r_kj, r_ik+1) = APz / AgeDifProp  ! scale relative to unrelated agedif distribution
      enddo
    enddo
    WHERE (AgePriorA > HUGE(0.0d0) .or. AgePriorA /= AgePriorA)  AgePriorA = 0.0d0
    ! 2nd check is for NaN (due to division by zero)
    
    deallocate(AgeDifCnt, AgeDifProp, BYP1, BYP2, APtmp, APz)
  endif
  
  end subroutine mk_APextra

  
end subroutine PrepAgeData

! ######################################################################

subroutine EstBYrange(A, k, MCI)  
use Global      
implicit none

integer, intent(IN) :: A, k
integer, intent(OUT) :: MCI(3)  ! mode - 95% lower - 95% upper
integer :: y, mx, CI(2)
double precision :: DBYP(nYears), cumProp, dd(nYears)

MCI = -999
if (A == 0) return
call getEstBY(A,k, 5, DBYP) ! all contributions from relatives
DBYP = 10**DBYP
mx = MAXLOC(DBYP, DIM=1)
cumProp = DBYP(mx)
CI = mx
do y=1, nYears 
  if (cumProp > 0.95) then
    dd = 0D0
    WHERE (DBYP > 0.0D0)  dd = ABS(DBYP - DBYP(mx))
    if (ANY(dd > 0.001) .or. CI(1)==CI(2)) then  ! else: DBYP is flat between BYmin & BYmax
      MCI(1) = mx
    endif
    MCI(2:3) = CI
    WHERE(MCI>0)  MCI = MCI + BYzero
    exit
  endif        

  if (CI(1) > 1 .and. CI(2) < nYears) then
    if (DBYP(CI(1)-1) > DBYP(CI(2)+1)) then
      CI(1) = CI(1)-1
      cumProp = cumProp + DBYP(CI(1))
    else 
      CI(2) = CI(2)+1
      cumProp = cumProp + DBYP(CI(2))
    endif
  else if (CI(1) > 1) then
    CI(1) = CI(1)-1
    cumProp = cumProp + DBYP(CI(1))
  else if (CI(2) < nYears) then
    CI(2) = CI(2)+1
    cumProp = cumProp + DBYP(CI(2))
  endif
enddo

end subroutine EstBYrange

! #####################################################################

subroutine getGenerations(Gen)
use Global
implicit none

integer, intent(OUT) :: Gen(nInd)
integer :: i, GenPar(2, nInd), g, Prevgen(nInd), m

Gen = -9
GenPar = 0   ! generation number of dam, sire
Prevgen = 0   ! ids of individuals in generations up to g-1

do i=1,nInd
  if (Parent(i,1)==0 .and. Parent(i,2)==0) then
    Gen(i) = 0  ! founder
    Prevgen(i) = i
  endif
enddo

do g = 0, 1000
  do i=1, nInd
    if (Gen(i) >= 0)  cycle
    do m=1,2
      if (Parent(i,m)==0 .or. GenPar(m,i) > 0)  cycle
      if (ANY(Prevgen == Parent(i,m))) then
        GenPar(m,i) = g
      endif
    enddo
    
    if (ALL(GenPar(:,i) <= g)) then  ! including Parent(i,m)==0 --> GenPar(m,i)==0
      Gen(i) = g+1
      PrevGen(i) = i
    endif    
  enddo
  if (.not. ANY(Gen < 0))  exit   ! all done
enddo

end subroutine getGenerations

! #####################################################################

subroutine deallocall
use Global
use OHfun
implicit none

! same order as module Global; allocated in AllocArrays
if (allocated(ToCheck)) deallocate(ToCheck)
if (allocated(SelfedIndiv)) deallocate(SelfedIndiv)
if (allocated(IsNewSibship)) deallocate(IsNewSibship)
if (allocated(mtDif)) deallocate(mtDif)

if (allocated(Sex)) deallocate(Sex)
if (allocated(BY)) deallocate(BY)
if (allocated(nFS)) deallocate(nFS)
if (allocated(Mate))  deallocate(Mate)
if (allocated(YearLast)) deallocate(YearLast)                                      

if (allocated(Genos)) deallocate(Genos)
if (allocated(Parent)) deallocate(Parent)
if (allocated(nS)) deallocate(nS)
if (allocated(FSID)) deallocate(FSID)
if (allocated(DumMate))  deallocate(DumMate) 
if (allocated(DumClone))  deallocate(DumClone) 

if (allocated(SibID)) deallocate(SibID)
if (allocated(GpID)) deallocate(GpID)

if (allocated(Lind)) deallocate(Lind)
if (allocated(AHWE)) deallocate(AHWE)
if (allocated(OHWE)) deallocate(OHWE) 
if (allocated(CLL)) deallocate(CLL)

if (allocated(AKAP)) deallocate(AKAP)
if (allocated(OKAP)) deallocate(OKAP)
if (allocated(OKOP)) deallocate(OKOP)
if (allocated(PPO)) deallocate(PPO)
if (allocated(PHS)) deallocate(PHS)
if (allocated(PFS)) deallocate(PFS)
if (allocated(LindX)) deallocate(LindX)
if (allocated(IndBY)) deallocate(IndBY)
if (allocated(AgePriorA)) deallocate(AgePriorA)

if (allocated(DumP)) deallocate(DumP)     
if (allocated(DumBY)) deallocate(DumBY)    
if (allocated(XPr)) deallocate(XPr)

! module OHfun, clean up private arrays
call dealloc_QPOsparse()
                                  
end subroutine deallocall

! ####################################################################

! end.                                                      