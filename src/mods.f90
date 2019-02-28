! vim:fdm=marker
module mostCommon
integer :: np
integer, dimension(:), allocatable :: typ, NatomsPerType
real(8), dimension(:,:), allocatable :: ipos, pos, vel, acc, force, ftmp
real(8), dimension(:), allocatable :: pene
real(8) :: dt, potential, kinetic, virial(6), e0
real(8), dimension(3) :: box, box0, initBox !for init
end module mostCommon
!*****************************************************************************80
module groupStuff
integer :: Ngrp, grpBuff, Nhold, holdBuff
integer, dimension(:), allocatable :: grp, hold
character(8), dimension(:), allocatable :: groupName, holdName
character(2), dimension(:), allocatable :: bcTyp
real(8), dimension(:,:), allocatable :: bcParm
logical, dimension(:,:), allocatable :: bcDir
contains
function groupID( t )
!{{{
! return the group ID given a group Name
implicit none
integer  :: groupID
character(*), intent(in) :: t
integer :: i
groupID = 0  !default
do i = 1, Ngrp
  if (trim(t).eq.groupName(i)) then
    groupID = i
    Return
  endif
enddo
end function groupID !}}}
function holdID( t )
!{{{
! return the group ID given a group Name
implicit none
integer  :: holdID
character(*), intent(in) :: t
integer :: i
holdID = 0  !default
do i = 1, Nhold
  if (trim(t).eq.holdName(i)) then
    holdID = i
    Return
  endif
enddo
end function holdID !}}}
end module groupStuff
!*****************************************************************************80
module interpot
real(8) :: rcut, rcut2, PTrMin, srcut(3)
integer :: Ntypes, NPairs, PTsize
integer, dimension(:,:), allocatable :: PairMat
real(8), dimension(:), allocatable :: mass, charge, PairChargeProd, PTdr2, PTrc2
real(8), dimension(:,:), allocatable :: PairParams, PPairTable, FPairTable
character(8), dimension(:), allocatable :: typeName, PairName
logical :: potOuts, asympots, springpots
real(8), parameter :: PI = 3.141592653589793D+00
contains
function typeID( t )
! {{{
! returns the ID of the given type ! a negative given type value causes a negative ID
implicit none
integer  :: typeID
character(8), intent(in) :: t
integer :: i
typeID = 0  !default
if (trim(t).eq."test")   Return
do i = 1, Ntypes
  if (trim(t).eq.trim(typeName(i))) then
    typeID = i
    Return
  endif
enddo
end function typeID !}}}
end module interpot
!*****************************************************************************80
module control
character(len=80) :: Model
character(10) :: ModelInType, ModelOutType
integer(4) :: step_num, thread_num, proc_num, seed, ConfOuts, OutSteps, step
integer(4) :: KEfreq, Pfreq, Efreq
real(8) :: wtime, KEtarget, Ptarget(3), Etarget, KEparm, Pparm, Eparm, boxVel(3)
logical :: isokinetic, isobaric, isoenergetic, boxLoad
end module control
!*****************************************************************************80
module listData
logical :: ReNei
real(8) :: skin, skinsq
integer :: NCELL, MCELL(3), NeiMax, NlistUpdates, ReNeiFreq
integer, allocatable, dimension(:) :: List, HEAD
integer, allocatable, dimension(:,:) :: VNL
real(8), allocatable, dimension(:,:) :: dispList
end module listData

