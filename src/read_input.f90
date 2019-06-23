! vim:fdm=marker
subroutine read_input
! {{{
! read input file directives
use mostCommon
use control
use domType
use interpot
use listData
use groupStuff
implicit none
integer :: ln, wc, i, j, p1, p2
character(len=80) :: tmp, word, var
real(4) :: r4tmp
integer, allocatable, dimension(:) :: seeda

!!! Define Default values for input parameters before allocations
np = 1000 !number of particles
dt = 0.0001d0
step = 0
!!! Control Defaults
box(1:3) = 10.d0
step_num = 400;       OutSteps = 40;      ConfOuts = 0
seed = 123456789
Model='NA';           ModelInType='NA';   ModelOutType='NA'
Ntypes=0
rcut = PI / 2.0D+00 !default for the sine potentail
PTsize = 1000;        PTrMin = 0.25d0 !table size and minimum particle separation
NeiMax = 0
potOuts = .false. 
skin = -1.d0
ReNei = .true.;       ReNeiFreq = 1 !frequency to check for reneighbouring
isokinetic = .false.; isobaric = .false.; isoenergetic = .false.
KEtarget = 0.0;       Ptarget(1:3) = 0.0; Etarget = 0.0
KEparm = 0.0;         Pparm = 0.0;        Eparm = 0.0
KEfreq = 1;           Pfreq = 1;          Efreq = 1
! group and boundary conditions
boxLoad = .false.;    boxVel(1:3) = 0.d0
Ngrp = 0;             grpBuff = 100
Nhold = 0;            holdBuff = 100
asympots = .false.;   springpots = .false.
!!!
! Read directives from Input file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(6,'(A)') "#  type 'START' to begin specifying directives"
write(6,'(A)') "#  type 'help' for a list of directives."; call flush(6)
lop: do ! infinite loop to find the beginning of the directives
 read(5,'(A)',end=200,err=800) tmp; ln=ln+1
 word=s_get_word(1,tmp)
 select case (trim(word))
  case("START"); exit lop
  case("help");  call help
  case("END","DONE","EXIT","quit","abort"); STOP
 end select
enddo lop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!start readin the input file directives
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ln=0
do
 read(5,'(A)',end=200,err=800) tmp; ln=ln+1
 call left_of("#",tmp) ! strip any comments
 word = s_get_word(1,tmp) ! get the first word
 wc   = s_word_count(tmp)   ! word count
 if (wc.lt.1) cycle !no words
 select case (word)
  case ("help"); call help
! Particles
  case ("BoxSize"); ! assuming 'nd' was defined
          var=trim(s_get_word(2,tmp)); read(var,*) box(1)
          var=trim(s_get_word(3,tmp)); read(var,*) box(2)
          var=trim(s_get_word(4,tmp)); read(var,*) box(3)
  case ("GroupBufferSize");  var=trim(s_get_word(2,tmp)); read(var,*) grpBuff
  case ("HoldBufferSize");   var=trim(s_get_word(2,tmp)); read(var,*) holdBuff
  case ("Model"); Model=trim(s_get_word(2,tmp))
       if (wc.ge.3)  ModelInType=trim(s_get_word(3,tmp))
       call read_model
  case ("Seed"); var=trim(s_get_word(2,tmp)); read(var,*) seed
       call random_seed(size=p1); allocate( seeda(1:p1) ); seeda(1) = seed
       call random_seed(PUT=seeda); deallocate(seeda) 
  case ("NTotal"); var=trim(s_get_word(2,tmp)); read(var,*) np
       if (Ntypes.gt.1) then !list number of atoms of each type after total number of atoms
         do i=1, Ntypes
           var=trim(s_get_word(2+i,tmp)); read(var,*) NatomsPerType(i)
         enddo
       endif
       call read_model ! still need to setup model stuff
  case ("Velocity"); var=trim(s_get_word(2,tmp)); read(var,*) r4tmp
       if (wc.ge.3) then;  var=trim(s_get_word(3,tmp))
       call normalVelocityDistribution( r4tmp, var ); else
       call normalVelocityDistribution( r4tmp, "N/A" ); endif
! Neighbour List stuff
  case ("ReNeibour","ReNeighbor"); var=trim(s_get_word(2,tmp)); read(var,*) ReNei
       if (wc.ge.3) then; var=trim(s_get_word(3,tmp)); read(var,*) ReNeiFreq; endif
  case ("Skin");   var=trim(s_get_word(2,tmp)); read(var,*) skin
  case ("NeiMax"); var=trim(s_get_word(2,tmp)); read(var,*) NeiMax
! Atom types and interactions
  case ("Types"); var=trim(s_get_word(2,tmp)); read(var,*) NTypes
       allocate( mass(Ntypes), typeName(Ntypes), charge(Ntypes), NatomsPerType(0:Ntypes) )
       NatomsPerType(:) = 0; charge(:) = 0.0
       do i=1, Ntypes
          read(5,'(A)',end=200,err=800) tmp; ln=ln+1
          typeName(i) = trim(s_get_word(1,tmp))
          var=trim(s_get_word(2,tmp)); read(var,*) mass(i)
          if (s_word_count(tmp).gt.2) then 
            var=trim(s_get_word(3,tmp)); read(var,*) charge(i)
          endif
       enddo
  case ("WritePotentialTable"); potOuts=.true.
  case ("Pairs"); var=trim(s_get_word(2,tmp)); read(var,*) NPairs
       if (wc.ge.3) then; var=trim(s_get_word(3,tmp)); read(var,*) PTrMin; endif !{{{
       if (wc.ge.4) then; var=trim(s_get_word(4,tmp)); read(var,*) PTsize; endif
       if (NTypes.eq.0) then
          write(0,*) "ERROR: must first define particle Types before Pairs"
          call flush(0); STOP; endif
       if (allocated(PairMat)) then !these are already allocated, unallocate then and remake them again
         deallocate( PairMat, PairName, PairParams, PairChargeProd, PTdr2, PTrc2 )
       endif
       allocate( PairMat(NTypes,NTypes), PairName(0:NPairs), PairParams(4,NPairs), &
         &  PairChargeProd(NPairs), PTdr2(NPairs), PTrc2(NPairs) )
       PairMat=0; PairName="N/A"; PairParams=0.d0; PairChargeProd=0.d0; PTdr2=0.d0
       PTrc2 = 0.d0
       PairParams(1,:) = rcut !set default value
       do i=1, NPairs
          read(5,'(A)',end=200,err=800) tmp; ln=ln+1
          call left_of("#",tmp)
          var=trim(s_get_word(1,tmp)); p1 = (typeID(var))
          var=trim(s_get_word(2,tmp)); p2 = (typeID(var))
          var=s_get_word(3,tmp)
          PairMat(p1,p2) = i !tabulate the interaction ID
          if (trim(var).eq."sym") then
                PairMat(p2,p1) = i !tabulate the interaction ID
          else; asympots=.true.; endif !flag to say that an asymmetric interaction exists
          PairChargeProd(i) = charge(p1)*charge(p2)
          PairName(i) = trim(s_get_word(4,tmp))
          if (trim(PairName(i)).eq."spring") springpots = .true. !flag for strain calc
          do j=5, s_word_count(tmp)
             var=trim(s_get_word(j,tmp)); read(var,*) PairParams(j-4,i)
          enddo
       enddo
       write(6,'(A)') '# pair potential type matrix'
       write(6,*) "#     ",(typeName(j),j=1,NTypes)
       do i = 1, NTypes
          write(6,*) "# "//typeName(i), (PairName(PairMat(i,j)),j=1,NTypes)
       enddo !}}}
       call buildPairTables
! System Control
  case ("Isokinetic"); isokinetic=.true.
       var=trim(s_get_word(2,tmp)); read(var,*) KEtarget
       var=trim(s_get_word(3,tmp)); read(var,*) KEparm
       if (wc.ge.4) then; var=trim(s_get_word(4,tmp)); read(var,*) KEfreq; endif
  case ("Anisokinetic"); isokinetic=.false.
  case ("Isoenergetic"); isoenergetic=.true.
       !var=trim(s_get_word(2,tmp)); read(var,*) Etarget
       var=trim(s_get_word(2,tmp)); read(var,*) Eparm
       if (wc.ge.3) then; var=trim(s_get_word(3,tmp)); read(var,*) Efreq; endif
  case ("Anisoenergetic"); isoenergetic=.false.
  case ("Isobaric"); isobaric=.true.
       var=trim(s_get_word(2,tmp)); read(var,*) Ptarget(1)
       var=trim(s_get_word(3,tmp)); read(var,*) Ptarget(2)
       var=trim(s_get_word(4,tmp)); read(var,*) Ptarget(3)
       var=trim(s_get_word(5,tmp)); read(var,*) Pparm
       if (wc.ge.6) then; var=trim(s_get_word(6,tmp)); read(var,*) Pfreq; endif
  case ("Anisobaric"); isobaric=.false.
  case ("InitPos"); if (springpots)  call initPositionSet
  case ("BoxVel"); boxLoad = .true.
       var=trim(s_get_word(2,tmp)); read(var,*) boxVel(1)
       var=trim(s_get_word(3,tmp)); read(var,*) boxVel(2)
       var=trim(s_get_word(4,tmp)); read(var,*) boxVel(3)
  case ("StopBoxVel"); boxLoad = .false.
  case ("Group") !group NAME [type|dom] parameters
       Ngrp = Ngrp + 1 !add another group to temporal system
       groupName(Ngrp) = trim(s_get_word(2,tmp))
       call setGroup( 'g', s_get_word(3,tmp), s_right_of_word(4,tmp) )
  case ("ClearGroups")
       Ngrp = 0; call setGroup( 'g', 'all', 'all' )
  case ("Hold") !hold NAME [type|dom] parameters
       Nhold = Nhold + 1 !add another HoldGroup to temporal system
       holdName(Nhold) = trim(s_get_word(2,tmp))
       call setGroup( 'h', s_get_word(3,tmp), s_right_of_word(4,tmp) )
  case ("ClearHolds")
       Nhold = 0; call setGroup( 'h', 'all', 'all' )
  case ("ChangeType"); call setType( s_get_word(2,tmp), s_get_word(3,tmp) )
  case ("Disp"); i = holdID(s_get_word(2,tmp)); bcTyp(i) = 'd'
       var=trim(s_get_word(3,tmp)); if (trim(var).ne."F") bcDir(1,i)=.true.
         read(var,*) bcParm(1,i)
       var=trim(s_get_word(4,tmp)); if (trim(var).ne."F") bcDir(2,i)=.true.
         read(var,*) bcParm(2,i)
       var=trim(s_get_word(5,tmp)); if (trim(var).ne."F") bcDir(3,i)=.true.
         read(var,*) bcParm(3,i)
  ! outputs
  case ("Exec"); call system(trim(tmp(6:80)))
  case ("Print"); write(6,'(A)') "# "//trim(tmp(6:80))
  case ("OutFreq"); var=trim(s_get_word(2,tmp)); read(var,*) OutSteps
  case ("StepNumber"); var=trim(s_get_word(2,tmp)); read(var,*) step
  case ("ConfOutFreq"); var=trim(s_get_word(2,tmp)); read(var,*) ConfOuts
       if (wc.ge.3)  ModelOutType=trim(s_get_word(3,tmp))
  case ("Run"); var=trim(s_get_word(2,tmp)); read(var,*) step_num
       if (wc.ge.3) then; var=trim(s_get_word(3,tmp)); read(var,*) dt; endif
       call dynamics
  case ("END","DONE","end","done"); exit
  case ("EXIT","exit","quit","abort"); STOP
  case DEFAULT
       write(0,'(A,i4)') "#  read_input: ERROR: Unknown directive: "//trim(word)//" line:",ln
       call flush(0)
 end select
enddo

Return
!!!!!!!!!!!!!!!!!!!!!!
200 continue
write(0,'(A)') "#  read_input: ERROR: premature end of input file. No END directive found"
call flush(0); STOP

800 continue
write(0,'(A,i4,A)') "#  read_input: ERROR: on line:",ln," of input file"
call flush(0); STOP

end subroutine read_input !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine help
! {{{
! Help message output
implicit none
write(6,'(A)') "Manditory arguments in [], optional ones in {}"
write(6,'(A)') "If a directive was typed incorrectly, just type it again."
write(6,'(A)') "Directive      Description"
write(6,'(A)') "========================================================="
write(6,'(A)') " help               This output."
write(6,'(A)') " Model [FILE]       Specify the model file with the particle position data. "
write(6,'(A)') "                     if not specified, will place random particles in box."
write(6,'(A)') " Seed [i]           seed for random number generator"
write(6,'(A)') " NTotal [i]         Total number of particles, if no Model file."
write(6,'(A)') " BoxSize [r][r]{r}  Total size of simulation box in each direction"
write(6,'(A)') " Velocity [r] {i}   set initial velocity distribution to target kinetic energy"
write(6,'(A)') "                     density [r] KE/N for group {i}, set this to 3/2*k*T,"
write(6,'(A)') "                     convert the output KE to T=2*KE/(3*N*k)"
write(6,'(A)') " ReNeighbor [b] {i}  Reneighbor, ([b]=T/F), DEFAULT=T, check frequency {i}=1"
write(6,'(A)') " NeiMax [i]         Maximum number of neighbours, default=90."
write(6,'(A)') " Types [i]          Number of particle types"
write(6,'(A)') "  [c4] [r] [r]      List of each particle type Name, Mass and charge"
write(6,'(A)') " Pairs [i] {r} {i}  Number of pair potential interactions. rMin, Table size."
write(6,'(A)') "  [c4] [c4] [sym|asym] [c4] {rN}   List of interactions between particle"
write(6,'(A)') "                     types, if they're [i-j]=[j-i] or not, potential type,"
write(6,'(A)') "                     list of N potential parameters, dependent on the type."
write(6,'(A)') " Isokinetic [r] [r] {i} Set a constant kinetic energy density KE/N [r] with"
write(6,'(A)') "                      parm [r] every [i] steps"
write(6,'(A)') " Isoenergetic [r] {i} Maintian initial total energy with parm [r] every [i] steps"
write(6,'(A)') " Isobaric [r(3)] [r] {i} Set a constant pressure vector [r(3)] with parm [r] every [i] steps"
write(6,'(A)') " Anisokinetic       Unset the kinetistat"
write(6,'(A)') " Anisoenergetic     Unset the energetistat"
write(6,'(A)') " Anisobaric         Unset the barostat"
write(6,'(A)') " Group [c8] [type|dom] parameters   specify a set of atoms as in one group"
write(6,'(A)') "                     either by atom type or geometric domain"
write(6,'(A)') " Hold [c8] [type|dom] parameters   specify a set of atoms as in one HoldGroup"
write(6,'(A)') "                     either by atom type or geometric domain"
write(6,'(A)') " ChangeType [c8] [c4] Set the atoms in group [c8] to the atom type [c4]"
write(6,'(A)') " Disp [c8] [r(3)|l(3)] Apply a displacement vector to a HoldGroup [c8]"
write(6,'(A)') " BoxVel [r(3)]      Change box dimensions at a constant rate"
write(6,'(A)') " StopBoxVel         stop changing the box by unsetting BoxVel"
write(6,'(A)') " ConfOutFreq [i]    Number of timesteps between configuration outputs"
write(6,'(A)') " StepNumber [i]     Step number to begin a run with, useful when restarting"
write(6,'(A)') " Run [i] {r}        Run for [i] steps with [r] deltat"
write(6,'(A)') " Exec  [s]          Executes [s] as a system command"
write(6,'(A)') " END or DONE        Marks the end of the directive list, then begins run."
write(6,'(A)') " EXIT or exit or quit or abort    Will Terminate the program"
write(6,'(A)') ""
write(6,'(A)') "Pair potential types and parameters, last parameter is cutoff"
write(6,'(A)') " sine               rcut=PI/2 [Default]"
write(6,'(A)') " soft rc A B        soft cosine useful for overlapping atoms"
write(6,'(A)') " spring rc D0 sc Vol  spring potential using strain based on initial configuration"
write(6,'(A)') " morse rc D0 alpha r0  Morse pair potential"
write(6,'(A)') " buck  rc A C rho   Buckingham pair potential"
call flush(6)
end subroutine help !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_model
! {{{
!read coordinates and types from a file
use mostCommon
use control
use interpot
use groupStuff
implicit none
integer :: err, ln, i, id, ty, j
character(8) :: Ctyp
real(4) :: lm, lx, x(3)

if (trim(Model).eq."NA") then
   write(0,'(A)') "#  WARNING: read_model: No model file specified. "
   write(0,'(A)') "#    continuing assuming random particles"
   ! allocate all the properties
   allocate( typ(np), pos(3,np), vel(3,np), acc(3,np), force(3,np), pene(np), &
     &  grp(np), groupName(grpBuff), hold(np), holdName(holdBuff), &
     &  bcTyp(holdBuff), bcParm(3,holdBuff), bcDir(3,holdBuff) )
   typ(1:np) = 1; grp(1:np) = 0; hold(1:np) = 0; bcTyp(:) = 'X'
   write(6,'(a,i8)') '#  Number of particles: ', np
   if (Ntypes.gt.1) then
     if (sum(NatomsPerType).gt.np) then
       write(0,*) 'ERROR: atoms per type exceeds total number of atoms'
       STOP
     endif
     do i=1, Ntypes
       write(6,'(a,2i8)') '#  Number of particles of type ', i, NatomsPerType(i)
       do j=sum(NatomsPerType(0:i-1))+1, sum(NatomsPerType(0:i-1))+NatomsPerType(i)
         typ(j) = i
       enddo
     enddo
   endif
   call randomPositions
   call flush(0); Return
endif

!########################################
! open the model file that contains all of the particle coordinates, etc.
open(10,file=trim(Model),action="READ",iostat=err)
! check for file access errors
if (err.ne.0) then
   write(0,*) "read_model: ERROR: opening model file: "//trim(Model)//" iostat:",err
   call flush(0); STOP
endif

ln=1
if (trim(ModelInType).eq.'NA') ModelInType='xyz'
write(6,'(A)') "# reading model input file:"//trim(Model)//" format: "//trim(ModelInType)
call flush(6)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read file header information
select case (trim(ModelInType))
  case("xyz","XYZ")
    read(10,*,err=800) np; ln=ln+1
    read(10,*) ! read column Properties for extended XYZ format
  case("dump","trj")
    read(10,*); ln=ln+1 !ITEM: TIMESTEP
    read(10,*); ln=ln+1 !0 
    read(10,*); ln=ln+1 !ITEM: NUMBER OF ATOMS
    read(10,*,err=800) np; ln=ln+1
    read(10,*); ln=ln+1 !ITEM: BOX BOUNDS pp pp pp
    read(10,*) lm, lx; box(1)=dble(lx-lm); ln = ln+1
    read(10,*) lm, lx; box(2)=dble(lx-lm); ln = ln+1
    read(10,*) lm, lx; box(3)=dble(lx-lm); ln = ln+1
    read(10,*) ! ITEM: ATOMS id type xs ys zs #assuming this pattern
  case default; write(0,*) "ERROR, input Model Type: "//trim(ModelInType) &
    //"not a valid input model format"
    STOP
end select

! allocate all the properties
allocate( typ(np), pos(3,np), vel(3,np), acc(3,np), force(3,np), pene(np), &
  &  grp(np), groupName(grpBuff), hold(np), holdName(holdBuff), &
  &  bcTyp(holdBuff), bcParm(3,holdBuff), bcDir(3,holdBuff) )
grp(1:np) = 0; hold(1:np) = 0;  bcTyp(:) = 'X'
 bcParm(:,:) = 0.d0; bcDir(:,:) = .false.
vel(1:3,1:np) = 0.d0
acc(1:3,1:np) = 0.d0

! case for the Body of file
select case (trim(ModelInType))
  case("xyz","XYZ")
    i=0
    do i=1, np
       ln=ln+1
       !    type of node, position, list of intensive property values
       read(10,*,err=800) Ctyp, pos(1:3,i)
       typ(i) = typeID(Ctyp) !convert string type into internal typeID
       if (typ(i).eq.0) then
         write(0,*) "ERROR: no atom type defined for atom",i, "attempted type name:"//Ctyp
         write(0,*) "ERROR:  perhaps the type definitions were not defined yet??"
         STOP
       endif
       pos(1:3,i) = pos(1:3,i)/box(1:3) !positions are real not scaled
    enddo
  case("dump","trj")
    i=0
    do i=1, np
       ln=ln+1
       !    type of node, position, list of intensive property values
       read(10,*,err=800) id, ty, x
       typ(id) = ty
       pos(1:3,id) = dble(x(1:3)-0.5) !positions should be scaled in input file
    enddo
end select
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(6,'(A,i8,A)') "#  read_model: successfully read ",np," model nodes."

Return

800 continue
write(0,*) "read_model: ERROR: on line: ",ln
call flush(0); STOP
end subroutine read_model !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_model( FD )
! {{{
! write an XYZ file to file descriptor FD
use mostCommon
use interpot
use control
implicit none
integer :: FD, i


select case (trim(ModelOutType))
  case("XYZ",'xyz')
    write(FD,'(i8)') np
    write(FD,'(A)')
    do i = 1, np
      write(FD,'(A4,x,7EN14.4)') typeName(typ(i)), pos(1:3,i), vel(1:3,i), pene(i)
    enddo
  case("dump","trj")
    write(FD,'(A)') 'ITEM: TIMESTEP'
    write(FD,'(i7)') step
    write(FD,'(A)') 'ITEM: NUMBER OF ATOMS'
    write(FD,'(i7)') np
    write(FD,'(A)') 'ITEM: BOX BOUNDS pp pp pp'
    write(FD,'(2f12.4)') 0.0, box(1)
    write(FD,'(2f12.4)') 0.0, box(2)
    write(FD,'(2f12.4)') 0.0, box(3)
    write(FD,'(A)') 'ITEM: ATOMS id type xs ys zs'
    do i = 1, np
      write(FD,'(i7,x,i2,3f9.5)') i, typ(i), real(pos(1:3,i))+0.0
    enddo
  case default
    write(0,*) "ERROR: model out type: "//trim(ModelOutType)//" not supported"
    STOP
end select

call flush(FD)

end subroutine !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_thermo
!{{{
! write some periodic system statistics
  use control, only : step
  use mostCommon, only : potential, kinetic, virial, box, e0
  implicit none
  
  write(6, '(x,i8,x,g14.6,x,g14.6,x,g14.6,x,g14.6,x,g14.6)' ) &
    step, potential, kinetic, ( potential + kinetic - e0 ) / e0, &
      &  (virial(1)+virial(2)+virial(3))/(box(1)*box(2)*box(3)*3.0), &
      &  box(1)*box(2)*box(3)
  call flush(6)

end !}}}
!!!!!!!!!!!!!!!!!!!!!!!!
subroutine timestamp ( )
!{{{**************************************************************************80
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!    31 May 2001   9:45:54.872
!
!  Licensing:
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!    18 Feb 2018 by Ross J. Stewart to remove AM/PM distinction
!
!  Author:
!    John Burkardt
!
!  Parameters:
!    None
!
  implicit none
  integer ( kind = 4 ) d, h, m, mm, n, s, values(8), y
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)

  call date_and_time ( values = values )

  ! DATE: YYYY MM DD
  y = values(1);   m = values(2);   d = values(3)
  ! TIME: HH mm ssss
  h = values(5);   n = values(6);   s = values(7);   mm = values(8)

  write ( 6, '(a,i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3)' ) &
    '# ',d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm

  return
end !}}}

