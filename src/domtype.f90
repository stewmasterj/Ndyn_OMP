! vim:fdm=marker
module domType
!
!  Defining geometric Domains can be complicated, this module consolodates it
!  into simple to use subroutines and functions.
!
!  Friday, June 15, 2012
!  Ross J. Stewart
!
! derived data type: domain to be able to refer to domains as objects
type domain
    character(3)          :: ct
    logical               :: inside
    integer               :: axis
    real(4), dimension(8) :: parm
end type domain
! define the null domain
type (domain), parameter  :: nullDom=domain("   ",.true.,0, &
                                (/0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/) )
real(4), parameter :: small=tiny(1.0)

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine readDomain( s, dom )
!
!  Extracts domain information from a string and put the data into a domain
!  object for easy reference in the rest of the code.
! {{{
implicit none
character(80), intent(in)    :: s
character(80)                ::  word
type (domain), intent(inout) :: dom
integer                      :: i, wc

call left_of("!",s) ! get rid of comment portion of string
dom%ct=trim(s_get_word( 1, s )) !set the type of domain
wc=s_word_count( s )

if     ( dom%ct .eq. "rec" ) then
!rectangular, inside the domain, xmin xmax, ymin ymax, zmin zmax
 word=s_get_word( 2, s )
 read(word,*) dom%inside
 do i=1,6
   word=s_get_word( i+2, s );    read(word,*) dom%parm(i)
 enddo
elseif ( dom%ct .eq. "per" ) then
!perimeter, inside the domain, xmin xmax, ymin ymax, zmin zmax, only go
! inside domain this much
 word=s_get_word( 2, s );   read(word,*) dom%inside
 do i=1,7
   word=s_get_word( i+2, s );   read(word,*) dom%parm(i)
 enddo
elseif ( dom%ct .eq. "sph" ) then
!sphere, inside, x y z, radius
 word=s_get_word( 2, s );  read(word,*) dom%inside
 do i=1,4
   word=s_get_word( i+2, s );    read(word,*) dom%parm(i)
 enddo
elseif ( dom%ct .eq. "ell" ) then
!circle, inside, axial direction i.e. 2=Y, x z, a b
 word=s_get_word( 2, s );  read(word,*) dom%inside
 word=s_get_word( 3, s )
 read(word,*) dom%axis
 do i=1,4
   word=s_get_word( i+3, s );    read(word,*) dom%parm(i)
 enddo
elseif ( dom%ct .eq. "tub" ) then
!tube, inside, alial direction, 2 of x y or z, outer a b, inner a b, zmin, zmax
 word=s_get_word( 2, s );  read(word,*) dom%inside
 word=s_get_word( 3, s )
 read(word,*) dom%axis
 do i=1,8
   word=s_get_word( i+3, s );    read(word,*) dom%parm(i)
 enddo
elseif ( dom%ct .eq. "wed" ) then
!wedge, inside, alial direction, 2 of x y or z, radius, largest Angle, angle
!between
 word=s_get_word( 2, s );  read(word,*) dom%inside
 word=s_get_word( 3, s )
 read(word,*) dom%axis
 do i=1,5
   word=s_get_word( i+3, s );    read(word,*) dom%parm(i)
 enddo
endif

end subroutine readDomain !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine writeDomain( dom )
!
!  Writes the specified domain to standard error
!  Useful for verifying wether the doamin was read correctly.
! {{{
implicit none
character(30)             ::  frm
type (domain), intent(in) :: dom

if     ( dom%ct .eq. "rec" ) then
!rectangular, inside the domain, xmin xmax, ymin ymax, zmin zmax
 frm='(a3,x,L,x,6f10.3)'
 write(0,frm) dom%ct,dom%inside,dom%parm(1:6)
elseif ( dom%ct .eq. "per" ) then
!perimeter, inside the domain, xmin xmax, ymin ymax, zmin zmax, only go
! inside domain this much
 frm='(a3,x,L,x,7f10.3)'
 write(0,frm) dom%ct,dom%inside,dom%parm
elseif ( dom%ct .eq. "sph" ) then
!sphere, inside, x y z, radius
 frm='(a3,x,L,x,4f10.3)'
 write(0,frm) dom%ct,dom%inside,dom%parm(1:4)
elseif ( dom%ct .eq. "ell" ) then
!circle, inside, axial direction i.e. 2=Y, x z, a b
 frm='(a3,x,L,x,i2,x,4f10.3)'
 write(0,frm) dom%ct,dom%inside,dom%axis,dom%parm(1:4)
elseif ( dom%ct .eq. "tub" ) then
!tube, inside, alial direction, 2 of x y or z, outer a b, inner a b
 frm='(a3,x,L,x,i2,x,8f10.3)'
 write(0,frm) dom%ct,dom%inside,dom%axis,dom%parm(1:8)
elseif ( dom%ct .eq. "wed" ) then
!wedge, inside, alial direction, 2 of x y or z, radius, largest Angle, angle
!between
 frm='(a3,x,L,x,i2,x,5f8.3)'
 write(0,frm) dom%ct,dom%inside,dom%axis,dom%parm(1:5)
endif

end subroutine writeDomain !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function s_word_count ( s )
!*****************************************************************************80
!! S_WORD_COUNT counts the number of "words" in a string.
! {{{
!  Licensing:
!    This code is distributed under the GNU LGPL license.
!  Modified:
!    14 April 1999, 30 May 2012
!  Author:
!    John Burkardt
!  Made into a function by Ross J. Stewart
!  Parameters:
!    Input, character ( len = * ) S, the string to be examined.
!    Output, integer ( kind = 4 ) NWORD, the number of "words" in the string.
!    Words are presumed to be separated by one or more blanks.
implicit none

logical blank
integer ( kind = 4 ) i
integer ( kind = 4 ) lens
integer ( kind = 4 ) s_word_count
character ( len = * ) s
character ( len = 12 ) :: delimiters

delimiters=" ,:()[]{}'"//char(09)//char(34)
s_word_count = 0
lens = len ( s )

if ( lens <= 0 ) Return

blank = .true.

do i = 1, lens
    if ( index(delimiters,s(i:i)).ne.0 ) then
      blank = .true.
    else if ( blank ) then
      s_word_count = s_word_count + 1
      blank = .false.
    end if
enddo
end function s_word_count !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function s_get_word( n, s )
!
! get a word from a string, s. the word at position, n
! where n corresponds to the nth word in the sentence, s
! {{{
! Ross J. Stewart, 01 June 2012
implicit none
logical blank
integer ( kind = 4 ) i, n
integer ( kind = 4 ) lens
integer ( kind = 4 ) s_word_count
character ( len = * ) s
character ( len = 80 )  s_get_word
character ( len = 12 ) :: delimiters

! delimiters are what separates words from eachother
delimiters=" ,:()[]{}'"//char(09)//char(34)
s_word_count = 0
lens = len ( s )
s_get_word = ""

if ( lens <= 0 ) Return

blank = .true.

do i = 1, lens
   if ( index(delimiters,s(i:i)).ne.0 ) then
      blank = .true.
   elseif ( blank ) then
      s_word_count = s_word_count + 1
      blank = .false.
   endif
   if (.not.blank.and.n.eq.s_word_count) then
      s_get_word=trim(s_get_word)//s(i:i)
   endif
enddo

end function !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function s_right_of_word( n, s )
!
! get all words right of the word at position, n
! where n corresponds to the nth word in the sentence, s
! {{{
! Ross J. Stewart, 16 Feb 2018
implicit none
logical blank
integer ( kind = 4 ) i, n
integer ( kind = 4 ) lens
integer ( kind = 4 ) s_word_count
character ( len = * ) s
character ( len = 80 )  s_right_of_word
character ( len = 12 ) :: delimiters

! delimiters are what separates words from eachother
delimiters=" ,:()[]{}'"//char(09)//char(34)
s_word_count = 0
lens = len ( s )
s_right_of_word = ""

if ( lens <= 0 ) Return

blank = .true.

do i = 1, lens
   if ( index(delimiters,s(i:i)).ne.0 ) then
      blank = .true.
   elseif ( blank ) then
      s_word_count = s_word_count + 1
      blank = .false.
   endif
   if (.not.blank.and.n.eq.s_word_count) then
      s_right_of_word=trim(s_right_of_word)//s(i:lens)
      Return
   endif
enddo

end function !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine left_of( w , s )
!
! Returns the string, s, as everything left of the character or word, w
! used for extracting useful information from comments denotes by the character,
! w. Typically a "#" or "!"
! {{{
! Ross J. Stewart, 31 May 2012
implicit none
logical comment
integer ( kind = 4 ) i
integer ( kind = 4 ) lens, lenw
character ( len = * ) s
character ( len = * ) w

lens = len ( s )
lenw = len ( w )

if ( lens <= 0 .or. lenw .le. 0 .or. index(s,w).eq.0 ) Return

comment = .false.

do i = 1, lens
   if ( index(s(i:i+lenw-1),w).ne.0 .or. comment) then
      comment=.true.
      s(i:i)=" " ! overwrite the rest with blanks
   endif
enddo

end subroutine left_of !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
function inDomain( pos, dom)
!
! Function to determine if an xyz position is in the domain,
! dom, the function returns a boolean value.
! This is useful for determinign which particles should be clamped and which to
! be monitored for stresses, etc. 
! {{{
! Ross Feb 2012 V16
implicit none
logical                           :: inDomain
integer                           :: i
logical, dimension(3)             :: good
real(4), dimension(3), intent(in) :: pos
real(4), dimension(3)             :: dmin, dmax
integer, dimension(3)             ::  dir
type (domain), intent(in)         :: dom
real(4)      :: depth, dist, slope, re, ro

inDomain=.false.
slope=0.0; dist=0.0; re=0.0

if     ( dom%ct .eq. "rec" .or. dom%ct .eq. "per" ) then
!rectangular, inside the domain, xmin xmax, ymin ymax, zmin zmax
 good=.false.
 dmin(1)=dom%parm(1); dmax(1)=dom%parm(2)
 dmin(2)=dom%parm(3); dmax(2)=dom%parm(4)
 dmin(3)=dom%parm(5); dmax(3)=dom%parm(6)
 dir=1 ! set dimension switches for a closed domain

 if ((abs(dmax(1))+abs(dmin(1))).lt.small) dir(1)= 0   ! entire domain [-INF,INF]
 if ((abs(dmax(2))+abs(dmin(2))).lt.small) dir(2)= 0   ! entire domain [-INF,INF]
 if ((abs(dmax(3))+abs(dmin(3))).lt.small) dir(3)= 0   ! entire domain [-INF,INF]

 do i=1,3 !loop through each dimension
  if     (dir(i).eq. 1) then ! closed domain
        if (pos(i).gt.dmin(i) .and. pos(i).le.dmax(i)) then
                               good(i)=.true.
        else;                  good(i)=.false.; endif
  elseif (dir(i).eq. 0) then ! entire domain [-INF,INF]
                               good(i)=.true.
  else
     write(6,*) "inDomain: WARNING: Something weird happend in domain, pos:", &
           &  pos,"good?",good, "dir:",dir; endif
 enddo

 if (good(1).and.good(2).and.good(3)) then; inDomain=.true.
 else;                                      inDomain=.false.; endif

 if ( dom%ct .eq. "per" ) then
 depth=dom%parm(7)
  do i=1,3 !loop through each dimension
   if     (dir(i).eq. 1) then ! closed domain
         if (pos(i).gt.dmin(i)+depth .and. pos(i).le.dmax(i)-depth) then
                                good(i)=.true.
         else;                  good(i)=.false.; endif
   elseif (dir(i).eq. 0) then ! entire domain [-INF,INF]
                                good(i)=.true.
   else
      write(6,*) "inDomain: WARNING: Something weird happend in domain, pos:", &
            &   pos,"good?",good, "dir:",dir; endif
  enddo

  if (good(1).and.good(2).and.good(3)) inDomain=.false.
 endif

 if (.not.dom%inside) inDomain=.not.inDomain ! invert result

elseif ( dom%ct .eq. "sph" ) then
!sphere, inside, x y z, radius
 if (sum((pos-dom%parm(1:3))**2.) .le. dom%parm(4)**2.) then
        inDomain=.true.; endif

 if (.not.dom%inside) inDomain=.not.inDomain ! invert result

elseif ( dom%ct .eq. "ell" ) then
!ellipse, inside, axial direction i.e. 2=Y, x z, a b
 if     (dom%axis.eq.1) then
   dist=((pos(2)-dom%parm(1))**2.+(pos(3)-dom%parm(2))**2.)
   slope=( (pos(2)-dom%parm(1)) + (pos(3)-dom%parm(2)) )/ &
         ( (pos(2)-dom%parm(1)) - (pos(3)-dom%parm(2)) )
   if (abs((pos(2)-dom%parm(1)) - (pos(3)-dom%parm(2))).lt.small) slope=0.0
 elseif (dom%axis.eq.2) then
   dist=((pos(1)-dom%parm(1))**2.+(pos(3)-dom%parm(2))**2.)
   slope=( (pos(1)-dom%parm(1)) + (pos(3)-dom%parm(2)) )/ &
         ( (pos(1)-dom%parm(1)) - (pos(3)-dom%parm(2)) )
   if (abs((pos(1)-dom%parm(1)) - (pos(3)-dom%parm(2))).lt.small) slope=0.0
 elseif (dom%axis.eq.3) then
   dist=((pos(1)-dom%parm(1))**2.+(pos(2)-dom%parm(2))**2.)
   slope=( (pos(1)-dom%parm(1)) + (pos(2)-dom%parm(2)) )/ &
         ( (pos(1)-dom%parm(1)) - (pos(2)-dom%parm(2)) )
   if (abs((pos(1)-dom%parm(1)) - (pos(2)-dom%parm(2))).lt.small) slope=0.0
 endif
 re=dom%parm(3)*dom%parm(4)*sqrt(slope**2.+1.)/ &
      sqrt(dom%parm(3)**2.*slope**2.+dom%parm(4)**2.)
 if (dist.le.re**2.) inDomain=.true.

 if (.not.dom%inside) inDomain=.not.inDomain ! invert result

elseif ( dom%ct .eq. "tub" ) then
!tube, inside, alial direction, 2 of x y or z, outer a b, inner a b, zmin, zmax
 good(1)=.false.
 if     (dom%axis.eq.1) then
   dist=((pos(2)-dom%parm(1))**2.+(pos(3)-dom%parm(2))**2.)
   slope=( (pos(2)-dom%parm(1)) + (pos(3)-dom%parm(2)) )/ &
         ( (pos(2)-dom%parm(1)) - (pos(3)-dom%parm(2)) )
   if (abs((pos(2)-dom%parm(1)) - (pos(3)-dom%parm(2))).lt.small) slope=0.0
 elseif (dom%axis.eq.2) then
   dist=((pos(1)-dom%parm(1))**2.+(pos(3)-dom%parm(2))**2.)
   slope=( (pos(1)-dom%parm(1)) + (pos(3)-dom%parm(2)) )/ &
         ( (pos(1)-dom%parm(1)) - (pos(3)-dom%parm(2)) )
   if (abs((pos(1)-dom%parm(1)) - (pos(3)-dom%parm(2))).lt.small) slope=0.0
 elseif (dom%axis.eq.3) then
   dist=((pos(1)-dom%parm(1))**2.+(pos(2)-dom%parm(2))**2.)
   slope=( (pos(1)-dom%parm(1)) + (pos(2)-dom%parm(2)) )/ &
         ( (pos(1)-dom%parm(1)) - (pos(2)-dom%parm(2)) )
   if (abs((pos(1)-dom%parm(1)) - (pos(2)-dom%parm(2))).lt.small) slope=0.0
 endif
 re=dom%parm(3)*dom%parm(4)*sqrt(slope**2.+1.)/ &
      sqrt(dom%parm(3)**2.*slope**2.+dom%parm(4)**2.)
 ro=dom%parm(5)*dom%parm(6)*sqrt(slope**2.+1.)/ &
      sqrt(dom%parm(5)**2.*slope**2.+dom%parm(6)**2.)
 if (pos(dom%axis).gt.dom%parm(7) .and. pos(dom%axis).lt.dom%parm(8)) good(1)=.true.
 if (dist.le.re**2. .and. dist.gt.ro**2. .and. good(1)) inDomain=.true.

 if (.not.dom%inside) inDomain=.not.inDomain ! invert result

elseif ( dom%ct .eq. "wed" ) then
!wedge, inside, alial direction, 2 of x y or z, radius, largest Angle, angle
!between
 if     (dom%axis.eq.1) then ! using yz plane
   dist=((pos(2)-dom%parm(2))**2.+(pos(3)-dom%parm(3))**2.)
   if (abs((pos(2)-dom%parm(2)) - (pos(3)-dom%parm(3))).lt.small) then
     re=0.0
   else
     re=fullangle((/(pos(2)-dom%parm(2)),(pos(3)-dom%parm(3))/)) ! particle angle
   endif
 elseif (dom%axis.eq.2) then ! using xz plane
   dist=((pos(1)-dom%parm(1))**2.+(pos(3)-dom%parm(3))**2.)
   if (abs((pos(1)-dom%parm(1)) - (pos(3)-dom%parm(3))).lt.small) then
     re=0.0
   else
     re=fullangle((/(pos(1)-dom%parm(1)),(pos(3)-dom%parm(3))/)) ! particle angle
   endif
 elseif (dom%axis.eq.3) then ! using xy plane
   dist=((pos(1)-dom%parm(1))**2.+(pos(2)-dom%parm(2))**2.)
   if (abs((pos(1)-dom%parm(1)) - (pos(2)-dom%parm(2))).lt.small) then
     re=0.0
   else
     re=fullangle((/(pos(1)-dom%parm(1)),(pos(2)-dom%parm(2))/)) ! particle angle
   endif
 endif
 if (dist.le.dom%parm(3)**2 .and. re.gt.dom%parm(4)-dom%parm(5) .and. &
     re.le.dom%parm(4)) inDomain=.true.
 if (.not.dom%inside) inDomain=.not.inDomain ! invert result

endif
end function inDomain !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
function DomainVolume( dom, box )
!
! Returns the volume of domain, dom, within the global box.
! USed in the local stress caluclation.
! {{{
implicit none
real(4)                   :: DomainVolume
integer                   :: i
type (domain), intent(in) :: dom
real(4), dimension(3), intent(in) :: box
real(4), dimension(3) :: dmin, dmax, lbox, sbox
integer, dimension(3) :: dir

DomainVolume=0.0
dmin(1)=dom%parm(1); dmax(1)=dom%parm(2)
dmin(2)=dom%parm(3); dmax(2)=dom%parm(4)
dmin(3)=dom%parm(5); dmax(3)=dom%parm(6)
dir=1 ! set dimension switches for a closed domain

if ((abs(dmax(1))+abs(dmin(1))).lt.small) dir(1)= 0   ! entire domain [-INF,INF]
if ((abs(dmax(2))+abs(dmin(2))).lt.small) dir(2)= 0   ! entire domain [-INF,INF]
if ((abs(dmax(3))+abs(dmin(3))).lt.small) dir(3)= 0   ! entire domain [-INF,INF]

lbox=0.0
if     ( dom%ct .eq. "rec" ) then
!rectangular, inside the domain, xmin xmax, ymin ymax, zmin zmax
 do i=1,3 !loop through each dimension
  if     (dir(i).eq. 1) then ! closed domain
        lbox(i)=dmax(i)-dmin(i)
  elseif (dir(i).eq. 0) then ! entire domain [-INF,INF]
        lbox(i)=box(i)
  else; write(6,*) "DomainVolume: WARNING: cannot find volume of domain:"
        call writeDomain(dom)
        write(6,*) "     in model size:", box
  endif
 enddo
 DomainVolume=product(lbox)
!perimeter, inside the domain, xmin xmax, ymin ymax, zmin zmax, only go
! inside domain this much
 do i=1,3 !loop through each dimension
  if     (dir(i).eq. 1) then ! closed domain
        lbox(i)=dmax(i)-dmin(i)
        sbox(i)=(dmax(i)-dom%parm(7))-(dmin(i)+dom%parm(7))
  elseif (dir(i).eq. 0) then ! entire domain [-INF,INF]
        lbox(i)=box(i)
        sbox(i)=box(i)
  else; write(6,*) "DomainVolume: WARNING: cannot find volume of domain:"
        call writeDomain(dom)
        write(6,*) "     in model size:", box
  endif
 enddo
 DomainVolume=product(lbox)-product(sbox)

elseif ( dom%ct .eq. "sph" ) then
!sphere, inside, x y z, radius
 DomainVolume=4./3.*3.1415926535*(dom%parm(4))**3.

elseif ( dom%ct .eq. "ell" ) then
!circle, inside, axial direction i.e. 2=Y, x z, a b
 DomainVolume=3.1415926535*dom%parm(3)*dom%parm(4)*box(dom%axis)

elseif ( dom%ct .eq. "tub" ) then
!tube, inside, alial direction, 2 of x y or z, outer a b, inner a b
 DomainVolume=3.1415926535*(dom%parm(3)*dom%parm(4)- &
                            dom%parm(5)*dom%parm(6))* & !box(dom%axis)
                           (dom%parm(8)-dom%parm(7))
elseif ( dom%ct .eq. "wed" ) then
!wedge, inside, alial direction, 2 of x y or z, radius, largest Angle, angle
!between
 DomainVolume=(dom%parm(4)-dom%parm(5))*dom%parm(3)**2*box(3)/2.0
endif

if (.not.dom%inside) DomainVolume=product(box)-DomainVolume

end function DomainVolume !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
function DomainCenter( dom )
!
! Returns the center position of the domain, dom
! USed to determine the farthest particles within a domain; part of the ADD
! balancing procedure.
! {{{
implicit none
type (domain), intent(in) :: dom
real(4)                   :: r, a  ! for wedge
real(4), dimension(3)     :: DomainCenter

DomainCenter=0.0
if     ( dom%ct .eq. "rec" ) then
!rectangular, inside the domain, xmin xmax, ymin ymax, zmin zmax
 DomainCenter(1)=(dom%parm(2)+dom%parm(1))/2.
 DomainCenter(2)=(dom%parm(4)+dom%parm(3))/2.
 DomainCenter(3)=(dom%parm(6)+dom%parm(5))/2.
elseif ( dom%ct .eq. "per" ) then
!perimeter, inside the domain, xmin xmax, ymin ymax, zmin zmax, only go
! inside domain this much
 DomainCenter(1)=(dom%parm(2)+dom%parm(1))/2.
 DomainCenter(2)=(dom%parm(4)+dom%parm(3))/2.
 DomainCenter(3)=(dom%parm(6)+dom%parm(5))/2.
elseif ( dom%ct .eq. "sph" ) then
!sphere, inside, x y z, radius
 DomainCenter=dom%parm(1:3)
elseif ( dom%ct .eq. "ell" ) then
!circle, inside, axial direction i.e. 2=Y, x z, a b
 if(dom%axis.eq.1)  DomainCenter(2:3)=dom%parm(1:2)
 if(dom%axis.eq.2) then
                    DomainCenter(1)=dom%parm(1)
                    DomainCenter(3)=dom%parm(2)
 endif
 if(dom%axis.eq.3)  DomainCenter(1:2)=dom%parm(1:2)
elseif ( dom%ct .eq. "tub" ) then
!tube, inside, alial direction, 2 of x y or z, outer a b, inner a b
 if(dom%axis.eq.1) then
                    DomainCenter(1)  =(dom%parm(7)+dom%parm(8))/2.0
                    DomainCenter(2:3)=dom%parm(1:2)
 endif
 if(dom%axis.eq.2) then
                    DomainCenter(1)=dom%parm(1)
                    DomainCenter(3)=dom%parm(2)
                    DomainCenter(2)=(dom%parm(7)+dom%parm(8))/2.0
 endif
 if(dom%axis.eq.3) then
                    DomainCenter(1:2)=dom%parm(1:2)
                    DomainCenter(3)  =(dom%parm(7)+dom%parm(8))/2.0
 endif
elseif ( dom%ct .eq. "wed" ) then
!wedge, inside, alial direction, 2 of x y or z, radius, largest Angle, angle
!between
 r=2.0/3.0*dom%parm(3)*sin(dom%parm(5))/dom%parm(5)
 a=dom%parm(4)-dom%parm(5)/2.0
 if(dom%axis.eq.1)  DomainCenter(2:3)=(/r*cos(a),r*sin(a)/)
 if(dom%axis.eq.2)  then
                    DomainCenter(3)=r*sin(a)
                    DomainCenter(3)=r*sin(a)
 endif
 if(dom%axis.eq.3)  DomainCenter(1:2)=(/r*cos(a),r*sin(a)/)
endif

end function DomainCenter !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
function fullangle(v)
!
! Returns the angle of a vector, v from the x axis.
! used for the wedge domain.
! {{{
implicit none
real(4) :: at
real(4) :: fullangle
real(4), dimension(2), intent(in) :: v
real(4), parameter :: pi=3.1415926535

fullangle=0.0

if(abs(v(1)).lt.small.and.abs(v(2)).lt.small)then
   at=pi/2.0
else
   at=atan2(v(2),v(1))
endif

fullangle=at
if (at.lt.0.0) fullangle=at+2*pi

endfunction fullangle !}}}

end module domType

