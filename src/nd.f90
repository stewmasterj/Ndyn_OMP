! vim: fdm=marker
!*****************************************************************************80
program main  
!{{{
!{{{**************************************************************************80
!
!! MAIN is the main program for Ndyn_OPENMP.
!
!  Discussion:
!    Ndyn_OPENMP implements a molecular dynamics simulation.
!    The program uses OpenMP directives to allow parallel computation.
!    The velocity Verlet time integration scheme is used. 
!    The particles interact with central pair potentials.
!
!  Licensing:
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!    30 July 2009
!    30 March 2018
!
!  Author:
!    Original FORTRAN90 version by Bill Magro.
!    Next FORTRAN90 version by John Burkardt.
!    This FORTRAN90 version by Ross J. Stewart.
!
!}}}**************************************************************************80
  use omp_lib
  use control
  implicit none
  real(8) :: tottime

  tottime = omp_get_wtime ( )
  call timestamp ( )

  proc_num = omp_get_num_procs ( )
  thread_num = omp_get_max_threads ( )

  write(6,'(a)') 
  write(6,'(a)') '#  N D Y N'
  write(6,'(a)') '#  FORTRAN90/OpenMP version 0.10.18'
  write(6,'(a)') 
  write(6,'(a)') '#  A newtonian dynamics program.'
  write(6,'(a)') 
  write(6,'(a)') '#  written by Ross J. Stewart'
  write(6,'(a)') 
  write(6,'(a,i3)') '#  Number of processors available: ', proc_num
  write(6,'(a,i3)') '#  Number of threads specified:    ', thread_num
  write(6,'(a)') 

  call read_input

  call timestamp ( )
  tottime = omp_get_wtime ( ) - tottime
  write(6, '(a,g14.6,a)' ) '# Total elapsted time: ',tottime, ' seconds'
  write(6, '(a)' ) ''
  write(6, '(a)' ) '#  Ndyn:  Normal end of execution.'

end !}}}
subroutine dynamics
!{{{
! integrate the equations of motion
  use omp_lib
  use mostCommon
  use control
  use interpot
  use listData
  implicit none
  integer(4) :: tt
  logical :: calcVirial, calcPot, calcKE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
! Initializing cell linked list parameters
  if (skin < 0.d0) then
    skin = rcut*1.2d0
  endif
  MCELL(1) = int(box(1)/(rcut+skin)) !number of cells in each direction
  MCELL(2) = int(box(2)/(rcut+skin)) !number of cells in each direction
  MCELL(3) = int(box(3)/(rcut+skin)) !number of cells in each direction
  NCELL = MCELL(1)*MCELL(2)*MCELL(3)
  write(6,'(A,3i4)') "# Number of link cell: ",MCELL
  if (NeiMax.eq.0) then
    ! approximate max neighbours +10% in cutoff radius
    NeiMax = int(4./3.*PI*((rcut+skin)*1.1)**3*np/(box(1)*box(2)*box(3)))
  endif
  write(6,'(A,i4)') "# Estimating max number of neighbours: ",NeiMax
  if (allocated(HEAD)) then
    deallocate( HEAD, VNL, List, dispList, ftmp )
  endif
  allocate( HEAD(NCELL), VNL(1:NeiMax+1,1:np), List(np), dispList(3,np), &
     &  ftmp(1:3,thread_num*np) )
  NlistUpdates = 0 !initialize initial times the lists are updated
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
  ! test that things exist
  if (.not.allocated(mass)) then
     allocate( mass(Ntypes), typeName(Ntypes) )
     mass(1) = 1.d0; typeName(1) = "test"
  endif
  if (ConfOuts.eq.0) ConfOuts = step_num
  if (trim(ModelOutType).eq.'NA') ModelOutType='xyz'
  if (springpots.and..not.allocated( ipos )) call initPositionSet
  skinsq = skin*skin

  call buildCellList( .true. ) ! output stats about neighbour lists
  ! time of ending "overhead"
  call timestamp ( )
!
  write(6,'(a)') ''
  write(6,'(a)') '#  Computing initial forces and energies.'
  write(6,'(a,i8)') '#  Number of time steps: ', step_num
  write(6,'(a,g14.6)') '#  Time step increment: ', dt
  write(6,'(a,g14.6)') '#  Cutoff radius: ', rcut
  write(6,'(a,g14.6)') '#  Cutoff skin: ', skin
  write(6,'(a)') '# Model Output format: '//trim(ModelOutType)
  write(6,'(a)') ''
  write(6,'(a)') '#  Step   Potential     Kinetic        EnergyError    pressure       volume'
  call flush(6)

  box0(1) = box(1)
  box0(2) = box(2)
  box0(3) = box(3)

  ! virial pressure has consistent units, usually eV/Ang^3
  !   multiply by 160.21766208 to get GPa.

  wtime = omp_get_wtime ( )
  open(11,file='out.'//trim(ModelOutType),position="append")
  call write_model(11)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
!  This is the main time stepping loop:
!    Compute forces and energies,
!    Update positions, velocities, accelerations.
  do tt = 1, step_num
    step = step + 1
    if (mod(step,OutSteps) == 0 .or. tt == 1) then ! for some delicious output data
      calcVirial = .true. ; calcPot = .true. ; calcKE = .true.
    endif
    !
    ! set flags for ensemble required data
    if (isokinetic) then;   if (mod(step,KEfreq) == 0) calcKE = .true.;     endif
    if (isobaric) then;     if (mod(step,Pfreq) == 0)  calcVirial = .true.; endif
    if (isoenergetic) then; if (mod(step,Efreq) == 0) then
      calcPot = .true.; calcKE = .true. 
    endif; endif
    !
    ! check to see if we need to update neighbour lists
    if (ReNei) then;    if (mod(step,ReNeiFreq) == 0) call CheckListUpdate; endif 
    !
    ! Compute forces etc
    call integrate1
    if     (asympots.and.springpots) then
      call compute3(calcVirial,calcPot)
    elseif (.not.asympots.and.springpots) then
      call compute2(calcVirial,calcPot)
    elseif (asympots.and..not.springpots) then
      call compute1(calcVirial,calcPot) ! compute forces, virial, potential
    else
      call compute0(calcVirial,calcPot) ! compute forces, virial, potential
    endif
    call integrate2(calcKE)
    ! after 'integrate2' we have x(t+dt) and v(t+dt) and a(t+dt)
    !
    if (tt == 1) e0 = potential + kinetic
    if (mod(step,OutSteps) == 0 .or. tt == 1) call write_thermo
    ! unset ensemble flags
    calcVirial = .false.; calcPot = .false.; calcKE = .false.
    
    if (mod(step,ConfOuts) == 0) call write_model(11)

    ! ensemble modifications
    if (isokinetic) then;   if (mod(step,KEfreq) == 0) call mkIsoKE; endif
    if (isobaric) then;     if (mod(step,Pfreq) == 0)  call mkIsoP;  endif
    if (isoenergetic) then; if (mod(step,Efreq) == 0)  call mkIsoE;  endif
    if (boxLoad) then;      box = box + boxVel*dt; endif

  end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80

  call flush(0);   call flush(6);   call flush(11)
  close(11)

  wtime = omp_get_wtime ( ) - wtime
  write(6, '(a,i8)' ) '# Number of neighbour list updates:',NlistUpdates
  write(6, '(a)' ) ''
  write(6, '(a)' ) '#  Elapsed time for main dynamics:'
  write(6, '(a,g14.6,a)' ) '#  ',wtime, ' seconds'
!
!  Terminate.
!
  call timestamp ( )

  deallocate( HEAD, VNL, List, dispList, ftmp )

  Return
end !}}}
subroutine compute0( CV, CU)
!{{{**************************************************************************80
!! COMPUTE computes the forces and energies.
!
!  Discussion:
!    The computation of forces and energies is fully parallel.
!    The potential function V(X) is a harmonic well which smoothly
!    saturates to a maximum value at PI/2:
!      v(x) = ( sin ( min ( x, PI2 ) ) )^2
!    The derivative of the potential is:
!      dv(x) = 2.0D+00 * sin ( min ( x, PI2 ) ) * cos ( min ( x, PI2 ) )
!            = sin ( 2.0 * min ( x, PI2 ) )
!
!  Licensing:
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!    Feb 2018
!
!  Author:
!    Original FORTRAN90 version by Bill Magro and later by John Burkardt.
!    This FORTRAN90 version by Ross J. Stewart.
!
!  Parameters:
!    Input, integer ( kind = 4 ) NP, the number of particles.
!    Input, integer ( kind = 4 ) ND, the number of spatial dimensions.
!    Input, real ( kind = 8 ) POS(ND,NP), the position of each particle.
!    Input, real ( kind = 8 ) VEL(ND,NP), the velocity of each particle.
!    Input, real ( kind = 8 ) MASS, the mass of each particle.
!    Output, real ( kind = 8 ) F(ND,NP), the forces.
!    Output, real ( kind = 8 ) POT, the total potential energy.
!    Output, real ( kind = 8 ) KIN, the total kinetic energy.
!
  use omp_lib
  use mostCommon
  use interpot
  use listData
  use control
  implicit none

  real(8) :: d2, rij(3), dp, df, kt, tmpi(3), sr(3), weight, ptmp(8), vtmp(6,8)
  real(4) :: rk, vPTdr2(NPairs)
  integer(4) :: i, jn, j, k, typi, PTij, istart, iend, tid, tido
  logical, intent(in) :: CV, CU

  if (CV) then
    virial(1:6) = 0.d0
    vtmp(1:6,1:8) = 0.d0
  endif
  if (CU) then
    potential = 0.0D+00
    pene(1:np) = 0.0d+00
    ptmp(1:8) = 0.d0
  endif
  vPTdr2(1:NPairs) = real(1.d0/PTdr2(1:NPairs))

!$omp parallel workshare &
!$omp shared ( force, np )
  force(1:3,1:np) = 0.0D+00
  ftmp(1:3,1:np*thread_num) = 0.0D+00
!$omp end parallel workshare

!$omp parallel &
!$omp shared ( force, np, pos, vel, mass, CV, vPTdr2, PTrMin, potential) &
!$omp shared ( virial, pene, ftmp, vtmp, ptmp, CU ) &
!$omp private ( d2, i, j, k, jn, rij, dp, df, kt, sr, weight, rk, tmpi) &
!$omp private ( typi, tid, istart, iend, PTij, tido )
tid=omp_get_thread_num()

istart=int(np/thread_num)*tid + 1
iend  =int(np/thread_num)*(tid+1)
if (iend<np.and.tid.eq.thread_num-1) iend=np

  do i = istart, iend !for each particle on this thread...
    tmpi(1) = pos(1,i)
    tmpi(2) = pos(2,i)
    tmpi(3) = pos(3,i)
    typi = typ(i)
    tido = i+tid*np
    do jn = 2, VNL(1,i) ! i particles in cell ic
      j = VNL(jn,i) !first particle ID in cell jc
      PTij = PairMat(typi,typ(j))
      ! calculate distance i-j
      sr(1) = tmpi(1) - pos(1,j)
      sr(2) = tmpi(2) - pos(2,j)
      sr(3) = tmpi(3) - pos(3,j)
      sr(1) = sr(1) - anint(sr(1))
      sr(2) = sr(2) - anint(sr(2))
      sr(3) = sr(3) - anint(sr(3))
      rij(1) = box(1)*sr(1) !real distance units
      rij(2) = box(2)*sr(2) !real distance units
      rij(3) = box(3)*sr(3) !real distance units
      d2 = rij(1)*rij(1)
      d2 = d2 + rij(2)*rij(2)
      d2 = d2 + rij(3)*rij(3)

      if (d2 > PTrc2(PTij) ) cycle
!  Attribute half of the potential energy to particle I and J.
      ! table interpolation
      rk = real(d2 - PTrMin)*vPTdr2(PTij) + 1.0 ! "continuous" index in table
      k = int(rk)                     ! discretized index
      if (k < 1) then
         k = 1                ! unlikely but just to protect
      !   write(0,*) "ERROR: WTF is k<1?",step,d2,i,j,pos(:,i),pos(:,j), rij
      endif
      weight = dble(rk) - dble(k)                 ! fractional part, in [0,1]
      ! do  linear
      df = weight* FPairTable(k+1,PTij) + &
         &  (1.d0-weight)* FPairTable(k,PTij)

      df = df / d2
      ftmp(1,tido) = ftmp(1,tido) + rij(1)*df
      ftmp(2,tido) = ftmp(2,tido) + rij(2)*df
      ftmp(3,tido) = ftmp(3,tido) + rij(3)*df
      ftmp(1,j+tid*np) = ftmp(1,j+tid*np) - rij(1)*df
      ftmp(2,j+tid*np) = ftmp(2,j+tid*np) - rij(2)*df
      ftmp(3,j+tid*np) = ftmp(3,j+tid*np) - rij(3)*df

      if (CU) then
        dp = weight* PPairTable(k+1,PTij) + &
           &  (1.d0-weight)* PPairTable(k,PTij)
        ptmp(tid+1) = ptmp(tid+1) + dp
      endif
      if (CV) then
        vtmp(1,tid+1) = vtmp(1,tid+1) + rij(1)*rij(1)*df
        vtmp(2,tid+1) = vtmp(2,tid+1) + rij(2)*rij(2)*df
        vtmp(3,tid+1) = vtmp(3,tid+1) + rij(3)*rij(3)*df
        vtmp(4,tid+1) = vtmp(4,tid+1) + rij(1)*rij(2)*df
        vtmp(5,tid+1) = vtmp(5,tid+1) + rij(1)*rij(3)*df
        vtmp(6,tid+1) = vtmp(6,tid+1) + rij(2)*rij(3)*df
      endif
    end do

  end do

!$OMP barrier
! manual force array reduction using each thread
  do j = 0, thread_num-1
    do i = istart, iend
      force(1,i) = force(1,i) + ftmp(1,i+j*np)
      force(2,i) = force(2,i) + ftmp(2,i+j*np)
      force(3,i) = force(3,i) + ftmp(3,i+j*np)
    enddo
  enddo

  if (CV) then
!$OMP do reduction(+:virial)
    do j = 0, thread_num-1
      virial(1:6) = virial(1:6) + vtmp(1:6,j+1)
    enddo
!$omp end do
  endif
  if (CU) then
!$OMP do reduction(+:potential)
    do j = 0, thread_num-1
      potential = potential + ptmp(j+1)
    enddo
!$omp end do
  endif

!$omp end parallel
  
  return
end !}}}
subroutine compute1( CV, CU)
!{{{**************************************************************************80
!! COMPUTE computes the forces and energies.
!
!  Discussion:
!    The computation of forces and energies is fully parallel.
!    The potential function V(X) is a harmonic well which smoothly
!    saturates to a maximum value at PI/2:
!      v(x) = ( sin ( min ( x, PI2 ) ) )^2
!    The derivative of the potential is:
!      dv(x) = 2.0D+00 * sin ( min ( x, PI2 ) ) * cos ( min ( x, PI2 ) )
!            = sin ( 2.0 * min ( x, PI2 ) )
!
!  Licensing:
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!    Feb 2018
!
!  Author:
!    Original FORTRAN90 version by Bill Magro and later by John Burkardt.
!    This FORTRAN90 version by Ross J. Stewart.
!
!  Parameters:
!    Input, integer ( kind = 4 ) NP, the number of particles.
!    Input, integer ( kind = 4 ) ND, the number of spatial dimensions.
!    Input, real ( kind = 8 ) POS(ND,NP), the position of each particle.
!    Input, real ( kind = 8 ) VEL(ND,NP), the velocity of each particle.
!    Input, real ( kind = 8 ) MASS, the mass of each particle.
!    Output, real ( kind = 8 ) F(ND,NP), the forces.
!    Output, real ( kind = 8 ) POT, the total potential energy.
!    Output, real ( kind = 8 ) KIN, the total kinetic energy.
!
  use omp_lib
  use mostCommon
  use interpot
  use listData
  use control
  implicit none

  real(8) :: d2, rij(3), dp, df, kt, tmpi(3), sr(3), weight, ptmp(8), vtmp(6,8)
  real(4) :: rk, vPTdr2(NPairs)
  integer(4) :: i, jn, j, k, typi, PTij, PTji, istart, iend, tid, tido, mxPTij
  logical, intent(in) :: CV, CU

  if (CV) then
    virial(1:6) = 0.d0
    vtmp(1:6,1:8) = 0.d0
  endif
  if (CU) then
    potential = 0.0D+00
    pene(1:np) = 0.0d+00
    ptmp(1:8) = 0.d0
  endif
  vPTdr2(1:NPairs) = real(1.d0/PTdr2(1:NPairs))

!$omp parallel workshare &
!$omp shared ( force, np )
  force(1:3,1:np) = 0.0D+00
  ftmp(1:3,1:np*thread_num) = 0.0D+00
!$omp end parallel workshare

!$omp parallel &
!$omp shared ( force, np, pos, vel, mass, CV, vPTdr2, PTrMin, potential) &
!$omp shared ( virial, pene, ftmp, vtmp, ptmp, CU ) &
!$omp private ( d2, i, j, k, jn, rij, dp, df, kt, sr, weight, rk, tmpi) &
!$omp private ( typi, tid, tido, istart, iend, PTij, PTji, mxPTij )
tid=omp_get_thread_num()

istart=int(np/thread_num)*tid + 1
iend  =int(np/thread_num)*(tid+1)
if (iend<np.and.tid.eq.thread_num-1) iend=np

  do i = istart, iend !for each particle on this thread...
    tmpi(1) = pos(1,i)
    tmpi(2) = pos(2,i)
    tmpi(3) = pos(3,i)
    typi = typ(i)
    tido = i+tid*np
    do jn = 2, VNL(1,i) ! i particles in cell ic
      j = VNL(jn,i) !first particle ID in cell jc
      PTij = PairMat(typi,(typ(j)))
      PTji = PairMat(typ(j),typi)
      mxPTij = max(PTij,PTji)
      ! calculate distance i-j
      sr(1) = tmpi(1) - pos(1,j)
      sr(2) = tmpi(2) - pos(2,j)
      sr(3) = tmpi(3) - pos(3,j)
      sr(1) = sr(1) - anint(sr(1))
      sr(2) = sr(2) - anint(sr(2))
      sr(3) = sr(3) - anint(sr(3))
      rij(1) = box(1)*sr(1) !real distance units
      rij(2) = box(2)*sr(2) !real distance units
      rij(3) = box(3)*sr(3) !real distance units
      d2 = rij(1)*rij(1)
      d2 = d2 + rij(2)*rij(2)
      d2 = d2 + rij(3)*rij(3)

      if (d2 > PTrc2(mxPTij) ) cycle
!  Attribute half of the potential energy to particle I and J.
      ! table interpolation
      rk = real(d2 - PTrMin)*vPTdr2(mxPTij) + 1.0 ! "continuous" index in table
      k = int(rk)                     ! discretized index
!      if (k < 1) k = 1                ! unlikely but just to protect
      weight = dble(rk) - dble(k)                 ! fractional part, in [0,1]
      ! do  linear
      df = weight* FPairTable(k+1,mxPTij) + &
         &  (1.d0-weight)* FPairTable(k,mxPTij)

      df = df / d2
      if (PTji.gt.0) then
        ftmp(1,tido) = ftmp(1,tido) + rij(1)*df
        ftmp(2,tido) = ftmp(2,tido) + rij(2)*df
        ftmp(3,tido) = ftmp(3,tido) + rij(3)*df
      endif
      if (PTij.gt.0) then
        ftmp(1,j+tid*np) = ftmp(1,j+tid*np) - rij(1)*df
        ftmp(2,j+tid*np) = ftmp(2,j+tid*np) - rij(2)*df
        ftmp(3,j+tid*np) = ftmp(3,j+tid*np) - rij(3)*df
      endif

      if (CU.and.PTij.eq.PTji) then !only count energy if force balanced
        dp = weight* PPairTable(k+1,mxPTij) + &
           &  (1.d0-weight)* PPairTable(k,mxPTij)
        ptmp(tid+1) = ptmp(tid+1) + dp
      endif
      if (CV.and.PTij.eq.PTji) then 
        vtmp(1,tid+1) = vtmp(1,tid+1) + rij(1)*rij(1)*df
        vtmp(2,tid+1) = vtmp(2,tid+1) + rij(2)*rij(2)*df
        vtmp(3,tid+1) = vtmp(3,tid+1) + rij(3)*rij(3)*df
        vtmp(4,tid+1) = vtmp(4,tid+1) + rij(1)*rij(2)*df
        vtmp(5,tid+1) = vtmp(5,tid+1) + rij(1)*rij(3)*df
        vtmp(6,tid+1) = vtmp(6,tid+1) + rij(2)*rij(3)*df
      endif
    end do

  end do

!$OMP barrier
! manual force array reduction using each thread
  do j = 0, thread_num-1
    do i = istart, iend
      force(1,i) = force(1,i) + ftmp(1,i+j*np)
      force(2,i) = force(2,i) + ftmp(2,i+j*np)
      force(3,i) = force(3,i) + ftmp(3,i+j*np)
    enddo
  enddo

  if (CV) then
!$OMP do reduction(+:virial)
    do j = 0, thread_num-1
      virial(1:6) = virial(1:6) + vtmp(1:6,j+1)
    enddo
!$omp end do
  endif
  if (CU) then
!$OMP do reduction(+:potential)
    do j = 0, thread_num-1
      potential = potential + ptmp(j+1)
    enddo
!$omp end do
  endif

!$omp end parallel
  
  return
end !}}}
subroutine compute2( CV, CU)
!{{{**************************************************************************80
!! COMPUTE computes the forces and energies.
!
!  Discussion:
!    The computation of forces and energies is fully parallel.
!    The potential function V(X) is a harmonic well which smoothly
!    saturates to a maximum value at PI/2:
!      v(x) = ( sin ( min ( x, PI2 ) ) )^2
!    The derivative of the potential is:
!      dv(x) = 2.0D+00 * sin ( min ( x, PI2 ) ) * cos ( min ( x, PI2 ) )
!            = sin ( 2.0 * min ( x, PI2 ) )
!
!  Licensing:
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!    Feb 2018
!
!  Author:
!    Original FORTRAN90 version by Bill Magro and later by John Burkardt.
!    This FORTRAN90 version by Ross J. Stewart.
!
!  Parameters:
!    Input, integer ( kind = 4 ) NP, the number of particles.
!    Input, integer ( kind = 4 ) ND, the number of spatial dimensions.
!    Input, real ( kind = 8 ) POS(ND,NP), the position of each particle.
!    Input, real ( kind = 8 ) VEL(ND,NP), the velocity of each particle.
!    Input, real ( kind = 8 ) MASS, the mass of each particle.
!    Output, real ( kind = 8 ) F(ND,NP), the forces.
!    Output, real ( kind = 8 ) POT, the total potential energy.
!    Output, real ( kind = 8 ) KIN, the total kinetic energy.
!
  use omp_lib
  use mostCommon
  use interpot
  use listData
  use control
  implicit none

  real(8) :: d2, rij(3), dp, df, kt, tmpi(3), sr(3), weight, ptmp(8), vtmp(6,8)
  real(8) :: rij0(3), d20, sij, rc2, Vi, r, r0
  real(4) :: rk, vPTdr2(NPairs)
  integer(4) :: i, jn, j, k, typi, PTij, istart, iend, tid, tido
  logical, intent(in) :: CV, CU

  if (CV) then
    virial(1:6) = 0.d0
    vtmp(1:6,1:8) = 0.d0
  endif
  if (CU) then
    potential = 0.0D+00
    pene(1:np) = 0.0d+00
    ptmp(1:8) = 0.d0
  endif
  vPTdr2(1:NPairs) = real(1.d0/PTdr2(1:NPairs))

!$omp parallel workshare &
!$omp shared ( force, np )
  force(1:3,1:np) = 0.0D+00
  ftmp(1:3,1:np*thread_num) = 0.0D+00
!$omp end parallel workshare

!$omp parallel &
!$omp shared ( force, np, pos, vel, mass, CV, vPTdr2, PTrMin, potential, ipos) &
!$omp shared ( virial, pene, ftmp, vtmp, ptmp, CU, initBox, PairName, pairParams ) &
!$omp private ( d2, d20, i, j, k, jn, rij, rij0, dp, df, kt, sr, weight, rk, tmpi) &
!$omp private ( typi, tid, tido, istart, iend, PTij, sij, Vi, rc2, r, r0 )
tid=omp_get_thread_num()

istart=int(np/thread_num)*tid + 1
iend  =int(np/thread_num)*(tid+1)
if (iend<np.and.tid.eq.thread_num-1) iend=np

  do i = istart, iend !for each particle on this thread...
    tmpi(1) = pos(1,i)
    tmpi(2) = pos(2,i)
    tmpi(3) = pos(3,i)
    typi = typ(i)
    tido = i+tid*np
    do jn = 2, VNL(1,i) ! i particles in cell ic
      j = VNL(jn,i) !first particle ID in cell jc
      PTij = PairMat(typi,typ(j))
      ! calculate distance i-j
      sr(1) = tmpi(1) - pos(1,j)
      sr(2) = tmpi(2) - pos(2,j)
      sr(3) = tmpi(3) - pos(3,j)
      sr(1) = sr(1) - anint(sr(1))
      sr(2) = sr(2) - anint(sr(2))
      sr(3) = sr(3) - anint(sr(3))
      rij(1) = box(1)*sr(1) !real distance units
      rij(2) = box(2)*sr(2) !real distance units
      rij(3) = box(3)*sr(3) !real distance units
      d2 = rij(1)*rij(1)
      d2 = d2 + rij(2)*rij(2)
      d2 = d2 + rij(3)*rij(3)
      if (trim(PairName(PTij)).eq."spring") then
      !  ! distance out of cutoff radius: cycle
      !  if (d2 > PTrc2(PTij) ) cycle
        ! calculate squared pair strain
        sr(1) = ipos(1,i) - ipos(1,j)
        sr(2) = ipos(2,i) - ipos(2,j)
        sr(3) = ipos(3,i) - ipos(3,j)
        sr(1) = sr(1) - anint(sr(1))
        sr(2) = sr(2) - anint(sr(2))
        sr(3) = sr(3) - anint(sr(3))
        rij0(1) = initBox(1)*sr(1) !real distance units
        rij0(2) = initBox(2)*sr(2) !real distance units
        rij0(3) = initBox(3)*sr(3) !real distance units
        d20 = rij0(1)*rij0(1)
        d20 = d20 + rij0(2)*rij0(2)
        d20 = d20 + rij0(3)*rij0(3)
        !sij = d2/d20 
        !rc2 = (PairParams(3,PTij)+1.d0)
        !rc2 = rc2*rc2
        ! strain beyond limit: cycle. so no force or energy interaction
        !if (sij > rc2) cycle
        r   = dsqrt(d2)  !current distance
        r0  = dsqrt(d20) !initial distance
        sij = r/r0 -1.d0 !actual strain
        if (sij > PairParams(3,PTij)) cycle
        Vi  = PairParams(4,PTij)
        ! f = C*s*Vi*Vj
        df = PairParams(2,PTij)*sij*Vi*Vi
        if (CU) then
          ! e = 0.5*(C*s*Vi*Vj)*s*r
          dp = 0.5d0*df*sij*r0
          ptmp(tid+1) = ptmp(tid+1) + dp
        endif
        df = -df / r
      else
        ! distance out of cutoff radius: cycle
        if (d2 > PTrc2(PTij) ) cycle
        ! table interpolation
        rk = real(d2 - PTrMin)*vPTdr2(PTij) + 1.0 ! "continuous" index in table
        k = int(rk)                     ! discretized index
!        if (k < 1) k = 1                ! unlikely but just to protect
        weight = dble(rk) - dble(k)                 ! fractional part, in [0,1]
        ! do  linear
        df = weight* FPairTable(k+1,PTij) + &
           &  (1.d0-weight)* FPairTable(k,PTij)
        if (CU) then
          dp = weight* PPairTable(k+1,PTij) + &
             &  (1.d0-weight)* PPairTable(k,PTij)
          ptmp(tid+1) = ptmp(tid+1) + dp
        endif
        df = df / d2
      endif

      ftmp(1,tido) = ftmp(1,tido) + rij(1)*df
      ftmp(2,tido) = ftmp(2,tido) + rij(2)*df
      ftmp(3,tido) = ftmp(3,tido) + rij(3)*df
      ftmp(1,j+tid*np) = ftmp(1,j+tid*np) - rij(1)*df
      ftmp(2,j+tid*np) = ftmp(2,j+tid*np) - rij(2)*df
      ftmp(3,j+tid*np) = ftmp(3,j+tid*np) - rij(3)*df

      if (CV) then
        vtmp(1,tid+1) = vtmp(1,tid+1) + rij(1)*rij(1)*df
        vtmp(2,tid+1) = vtmp(2,tid+1) + rij(2)*rij(2)*df
        vtmp(3,tid+1) = vtmp(3,tid+1) + rij(3)*rij(3)*df
        vtmp(4,tid+1) = vtmp(4,tid+1) + rij(1)*rij(2)*df
        vtmp(5,tid+1) = vtmp(5,tid+1) + rij(1)*rij(3)*df
        vtmp(6,tid+1) = vtmp(6,tid+1) + rij(2)*rij(3)*df
      endif
    end do

  end do

!$OMP barrier
! manual force array reduction using each thread
  do j = 0, thread_num-1
    do i = istart, iend
      force(1,i) = force(1,i) + ftmp(1,i+j*np)
      force(2,i) = force(2,i) + ftmp(2,i+j*np)
      force(3,i) = force(3,i) + ftmp(3,i+j*np)
    enddo
  enddo

  if (CV) then
!$OMP do reduction(+:virial)
    do j = 0, thread_num-1
      virial(1:6) = virial(1:6) + vtmp(1:6,j+1)
    enddo
!$omp end do
  endif
  if (CU) then
!$OMP do reduction(+:potential)
    do j = 0, thread_num-1
      potential = potential + ptmp(j+1)
    enddo
!$omp end do
  endif

!$omp end parallel
  
  return
end !}}}
subroutine compute3( CV, CU)
!{{{**************************************************************************80
!! COMPUTE computes the forces and energies.
!
!  Discussion:
!    The computation of forces and energies is fully parallel.
!    The potential function V(X) is a harmonic well which smoothly
!    saturates to a maximum value at PI/2:
!      v(x) = ( sin ( min ( x, PI2 ) ) )^2
!    The derivative of the potential is:
!      dv(x) = 2.0D+00 * sin ( min ( x, PI2 ) ) * cos ( min ( x, PI2 ) )
!            = sin ( 2.0 * min ( x, PI2 ) )
!
!  Licensing:
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!    Feb 2018
!
!  Author:
!    Original FORTRAN90 version by Bill Magro and later by John Burkardt.
!    This FORTRAN90 version by Ross J. Stewart.
!
!  Parameters:
!    Input, integer ( kind = 4 ) NP, the number of particles.
!    Input, integer ( kind = 4 ) ND, the number of spatial dimensions.
!    Input, real ( kind = 8 ) POS(ND,NP), the position of each particle.
!    Input, real ( kind = 8 ) VEL(ND,NP), the velocity of each particle.
!    Input, real ( kind = 8 ) MASS, the mass of each particle.
!    Output, real ( kind = 8 ) F(ND,NP), the forces.
!    Output, real ( kind = 8 ) POT, the total potential energy.
!    Output, real ( kind = 8 ) KIN, the total kinetic energy.
!
  use omp_lib
  use mostCommon
  use interpot
  use listData
  use control
  implicit none

  real(8) :: d2, rij(3), dp, df, kt, tmpi(3), sr(3), weight, ptmp(8), vtmp(6,8)
  real(8) :: rij0(3), d20, sij, rc2, Vi, r, r0
  real(4) :: rk, vPTdr2(NPairs)
  integer(4) :: i, jn, j, k, typi, PTij, PTji, mxPTij, istart, iend, tid, tido
  logical, intent(in) :: CV, CU

  if (CV) then
    virial(1:6) = 0.d0
    vtmp(1:6,1:8) = 0.d0
  endif
  if (CU) then
    potential = 0.0D+00
    pene(1:np) = 0.0d+00
    ptmp(1:8) = 0.d0
  endif
  vPTdr2(1:NPairs) = real(1.d0/PTdr2(1:NPairs))

!$omp parallel workshare &
!$omp shared ( force, np )
  force(1:3,1:np) = 0.0D+00
  ftmp(1:3,1:np*thread_num) = 0.0D+00
!$omp end parallel workshare

!$omp parallel &
!$omp shared ( force, np, pos, vel, mass, CV, vPTdr2, PTrMin, potential, ipos) &
!$omp shared ( virial, pene, ftmp, vtmp, ptmp, CU, initBox, PairName, pairParams ) &
!$omp private ( d2, d20, i, j, k, jn, rij, rij0, dp, df, kt, sr, weight, rk, tmpi) &
!$omp private ( typi, tid, tido, istart, iend, PTij, PTji, mxPTij, sij, Vi, rc2, r, r0 )
tid=omp_get_thread_num()

istart=int(np/thread_num)*tid + 1
iend  =int(np/thread_num)*(tid+1)
if (iend<np.and.tid.eq.thread_num-1) iend=np

  do i = istart, iend !for each particle on this thread...
    tmpi(1) = pos(1,i)
    tmpi(2) = pos(2,i)
    tmpi(3) = pos(3,i)
    typi = typ(i)
    tido = i+tid*np
    do jn = 2, VNL(1,i) ! i particles in cell ic
      j = VNL(jn,i) !first particle ID in cell jc
      PTij = PairMat(typi,typ(j))
      PTji = PairMat(typ(j),typi)
      mxPTij = max(PTij,PTji)
      ! calculate distance i-j
      sr(1) = tmpi(1) - pos(1,j)
      sr(2) = tmpi(2) - pos(2,j)
      sr(3) = tmpi(3) - pos(3,j)
      sr(1) = sr(1) - anint(sr(1))
      sr(2) = sr(2) - anint(sr(2))
      sr(3) = sr(3) - anint(sr(3))
      rij(1) = box(1)*sr(1) !real distance units
      rij(2) = box(2)*sr(2) !real distance units
      rij(3) = box(3)*sr(3) !real distance units
      d2 = rij(1)*rij(1)
      d2 = d2 + rij(2)*rij(2)
      d2 = d2 + rij(3)*rij(3)
      if (trim(PairName(mxPTij)).eq."spring") then
      !  ! distance out of cutoff radius: cycle
      !  if (d2 > PTrc2(PTij) ) cycle
        ! calculate squared pair strain
        sr(1) = ipos(1,i) - ipos(1,j)
        sr(2) = ipos(2,i) - ipos(2,j)
        sr(3) = ipos(3,i) - ipos(3,j)
        sr(1) = sr(1) - anint(sr(1))
        sr(2) = sr(2) - anint(sr(2))
        sr(3) = sr(3) - anint(sr(3))
        rij0(1) = initBox(1)*sr(1) !real distance units
        rij0(2) = initBox(2)*sr(2) !real distance units
        rij0(3) = initBox(3)*sr(3) !real distance units
        d20 = rij0(1)*rij0(1)
        d20 = d20 + rij0(2)*rij0(2)
        d20 = d20 + rij0(3)*rij0(3)
        !sij = d2/d20 
        !rc2 = (PairParams(3,PTij)+1.d0)
        !rc2 = rc2*rc2
        ! strain beyond limit: cycle. so no force or energy interaction
        !if (sij > rc2) cycle
        r   = dsqrt(d2)  !current distance
        r0  = dsqrt(d20) !initial distance
        sij = r/r0 -1.d0 !actual strain
        if (sij > PairParams(3,mxPTij)) cycle
        Vi  = PairParams(4,mxPTij)
        ! f = C*s*Vi*Vj
        df = PairParams(2,mxPTij)*sij*Vi*Vi
        if (CU.and.PTji.eq.PTij) then
          ! e = 0.5*(C*s*Vi*Vj)*s*r
          dp = 0.5d0*df*sij*r0
          ptmp(tid+1) = ptmp(tid+1) + dp
        endif
        df = -df / r
      else
        ! distance out of cutoff radius: cycle
        if (d2 > PTrc2(mxPTij) ) cycle
        ! table interpolation
        rk = real(d2 - PTrMin)*vPTdr2(mxPTij) + 1.0 ! "continuous" index in table
        k = int(rk)                     ! discretized index
        if (k < 1) then
          k = 1                ! unlikely but just to protect
          write(0,*) "ERROR: WTF is k<1?",step,d2,i,j,pos(:,i),pos(:,j), rij, mxPTij, typ(i), typ(j), mass(typ(j))
        endif
        weight = dble(rk) - dble(k)                 ! fractional part, in [0,1]
        ! do  linear
        df = weight* FPairTable(k+1,mxPTij) + &
           &  (1.d0-weight)* FPairTable(k,mxPTij)
        if (CU.and.PTij.eq.PTji) then
          dp = weight* PPairTable(k+1,mxPTij) + &
             &  (1.d0-weight)* PPairTable(k,mxPTij)
          ptmp(tid+1) = ptmp(tid+1) + dp
        endif
        df = df / d2
      endif

      if (PTji.gt.0) then
        ftmp(1,tido) = ftmp(1,tido) + rij(1)*df
        ftmp(2,tido) = ftmp(2,tido) + rij(2)*df
        ftmp(3,tido) = ftmp(3,tido) + rij(3)*df
      endif
      if (PTij.gt.0) then
        ftmp(1,j+tid*np) = ftmp(1,j+tid*np) - rij(1)*df
        ftmp(2,j+tid*np) = ftmp(2,j+tid*np) - rij(2)*df
        ftmp(3,j+tid*np) = ftmp(3,j+tid*np) - rij(3)*df
      endif

      if (CV.and.PTij.eq.PTji) then
        vtmp(1,tid+1) = vtmp(1,tid+1) + rij(1)*rij(1)*df
        vtmp(2,tid+1) = vtmp(2,tid+1) + rij(2)*rij(2)*df
        vtmp(3,tid+1) = vtmp(3,tid+1) + rij(3)*rij(3)*df
        vtmp(4,tid+1) = vtmp(4,tid+1) + rij(1)*rij(2)*df
        vtmp(5,tid+1) = vtmp(5,tid+1) + rij(1)*rij(3)*df
        vtmp(6,tid+1) = vtmp(6,tid+1) + rij(2)*rij(3)*df
      endif
    end do

  end do

!$OMP barrier
! manual force array reduction using each thread
  do j = 0, thread_num-1
    do i = istart, iend
      force(1,i) = force(1,i) + ftmp(1,i+j*np)
      force(2,i) = force(2,i) + ftmp(2,i+j*np)
      force(3,i) = force(3,i) + ftmp(3,i+j*np)
    enddo
  enddo

  if (CV) then
!$OMP do reduction(+:virial)
    do j = 0, thread_num-1
      virial(1:6) = virial(1:6) + vtmp(1:6,j+1)
    enddo
!$omp end do
  endif
  if (CU) then
!$OMP do reduction(+:potential)
    do j = 0, thread_num-1
      potential = potential + ptmp(j+1)
    enddo
!$omp end do
  endif

!$omp end parallel
  
  return
end !}}}
subroutine randomPositions
!{{{**************************************************************************80
!
!! INITIALIZE initializes the positions, velocities, and accelerations.
!
!  Licensing:
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!    21 November 2007
!
!  Author:
!    Original FORTRAN90 version by Bill Magro.
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!    Input, integer ( kind = 4 ) NP, the number of particles.
!    Input, integer ( kind = 4 ) ND, the number of spatial dimensions.
!    Input, real ( kind = 8 ) BOX(ND), specifies the maximum position
!    of particles in each dimension.
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!    Output, real ( kind = 8 ) POS(ND,NP), the position of each particle.
!    Output, real ( kind = 8 ) VEL(ND,NP), the velocity of each particle.
!    Output, real ( kind = 8 ) ACC(ND,NP), the acceleration of each particle.
!
  use mostCommon
  implicit none
  integer ( kind = 8 ) j
  write(6,'(A,3f10.4)') "# using system box: ",box(1:3) 
!
!  Start by setting the positions to random numbers between 0 and 1.
!
    call random_number ( harvest = pos(1:3,1:np) )
!
!  Use these random values as scale factors to pick random locations
!  inside the box.
!
!$omp parallel do &
!$omp shared ( box, np, pos ) &
!$omp private ( j )
    do j = 1, np
      !pos(1:nd,j) = box(1:nd) * pos(1:nd,j)
      pos(1,j) = pos(1,j) - 0.5d0 !to range from [-0.5:0.5]
      pos(2,j) = pos(2,j) - 0.5d0 !to range from [-0.5:0.5]
      pos(3,j) = pos(3,j) - 0.5d0 !to range from [-0.5:0.5]
    end do
!$omp end parallel do

end !}}}
subroutine integrate1
!{{{**************************************************************************80
!! UPDATE updates positions, velocities and accelerations.
!
!  Discussion:
!    The time integration is fully parallel.
!    A velocity Verlet algorithm is used for the updating.
!    x(t+dt) = x(t) + v(t) * dt + 0.5 * a(t) * dt * dt
!    v(t+dt) = v(t) + 0.5 * ( a(t) + a(t+dt) ) * dt
!    a(t+dt) = f(t) / m
!
!  Licensing:
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!    21 November 2007
!
!  Author:
!    Original FORTRAN90 version by Bill Magro.
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!    Input, integer ( kind = 4 ) NP, the number of particles.
!    Input, integer ( kind = 4 ) ND, the number of spatial dimensions.
!    Input/output, real ( kind = 8 ) POS(ND,NP), the position of each particle.
!    Input/output, real ( kind = 8 ) VEL(ND,NP), the velocity of each particle.
!    Input, real ( kind = 8 ) F(ND,NP), the force on each particle.
!    Input/output, real ( kind = 8 ) ACC(ND,NP), the acceleration of each
!    particle.
!    Input, real ( kind = 8 ) MASS, the mass of each particle.
!    Input, real ( kind = 8 ) DT, the time step.
  use mostCommon
  use interpot, only : mass !for atom masses
  use listData, only : dispList
  use groupStuff
  implicit none
  integer :: j, i
  real(8) :: li(3), du(3)

  li(1) = 1.d0/box(1)
  li(2) = 1.d0/box(2)
  li(3) = 1.d0/box(3)

IF (Nhold.eq.0) then ! if no domain being held (runs faster than checking)

!$omp parallel do &
!$omp shared ( dt, np, pos, vel, mass, li, dispList ) &
!$omp private ( j, du )
  do j = 1, np
    ! calculate discrete displacement interval (du)
    du(1) = li(1)*(vel(1,j) * dt + 0.5D+00 * acc(1,j) * dt * dt)
    du(2) = li(2)*(vel(2,j) * dt + 0.5D+00 * acc(2,j) * dt * dt)
    du(3) = li(3)*(vel(3,j) * dt + 0.5D+00 * acc(3,j) * dt * dt)
    ! accumulate displacements to check for neighbour list updates
    dispList(1,j) = dispList(1,j) + du(1)
    dispList(2,j) = dispList(2,j) + du(2)
    dispList(3,j) = dispList(3,j) + du(3)
    ! update positions
    pos(1,j) = pos(1,j) + du(1)
    pos(2,j) = pos(2,j) + du(2)
    pos(3,j) = pos(3,j) + du(3)
  end do
!$omp end parallel do

ELSE

!$omp parallel do &
!$omp shared ( dt, np, pos, vel, mass, li, dispList, hold, bcTyp, bcDir, bcParm ) &
!$omp private ( j, i, du )
  do j = 1, np
    ! calculate discrete displacement interval (du)
    du(1) = li(1)*(vel(1,j) * dt + 0.5D+00 * acc(1,j) * dt * dt)
    du(2) = li(2)*(vel(2,j) * dt + 0.5D+00 * acc(2,j) * dt * dt)
    du(3) = li(3)*(vel(3,j) * dt + 0.5D+00 * acc(3,j) * dt * dt)
    i = hold(j)
    if (i.ne.0) then !this has some BC
      if (bcTyp(i).eq.'d') then ! displacement conditions
        if (bcDir(1,i)) du(1) = bcParm(1,i) ! which directions are constrained?
        if (bcDir(2,i)) du(2) = bcParm(2,i)
        if (bcDir(3,i)) du(3) = bcParm(3,i)
      endif
    endif
    ! accumulate displacements to check for neighbour list updates
    dispList(1,j) = dispList(1,j) + du(1)
    dispList(2,j) = dispList(2,j) + du(2)
    dispList(3,j) = dispList(3,j) + du(3)
    ! update positions
    pos(1,j) = pos(1,j) + du(1)
    pos(2,j) = pos(2,j) + du(2)
    pos(3,j) = pos(3,j) + du(3)
  end do
!$omp end parallel do

ENDIF

  return
end !}}}
subroutine integrate2( CK )
! {{{
  use mostCommon
  use interpot, only : mass
  implicit none
  integer :: j
  real(8) :: rmass, kt
  logical :: CK
  
  kinetic = 0.d0
!$omp parallel do &
!$omp shared ( acc, dt, force, np, vel, mass, CK, typ ) &
!$omp private ( j, rmass, kt ) &
!$omp reduction (+:kinetic)
  do j = 1, np
    rmass = 1.0d0/mass(typ(j))
    ! update velocities (using force/mass as the 'old' acceleration to average)
    vel(1,j) = vel(1,j) + 0.5D+00 * dt * ( force(1,j) * rmass + acc(1,j) )
    vel(2,j) = vel(2,j) + 0.5D+00 * dt * ( force(2,j) * rmass + acc(2,j) )
    vel(3,j) = vel(3,j) + 0.5D+00 * dt * ( force(3,j) * rmass + acc(3,j) )
    ! update accelerations
    acc(1,j) = force(1,j) * rmass
    acc(2,j) = force(2,j) * rmass
    acc(3,j) = force(3,j) * rmass
    if (CK) then
!     Compute the kinetic energy when we only loop over i once
      kt =      vel(1,j)*vel(1,j)
      kt = kt + vel(2,j)*vel(2,j) 
      kt = kt + vel(3,j)*vel(3,j) 
      kt = kt * mass(typ(j))
      kinetic = kinetic + kt
!if (kt.ne.kt) then
!write(0,*) "ERROR: kinetics:",mass(typ(j)), typ(j), vel(:,j), rmass, force(:,j)
!endif
    endif
  end do
!$omp end parallel do
  if (CK) kinetic = kinetic * 0.5D+00

  return
end !}}}
subroutine mkIsoKE
!{{{
! scale velocities using berendsen "thermo"stat for constant kinetic energy
  use mostCommon
  use control
  implicit none
  real(8) :: lam
  integer :: j !, gt
!  character(80) :: gN
!  logical :: ss

  lam = 1.d0 + dt/KEparm * (KEtarget*dble(np)/kinetic - 1.d0)
  lam = dsqrt(lam)

! using a system subset?
!if (trim(gN).ne."N/A") then
!   ss=.true.
!   do i = 1, Ngrp
!     if (trim(groupname(i)).eq.trim(gN)) then
!       gt = i; endif; enddo
!else; ss = .false.; endif
  
!$omp parallel do &
!$omp shared ( vel, lam, np ) &
!$omp private ( j )
  do j = 1, np
!    if (ss) then !use a subset of particle groups?
!      if (grp(j).ne.gt) cycle
!    endif
    vel(1,j) = vel(1,j)*lam
    vel(2,j) = vel(2,j)*lam
    vel(3,j) = vel(3,j)*lam
  enddo
!$omp end parallel do

  return
end !}}}
subroutine mkIsoE
!{{{
! scale velocities using berendsen "thermo"stat for constant total energy
  use mostCommon
  use control
  implicit none
  real(8) :: lam
  integer :: j

  lam = 1.d0 - dt/Eparm * (e0/(potential+kinetic) - 1.d0)
  !lam = dsqrt(lam)
  
!$omp parallel do &
!$omp shared ( vel, lam, np ) &
!$omp private ( j )
  do j = 1, np
    vel(1,j) = vel(1,j)*lam
    vel(2,j) = vel(2,j)*lam
    vel(3,j) = vel(3,j)*lam
  enddo
!$omp end parallel do

  return
end !}}}
subroutine mkIsoP
!{{{
! scale boxsize using berendsen barostat for constant pressure
  use mostCommon
  use control
  implicit none
  real(8) :: lam(6), press(6), vol

  vol = box(1)*box(2)*box(3)

  ! if this uses units of eV and Ang, then multiply this by 160.21773 to get GPa
  press(1) = virial(1)/vol
  press(2) = virial(2)/vol
  press(3) = virial(3)/vol
  press(4) = virial(4)/vol
  press(5) = virial(5)/vol
  press(6) = virial(6)/vol

  lam(1) = 1.d0 - dt/Pparm * (Ptarget(1) - press(1))
  lam(2) = 1.d0 - dt/Pparm * (Ptarget(2) - press(2))
  lam(3) = 1.d0 - dt/Pparm * (Ptarget(3) - press(3))
  lam(4) = dt/Pparm * press(4) !because Ptarget(4:6)=0 due to orthogonal cell shape
  lam(5) = dt/Pparm * press(5)
  lam(6) = dt/Pparm * press(6)
  
  box(1) = box(1)*lam(1) + box(2)*lam(6) + box(3)*lam(5)  
  box(2) = box(1)*lam(6) + box(2)*lam(2) + box(3)*lam(4)
  box(3) = box(1)*lam(5) + box(2)*lam(4) + box(3)*lam(3)

  return
end !}}}
subroutine buildCellList( outs )
!{{{ 
! Bin each particle into a cell, link them together via a List that points to the next
! particle in the cell.
! Head(cell) gives the ID of the first particle in 'cell'. List(Head(cell)) gives the
! particle ID of the next particle in that cell. When it is '0' that's the end
  use mostCommon
  use interpot
  use listData
  implicit none
  integer :: i
  integer(4) :: IXYZ(3), ICELL
  real(8) :: sr(3)
  logical :: outs

!$omp parallel workshare &
!$omp shared ( HEAD, List, NCELL, np )
  HEAD(1:NCELL) = 0
  List(1:np) = 0
!$omp end parallel workshare

!omp parallel do &
!omp shared ( pos, MCELL, np, NCELL) &
!omp private ( i, IXYZ, ICELL, sr )
  DO i = 1, np
      ! pos(j,i) is in [-.5:.5] so +0.5 to be in [0:1]
    sr(1) = pos(1,i)-anint(pos(1,i))
    sr(2) = pos(2,i)-anint(pos(2,i))
    sr(3) = pos(3,i)-anint(pos(3,i))
    IXYZ(1) = Int( (sr(1)+0.4999d0)*dble(MCELL(1)) )
    IXYZ(2) = Int( (sr(2)+0.4999d0)*dble(MCELL(2)) )
    IXYZ(3) = Int( (sr(3)+0.4999d0)*dble(MCELL(3)) )
    ICELL   = 1+ IXYZ(1)+ IXYZ(2)*MCELL(1)+ IXYZ(3)*MCELL(1)*MCELL(2)
    !if(ICELL.lt.1 .or. ICELL.gt.NCELL) then
    ! write(0,*) "ERROR: CLL: ICELL problem, particle Outside of simulation supercell",&
    !         &   ICELL, " of ( 1,",NCELL,")", IXYZ, MCELL
    ! write(0,*) "pos", pos(1:nd,i), "particle:",i
    ! call flush(0); STOP
    !endif
    List(i) = HEAD(ICELL)
    HEAD(ICELL) = i
  ENDDO
!omp end parallel do

call buildVerletList( outs )

end subroutine buildCellList !}}}
subroutine getJcellNP( IXYZ, k, jc )
! {{{ 
! Returns JCELL ID given an ICELL XYZ vector and neighbour cell ID, k
use listData
implicit none
integer :: IXYZ(3), k, jc, t(3)
real(4) :: tmp
select case (k)
  case(1); t=(/0,0,0/)
  case(2); t=(/0,0,1/)
  case(3); t=(/0,1,0/)
  case(4); t=(/0,1,-1/) !E3
  case(5); t=(/0,1,1/) ! E3
  case(6); t=(/1,0,0/)
  case(7); t=(/1,0,-1/) ! A3
  case(8); t=(/1,0,1/) !  A3
  case(9); t=(/1,1,0/) ! D2
  case(10); t=(/1,1,-1/) ! B3 F2
  case(11); t=(/1,1,1/) !  B3 G2
  case(12); t=(/1,-1,0/) ! D2
  case(13); t=(/1,-1,-1/) ! C3 F2
  case(14); t=(/1,-1,1/) !  C3 G2
  case default; t=(/0,0,0/)
end select
!if (MCELL(2).eq.2) then
!  if (t(2).eq.-1) then; jc = -1; RETURN; endif
!endif
!if (MCELL(3).eq.2) then
!  if (t(3).eq.-1) then; jc = -1; RETURN; endif
!endif
t(1) = IXYZ(1) + t(1)
t(2) = IXYZ(2) + t(2)
t(3) = IXYZ(3) + t(3)
! check for cells out of box, and fold them in.
tmp = (real(t(1)+1)/real(MCELL(1)+1)-0.5)
t(1) = t(1) -MCELL(1)*nint(tmp)
tmp = (real(t(2)+1)/real(MCELL(2)+1)-0.5)
t(2) = t(2) -MCELL(2)*nint(tmp)
tmp = (real(t(3)+1)/real(MCELL(3)+1)-0.5)
t(3) = t(3) -MCELL(3)*nint(tmp)
!return JCELL
jc = 1+ t(1)+ t(2)*MCELL(1)+ t(3)*MCELL(1)*MCELL(2)
end subroutine !}}}
subroutine getCellVector( I, VEC, fold )
! {{{
! return a cell vector given the cell ID
use listData
implicit none
integer, intent(in) :: I
integer, intent(out) ::  VEC(1:3)
logical, intent(in) :: fold
real(4) :: tmp
    VEC = 0
    VEC(1) = mod(I-1,MCELL(1))
    VEC(3) = int((I-1)/(MCELL(1)*MCELL(2)))
    !VEC(2) = int((I-VEC(3)*MCELL(1)*MCELL(3)-1-VEC(1))/MCELL(1))
    VEC(2) = int(mod(I-1,MCELL(1)*MCELL(2))/MCELL(1))
if (fold) then
  ! check for cells out of box, and fold them in.
  tmp = (real(VEC(1)+1)/real(MCELL(1)+1)-0.5)
  VEC(1) = VEC(1) -MCELL(1)*nint(tmp)
  tmp = (real(VEC(2)+1)/real(MCELL(2)+1)-0.5)
  VEC(2) = VEC(2) -MCELL(2)*nint(tmp)
  tmp = (real(VEC(3)+1)/real(MCELL(3)+1)-0.5)
  VEC(3) = VEC(3) -MCELL(3)*nint(tmp)
endif
end subroutine !}}}
subroutine buildVerletList( outs )
! {{{
! build a list per atom i of all j neighbours
use mostCommon
use listData
use interpot
implicit none
integer :: nj, dreiPnd, IXYZ(3), ic, jc, i, j, k, typi, PTij, PTji
real(8) :: d2, rij(3), tmpi(3), sr(3), rcuts2
logical :: outs

  dreiPnd = 14 !(3**3-1)/2+1
!omp parallel do &
!omp shared ( np, pos, dreiPnd, list, dispList, skin, VNL ) &
!omp private ( d2, i, j, k, rij, IXYZ, sr, tmpi, typi, rcuts2, PTij, PTji, nj, jc )
  do ic = 1, NCELL !for each cell in system, do...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first i-j loop is within the same cell, ic=jc
    i = HEAD(ic) !first particle ID in cell ic
    ! first neighbour cell is same cell.
    jc = ic
    do while (i > 0) ! i particles in cell ic
      typi = typ(i)
      VNL(1:NeiMax,i) = 0
      nj = 1
      tmpi(1) = pos(1,i)
      tmpi(2) = pos(2,i)
      tmpi(3) = pos(3,i)
      j = HEAD(jc) !first particle ID in cell jc
      do while (j > 0)
        if (i >= j) then; j = List(j); cycle; endif
        PTij = PairMat( typi, typ(j))
        PTji = PairMat( typ(j), typi)
        if (max(PTij,PTji).le.0) then; j = List(j); cycle; endif ! no interaction potential 
        ! calculate distance i-j
        sr(1) = tmpi(1) - pos(1,j)
        sr(2) = tmpi(2) - pos(2,j)
        sr(3) = tmpi(3) - pos(3,j)
        sr(1) = sr(1)-anint(sr(1))
        sr(2) = sr(2)-anint(sr(2))
        sr(3) = sr(3)-anint(sr(3))
        rij(1) = box(1)*sr(1) !real distance units
        rij(2) = box(2)*sr(2) !real distance units
        rij(3) = box(3)*sr(3) !real distance units
        d2 = rij(1)*rij(1)
        d2 = d2 + rij(2)*rij(2)
        d2 = d2 + rij(3)*rij(3)
        rcuts2 = skin + PairParams(1, max(PTij,PTji) )
        if (d2 < PTrMin) then; j=List(j); cycle; endif !neighbour on top of this
        if (d2 > rcuts2*rcuts2) then; j=List(j); cycle; endif
        if (nj>NeiMax) STOP 'too many neighbours in list'
        VNL(nj+1,i) = j
        nj = nj + 1
        j = List(j) !get next j particle in cell jc
      end do
      VNL(1,i) = nj !this is actually one greater than number of neighbours
      dispList(1:3,i) = 0.d0 !reset i's displacement since last list build
      i = List(i) ! get next i particle in cell ic
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! second i-j loop with various different jcells.
    ! the neighbouring cells
    call getCellVector( ic, IXYZ, .false. )
    do k = 2, dreiPnd ! neighbour cell search pattern ID
      call getJcellNP( IXYZ, k, jc )
      if (jc.eq.-1) cycle 
      i = HEAD(ic) !first particle ID in cell ic
      do while (i > 0) ! i particles in cell ic
        typi = typ(i)
        nj = VNL(1,i) 
        tmpi(1) = pos(1,i)
        tmpi(2) = pos(2,i)
        tmpi(3) = pos(3,i)
        j = HEAD(jc) !first particle ID in cell jc
        do while (j > 0)
          PTij = PairMat( typi, typ(j))
          PTji = PairMat( typ(j), typi)
          if (max(PTij,PTji).le.0) then; j = List(j); cycle; endif ! no interaction potential 
          ! calculate distance i-j
          sr(1) = tmpi(1) - pos(1,j)
          sr(2) = tmpi(2) - pos(2,j)
          sr(3) = tmpi(3) - pos(3,j)
          sr(1) = sr(1)-anint(sr(1))
          sr(2) = sr(2)-anint(sr(2))
          sr(3) = sr(3)-anint(sr(3))
          rij(1) = box(1)*sr(1) !real distance units
          rij(2) = box(2)*sr(2) !real distance units
          rij(3) = box(3)*sr(3) !real distance units
          d2 = rij(1)*rij(1)
          d2 = d2 + rij(2)*rij(2)
          d2 = d2 + rij(3)*rij(3)
          rcuts2 = skin + PairParams(1, max(PTij,PTji) )
          if (d2 > rcuts2*rcuts2) then; j=List(j); cycle; endif
          if (nj>NeiMax) STOP 'too many neighbours in list'
          VNL(nj+1,i) = j
          nj = nj + 1
          j = List(j) !get next j particle in cell jc
        end do
        VNL(1,i) = nj
        dispList(1:3,i) = 0.d0 !reset i's displacement since last list build
        i = List(i) ! get next i particle in cell ic
      end do
    end do
  end do
!omp end parallel do

if (outs) then
  write(6,'(A,2i4)') "# min and max number of neighbours+skin:", &
    & minval(VNL(1,:)), maxval(VNL(1,:))
  write(6,'(A,3f10.6)') "# position of atom with least neighbours: ",pos(:,minloc(VNL(1,:)))
  write(6,'(A,3f10.6)') "# position of atom with most neighbours:  ",pos(:,maxloc(VNL(1,:)))
!open(32,file='derp')
!do i = 1, np
!write(32,*) pos(:,i)*box(:), VNL(1,i)
!enddo
!close(32)
endif

end !}}}
subroutine CheckListUpdate
! {{{
! add the two largest displacements together since the last neighbour list update
! and if it's larger than the SKIN then call a list rebuild.
use omp_lib
use mostCommon, only : np
use listData
use control
implicit none
real(8) :: disp
real(8), dimension(2*thread_num) :: displ
integer :: istart, iend, tid, i, t1, t2, dl(1)

!$omp parallel &
!$omp shared(displ, thread_num, np) &
!$omp private(tid, istart, iend, disp, t1, t2)
tid = omp_get_thread_num()
t1 = 2*(tid+1)-1
t2 = 2*(tid+1)
displ(t1) = 0.d0
displ(t2) = 0.d0

istart=int(np/thread_num)*tid + 1
iend  =int(np/thread_num)*(tid+1)
if (iend<np.and.tid.eq.thread_num-1) iend=np

do i = istart, iend
  disp =        dispList(1,i)*dispList(1,i)
  disp = disp + dispList(2,i)*dispList(2,i)
  disp = disp + dispList(3,i)*dispList(3,i)
  if (disp >= displ(t1)) then
    displ(t2) = displ(t1)  ! push current 1st into 2nd place
    displ(t1) = disp       ! and put this one into current 1st
  elseif (disp >= displ(t2)) then
    displ(t2) = disp
  endif
enddo
!$omp end parallel
!$omp barrier

! just find the top two max values to sum
disp = maxval(displ(1:2*thread_num))  ! max value
dl = maxloc(displ(1:2*thread_num))    ! max location
displ(dl(1)) = 0.d0                      ! zero this to check 2nd max
disp = disp + maxval(displ(1:2*thread_num))  !add the 2nd max to the first

if (disp > skinsq) then
  if (disp > 1.0)  STOP "Particle outside of box"
  NlistUpdates = NlistUpdates + 1 !count how many times this list updates
  call buildCellList( .false. ) !don't output stats during run
endif

end !}}}
subroutine normalVelocityDistribution( kt, gN )
! {{{
! create a velocity distribution that is normally distributed
! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.
use omp_lib
use mostCommon
use control
use interpot, only : mass
use groupStuff
implicit none
integer :: i
integer :: istart, iend, tid, gt
real(4) :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472 
real(4) :: r1 = 0.27597, r2 = 0.27846, half = 0.5
real(4) :: u, v, x, y, q, nn(3), coeff, kt
real(8) :: vtot(3,thread_num), vsum(3)
logical :: ss
character(80) :: gN

! using a system subset?
if (trim(gN).ne."N/A") then
   ss=.true.
   gt = groupID(gN) 
else; ss = .false.; endif

vtot(1:3,1:thread_num) = 0.d0

!$omp parallel &
!$omp shared (thread_num, np, s, t, a, b, r1, r2, half, kt, vel, gt, grp, vtot) &
!$omp private (tid, istart, iend, u, v, y, q, coeff) 
tid = omp_get_thread_num() + 1

istart=int(np/thread_num)*(tid-1) + 1
iend  =int(np/thread_num)*tid
if (iend<np.and.tid.eq.thread_num) iend=np

do i = istart, iend
  if (ss) then !use a subset of particle groups?
    if (grp(i).ne.gt) cycle
  endif
  do 
    call random_number(u)
    call random_number(v)
    v = 1.7156 * (v - half)
    ! evaluate the quadratic form
    x = u - s
    y = abs(v) - t
    q = x*x + y*(a*y - b*x)
    ! accept P if inside inner ellipse
    if (q < r1) exit
    ! reject P if outside outer ellipse
    if (q > r2) cycle
    ! reject P if outside acceptance region
    if (v*v < -4.0*log(u)*u*u) exit
  enddo
  nn(1) = v/u
  do 
    call random_number(u)
    call random_number(v)
    v = 1.7156 * (v - half)
    ! evaluate the quadratic form
    x = u - s
    y = abs(v) - t
    q = x*x + y*(a*y - b*x)
    ! accept P if inside inner ellipse
    if (q < r1) exit
    ! reject P if outside outer ellipse
    if (q > r2) cycle
    ! reject P if outside acceptance region
    if (v*v < -4.0*log(u)*u*u) exit
  enddo
  nn(2) = v/u
  do 
    call random_number(u)
    call random_number(v)
    v = 1.7156 * (v - half)
    ! evaluate the quadratic form
    x = u - s
    y = abs(v) - t
    q = x*x + y*(a*y - b*x)
    ! accept P if inside inner ellipse
    if (q < r1) exit
    ! reject P if outside outer ellipse
    if (q > r2) cycle
    ! reject P if outside acceptance region
    if (v*v < -4.0*log(u)*u*u) exit
  enddo
  nn(3) = v/u

  coeff = sqrt(2.0*kt/(3.0*real(mass(typ(i)))))
  vel(1,i) = dble(nn(1)*coeff)
  vel(2,i) = dble(nn(2)*coeff)
  vel(3,i) = dble(nn(3)*coeff)
  vtot(1,tid) = vtot(1,tid) + vel(1,i)
  vtot(2,tid) = vtot(2,tid) + vel(2,i)
  vtot(3,tid) = vtot(3,tid) + vel(3,i)
enddo
!$omp end parallel

vsum(1) = sum(vtot(1,1:thread_num))/dble(np)
vsum(2) = sum(vtot(2,1:thread_num))/dble(np)
vsum(3) = sum(vtot(3,1:thread_num))/dble(np)

! subtract off the average drift velocity
!
!$omp parallel do &
!$omp shared ( np, vel, vsum ) &
!$omp private ( i )
do i = 1, np
  vel(1,i) = vel(1,i) - vsum(1)
  vel(2,i) = vel(2,i) - vsum(2)
  vel(3,i) = vel(3,i) - vsum(3)
enddo
!$omp end parallel do

end !}}}
subroutine setGroup( g, t, dp )
!{{{
! set group ID to either an entire atom type or geometric domain
! t   is the type of group setting [typ|dom]
! dp  is a string of parameters for these types, either a type ID or a domain
  use domType
  use mostCommon, only : np, typ, pos
  use groupStuff
  implicit none
  integer :: i, cht
  character :: g
  character(80) :: t, dp
  type (domain) :: dom
  
  call left_of("#",dp)
  if (t(1:3).eq."typ") then
    read(dp,*) cht
    if ( g == 'g' ) then
!$omp parallel do &
!$omp shared ( cht, grp, Ngrp, typ, np) &
!$omp private ( i )
      do i = 1, np
        if (typ(i).eq.cht) grp(i) = Ngrp
      enddo
!$omp end parallel do
    elseif ( g == 'h') then
!$omp parallel do &
!$omp shared ( cht, hold, Nhold, typ, np) &
!$omp private ( i )
      do i = 1, np
        if (typ(i).eq.cht) hold(i) = Nhold
      enddo
!$omp end parallel do
    endif
  elseif (t(1:3).eq."dom") then
    call readDomain( dp, dom )
    if ( g == 'g' ) then
!$omp parallel do &
!$omp shared ( grp, Ngrp, dom, np, pos ) &
!$omp private ( i )
      do i = 1, np
        if (inDomain( real(pos(1:3,i)), dom )) grp(i) = Ngrp
      enddo
!$omp end parallel do
    elseif ( g == 'h' ) then
!$omp parallel do &
!$omp shared ( hold, Nhold, dom, np, pos ) &
!$omp private ( i )
      do i = 1, np
        if (inDomain( real(pos(1:3,i)), dom )) hold(i) = Nhold
      enddo
!$omp end parallel do
    endif
  elseif (t(1:3).eq."all") then ! clear all groups or holds
    if ( g == 'g' ) then
!$omp parallel do &
!$omp shared ( grp, np ) &
!$omp private ( i )
      do i = 1, np
        grp(i) = 0
      enddo
!$omp end parallel do
    elseif ( g == 'h' ) then
!$omp parallel do &
!$omp shared ( hold, np ) &
!$omp private ( i )
      do i = 1, np
        hold(i) = 0
      enddo
!$omp end parallel do
    endif

  endif

end subroutine !}}}
subroutine setType( t, dp )
!{{{
! set atoms in group ID to a different atom type
! t   is the group name
! dp  is a string of the new atom type
  use mostCommon, only : typ, np
  use groupStuff
  use interpot
  implicit none
  integer :: i, gID, tID
  character(80) :: t, dp
  character(8) :: cID

  ! find the group ID based on the group Name 
  do i = 1, Ngrp
    if (trim(t).eq.groupName(i)) then
      gID = i; exit
    endif
  enddo
  cID = trim(dp)
  tID = typeID( cID )
 
!$omp parallel do &
!$omp shared ( gID, grp, typ, np, tID) &
!$omp private ( i )
  do i = 1, np
    if (grp(i).eq.gID) typ(i) = tID
  enddo
!$omp end parallel do

end subroutine !}}}
subroutine initPositionSet
!{{{
! save the current system positions to be used as initial positions used for spring potentials
  use mostCommon
  use interpot, only : springpots
  implicit none
  integer :: i

  if (.not.springpots) then
    write(6,'(A)') "# WARNING: InitPos not applicable without a spring potential"
    Return
  endif

  if (.not.allocated(ipos)) then
    allocate( ipos(1:3,1:np) )
    write(6,'(A)') '#  Setting current configuration positions as initial positions for'
    write(6,'(A)') '#   spring potential. Use directive InitPos to handle this manually.'
  endif
  initBox(1) = box(1)
  initBox(2) = box(2)
  initBox(3) = box(3)

!$omp parallel do &
!$omp shared (pos, ipos) &
!$omp private (i)
  do i = 1, np
    ipos(1,i) = pos(1,i)
    ipos(2,i) = pos(2,i)
    ipos(3,i) = pos(3,i)
  enddo
!$omp end parallel do

end subroutine !}}}
