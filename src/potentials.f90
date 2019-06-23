! vim:fdm=marker
subroutine buildPairTables
!{{{
! create table of pair potentials and forces
  use interpot
  use mostCommon, only : box
  implicit none
  integer :: i, j
  real(8) :: s, r, alpha, vcon, fcon

  alpha = 0.2d0 ! for the damped shifted coulomb interaction

  if (allocated(PPairTable)) then !reallocate and rebuild tables
    deallocate( PPairTable, FPairTable )
  endif
  allocate( PPairTable(PTsize+1,NPairs), FPairTable(PTsize+1,NPairs) )

  do j = 1, NPairs
    rcut = PairParams(1,j)
    PTrc2(j) = rcut*rcut
    PTdr2(j) = (PTrc2(j)-PTrMin)/(PTsize-1)
    write(6,'(a,i3,2f12.6)') "# pair,rcut2,dr2",j,PTrc2(j),PTdr2(j)
    ! set initial parameters for the damped shited coulomb potential
    vcon=  erfc( alpha*   rcut)   / rcut
    fcon=( erfc( alpha*   rcut)   /PTrc2(j) + 2.*alpha/dsqrt(PI) &
           *exp(-alpha*alpha*PTrc2(j))/ rcut )
    do i = 1, PTsize+1
      s = PTrMin + dble(i-1)*PTdr2(j)
      r = dsqrt(s) 
      ! Force is multiplied by 'r' so that we can do Rij(:)/dist^2
      PPairTable(i,j) =    pot(PairName(j), r, PairParams(1:4,j), PairChargeProd(j) )
      FPairTable(i,j) = r*dpot(PairName(j), r, PairParams(1:4,j), PairChargeProd(j) )
    enddo
  enddo
  
  rcut = maxval(PairParams(1,1:NPairs)) !set global cutoff as maximum
  rcut2 = rcut*rcut
  srcut(1:3) = rcut/box(1:3) !normalized cutoff radius

  write(6,'(A)') "# Built potential lists"
  write(6,'(A,g11.4)') "# Global maximum cutoff radius: ",rcut
  call flush(6)

  if (potOuts) then
    open(37,file="pottab.dat")
    open(38,file="forcetab.dat")
    do i = 1, PTsize
      ! fraction of cutoff radius, pair values
      write(37,*) sqrt(real(i)/real(PTsize)), PPairTable(i,1:NPairs)
      write(38,*) sqrt(real(i)/real(PTsize)), FPairTable(i,1:NPairs)
    enddo
    close(37)
    close(38)
  endif

!}}}
contains 
!!!!!!!!!!!!!!
function pot(typ, r, parm, qprd )
!{{{
!  This function returns the potential energy for a given inter-atomic distance
!  and a given potential form.
!  Available Potentials are:
!     Morse, Morse with quadratic repulsion, and/or coulomb
!     Buckingham, Buckingham with coulomb
!     Lennard-Jones, Lennard-Jones with coulomb
!     Johnson
!
implicit none
real(8) :: pot
character(8) :: typ
real(8), parameter :: ke=1.d0 !14.399637 !Ang*eV/e^2
real(8) :: r
real(8) :: alfa, tmp, D, r0, rc, alpha0  !morse
real(8) :: qprd         !coulomb charge product
real(8) :: A, B, rho, C    !buckingham
real(8) :: eps, sig     ! LJ
real(8), dimension(4) :: parm
pot=0.0

select case (trim(typ))
  case("sine")
    pot= dsin(r)*dsin(r)-1.d0
  case("cosine")
    pot= dcos(r)*dcos(r)-1.d0
  case("soft")
    pot= parm(2)*(1.d0+dcos(PI*r/parm(3)))
  case("spring")
    rc = parm(1) !cutoff radius
    D  = parm(2) !well depth
    C  = parm(3) !max strain
    r0 = r/rc*C  !rescale the input r which ranges to rc, to range to C
    pot= D*( r0*r0 - C*C )
  case("morse","Morse","morsec","morseC","MorseC")
    !!!!!! Morse !!!!!!!!!!!!!!!
        rc    =parm(1)
        D     =parm(2)  !eV
        r0    =parm(4)  !Ang
        alpha0=parm(3)  !Ang^-1
        tmp   =(r-r0) !/r0
        alfa  =alpha0 !*r0
     if (r0.le.0.0) tmp=0.0 !catch NaN

                  pot= exp(-Alfa*tmp) 
                  pot = D*(pot*pot - 2.d0*pot)
                  !pot= D*(exp(-2.*Alfa*tmp) -2.*exp(-Alfa*tmp))
   ! addition of coulomb potential
   if (typ.eq."MorseC".or.typ.eq."morseC") then
     pot= pot + ke*qprd*(erfc(alpha*r)/r -vcon +fcon*(r-rc))
   endif

  case("Buck","BuckC","buck","buckC","buckc")
    !!!!!! Buckingham !!!!!!!!!!!!!!!
        rc    =parm(1)
        A     =parm(2)  !eV
        rho   =parm(3)  !Ang
        C     =parm(4)  !eV*Ang^6

                  pot= A*exp(-r/rho) -C/(r**6)
     if (rho.le.0.0) pot=A -C/(r**6) !catch NaN
   ! addition of coulomb potential
   if (typ.eq."BuckC") pot= pot + ke*qprd*(erfc(alpha*r)/r -vcon +fcon*(r-rc))

  case("LJ","LJC")
    !!!!!! Lennard-Jones !!!!!!!!!!!!!!!tt
        rc    =parm(1)
        eps   =parm(2)  !eV
        sig   =parm(3)  !Ang
        tmp   =sig/r

                  pot= 4.*eps*(tmp**12-tmp**6)
   ! addition of coulomb potential
   if (trim(typ).eq."LJC") pot= pot + ke*qprd*(erfc(alpha*r)/r -vcon +fcon*(r-rc))

  case("John")
    !!!!!! Johnson !!!!!!!!!!!!!!!
        eps   =parm(2)  !eV
        r0    =parm(3)  !Ang
        rc    =parm(4)  !Ang
        A     =2.*r0**3/((rc-r0)**3)
        B     =3.*r0**2/((rc-r0)**2)
        tmp   =(r-rc)/r0
        if (r.lt.rc) then
                  pot= -eps*(A*tmp**3 + B*tmp**2)
        else
                  pot=0.0
        endif
  case default
    write(0,'(A)') "ERROR: No potential available for Potential Type: "//typ
    call flush(0); STOP
end select
end function pot !}}}
!!!!!!!!!!!!!!
function dpot(typ, r, parm, qprd)
!{{{
!  This function returns the derivative of the potential energy for a given
!  inter-atomic distance and a given potential form.
!
implicit none
real(8) :: dpot
character(8) :: typ
real(8), parameter :: ke=1.d0 !14.399637d0 !Ang*eV/e^2
real(8) :: r
real(8) :: alfa, tmp, D, r0, rc, alpha0  !morse
real(8) :: qprd         !coulomb charge product
real(8) :: A, B, rho, C    !buckingham
real(8) :: eps, sig     ! LJ
real(8), dimension(4) :: parm
dpot=0.0

select case (trim(typ))
  case("sine")
    dpot= dsin(2.d0*r)
  case("cosine")
    dpot= -dsin(2.d0*r)
  case("soft")
    dpot= PI*parm(1)*dsin(PI*r/parm(2))/parm(2)
  case("spring") !linear force response
    rc = parm(1)
    D  = parm(2) !well depth
    C  = parm(3) !max strain
    r0 = r/rc*C  !rescale the input r which ranges to rc, to range to C
    dpot= -2.d0*D*r0 
  case("morse","Morse","morseC","morsec","MorseC")
    !!!!!! Morse derivative !!!!!!!!!!!!!!!
        D     =parm(2)  !eV
        r0    =parm(4)  !Ang
        alpha0=parm(3)  !Ang^-1
        tmp   =(r-r0) !/r0
        alfa  =alpha0 !*r0
     if (r0.le.0.0) then
         tmp=0.0
         r0=1.0 !catch NaN
     endif

                  !dpot= 2.*Alfa/r0*D*(exp(-2.*Alfa*tmp) -exp(-Alfa*tmp))
                  dpot = exp(-Alfa*tmp)
                  dpot = 2.*D*Alfa*(dpot*dpot - dpot)/r
                  !dpot= 2.*Alfa*D*(exp(-2.*Alfa*tmp) -exp(-Alfa*tmp))
   ! addition of coulomb potential
     if (typ.eq."MorseC".or.typ.eq."morseC") then
      dpot= dpot + ke*qprd*(erfc(alpha*r)/(r**2) +2.*alpha/ &
                 &      dsqrt(PI)*exp(-alpha**2*r**2)/r -fcon)
     endif

  case("Buck","BuckC","buck","buckC","buckc")
    !!!!!! Buckingham !!!!!!!!!!!!!!!
        A     =parm(2)  !eV
        rho   =parm(3)  !Ang
        C     =parm(4)  !eV*Ang^6

                  dpot= A*exp(-r/rho)/rho -6.*C/(r**7)
     if (rho.le.0.0) dpot=A -6.*C/(r**7) !catch NaN
   ! addition of coulomb potential
     if (typ.eq."BuckC") dpot= dpot + ke*qprd*(erfc(alpha*r)/(r**2) +2.*alpha/ &
             &   dsqrt(PI)*exp(-alpha**2*r**2)/r -fcon)

  case("LJ","LJC")
    !!!!!! Lennard-Jones !!!!!!!!!!!!!!!
        eps   =parm(2)  !eV
        sig   =parm(3)  !Ang
        tmp   =sig/r

                  dpot= 24.*eps*(2.*tmp**12/r-tmp**6/r)
   ! addition of coulomb potential
     if (trim(typ).eq."LJC") dpot= dpot + ke*qprd*(erfc(alpha*r)/(r**2) +2.*alpha/ &
            & dsqrt(PI)*exp(-alpha**2*r**2)/r -fcon)

  case("John")
    !!!!!! Johnson !!!!!!!!!!!!!!!
        eps   =parm(2)  !eV
        r0    =parm(3)  !Ang
        rc    =parm(4)  !Ang
        A     =2.*r0**3/((rc-r0)**3)
        B     =3.*r0**2/((rc-r0)**2)
        tmp   =(r-rc)/r0
        if (r.lt.rc) then
                  dpot= eps/r0*(3.*A*tmp**2 + 2.*B*tmp)
        else
                  dpot=0.0
        endif
  
  case default
    write(0,'(A)') "ERROR: no differentiated potential available for potential type: "//typ
    call flush(0); STOP
end select
end function dpot !}}}

end subroutine buildPairTables

