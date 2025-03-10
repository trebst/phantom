!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module cooling
!
! Gas cooling
!  Current options:
!     0 = none
!     1 = 'explicit' cooling [implicitly calculated, ironically] [default]
!         (not sure what this actually means, but requires a lot of non-default inputs)
!     2 = Townsend (2009) cooling tables [implicitly calculated]
!     3 = Gammie cooling [explicitly calculated]
!     5 = Koyama & Inutuska (2002) [explicitly calculated]
!     6 = Koyama & Inutuska (2002) [implicitly calculated]
!
! :References:
!   Koyama & Inutsuka (2002), ApJL 564, 97-100
!   Vazquez-Semadeni, et.al (2007), ApJ 657, 870-883
!   Townsend (2009), ApJS 181, 391-397
!   Gail & Sedlmayr textbook Physics and chemistry of Circumstellar dust shells
!
! :Owner: Lionel Siess
!
! :Runtime parameters:
!   - C_cool         : *factor controlling cooling timestep*
!   - Tfloor         : *temperature floor (K); on if > 0*
!   - beta_cool      : *beta factor in Gammie (2001) cooling*
!   - bowen_Cprime   : *radiative cooling rate (g.s/cm³)*
!   - cooltable      : *data file containing cooling function*
!   - dust_collision : *dust collision [1=on/0=off]*
!   - excitation_HI  : *cooling via electron excitation of HI [1=on/0=off]*
!   - habund         : *Hydrogen abundance assumed in cooling function*
!   - icooling       : *cooling function (0=off, 1=explicit, 2=Townsend table, 3=Gammie, 5,6=KI02)*
!   - relax_bowen    : *Bowen (diffusive) relaxation [1=on/0=off]*
!   - relax_stefan   : *radiative relaxation [1=on/0=off]*
!   - temp_floor     : *Minimum allowed temperature in K for Townsend cooling table*
!
! :Dependencies: chem, cooling_molecular, datafiles, dim, eos, h2cooling,
!   infile_utils, io, options, part, physcon, timestep, units
!

 use options,  only:icooling
 use timestep, only:C_cool
 use physcon

 implicit none
 character(len=*), parameter :: label = 'cooling'

 public :: init_cooling,calc_cooling_rate,energ_cooling
 public :: write_options_cooling, read_options_cooling
 public :: find_in_table, implicit_cooling, exact_cooling
 logical, public :: cooling_implicit
 logical, public :: cooling_explicit
 real,    public :: bowen_Cprime = 3.000d-5
 real,    public :: GammaKI_cgs = 2.d-26 ! [erg/s] heating rate for Koyama & Inutuska cooling

 private
 integer, parameter :: nTg = 64
 integer, parameter :: maxt = 1000
 real,    parameter :: Tref = 1.d5, T_floor = 10.   ! required for exact_cooling
 integer :: nt
 real    :: temper(maxt),lambda(maxt),slope(maxt),yfunc(maxt),rhov4_KI02(2,maxt)
 real    :: beta_cool  = 3.
 real    :: habund     = 0.7
 real    :: temp_floor = 1.e4                       ! required for exact_cooling_table
 real    :: Tgrid(nTg)
 real    :: LambdaKI_coef,GammaKI
 real    :: KI02_rho_min_cgs = 1.0d-30  ! minimum density of the KI02 cooling curve
 real    :: KI02_rho_max_cgs = 1.0d-14  ! maximum density of the KI02 cooling curve
 real    :: KI02_rho_min,KI02_rho_max
 integer :: excitation_HI = 0, relax_Bowen = 0, dust_collision = 0, relax_Stefan = 0
 integer :: icool_radiation_H0 = 0, icool_relax_Bowen = 0, icool_dust_collision = 0, icool_relax_Stefan = 0
 character(len=120) :: cooltable = 'cooltable.dat'
 !--Minimum temperature (failsafe to prevent u < 0); optional for ALL cooling options
 real,    public :: Tfloor = 0. ! [K]; set in .in file.  On if Tfloor > 0.
 real,    public :: ufloor = 0. ! [code units]; set in init_cooling

contains

!-----------------------------------------------------------------------
!+
!  Initialise cooling
!+
!-----------------------------------------------------------------------
subroutine init_cooling(id,master,iprint,ierr)
 use dim,       only:maxvxyzu
 use units,     only:utime,umass,udist,unit_ergg
 use physcon,   only:mass_proton_cgs,kboltz
 use io,        only:fatal
 use eos,       only:gamma,gmw
 use part,      only:h2chemistry
 use h2cooling, only:init_h2cooling
 use chem,      only:init_chem
 use cooling_molecular, only:do_molecular_cooling,init_cooling_molec
 integer, intent(in)  :: id,master,iprint
 integer, intent(out) :: ierr

 if (h2chemistry) then
    if (id==master) write(iprint,*) 'initialising cooling function...'
    call init_chem()
    call init_h2cooling()
 else
    !you can't have cool_relaxation_Stefan and cool_relaxation_Bowen at the same time
    if (relax_Bowen == 1 .and. relax_Stefan == 1) then
       call fatal(label,'you can"t have bowen and stefan cooling at the same time')
    endif

#ifdef KROME
    !krome calculates its own cooling rate
    excitation_HI = 0
    dust_collision = 0
#else
    !if no cooling flag activated, disable cooling
    if (icooling == 1 .and. (excitation_HI+relax_Bowen+dust_collision+relax_Stefan) == 0 &
       .and. do_molecular_cooling) then
       icooling = 0
       return
    elseif (icooling == 1 .and. do_molecular_cooling) then
       ! Initialise cooling tables
       call init_cooling_molec
    endif
#endif

    !--initialise remaining variables
    if (icooling == 2) then
       call init_cooltable(ierr)
    elseif (icooling == 5 .or. icooling == 6) then
       LambdaKI_coef = GammaKI_cgs*umass*utime**3/(mass_proton_cgs**2 * udist**5)
       GammaKI       = GammaKI_cgs*utime**3/(mass_proton_cgs*udist**2)
       call init_hv4table(ierr)
       if (ierr > 0) call fatal('init_cooling','Failed to create KI02 cooling table')
    elseif (icooling > 0) then
       call set_Tgrid
    endif
 endif

 !--Determine if this is implicit or explicit cooling
 cooling_implicit = .false.
 cooling_explicit = .false.
 if (h2chemistry) then
    if (icooling > 0) cooling_implicit = .true.    ! cooling is calculated implicitly in step
 elseif (icooling > 0) then
    if (icooling == 3 .or. icooling == 5) then
       cooling_explicit = .true.                   ! cooling is calculated explicitly in force
    else
       cooling_implicit = .true.                   ! cooling is calculated implicitly in step
    endif
 endif

 !--calculate the energy floor in code units
 if (Tfloor > 0.) then
    if (gamma > 1.) then
       ufloor = kboltz*Tfloor/((gamma-1.)*gmw*mass_proton_cgs)/unit_ergg
    else
       ufloor = 3.0*kboltz*Tfloor/(2.0*gmw*mass_proton_cgs)/unit_ergg
    endif
    if (maxvxyzu < 4) ierr = 1
 else
    ufloor = 0.
 endif

end subroutine init_cooling

!-----------------------------------------------------------------------
!+
!  read cooling table from file and initialise arrays
!+
!-----------------------------------------------------------------------
subroutine init_cooltable(ierr)
 use io,        only:fatal
 use datafiles, only:find_phantom_datafile
 integer, intent(out) :: ierr
 integer, parameter :: iu = 127
 integer :: i
 character(len=120) :: filepath

 !
 ! read the cooling table from file
 !
 filepath=find_phantom_datafile(cooltable,'cooling')
 open(unit=iu,file=filepath,status='old',iostat=ierr)
 if (ierr /= 0) call fatal('cooling','error opening cooling table')
 i = 0
 do while(ierr==0 .and. i < maxt)
    i = i + 1
    read(iu,*,iostat=ierr) temper(i),lambda(i)
 enddo
 nt = i-1
 if (nt==maxt) call fatal('cooling','size of cooling table exceeds array size')
 if (nt < 2) call fatal('cooling','size of cooling table is too small',ival=nt,var='nt')
 !
 ! calculate the slope of the cooling function
 !
 do i=1,nt-1
    slope(i) = log(lambda(i+1)/lambda(i))/log(temper(i+1)/temper(i))
 enddo
 slope(nt) = slope(nt-1)

 !
 ! initialise the functions required for Townsend method
 !
 yfunc(nt) = 0.
 do i=nt-1,1,-1
    !Lionel Siess : I think there is an error. in yfunc slope(nt) should be replaced by lambda(nt)
    ! Eq A6
    if (abs(slope(i)-1.) < tiny(0.)) then
       !ori yfunc(i) = yfunc(i+1) - slope(nt)*temper(i)/(lambda(i)*temper(nt))*log(temper(i)/temper(i+1))
       yfunc(i) = yfunc(i+1) - lambda(nt)*temper(i)/(lambda(i)*temper(nt))*log(temper(i)/temper(i+1))
    else
       !ori yfunc(i) = yfunc(i+1) - slope(nt)*temper(i)/((1. - slope(i))*lambda(i)*temper(nt))&
       yfunc(i) = yfunc(i+1) - lambda(nt)*temper(i)/((1. - slope(i))*lambda(i)*temper(nt))&
                 *(1.- (temper(i)/temper(i+1))**(slope(i) - 1.))
    endif
 enddo
end subroutine init_cooltable

!-----------------------------------------------------------------------
!+
!  create a h-v4 table based upon the cooling curve of KI02
!+
!-----------------------------------------------------------------------
subroutine init_hv4table(ierr)
 use part,    only:hrho,igas
 use physcon, only:mass_proton_cgs,kboltz
 use units,   only:unit_density,unit_velocity
 use eos,     only:gmw,gamma
 integer, intent(out) :: ierr
 integer              :: i,ctr
 real                 :: nrho0_min,nrho0_max,nrho,dnrho,dGammaKI,Lambda,dLambda
 real                 :: T,Tnew,Trat,fatT,faTdT
 logical              :: iterate
 logical              :: print_cc = .false. ! Print the cooling curve (for testing)

 !--Initialise densities
 KI02_rho_min = KI02_rho_min_cgs/unit_density
 KI02_rho_max = KI02_rho_max_cgs/unit_density
 nrho0_min    = KI02_rho_min_cgs/mass_proton_cgs
 nrho0_max    = KI02_rho_max_cgs/mass_proton_cgs
 dnrho        = (log10(nrho0_max) - log10(nrho0_min))/maxt
 !--Initialise additional variables
 dGammaKI     = 0.0
 ierr         = 0

 if (print_cc) open(unit=1031,file='coolingcurve.dat')

 !--Iterate (in cgs units)!
 T = 20000.
 do i = 1,maxt
    ctr     = 0
    iterate = .true.
    nrho    = 10**(log10(nrho0_min) + (i-1)*dnrho)
    do while ( iterate )
       Lambda  = 1.d7*exp(-1.184d5/(T+1.d3)) + 0.014*sqrt(T)*exp(-92./T) ! This is actually Lamda / Gamma
       dLambda = 0.007*exp(-92./T)*(T+184.)*T**(-1.5) + 1.184d12*exp(-1.184d5/(T+1.d3))*(T+1.d3)**(-2)
       fatT    =  Lambda*GammaKI_cgs*nrho -  GammaKI_cgs
       faTdT   = dLambda*GammaKI_cgs*nrho - dGammaKI
       Tnew    = abs(T - fatT/faTdT)
       Trat    = abs( 1.0 - T/Tnew )
       T       = Tnew
       ctr     = ctr + 1
       !--converged
       if (Trat < 1.0d-6) iterate = .false.
       !--failed to converge
       if (T < 0. .or. ctr > 2000) then
          iterate = .false.
          ierr    = 1
       endif
    enddo
    if (print_cc) write(1031,*) nrho,nrho*mass_proton_cgs,T,T*nrho,Lambda*GammaKI_cgs
    rhov4_KI02(1,i) = nrho
    rhov4_KI02(2,i) = T
 enddo
 if (print_cc) close(1031)

 !--Convert to useful values
 do i = 1,maxt
    rhov4_KI02(1,i) = rhov4_KI02(1,i)*mass_proton_cgs/unit_density                               ! number density (cm^-3) -> mass density (code units)
    rhov4_KI02(2,i) = kboltz*rhov4_KI02(2,i)/(gmw*mass_proton_cgs*(gamma-1.0))/unit_velocity**2  ! T -> internal energy (code units)
 enddo

end subroutine init_hv4table
!-----------------------------------------------------------------------
!
!  calculate cooling rates
!
!-----------------------------------------------------------------------
subroutine calc_cooling_rate(r, Q, dlnQ_dlnT, rho, T, Teq, mu, K2, kappa)
 use units,             only:unit_ergg,unit_density
 use cooling_molecular, only:do_molecular_cooling,calc_cool_molecular

 real, intent(in) :: rho, T, Teq               !rho in code units
 real, intent(in) :: r
 real, intent(in), optional :: mu, K2, kappa   !cgs
 real, intent(out) :: Q, dlnQ_dlnT             !code units
 real :: Q_cgs,Q_H0, Q_relax_Bowen, Q_col_dust, Q_relax_Stefan, Q_molec, rho_cgs
 real :: dlnQ_H0, dlnQ_relax_Bowen, dlnQ_col_dust, dlnQ_relax_Stefan, dlnQ_molec

 rho_cgs = rho*unit_density
 Q_H0              = 0.
 Q_relax_Bowen     = 0.
 Q_col_dust        = 0.
 Q_relax_Stefan    = 0.
 Q_molec           = 0.

 dlnQ_H0           = 0.
 dlnQ_relax_Bowen  = 0.
 dlnQ_col_dust     = 0.
 dlnQ_relax_Stefan = 0.
 dlnQ_molec        = 0.

 if (excitation_HI  == 1) call cooling_neutral_hydrogen(T, rho_cgs, Q_H0, dlnQ_H0)
 if (relax_Bowen    == 1) call cooling_Bowen_relaxation(T, Teq, rho_cgs, mu, Q_relax_Bowen, dlnQ_relax_Bowen)
 if (dust_collision == 1) call cooling_dust_collision(T, Teq, rho_cgs, K2, mu, Q_col_dust, dlnQ_col_dust)
 if (relax_Stefan   == 1) call cooling_radiative_relaxation(T, Teq, kappa, Q_relax_Stefan, dlnQ_relax_Stefan)
 if (do_molecular_cooling) call calc_cool_molecular(T, r, rho_cgs, Q_molec, dlnQ_molec)

 Q_cgs = Q_H0 + Q_relax_Bowen + Q_col_dust + Q_relax_Stefan + Q_molec
 if (Q_cgs == 0.) then
    dlnQ_dlnT = 0.
 else
    dlnQ_dlnT = (Q_H0*dlnQ_H0 + Q_relax_Bowen*dlnQ_relax_Bowen + Q_col_dust*dlnQ_col_dust&
   + Q_relax_Stefan*dlnQ_relax_Stefan + Q_molec*dlnQ_molec)/Q_cgs
 endif
 !limit exponent to prevent overflow
 dlnQ_dlnT = sign(min(50.,abs(dlnQ_dlnT)),dlnQ_dlnT)
 Q = Q_cgs/unit_ergg
end subroutine calc_cooling_rate


!-----------------------------------------------------------------------
!+
!  Bowen 1988 cooling prescription
!+
!-----------------------------------------------------------------------
subroutine cooling_Bowen_relaxation(T, Teq, rho, mu, Q, dlnQ_dlnT)
! all quantities in cgs
 use eos,     only:gamma
 use physcon, only:Rg
 real, intent(in) :: T, Teq, rho, mu
 real, intent(out) :: Q,dlnQ_dlnT

 Q         = Rg/((gamma-1.)*mu)*rho*(Teq-T)/bowen_Cprime
 dlnQ_dlnT = -T/(Teq-T+1.d-10)

end subroutine cooling_Bowen_relaxation

!-----------------------------------------------------------------------
!+
!  collisionnal cooling
!+
!-----------------------------------------------------------------------
subroutine cooling_dust_collision(T, Teq, rho, K2, mu, Q, dlnQ_dlnT)
! all quantities in cgs
 use physcon, only: kboltz, mass_proton_cgs, pi
 real, intent(in) :: T, Teq, rho, K2, mu
 real, intent(out) :: Q,dlnQ_dlnT

 real, parameter :: f = 0.15, a0 = 1.28e-8
 real :: A

 A = 2. * f * kboltz * a0**2/(mass_proton_cgs**2*mu) &
         * (1.05/1.54) * sqrt(2.*pi*kboltz/mass_proton_cgs) * 2.*K2 * rho
 Q = A * sqrt(T) * (Teq-T)
 if (Q  >  1.d6) then
    print *, f, kboltz, a0, mass_proton_cgs, mu
    print *, mu, K2, rho, T, Teq, A, Q
    stop 'cooling'
 else
    dlnQ_dlnT = 0.5+T/(Teq-T+1.d-10)
 endif
end subroutine cooling_dust_collision

!-----------------------------------------------------------------------
!+
!  Woitke (2006 A&A) cooling term
!+
!-----------------------------------------------------------------------
subroutine cooling_radiative_relaxation(T, Teq, kappa, Q, dlnQ_dlnT)
 use physcon, only: steboltz
 real, intent(in) :: T, Teq, kappa
 real, intent(out) :: Q,dlnQ_dlnT

 Q         = 4.*steboltz*(Teq**4-T**4)*kappa
 dlnQ_dlnT = -4.*T**4/(Teq**4-T**4+1.d-10)

end subroutine cooling_radiative_relaxation

!-----------------------------------------------------------------------
!+
!  Cooling due to electron excitation of neutral H (Spitzer 1978)
!+
!-----------------------------------------------------------------------
subroutine cooling_neutral_hydrogen(T, rho, Q, dlnQ_dlnT)
 use physcon, only: mass_proton_cgs, pi
 real, intent(in) :: T, rho
 real, intent(out) :: Q,dlnQ_dlnT

 real, parameter :: f = 1.0d0
 real :: eps_e

 if (T > 3000.) then
    eps_e = calc_eps_e(T)
    Q = -f*7.3d-19*eps_e*exp(-118400./T)*rho/(1.4*mass_proton_cgs)**2
    dlnQ_dlnT = 118400.d0/T+log(calc_eps_e(1.001*T)/eps_e)/log(1.001)
 else
    Q = 0.
    dlnQ_dlnT = 0.
 endif
end subroutine cooling_neutral_hydrogen

!-----------------------------------------------------------------------
!+
!  compute electron equilibrium abundance (Palla et al 1983)
!+
!-----------------------------------------------------------------------
real function calc_eps_e(T)
 real, intent(in) :: T
 real :: k1, k2, k3, k8, k9, p, q

 k1 = 1.88d-10 / T**6.44e-1
 k2 = 1.83d-18 * T
 k3 = 1.35d-9
 k8 = 5.80d-11 * sqrt(T) * exp(-1.58d5/T)
 k9 = 1.7d-4 * k8
 p = .5*k8/k9
 q = k1*(k2+k3)/(k3*k9)
 calc_eps_e = (p + sqrt(q+p**2))/q
end function calc_eps_e

!-----------------------------------------------------------------------
!+
!  Set Temperature grid
!+
!-----------------------------------------------------------------------

subroutine set_Tgrid
 integer :: i
 real :: dlnT
 dlnT = log(Tref)/(nTg-1)

 do i = 1,nTg
    Tgrid(i) = exp((i-1)*dlnT)
 enddo
end subroutine set_Tgrid

!-----------------------------------------------------------------------
!+
!   Gammie (2001) cooling
!+
!-----------------------------------------------------------------------
subroutine cooling_Gammie(xi,yi,zi,ui,dudti)
 real, intent(in)    :: ui,xi,yi,zi
 real, intent(inout) :: dudti

 real :: omegai,r2,tcool1

 r2     = xi*xi + yi*yi + zi*zi
 Omegai = r2**(-0.75)
 tcool1 = Omegai/beta_cool
 dudti  = dudti - ui*tcool1

end subroutine cooling_Gammie

!-----------------------------------------------------------------------
!+
!   Cooling rate as per Koyama & Inutuska (2002; eqns 4 & 5);
!   typos corrected as per Vazquez-Semadeni+ (2007)
!   This is for the explicit calculation
!   In equilibrium, n*LambdaKI = (rho/mp)*LambdaKI = GammaKI
!+
!-----------------------------------------------------------------------
subroutine cooling_KoyamaInutuska_explicit(rhoi,Tgas,dudti)
 real, intent(in)    :: rhoi,Tgas
 real, intent(inout) :: dudti
 real                :: LambdaKI

 ! Derivation to obtain correct units; used Koyama & Inutuska (2002) as the reference
 !LambdaKI = GammaKI_cgs * (1.d7*exp(-118400./(Tgas+1000))+0.014*sqrt(Tgas)*exp(-92./Tgas)) ! The cooling rate in erg cm^3/s = g cm^5/s^3
 !LambdaKI = LambdaKI/mass_proton_cgs**2                                                    ! units are now cm^5/(g s^3) ! since [u] = erg/g = cm^2/s^2
 !LambdaKI = LambdaKI*umass*utime**3/udist**5                                               ! convert to from cm^5/(g s^3) to code units
 !dudti    = dudti - LambdaKI*rhoi*fac                                                      ! multiply by rho (code) to get l^5/(m t^3) * m/l^3 = l^2/s^3 = [u]
 !
 !GammaKI = GammaKI_cgs                                                                     ! The heating rate in erg /s = g cm^2/s^3
 !GammaKI = GammaKI/mass_proton_cgs                                                         ! divide by proton mass.  Units are now g cm^2 / s^3 / g = cm^2/s^3
 !GammaKI = GammaKI*utime**3/udist**2                                                       ! convert from cm^2/s^3 to code units
 !dudti   = dudti + GammaKI                                                                 ! units and dependencies are correct

 LambdaKI = LambdaKI_coef*(1.d7*exp(-118400./(Tgas+1000.))+0.014*sqrt(Tgas)*exp(-92./Tgas))
 dudti    = dudti - LambdaKI*rhoi + GammaKI

end subroutine cooling_KoyamaInutuska_explicit

!-----------------------------------------------------------------------
!+
!   Cooling rate as per Koyama & Inutuska (2002; eqns 4 & 5);
!   typos corrected as per Vazquez-Semadeni+ (2007)
!   This is the implicit method given by (5)-(6) in Vazquez-Semadeni+ (2007)
!+
!-----------------------------------------------------------------------
subroutine cooling_KoyamaInutuska_implicit(eni,rhoi,dt,dudti)
 use eos, only: gamma,temperature_coef,gmw
 real, intent(in)    :: rhoi,eni,dt
 real, intent(out)   :: dudti
 integer             :: i,j,jm1
 real                :: ponrhoi,tempi,eni_equil,eni_final,deni,tau1,LambdaKI

 !--Determine the indicies surrounding the input h
 i = minloc(abs(rhov4_KI02(1,1:maxt)-rhoi), 1)
 if (i==1) then
    !print*, 'min density too large! extrapolating using two smallest densities'
    j = 2
 elseif (i==maxt) then
    !print*, 'max density too small! extrapolating using two largest densities'
    j = maxt
 elseif (rhov4_KI02(1,i-1) <= rhoi .and. rhoi <= rhov4_KI02(1,i  )) then
    j = i
 elseif (rhov4_KI02(1,i  ) <= rhoi .and. rhoi <= rhov4_KI02(1,i+1)) then
    j = i+1
 else
    print*, rhoi,rhov4_KI02(1,i-1:i+1)
    print*, 'this should not happen'
    stop
 endif

 !--Calculate the equilibrium energy by linear interpolation
 jm1       = j - 1
 eni_equil = rhov4_KI02(2,j) + (rhov4_KI02(2,jm1)-rhov4_KI02(2,j))/(rhov4_KI02(1,jm1)-rhov4_KI02(1,j))*(rhoi-rhov4_KI02(1,j))

 !--Determine the inverse time require to radiate/acquire excess/deficit energy & Update energy
 ponrhoi  = (gamma-1.)*eni
 tempi    = temperature_coef*gmw*ponrhoi
 LambdaKI = LambdaKI_coef*(1.d7*exp(-118400./(tempi+1000.))+0.014*sqrt(tempi)*exp(-92./tempi))
 dudti    = LambdaKI*rhoi - GammaKI
 deni     = eni - eni_equil

 if (abs(deni) > 0.) then
    ! in both limits, this will approach the correct value
    tau1      = abs(dudti/deni)
    eni_final = eni_equil + deni*exp(-dt*tau1)
    dudti     = -(eni - eni_final)/dt
 else
    ! in the unlikly chance deni = 0
    dudti = -dudti
 endif

end subroutine cooling_KoyamaInutuska_implicit

!-----------------------------------------------------------------------
!
!   explicit cooling
!
!-----------------------------------------------------------------------
subroutine explicit_cooling (xi,yi,zi,ui, dudt, rho, dt, Trad, mu_in, K2, kappa)
 use eos,     only:gamma,gmw
 use physcon, only:Rg
 use units,   only:unit_ergg
 real, intent(in) :: xi, yi, zi, ui, rho, dt, Trad !code units
 real, intent(in), optional :: mu_in, K2, kappa
 real, intent(out) :: dudt !code units
 real :: u,Q,dlnQ_dlnT,T,mu,T_on_u
 real :: r                         !in au

 r  = sqrt(xi*xi + yi*yi + zi*zi)

 if (.not.present(mu_in)) then
    mu = gmw
 else
    mu = mu_in
 endif
 T_on_u = (gamma-1.)*mu*unit_ergg/Rg
 T = T_on_u*ui
 call calc_cooling_rate(r, Q, dlnQ_dlnT, rho, T, Trad, mu, K2, kappa)
 if (-Q*dt  > ui) then   ! assume thermal equilibrium
    u = Trad/T_on_u
    dudt = (u-ui)/dt
 else
    dudt = Q
 endif
 !print *,T,Teq,T_on_u*u,'dT=',T_on_u*Q*dt,u,Q*dt

end subroutine explicit_cooling

!-----------------------------------------------------------------------
!
!   implicit cooling
!
!-----------------------------------------------------------------------
subroutine implicit_cooling (r, ui, dudt, rho, dt, Trad, mu_in, K2, kappa)
 use eos,     only:gamma,gmw
 use physcon, only:Rg
 use units,   only:unit_ergg
 real, intent(in) :: ui, rho, dt
 real, intent(in), optional :: Trad, mu_in, K2, kappa
 real, intent(out) :: dudt

 real, parameter :: tol = 1.d-4 ! to be adjusted
 integer, parameter :: iter_max = 200
 real :: u,Q,dlnQ_dlnT,T,mu,T_on_u,delta_u,term1,term2,term3
 real :: r    ! in au
 integer :: iter

 if (.not.present(mu_in)) then
    mu = gmw
 else
    mu = mu_in
 endif
 u = ui
 T_on_u = (gamma-1.)*mu*unit_ergg/Rg
 delta_u = 1.d-3
 iter = 0
 !The pdv_work also depends on the internal energy and could also be included
 !in this loop provided this contribution was not accounted for in Force.F90
 ! see PP flag : IMPLICIT COOLING - pb: we need div(v) and it is only real*4
 !term2 = 1.-(gamma-1.)*dt*divcurlv !pdv=(gamma-1.)*vxyzu(4,i)*divcurlv(1,i)*dt
 term2 = 1.
 term1 = u !initial internal energy without cooling contributions
 do while (abs(delta_u) > tol .and. iter < iter_max)
    T = u*T_on_u
    call calc_cooling_rate(r,Q,dlnQ_dlnT, rho, T, Trad, mu, K2, kappa)
    term3 = u*term2-Q*dt
    delta_u = (term1-term3)/(term2-Q*dlnQ_dlnT*dt/u)
    u = u+delta_u
    iter = iter + 1
 enddo
 dudt =(u-term1)/dt
 if (u < 0. .or. isnan(u)) then
    print *,u
    stop ' u<0'
 endif

end subroutine implicit_cooling

!-----------------------------------------------------------------------
!
!   this routine returns the effective cooling rate du/dt
!
!-----------------------------------------------------------------------
subroutine energ_cooling(xi,yi,zi,ui,dudt,rho,dt,Trad,mu_in,K2,kappa,Tgas)
 use io, only: fatal
 real, intent(in)           :: xi,yi,zi,ui,rho,dt         ! in code units
 real, intent(in), optional :: Tgas,Trad,mu_in,K2,kappa   ! in cgs units
 real, intent(inout)        :: dudt                       ! in code units
 real                       :: r

 r  = sqrt(xi*xi + yi*yi + zi*zi)

 select case (icooling)
 case(1)
    call explicit_cooling(xi,yi,zi,ui, dudt, rho, dt, Trad, mu_in, K2, kappa)
 case (3)
    call cooling_Gammie(xi,yi,zi,ui,dudt)
 case (2)
    call exact_cooling_table(ui,rho,dt,dudt)
 case (5)
    if (present(Tgas)) then
       call cooling_KoyamaInutuska_explicit(rho,Tgas,dudt)
    else
       call fatal('energ_cooling','Koyama & Inutuska cooling requires gas temperature')
    endif
 case (6)
    call cooling_KoyamaInutuska_implicit(ui,rho,dt,dudt)
 case default
    !call exact_cooling(r, u, dudt, rho, dt, Trad, mu_in, K2, kappa)
    !call implicit_cooling(u, dudt, rho, dt, Trad, mu_in, K2, kappa)
    if (present(Trad) .and. present(mu_in) .and. present(K2) .and. present(kappa)) then
       call explicit_cooling(xi,yi,zi,ui, dudt, rho, dt, Trad, mu_in, K2, kappa)
    else
       call fatal('energ_cooling','default requires optional arguments; change icooling or ask D Price or L Siess to patch')
    endif
 end select

end subroutine energ_cooling

!-----------------------------------------------------------------------
!
!   cooling using Townsend (2009), ApJS 181, 391-397 method with
!   analytical cooling rate prescriptions
!
!-----------------------------------------------------------------------
subroutine exact_cooling (r, u, dudt, rho, dt, Trad, mu_in, K2, kappa)
 use eos,     only:gamma,gmw
 use physcon, only:Rg
 use units,   only:unit_ergg
 real, intent(in) :: u, rho, dt, Trad
 real, intent(in), optional :: mu_in, K2, kappa
 real, intent(out) :: dudt

 real, parameter :: tol = 1.d-12
 real :: Qref,dlnQref_dlnT,Q,dlnQ_dlnT,Y,Yk,Yinv,Temp,dy,T,mu,T_on_u
 real :: r  ! in au
 integer :: k

 if (.not.present(mu_in)) then
    mu = gmw
 else
    mu = mu_in
 endif
 T_on_u = (gamma-1.)*mu*unit_ergg/Rg
 T = T_on_u*u

 if (T < T_floor) then
    Temp = T_floor
 elseif (T > Tref) then
    call calc_cooling_rate(r,Q, dlnQ_dlnT, rho, T, Trad, mu, K2, kappa)
    Temp = T+T_on_u*Q*dt
 else
    call calc_cooling_rate(r,Qref,dlnQref_dlnT, rho, Tref, Trad, mu, K2, kappa)
    Y = 0.
    k = nTg
    Q = Qref                  ! default value if Tgrid < T for all k
    dlnQ_dlnT = dlnQref_dlnT  ! default value if Tgrid < T for all k
    do while (Tgrid(k) > T)
       k = k-1
       call calc_cooling_rate(r,Q, dlnQ_dlnT, rho, Tgrid(k), Trad, mu, K2, kappa)
       ! eqs A6
       if (abs(dlnQ_dlnT-1.) < tol) then
          y = y - Qref*Tgrid(k)/(Q*Tref)*log(Tgrid(k)/Tgrid(k+1))
       else
          y = y - Qref*Tgrid(k)/(Q*Tref*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/Tgrid(k+1))**(dlnQ_dlnT-1.))
       endif
    enddo
    !eqs A5
    yk = y
    if (abs(dlnQ_dlnT-1.) < tol) then
       y = yk + Qref*Tgrid(k)/(Q*Tref)*log(Tgrid(k)/T)
    else
       y = yk + Qref*Tgrid(k)/((Q*Tref)*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/T)**(dlnQ_dlnT-1))
    endif
    !eq 26
    dy = Qref*dt*T_on_u/Tref
    y = y + dy
    !compute Yinv (eqs A7)
    if (abs(dlnQ_dlnT-1.) < tol) then
       Temp = max(Tgrid(k)*exp(-Q*Tref*(y-yk)/(Qref*Tgrid(k))),T_floor)
    else
       Yinv = 1.-(1.-dlnQ_dlnT)*Q*Tref/(Qref*Tgrid(k))*(y-yk)
       if (Yinv > 0.) then
          Temp = Tgrid(k)*(Yinv**(1./(1.-dlnQ_dlnT)))
       else
          Temp = T_floor
       endif
    endif
 endif

 dudt = (Temp-T)/T_on_u/dt
 !note that u = Temp/T_on_u

end subroutine exact_cooling

!-----------------------------------------------------------------------
!+
!  cooling using Townsend (2009) method with tabulated rate
!
!   Implements cooling defined using a tabulated cooling table
!   produced e.g. by CLOUDY.
!+
!-----------------------------------------------------------------------
subroutine exact_cooling_table(uu,rho,dt,dudt)
 use eos,     only:gamma,gmw
 use physcon, only:atomic_mass_unit,kboltz,Rg
 use units,   only:unit_density,unit_ergg,utime
 real, intent(in)  :: uu, rho,dt
 real, intent(out) :: dudt
 real    :: gam1,density_cgs,dt_cgs,amue,amuh,dtemp
 real    :: sloperef,slopek,temp,temp1,tref,yfunx,yinv0
 integer :: k

 gam1 = gamma - 1.
 temp = gam1*uu/Rg*gmw*unit_ergg

 tref     = temper(nt)
 sloperef = slope(nt)

 if (temp < temp_floor) then
    temp1 = temp_floor
 else
    amue = 2.*atomic_mass_unit/(1. + habund)
    amuh = atomic_mass_unit/habund
    density_cgs = rho*unit_density
    dt_cgs      = dt*utime

    !Lionel Siess : I think there is an error. in dtemp sloperef should be replaced by lambda(nt)
    !original dtemp = gam1*density_cgs*(atomic_mass_unit*gmw/(amue*amuh*kboltz))* &
    !     sloperef/tref*dt_cgs
    ! Eq 26
    dtemp = gam1*density_cgs*(atomic_mass_unit*gmw/(amue*amuh*kboltz))*lambda(nt)/tref*dt_cgs

    k = find_in_table(nt,temper,temp)
    slopek = slope(k)
    ! Eq A5
    if (abs(slopek - 1.) < tiny(0.)) then
       yfunx = yfunc(k) + lambda(nt)*temper(k)/(lambda(k)*temper(nt))*log(temper(k)/temp)
    else
       yfunx = yfunc(k) + lambda(nt)*temper(k)/(lambda(k)*temper(nt)*(1. - slopek)) &
                          *(1. - (temper(k)/temp)**(slopek-1.))
    endif
    yfunx = yfunx + dtemp
    ! Eq A7
    if (abs(slopek - 1.) < tiny(0.)) then
       temp1 = max(temper(k)*exp(-lambda(k)*temper(nt)/(lambda(nt)*temper(k))*(yfunx-yfunc(k))),temp_floor)
    else
       yinv0 = 1. - (1. - slopek)*lambda(k)*temper(nt)/(lambda(nt)*temper(k))*(yfunx-yfunc(k))
       if (yinv0 > 0.) then
          temp1 = max(temper(k)*yinv0**(1./(1. - slopek)),temp_floor)
       else
          temp1 = temp_floor
       endif
    endif
 endif

 dudt = (temp1 - temp)*Rg/(gam1*gmw*unit_ergg)/dt

end subroutine exact_cooling_table

!-----------------------------------------------------------------------
!+
!  utility to find the index of closest value in a table
!+
!-----------------------------------------------------------------------
pure integer function find_in_table(n,table,val) result(i)
 integer, intent(in) :: n
 real,    intent(in) :: table(n), val
 integer :: i0,i1

 i0 = 0
 i1 = n + 1
 do while (i1 - i0 > 1)
    i = (i0 + i1)/2
    if ((table(n) >= table(1)).eqv.(val >= table(i))) then
       i0 = i
    else
       i1 = i
    endif
 enddo
 if (abs(val-table(1)) < tiny(0.)) then
    i = 1
 elseif (abs(val-table(n)) < tiny(0.)) then
    i = n-1
 else
    i = i0
 endif

end function find_in_table

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_cooling(iunit)
 use infile_utils, only:write_inopt
 use h2cooling,    only:write_options_h2cooling
 use part,         only:h2chemistry
 use cooling_molecular, only: write_options_molecularcooling
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options controlling cooling'
 call write_inopt(C_cool,'C_cool','factor controlling cooling timestep',iunit)
 if (h2chemistry) then
    call write_inopt(icooling,'icooling','cooling function (0=off, 1=on)',iunit)
    if (icooling > 0) then
       call write_options_h2cooling(iunit)
    endif
 else
    call write_inopt(icooling,'icooling','cooling function (0=off, 1=explicit, 2=Townsend table, 3=Gammie, 5,6=KI02)',iunit)
    select case(icooling)
    case(1)
       call write_options_molecularcooling(iunit)
       call write_inopt(excitation_HI,'excitation_HI','cooling via electron excitation of HI [1=on/0=off]',iunit)
       call write_inopt(relax_bowen,'relax_bowen','Bowen (diffusive) relaxation [1=on/0=off]',iunit)
       call write_inopt(relax_stefan,'relax_stefan','radiative relaxation [1=on/0=off]',iunit)
       call write_inopt(dust_collision,'dust_collision','dust collision [1=on/0=off]',iunit)
       call write_inopt(bowen_Cprime,'bowen_Cprime','radiative cooling rate (g.s/cm³)',iunit)
    case(2)
       call write_inopt(cooltable,'cooltable','data file containing cooling function',iunit)
       call write_inopt(habund,'habund','Hydrogen abundance assumed in cooling function',iunit)
       call write_inopt(temp_floor,'temp_floor','Minimum allowed temperature in K for Townsend cooling table',iunit)
    case(3)
       call write_inopt(beta_cool,'beta_cool','beta factor in Gammie (2001) cooling',iunit)
    end select
 endif
 if (icooling > 0) call write_inopt(Tfloor,'Tfloor','temperature floor (K); on if > 0',iunit)

end subroutine write_options_cooling

!-----------------------------------------------------------------------
!+
!  reads options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_cooling(name,valstring,imatch,igotall,ierr)
 use part,              only:h2chemistry
 use h2cooling,         only:read_options_h2cooling
 use io,                only:fatal
 use cooling_molecular, only:read_options_molecular_cooling
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 logical :: igotallh2,igotallcf,igotallmol

 imatch  = .true.
 igotall = .false.  ! cooling options are compulsory
 igotallh2 = .true.
 igotallcf = .true.
 igotallmol = .true.
 select case(trim(name))
 case('icooling')
    read(valstring,*,iostat=ierr) icooling
    ngot = ngot + 1
 case('excitation_HI')
    read(valstring,*,iostat=ierr) excitation_HI
    ngot = ngot + 1
 case('relax_bowen')
    read(valstring,*,iostat=ierr) relax_bowen
    ngot = ngot + 1
 case('relax_stefan')
    read(valstring,*,iostat=ierr) relax_stefan
    ngot = ngot + 1
 case('dust_collision')
    read(valstring,*,iostat=ierr) dust_collision
    ngot = ngot + 1
 case('C_cool')
    read(valstring,*,iostat=ierr) C_cool
    ngot = ngot + 1
 case('cooltable')
    read(valstring,*,iostat=ierr) cooltable
    ngot = ngot + 1
 case('habund')
    read(valstring,*,iostat=ierr) habund
    ngot = ngot + 1
 case('temp_floor')
    read(valstring,*,iostat=ierr) temp_floor
    ngot = ngot + 1
 case('bowen_Cprime')
    read(valstring,*,iostat=ierr) bowen_Cprime
    ngot = ngot + 1
 case('beta_cool')
    read(valstring,*,iostat=ierr) beta_cool
    ngot = ngot + 1
    if (beta_cool < 1.) call fatal('read_options','beta_cool must be >= 1')
 case('Tfloor')
    ! not compulsory to read in
    read(valstring,*,iostat=ierr) Tfloor
 case default
    imatch = .false.
    if (h2chemistry) then
       call read_options_h2cooling(name,valstring,imatch,igotallh2,ierr)
    endif
 end select
 if (icooling > 0 ) call read_options_molecular_cooling(name,valstring,imatch,igotallmol,ierr)
 if (icooling == 3 .and. ngot >= 1) igotall = .true.
 if (icooling == 2 .and. ngot >= 3) igotall = .true.
 if (icooling == 1 .and. ngot >= 9) igotall = .true.
 if (igotallh2 .and. ngot >= 1) igotall = .true.

end subroutine read_options_cooling

end module cooling
