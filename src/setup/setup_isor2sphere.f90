!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setup
!
! This module sets up a sphere-in-a-box: a cold, dense sphere placed in
!   a warm medium; the two media are in pressure-equilibrium.
!   This currently works for gas-only and two-fluid dust.
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - ang_Bomega       : *Angle (degrees) between B and rotation axis*
!   - angvel           : *angular velocity in rad/s*
!   - cs_sphere_cgs    : *sound speed in sphere in cm/s*
!   - density_contrast : *density contrast in code units*
!   - dist_unit        : *distance unit (e.g. au)*
!   - dusttogas        : *dust-to-gas ratio*
!   - form_binary      : *the intent is to form a central binary*
!   - h_acc            : *accretion radius (code units)*
!   - h_soft_sinksink  : *sink-sink softening radius (code units)*
!   - icreate_sinks    : *1: create sinks.  0: do not create sinks*
!   - lbox             : *length of a box side in terms of spherical radii*
!   - mass_unit        : *mass unit (e.g. solarm)*
!   - masstoflux       : *mass-to-magnetic flux ratio in units of critical value*
!   - np               : *requested number of particles in sphere*
!   - pmass_dusttogas  : *dust-to-gas particle mass ratio*
!   - r_crit           : *critical radius (code units)*
!   - r_sphere         : *radius of sphere in code units*
!   - rho_pert_amp     : *amplitude of density perturbation*
!   - totmass_sphere   : *mass of sphere in code units*
!
! :Dependencies: boundary, centreofmass, dim, eos, eos_barotropic,
!   infile_utils, io, kernel, mpidomain, options, part, physcon, prompting,
!   ptmass, rho_profile, setup_params, spherical, timestep, unifdis, units
!
 use part,    only:periodic
 use dim,     only:use_dust,maxvxyzu,periodic
 use options, only:calc_erot
 implicit none
 public :: setpart

 private
 !--private module variables
 real :: xmini(3), xmaxi(3)
 real :: density_contrast,totmass_sphere,r_sphere,cs_sphere,cs_sphere_cgs
 real :: angvel,masstoflux,dusttogas,pmass_dusttogas,ang_Bomega
 real :: rho_pert_amp,lbox
 real :: r_crit_setup,h_acc_setup,h_soft_sinksink_setup
 real(kind=8)                 :: udist,umass
 integer                      :: np,icreate_sinks_setup
 logical                      :: binary,mu_not_B,cs_in_code
 character(len=20)            :: dist_unit,mass_unit
 character(len= 1), parameter :: labelx(3) = (/'x','y','z'/)

contains

!----------------------------------------------------------------
!+
!  setup for a sphere-in-a-box
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use physcon,      only:pi,solarm,hours,years,au
 use setup_params, only:rhozero,npart_total,rmax,ihavesetupB
 use io,           only:master,fatal
 use unifdis,      only:set_unifdis
 use spherical,    only:set_sphere, rho_func
 use boundary,     only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use prompting,    only:prompt
 use units,        only:set_units,select_unit,utime,unit_density,unit_velocity
 use eos,          only:polyk2,ieos,gmw
 use eos_barotropic, only:rhocrit0cgs
 use part,         only:igas,idust,set_particle_type
 use timestep,     only:dtmax,tmax,dtmax_dratio,dtmax_min
 use centreofmass, only:reset_centreofmass
 use options,      only:nfulldump,rhofinal_cgs
 use kernel,       only:hfact_default
 use mpidomain,    only:i_belong
 use ptmass,       only:icreate_sinks,r_crit,h_acc,h_soft_sinksink
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real(kind=8)       :: h_acc_in
 integer            :: i,nx,np_in,npartsphere,npmax,ierr
 real               :: totmass,vol_box,psep,psep_box
 real               :: vol_sphere,dens_sphere,dens_medium,cs_medium,angvel_code,przero
 real               :: totmass_box,t_ff,r2,area
 real               :: rxy2,rxyz2,phi,dphi,central_density,edge_density
 logical            :: iexist
 logical            :: make_sinks = .true.
 logical            :: turb
 character(len=100) :: filename
 character(len=40)  :: fmt
 character(len=10)  :: h_acc_char
 procedure(rho_func), pointer :: density_func ! pointer to created function

 npmax    = size(xyzh(1,:))
 filename = trim(fileprefix)//'.setup'
 print "(/,1x,63('-'),1(/,a),/,1x,63('-'),/)",&
   '  Sphere-in-box setup: Almost Archimedes'' greatest achievement.'
 inquire(file=filename,exist=iexist)
 if (iexist) then
    call read_setupfile(filename,ierr)
    np_in = np
    if (ierr /= 0) then
       if (id==master) call write_setupfile(filename)
       stop
    endif
 elseif (id==master) then
    print "(a,/)",trim(filename)//' not found: using interactive setup'
    dist_unit = '1.0d16cm'
    mass_unit = 'solarm'
    ierr = 1
    do while (ierr /= 0)
       call prompt('Enter mass unit (e.g. solarm,jupiterm,earthm)',mass_unit)
       call select_unit(mass_unit,umass,ierr)
       if (ierr /= 0) print "(a)",' ERROR: mass unit not recognised'
    enddo
    ierr = 1
    do while (ierr /= 0)
       call prompt('Enter distance unit (e.g. au,pc,kpc,0.1pc)',dist_unit)
       call select_unit(dist_unit,udist,ierr)
       if (ierr /= 0) print "(a)",' ERROR: length unit not recognised'
    enddo
    !
    ! units
    !
    call set_units(dist=udist,mass=umass,G=1.d0)
    !
    ! prompt user for settings
    !
    npmax = int(2.0/3.0*size(xyzh(1,:))) ! approx max number allowed in sphere given size(xyzh(1,:))
    if (npmax < 300000) then
       np = npmax
    elseif (npmax < 1000000) then
       np = 300000
    else
       np = 1000000
    endif
    
    call prompt('Enter the approximate number of particles in the sphere',np,0,npmax)
    np_in    = np
    r_sphere = 4.
    call prompt('Enter radius of sphere in units of '//dist_unit,r_sphere,0.)
    lbox     = 4.
    call prompt('Enter the box size in units of spherical radii: ',lbox,1.)
    do i=1,3
      ! note that these values will be saved to the .setup file rather than lbox so user can convert
      ! box to a rectangle if desired
      xmini(i) = -0.5*(lbox*r_sphere)
      xmaxi(i) = -xmini(i)
    enddo

   totmass_sphere = 1.0
   call prompt('Enter total mass in sphere in units of '//mass_unit,totmass_sphere,0.)

    density_contrast = 30.0
    call prompt('Enter density contrast between sphere and box ',density_contrast,1.)

    binary = .false.
    call prompt('Do you intend to form a binary system (i.e. add an m=2 perturbation)?',binary)

    if (binary) then
       cs_sphere_cgs = 18696.96 ! cm/s ~ 5K assuming mu = 2.31 & gamma = 5/3
    else
       cs_sphere_cgs = 21888.0  ! cm/s ~ 8K assuming mu = 2.31 & gamma = 5/3
    endif
    call prompt('Enter sound speed in sphere in units of cm/s',cs_sphere_cgs,0.)

    if (binary) then
       angvel = 1.006d-12
    else
       angvel = 1.77d-13
    endif
    call prompt('Enter angular rotation speed in rad/s ',angvel,0.)

    ! ask about sink particle details; these will not be saved to the .setup file since they exist in the .in file
    !
    call prompt('Do you wish to dynamically create sink particles? ',make_sinks)
    if (make_sinks) then
       if (binary) then
          h_acc_char  = '3.35au'
       else
          h_acc_char  = '1.0d-3'
       endif
       call prompt('Enter the accretion radius of the sink (with units; e.g. au,pc,kpc,0.1pc) ',h_acc_char)
       call select_unit(h_acc_char,h_acc_in,ierr)
       h_acc_setup = h_acc_in
       if (ierr==0 ) h_acc = h_acc/udist
       r_crit_setup        = 5.0*h_acc_setup
       icreate_sinks_setup = 1
       if (binary) h_soft_sinksink_setup = 0.4*h_acc_setup
    else
       icreate_sinks_setup = 0
    endif
    if (id==master) call write_setupfile(filename)
    stop 'please edit .setup file and rerun phantomsetup'
 else
    stop ! MPI, stop on other threads, interactive on master
 endif
 !
 ! units
 !
 call set_units(dist=udist,mass=umass,G=1.d0)
 !
 ! convert units of sound speed
 !
 if (cs_in_code) then
    cs_sphere_cgs = cs_sphere*unit_velocity
 else
    cs_sphere     = cs_sphere_cgs/unit_velocity
 endif
 !
 ! boundaries
 !
 call set_boundary(xmini(1),xmaxi(1),xmini(2),xmaxi(2),xmini(3),xmaxi(3))
 !
 ! general parameters
 !
 time        = 0.
 hfact       = hfact_default
 if (maxvxyzu >=4 ) then
    gamma    = 5./3.
 else
    gamma    = 1.
 endif
 rmax        = r_sphere
 angvel_code = angvel*utime
 vol_box     = dxbound*dybound*dzbound
 vol_sphere  = 4./3.*pi*r_sphere**3
 rhozero     = totmass_sphere / vol_sphere
 dens_sphere = rhozero
 dens_medium = dens_sphere/density_contrast
 cs_medium   = cs_sphere*sqrt(density_contrast)
 totmass_box = (vol_box - vol_sphere)*dens_medium
 totmass     = totmass_box + totmass_sphere
 t_ff        = sqrt(3.*pi/(32.*dens_sphere))
 !
 ! setup particles in the sphere; use this routine to get N_sphere as close to np as possible
 !
 density_func => r2_func
 call set_sphere('closepacked',id,master,0.,r_sphere,psep,hfact,npart,xyzh, &
               rhofunc=density_func,nptot=npart_total,&
               exactN=.true.,np_requested=np,mask=i_belong)

 npartsphere = npart
 if (np_in /= npartsphere) np = npartsphere
 !
 ! setup surrounding low density medium
 !
 !psep_box = psep*(density_contrast)**(1./3.)  ! calculate psep in box
 !call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep_box, &
 !                  hfact,npart,xyzh,periodic,rmin=r_sphere,nptot=npart_total,mask=i_belong,err=ierr)
 !print "(a,es10.3)",' Particle separation in low density medium = ',psep_box
 !print "(a,i10,a)",' added ',npart-npartsphere,' particles in low-density medium'
 !print*, ""
 !
 ! set particle properties
 !
 npartoftype(:)    = 0
 npartoftype(igas) = npart
 massoftype(igas)  = totmass/npart_total
 do i = 1,npartoftype(igas)
    call set_particle_type(i,igas)
 enddo
 !
 ! reset to centre of mass
 !
 call reset_centreofmass(npart,xyzh,vxyzu)
 !
 ! temperature set to give a pressure equilibrium
 !
 polyk  = cs_sphere**2
 polyk2 = cs_medium**2
 !
 !--Stretching the spatial distribution to perturb the density profile (for binary==.true. only)
 !
 if (binary) then
    do i = 1,npart
       rxy2  = xyzh(1,i)*xyzh(1,i) + xyzh(2,i)*xyzh(2,i)
       rxyz2 = rxy2 + xyzh(3,i)*xyzh(3,i)
       if (rxyz2 <= r_sphere**2) then
          phi       = atan(xyzh(2,i)/xyzh(1,i))
          if (xyzh(1,i) < 0.0) phi = phi + pi
          dphi      = 0.5*rho_pert_amp*sin(2.0*phi)
          phi       = phi - dphi
          xyzh(1,i) = sqrt(rxy2)*cos(phi)
          xyzh(2,i) = sqrt(rxy2)*sin(phi)
       endif
    enddo
 endif
 !
 ! velocity field corresponding to uniform rotation
 !
 do i=1,npart
   r2 = dot_product(xyzh(1:3,i),xyzh(1:3,i))
   if (r2 < r_sphere**2) then
      vxyzu(1,i) = -angvel_code*xyzh(2,i)
      vxyzu(2,i) = angvel_code*xyzh(1,i)
      vxyzu(3,i) = 0.
      if (maxvxyzu >= 4) vxyzu(4,i) = 1.5*polyk
   else
      vxyzu(1:3,i) = 0.
      if (maxvxyzu >= 4) vxyzu(4,i) = 1.5*polyk2
   endif
 enddo
 !
 ! set default runtime parameters if .in file does not exist
 !
 filename=trim(fileprefix)//'.in'
 inquire(file=filename,exist=iexist)
 dtmax = t_ff/100.  ! Since this variable can change, always reset it if running phantomsetup
 if (.not. iexist) then
    if (binary) then
       tmax      = 1.50*t_ff ! = 13.33 for default settings
    else
       tmax      = 1.21*t_ff ! = 10.75 for default settings
    endif
    ieos         = 1
    nfulldump    = 1
    calc_erot    = .true.
    dtmax_dratio = 1.258
    icreate_sinks   = icreate_sinks_setup
    r_crit          = r_crit_setup
    h_acc           = h_acc_setup
    h_soft_sinksink = h_soft_sinksink_setup
    if (icreate_sinks==1) then
       dtmax_min = dtmax/8.0
    else
       dtmax_min = 0.0
       rhofinal_cgs = 0.15
    endif
 endif
 !
 !--Summarise the sphere
 !
 print "(a,i10)",' Input npart_sphere = ',np
 print "(1x,50('-'))"
 print "(a)",'  Quantity         (code units)  (physical units)'
 print "(1x,50('-'))"
 fmt = "((a,1pg10.3,3x,1pg10.3),a)"
 print fmt,' Total mass       : ',totmass,totmass*umass,' g'
 print fmt,' Mass in sphere   : ',totmass_sphere,totmass_sphere*umass,' g'
 print fmt,' Radius of sphere : ',r_sphere,r_sphere*udist,' cm'
 print fmt,' Density sphere   : ',dens_sphere,dens_sphere*unit_density,' g/cm^3'
 print fmt,' Density medium   : ',dens_medium,dens_medium*unit_density,' g/cm^3'
 print fmt,' cs in sphere     : ',cs_sphere,cs_sphere_cgs,' cm/s'
 print fmt,' cs in medium     : ',cs_medium,cs_medium*unit_velocity,' cm/s'
 print fmt,' Free fall time   : ',t_ff,t_ff*utime/years,' yrs'
 print fmt,' Angular velocity : ',angvel_code,angvel,' rad/s'
 print fmt,' Omega*t_ff       : ',angvel_code*t_ff
 print "(1x,50('-'))"

end subroutine setpart

!----------------------------------------------------------------
!+
!  write parameters to setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only: write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter           :: iunit = 20
 integer                      :: i

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for sphere-in-box setup routines'
 write(iunit,"(/,a)") '# units'
 call write_inopt(dist_unit,'dist_unit','distance unit (e.g. au)',iunit)
 call write_inopt(mass_unit,'mass_unit','mass unit (e.g. solarm)',iunit)
 write(iunit,"(/,a)") '# resolution'
 call write_inopt(np,'np','requested number of particles in sphere',iunit)
 write(iunit,"(/,a)") '# options for box'
 do i=1,3
    call write_inopt(xmini(i),labelx(i)//'min',labelx(i)//' min',iunit)
    call write_inopt(xmaxi(i),labelx(i)//'max',labelx(i)//' max',iunit)
 enddo
 write(iunit,"(/,a)") '# intended result'
 call write_inopt(binary,'form_binary','the intent is to form a central binary',iunit)
 write(iunit,"(/,a)") '# options for sphere'
 !call write_inopt(BEsphere,'use_BE_sphere','centrally condense as a BE sphere',iunit)
 call write_inopt(r_sphere,'r_sphere','radius of sphere in code units',iunit)
 call write_inopt(totmass_sphere,'totmass_sphere','mass of sphere in code units',iunit)
 call write_inopt(density_contrast,'density_contrast','density contrast in code units',iunit)
 call write_inopt(cs_sphere_cgs,'cs_sphere_cgs','sound speed in sphere in cm/s',iunit)
 call write_inopt(angvel,'angvel','angular velocity in rad/s',iunit)
 write(iunit,"(/,a)") '# Sink properties (values in .in file, if present, will take precedence)'
 call write_inopt(icreate_sinks_setup,'icreate_sinks','1: create sinks.  0: do not create sinks',iunit)
 if (icreate_sinks_setup==1) then
    call write_inopt(h_acc_setup,'h_acc','accretion radius (code units)',iunit)
    call write_inopt(r_crit_setup,'r_crit','critical radius (code units)',iunit)
    if (binary) then
       call write_inopt(h_soft_sinksink_setup,'h_soft_sinksink','sink-sink softening radius (code units)',iunit)
    endif
 endif
 close(iunit)

end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  Read parameters from setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only: open_db_from_file,inopts,read_inopt,close_db
 use io,           only: error
 use units,        only: select_unit
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit = 21
 integer                       :: i,nerr,jerr,kerr
 type(inopts), allocatable     :: db(:)

 !--Read values
 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(mass_unit,'mass_unit',db,ierr)
 call read_inopt(dist_unit,'dist_unit',db,ierr)
 call read_inopt(binary,'form_binary',db,ierr)
 call read_inopt(np,'np',db,ierr)
 
 do i=1,3
    call read_inopt(xmini(i),labelx(i)//'min',db,ierr)
    call read_inopt(xmaxi(i),labelx(i)//'max',db,ierr)
 enddo
 call read_inopt(r_sphere,'r_sphere',db,ierr)
 call read_inopt(totmass_sphere,'totmass_sphere',db,ierr)
 lbox = -2.0*xmini(1)/r_sphere

 call read_inopt(density_contrast,'density_contrast',db,ierr)
 call read_inopt(cs_sphere,'cs_sphere',db,jerr)
 call read_inopt(cs_sphere_cgs,'cs_sphere_cgs',db,kerr)
 cs_in_code = .false.  ! for backwards compatibility
 if (jerr /= 0 .and. kerr == 0) then
    cs_in_code = .false.
 elseif (jerr == 0 .and. kerr /= 0) then
    cs_in_code = .true.
 else
    ierr = ierr + 1
 endif
 call read_inopt(angvel,'angvel',db,ierr)
 mu_not_B = .true.
 call read_inopt(icreate_sinks_setup,'icreate_sinks',db,ierr)
 if (icreate_sinks_setup==1) then
    call read_inopt(h_acc_setup,'h_acc',db,ierr)
    call read_inopt(r_crit_setup,'r_crit',db,ierr)
    if (binary) then
       call read_inopt(h_soft_sinksink_setup,'h_soft_sinksink',db,ierr)
    endif
 endif
 call close_db(db)
 !
 ! parse units
 !
 call select_unit(mass_unit,umass,nerr)
 if (nerr /= 0) then
    call error('setup_sphereinbox','mass unit not recognised')
    ierr = ierr + 1
 endif
 call select_unit(dist_unit,udist,nerr)
 if (nerr /= 0) then
    call error('setup_sphereinbox','length unit not recognised')
    ierr = ierr + 1
 endif

 if (ierr > 0) then
    print "(1x,a,i2,a)",'Setup_sphereinbox: ',nerr,' error(s) during read of setup file.  Re-writing.'
 endif

end subroutine read_setupfile

!----------------------------------------------------------------
!
! Spherical density profile as a function of radius
!
!----------------------------------------------------------------
real function r2_func(r)
   real, intent(in) :: r

   r2_func = 1./r**2

end function r2_func


end module setup
