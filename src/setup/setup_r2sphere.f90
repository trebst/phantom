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
!   - angvel           : *angular velocity in rad/s*
!   - dist_unit        : *distance unit (e.g. au)*
!   - dusttogas        : *dust-to-gas ratio*
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
 real :: density_contrast,totmass_sphere,r_sphere,cs_sphere,cs_sphere_cgs,Temperature,mu
 real :: angvel,masstoflux,dusttogas,pmass_dusttogas,ang_Bomega
 real :: rho_pert_amp,lbox
 real :: r_crit_setup,h_acc_setup,h_soft_sinksink_setup
 real(kind=8)                 :: udist,umass
 integer                      :: np,icreate_sinks_setup,ieos_in,ierr
 logical                      :: mu_not_B,cs_in_code
 logical                      :: turb = .false.
 character(len=20)            :: dist_unit,mass_unit
 character(len=10)  :: h_acc_char
 character(len= 1), parameter :: labelx(3) = (/'x','y','z'/)

contains

!----------------------------------------------------------------
!+
!  setup for a sphere-in-a-box
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use physcon,       only:pi,solarm,hours,years,au,kboltz,mass_proton_cgs
 use setup_params,  only:rhozero,npart_total,rmax,ihavesetupB
 use io,            only:master,fatal
 use unifdis,       only:set_unifdis
 use spherical,     only:set_sphere, rho_func
 use boundary,      only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use prompting,     only:prompt
 use units,         only:set_units,select_unit,utime,unit_density,unit_velocity
 use eos,           only:ieos,gmw
 use eos_barotropic,only:rhocrit0cgs
 use part,          only:igas,idust,set_particle_type
 use timestep,      only:dtmax,tmax,dtmax_dratio,dtmax_min
 use centreofmass,  only:reset_centreofmass
 use options,       only:nfulldump,rhofinal_cgs
 use kernel,        only:hfact_default
 use mpidomain,     only:i_belong
 use ptmass,        only:icreate_sinks,r_crit,h_acc,h_soft_sinksink
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real(kind=8)                 :: h_acc_in
 integer                      :: i,nx,np_in,npartsphere,npmax
 real                         :: vol_box,psep,psep_box,epotgrav
 real                         :: vol_sphere,dens_sphere,dens_medium,cs_medium,angvel_code,przero
 real                         :: t_ff,r2,area
 real                         :: rxy2,rxyz2,phi,dphi,central_density,edge_density
 logical                      :: inexists,setexists
 logical                      :: make_sinks = .true.
 character(len=120)           :: filename,filein,fileset
 character(len=40)            :: fmt
 character(len=16)            :: lattice
 procedure(rho_func), pointer :: density_func ! pointer to created function

 !--Check for existence of the .in and .setup files
 filein=trim(fileprefix)//'.in'
 inquire(file=filein,exist=inexists)
 fileset=trim(fileprefix)//'.setup'
 inquire(file=fileset,exist=setexists)

 npmax          = size(xyzh(1,:))
 np             = npmax
 Temperature    = 10.
 mu             = 2.381
 ieos_in        = 8
 r_sphere       = 0.05
 totmass_sphere = 2.
 h_acc_char  = '1.0d-3'

 !--Read values from .setup
 if (setexists) then
   call read_setupfile(fileset,ierr)
   if (ierr /= 0) then
      if (id==master) call write_setupfile(fileset)
      stop
   endif
   !--Prompt to get inputs and write to file
elseif (id==master) then
   print "(a,/)",trim(fileset)//' not found: using interactive setup'
   call get_input_from_prompts()

   call select_unit(h_acc_char,h_acc_in,ierr)
   h_acc_setup = h_acc_in
   if (ierr==0 ) h_acc = h_acc/udist
   r_crit_setup        = 5.0*h_acc_setup

   call write_setupfile(fileset)
endif
 !
 ! units
 !
 call set_units(dist=udist,mass=umass,G=1.d0)

 polyk         = kboltz*Temperature/(mu*mass_proton_cgs)*(utime/udist)**2
 vol_sphere    = 4./3.*pi*r_sphere**3
 rhozero       = totmass_sphere / vol_sphere
 dens_sphere   = rhozero
 hfact         = hfact_default
 t_ff          = sqrt(3.*pi/(32.*rhozero))  ! free-fall time (the characteristic timescale)
 epotgrav      = 3./5.*totmass_sphere**2/r_sphere      ! Gravitational potential energy
 lattice       = 'hcp'
 angvel_code = angvel*utime
 vol_box     = dxbound*dybound*dzbound
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
 
 angvel_code = angvel*utime
 vol_box     = dxbound*dybound*dzbound
 !
 ! setup particles in the sphere; use this routine to get N_sphere as close to np as possible
 !
 density_func => r2_func
 call set_sphere(trim(lattice),id,master,0.,r_sphere,psep,hfact,npart,xyzh, & !formerly lattice: 'closepacked'
               rhofunc=density_func,nptot=npart_total,&
               exactN=.true.,np_requested=np,mask=i_belong)

 npartsphere = npart
 if (np_in /= npartsphere) np = npartsphere
 !
 ! set particle properties
 !
 npartoftype(:)    = 0
 npartoftype(igas) = npart
 massoftype(igas)  = totmass_sphere/npart_total
 do i = 1,npartoftype(igas)
    call set_particle_type(i,igas)
 enddo
 !
 ! reset to centre of mass
 !
 call reset_centreofmass(npart,xyzh,vxyzu)
 !
 ! velocity field corresponding to uniform rotation
 !
 do i=1,npart
   if (turb) then

   else
      vxyzu(1,i) = -angvel_code*xyzh(2,i)
      vxyzu(2,i) =  angvel_code*xyzh(1,i)
      vxyzu(3,i) = 0.
      if (maxvxyzu >= 4) vxyzu(4,i) = 1.5*polyk
   endif
 enddo
 !
 ! set default runtime parameters if .in file does not exist
 !
 filename=trim(fileprefix)//'.in'
 dtmax = t_ff/100.  ! Since this variable can change, always reset it if running phantomsetup
 if (.not. inexists) then
    tmax            = 2.*t_ff 
    dtmax           = 0.01*t_ff
    ieos            = ieos_in
    !nfulldump       = 1
    calc_erot       = .true.
    dtmax_dratio    = 1.258
    icreate_sinks   = 1
    h_acc           = h_acc_setup
    r_crit          = r_crit_setup
    h_soft_sinksink = h_soft_sinksink_setup
 endif
 !
 !--Summarise the sphere
 !
 print "(a,i10)",' Input npart_sphere = ',np
 print "(1x,50('-'))"
 print "(a)",'  Quantity         (code units)  (physical units)'
 print "(1x,50('-'))"
 fmt = "((a,1pg10.3,3x,1pg10.3),a)"
 print fmt,' Total mass       : ',totmass_sphere,totmass_sphere*umass,' g'
 print fmt,' Mass in sphere   : ',totmass_sphere,totmass_sphere*umass,' g'
 print fmt,' Radius of sphere : ',r_sphere,r_sphere*udist,' cm'
 print fmt,' Density sphere   : ',dens_sphere,dens_sphere*unit_density,' g/cm^3'
 print fmt,' cs in sphere     : ',cs_sphere,cs_sphere_cgs,' cm/s'
 print fmt,' Free fall time   : ',t_ff,t_ff*utime/years,' yrs'
 print fmt,' Angular velocity : ',angvel_code,angvel,' rad/s'
 print fmt,' Omega*t_ff       : ',angvel_code*t_ff
 print "(1x,50('-'))"

end subroutine setpart

!----------------------------------------------------------------
!
!  Prompt user for inputs
!
!----------------------------------------------------------------
subroutine get_input_from_prompts()
   use prompting, only:prompt
   use units, only:select_unit
  
   call prompt('Enter the number of particles in the sphere',np,0,np)
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
   call prompt('Enter the mass of the cloud (in Msun)',totmass_sphere)
   call prompt('Enter the radius of the cloud (in pc)',r_sphere)
   call prompt('Enter the Temperature of the cloud (used for initial sound speed)',Temperature)
   call prompt('Enter the mean molecular mass (used for initial sound speed)',mu)
   if (maxvxyzu < 4) call prompt('Enter the EOS id (1: isothermal, 8: barotropic)',ieos_in)
   call prompt('Do you wish to have a turbulent velocity field?', turb)
   if (.not. turb) then
      call prompt('Enter angular rotation speed around z-axis in rad/s ',angvel,0.)
   endif
   call prompt('Enter the accretion radius of the sink (with units; e.g. au,pc,kpc,0.1pc) ',h_acc_char)
  end subroutine get_input_from_prompts
  !----------------------------------------------------------------
  !+
  !  write parameters to setup file
  !+
  !----------------------------------------------------------------
  subroutine write_setupfile(filename)
   use infile_utils, only: write_inopt
   character(len=*), intent(in) :: filename
   integer, parameter           :: iunit = 20
  
   print "(a)",' writing setup options file '//trim(filename)
   open(unit=iunit,file=filename,status='replace',form='formatted')
   write(iunit,"(a)") '# input file for r2-sphere setup routines'
   write(iunit,"(/,a)") '# units'
   call write_inopt(dist_unit,'dist_unit','distance unit (e.g. au)',iunit)
   call write_inopt(mass_unit,'mass_unit','mass unit (e.g. solarm)',iunit)
   write(iunit,"(/,a)") '# resolution'
   call write_inopt(np,'n_particles','number of particles in sphere',iunit)
   write(iunit,"(/,a)") '# options for sphere'
   call write_inopt(r_sphere,'r_sphere','radius of sphere in code units',iunit)
   call write_inopt(totmass_sphere,'totmass_sphere','mass of sphere in code units',iunit)
   if (.not. turb) call write_inopt(angvel,'angvel','angular velocity in rad/s',iunit)
   call write_inopt(turb,'use_turb','turbulent velocity field', iunit)
   write(iunit,"(/,a)") '# options required for initial sound speed'
   call write_inopt(Temperature,'Temperature','Temperature',iunit)
   call write_inopt(mu,'mu','mean molecular mass',iunit)
   write(iunit,"(/,a)") '# Sink properties (values in .in file, if present, will take precedence)'
   call write_inopt(h_acc_setup,'h_acc','accretion radius (code units)',iunit)
   call write_inopt(r_crit_setup,'r_crit','critical radius (code units)',iunit)
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
   integer                       :: nerr
   type(inopts), allocatable     :: db(:)
  
   print "(a)",' reading setup options from '//trim(filename)
   call open_db_from_file(db,filename,iunit,ierr)
   call read_inopt(np,'n_particles',db,ierr)
   call read_inopt(mass_unit,'mass_unit',db,ierr)
   call read_inopt(dist_unit,'dist_unit',db,ierr)
   call read_inopt(r_sphere,'r_sphere',db,ierr)
   call read_inopt(totmass_sphere,'totmass_sphere',db,ierr)

   call read_inopt(Temperature,'Temperature',db,ierr)
   call read_inopt(mu,'mu',db,ierr)
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
      print "(1x,a,i2,a)",'Setup_r2sphere: ',nerr,' error(s) during read of setup file.  Re-writing.'
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
