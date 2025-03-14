#--------------------------------------------------------------------------!
# The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
# Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
# See LICENCE file for usage and distribution conditions                   !
# http://users.monash.edu.au/~dprice/phantom                               !
#--------------------------------------------------------------------------!
#+
#  The Phantom Makefile
#
#  DESCRIPTION:
#   This is the main Makefile for all of the code and utilities
#   Compiler settings are grouped under the SYSTEM variable while
#   compile-time settings for different problems are grouped under
#   the SETUP variable
#
#  OWNER: Daniel Price
#
#  $Id: 2788b71b1c08e560e77dce9849c5cb24a668f4b9 $
#+
#--------------------------------------------------------------------------

.KEEP_STATE:

PHANTOM_VERSION_MAJOR=2022
PHANTOM_VERSION_MINOR=0
PHANTOM_VERSION_MICRO=1
VERSION=$(PHANTOM_VERSION_MAJOR).$(PHANTOM_VERSION_MINOR).$(PHANTOM_VERSION_MICRO)

KNOWN_SYSTEM=no
SHELL = /bin/bash
VPATH = "${RUNDIR}" ../src/main ../src/utils ../src/setup ../src/tests ../src/lib/NICIL/src
BINDIR= ../bin
UNAME=${shell uname}
#----------------------------------------------------------------
# Sensible defaults for phantom configuration
#----------------------------------------------------------------
CONFIG        = config.F90
SETUPFILE     = setup_unifdis.f90
MODFILE       = moddump_default.f90
ANALYSIS      = analysis_dtheader.f90
MULTIRUNFILE  = multirun.f90
LIVE_ANALYSIS = no
DEBUG         = no
#
# can comment out the following lines and instead set
# the parameters as environment variables
#
ifndef DOUBLEPRECISION
    DOUBLEPRECISION= yes
endif
ifndef EDITOR
    EDITOR= emacs
endif
ifndef OPENMP
    OPENMP= yes
endif
ifndef SPLASH_DIR
    SPLASH_DIR=${shell if [ -d $$HOME/splash ]; then echo $$HOME/splash; fi}
endif
#MPI= yes
#
# endian can be "BIG", "LITTLE" or anything else which has no effect
#
# ENDIAN= default
#
CC = gcc
CCFLAGS = -O5
LIBCXX = -lstdc++
#FPPFLAGS=
LDFLAGS=
SRCPHOTO=
LIBTOOL=libtool

#----------------------------------------------------------------
# here follows specific configuration options used
# for various types of simulations
#
# preprocessor options are as follows:
#
# -DPERIODIC            ! periodic boundaries
# -DIND_TIMESTEPS       ! individual particle timesteps
# -DSTS_TIMESTEPS       ! super-timestepping
# -DDISC_VISCOSITY      ! use artificial disc viscosity ( nu \propto alpha_sph cs h
#                       ! and calculated for both approaching and receding particles
# -DDRIVING             ! use turbulence driving
# -DMHD                 ! magnetic fields
# -DNONIDEALMHD         ! non-ideal magnetic fields including ionisation; uses NICIL
# -DPHOTO               ! turn on the photoevaporation
# -DLIGHTCURVE          ! lightcurve estimation

include Makefile_setups

ifeq ($(SYSTEM), cipcluster)
    include Makefile_defaults_gfortran
    QSYS = slurm
    QPROJECT='p01'
endif 

ifndef SETUPFILE
    SETUPFILE= setup_unifdis.f90
endif

ifndef SRCNIMHD
    SRCNIMHD = nicil.F90 nicil_supplement.F90
endif

ifndef SRCDUST
    SRCDUST = dust.F90 growth.F90
endif

ifdef SMOL
    SRCDUST+= growth_smol.f90
endif

#ifndef SRCGROWTH
#    SRCGROWTH = growth.F90
#endif

#---  live feedback from mcfost
ifeq ($(MCFOST), yes)
    ANALYSIS= analysis_mcfost.f90
    LIVE_ANALYSIS=yes
    ISOTHERMAL=no
    MCFOST=yes

    MCFOST_LIBS = $(MCFOST_INSTALL)/lib/$(FC)
    MCFOST_INCLUDE = $(MCFOST_INSTALL)/include

    FPPFLAGS+= -DMCFOST
    LDFLAGS+= -I$(MCFOST_INCLUDE) -I$(MCFOST_INCLUDE)/voro++ -I$(MCFOST_INCLUDE)/hdf5 -I$(MCFOST_INCLUDE)/$(FC) \
	-L$(MCFOST_DIR)/src -lmcfost -L$(MCFOST_LIBS) $(LIBCXX) -lcfitsio -lvoro++ -lsprng \
	-L$(HDF5ROOT)/lib/Intel -lhdf5_fortran -L$(HDF5ROOT)/lib -lhdf5 -lz -lxgboost -ldmlc -lrabit -lhdf5_fortran -lhdf5 -lz
endif

include Makefile_systems

FFLAGS+=-fPIC
# Set some default files if not defined above
ifdef MAXP
    FPPFLAGS += -DMAXP=${MAXP}
endif
ifdef MAXPTMASS
    FPPFLAGS += -DMAXPTMASS=${MAXPTMASS}
endif
ifdef MAXNEIGH
    FPPFLAGS += -DMAXNEIGH=${MAXNEIGH}
endif
ifdef NCELLSMAX
    FPPFLAGS += -DNCELLSMAX=${NCELLSMAX}
endif
ifdef STACKSIZE
    FPPFLAGS += -DSTACKSIZE=${STACKSIZE}
endif
# Set other optional flags depending on settings

ifeq ($(DEBUG), yes)
    FFLAGS += ${DEBUGFLAG}
    FFLAGS := $(FFLAGS:-O3=-O0)
    FFLAGS := $(FFLAGS:-ipo= )
endif

ifeq ($(ENDIAN), BIG)
    FFLAGS += ${ENDIANFLAGBIG}
endif

ifeq ($(ENDIAN), LITTLE)
    FFLAGS += ${ENDIANFLAGLITTLE}
endif

ifeq ($(OPENMP), yes)
    FFLAGS += ${OMPFLAGS}
endif

ifeq ($(PERIODIC), yes)
    FPPFLAGS += -DPERIODIC
endif

ifeq ($(GRAVITY), yes)
    FPPFLAGS += -DGRAVITY
endif

ifeq ($(ISOTHERMAL), yes)
    FPPFLAGS += -DISOTHERMAL
endif

ifeq ($(MHD), yes)
    FPPFLAGS += -DMHD
endif

ifeq ($(GR), yes)
    FPPFLAGS += -DGR
ifeq ($(METRIC), kerr)
    FPPFLAGS += -DKERR
endif
ifeq ($(ISENTROPIC), yes)
    FPPFLAGS += -DISENTROPIC
    FPPFLAGS += -DLIGHTCURVE
endif
    CONST_AV=yes
endif

ifeq ($(DUST), yes)
    FPPFLAGS += -DDUST
    ifndef KERNEL
        KERNEL=quintic
    endif
endif

ifdef MAXDUSTSMALL
    FPPFLAGS += -DMAXDUSTSMALL=${MAXDUSTSMALL}
endif
ifdef MAXDUSTLARGE
    FPPFLAGS += -DMAXDUSTLARGE=${MAXDUSTLARGE}
endif

ifeq ($(DUSTGROWTH), yes)
    FPPFLAGS += -DDUSTGROWTH
endif

ifeq ($(SMOL), yes)
    FPPFLAGS += -DSMOL
    LDFLAGS += -L$(SMOL_DIR) -lsmol
endif

ifeq ($(NONIDEALMHD), yes)
    FPPFLAGS += -DNONIDEALMHD
endif

ifeq ($(H2CHEM), yes)
    FPPFLAGS += -DH2CHEM
endif

ifeq ($(DISC_VISCOSITY), yes)
    FPPFLAGS += -DDISC_VISCOSITY
endif

ifeq ($(CONST_AV), yes)
    FPPFLAGS += -DCONST_AV
endif

ifeq ($(MORRIS_MONAGHAN), yes)
    FPPFLAGS += -DUSE_MORRIS_MONAGHAN
endif

ifeq ($(CONST_ARTRES), yes)
    FPPFLAGS += -DCONST_ARTRES
endif

ifeq ($(CURLV), yes)
    FPPFLAGS += -DCURLV
endif

ifeq ($(IND_TIMESTEPS), yes)
    FPPFLAGS += -DIND_TIMESTEPS
endif

ifeq ($(STS_TIMESTEPS), yes)
    FPPFLAGS += -DSTS_TIMESTEPS
endif

ifeq ($(DEBUG), yes)
    FFLAGS += -DDEBUG
endif

ifeq ($(RADIATION), yes)
    FPPFLAGS += -DRADIATION
endif

ifdef SRCTURB
    FPPFLAGS += -DDRIVING
endif

ifeq ($(KROME), krome)
    FPPFLAGS += -DKROME
ifeq ($(SYSTEM), ifort)
    LDFLAGS  += -mkl
else
    LDFLAGS  += -L/usr/lib/x86_64-linux-gnu -lmkl_core -lmkl_gnu_thread -lmkl_gf_lp64 -fopenmp
endif
endif

#
# kernel choice
#
ifndef SRCKERNEL
ifdef KERNEL
    SRCKERNEL= kernel_${KERNEL}.f90
else
    SRCKERNEL= kernel_cubic.f90
    KERNEL=cubic
endif
endif

#
# can turn particle injection off
# by setting INJECT_PARTICLES=no
# on command line. Otherwise on
# if injection module selected
#
ifeq ($(INJECT_PARTICLES), no)
    SRCINJECT=
else
ifdef SRCINJECT
    FPPFLAGS += -DINJECT_PARTICLES
endif
endif

ifdef LIGHTCURVE
    FPPFLAGS += -DLIGHTCURVE
endif

# do double precision flag last (append only to FFLAGS)

ZZFFLAGS := ${FFLAGS}
ifeq ($(DOUBLEPRECISION), yes)
    FFLAGS += ${DBLFLAG}
endif

ifeq ($(ANALYSISONLY), yes)
    FPPFLAGS += -DANALYSIS
endif

ifeq ($(LIVE_ANALYSIS), yes)
    FPPFLAGS += -DLIVE_ANALYSIS
    SRCAN = $(ANALYSIS)
else
    SRCAN=
endif

#
# MPI flavour (mostly openmpi these days)
#
ifeq ($(MPI), yes)
    FC= mpif90 `mpif90 --showme:compile`
    CC= mpicc `mpicc --showme:compile`
    USEMPI=yes
endif

ifeq ($(MPI), openmpi)
    FC= openmpif90 `openmpif90 --showme:compile`
    LDFLAGS+= `openmpif90 --showme:link`
    USEMPI=yes
endif

ifeq ($(MPI), zen)
    FC= mpif90
    LDFLAGS+= -lwmpi -lmpiif
    USEMPI=yes
endif

ifeq ($(MPI), psxe)
    FC= mpiifort
    LDFLAGS+= `mpiifort--showme:link`
    USEMPI=yes
endif

ifeq ($(MPI), mpifort)
    FC= mpifort
    USEMPI=yes
endif

ifeq ($(MPI), mpiifort)
    FC= mpiifort
    USEMPI=yes
endif

ifeq ($(MPI), intel)
    FC= mpif90
    USEMPI=yes
endif

ifeq ($(MPI), cray)
    FC= ftn
    CC= cc
    USEMPI=yes
endif

ifeq ($(USEMPI), yes)
    RUNMPI=$(MPIEXEC)
    FPPFLAGS += -DMPI
else
    RUNMPI=
endif

#
# HDF5 libraries (if required)
#
# Requires two directories:
#   - include for Fortran .mod files
#   - lib for the shared library .so files
#
# Often both directories are under one root,
# e.g. HDF5ROOT= /usr/local/opt/hdf5
# In this case just set HDF5ROOT for your machine.
#
# However, sometimes these directories are separate,
# then you must set both HDF5INCLUDE and HDF5LIB.
#
ifeq ($(HDF5), yes)
ifeq (X$(HDF5ROOT), X)
    HDF5ROOT= /usr/local/opt/hdf5
endif
ifeq (X$(HDF5INCLUDE), X)
    HDF5INCLUDE= $(HDF5ROOT)/include
endif
ifeq (X$(HDF5LIB), X)
    HDF5LIB= $(HDF5ROOT)/lib
endif
    FFLAGS+= -I$(HDF5INCLUDE)
    CCFLAGS+= -I$(HDF5INCLUDE)
    LDFLAGS+= -L$(HDF5LIB) -lhdf5 -lhdf5_fortran
    FPPFLAGS+= -DHDF5
endif

IDFLAGS=$(FPPFLAGS)
ifeq ($(DEBUG), yes)
    IDFLAGS += -DDEBUG
endif
#
# select domain decomposition type
#
DOMAIN= mpi_domain.F90
OBJDIR=obj

# define the implicit rule to make a .o file from a .f90 file

.SUFFIXES:
.SUFFIXES: .o .f90 .F90 .c .f

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

%.o : %.F90
	$(FC) -c $(FFLAGS) ${FPP_PREFIX} $(FPPFLAGS) $< -o $@

%.o : %.c
	$(CC) -c $(CCFLAGS) $< -o $@

%.o : %.f
	$(FC) -c $(FFLAGS) $< -o $@

#
# external forces
#
ifeq (X$(SRCPOTS), X)
SRCPOTS= extern_corotate.f90 \
         extern_binary.f90 \
         extern_spiral.f90 \
         extern_lensethirring.f90 \
         extern_gnewton.F90 \
         lumin_nsdisc.F90 extern_prdrag.F90 \
         extern_Bfield.f90 \
         extern_densprofile.f90 \
         extern_staticsine.f90 \
         extern_gwinspiral.f90 \
         externalforces.F90
endif
ifeq (X$(SRCPOT), X)
SRCPOT=${SRCPOTS}
endif
#
# metrics for GR
#
ifeq ($(GR),yes)
    SRCPOT=extern_gr.F90 $(SRCPOTS:externalforces.F90=externalforces_gr.F90)
endif
ifdef METRIC
    SRCMETRIC= metric_${METRIC}.f90
else
    SRCMETRIC= metric_minkowski.f90
endif
SRCGR=inverse4x4.f90 $(SRCMETRIC) metric_tools.f90 utils_gr.f90
#
# chemistry
#
SRCCHEM= fs_data.f90 mol_data.f90 utils_spline.f90 h2cooling.f90 cooling_molecular.f90 h2chem.f90 cooling.F90
#
# equations of state
#
SRCMESA= eos_mesa_microphysics.f90 eos_mesa.f90
SRCEOS = eos_barotropic.f90 eos_piecewise.f90 ${SRCMESA} eos_shen.f90 eos_helmholtz.f90 eos_idealplusrad.f90 ionization.F90 eos_gasradrec.f90 eos.F90

ifeq ($(HDF5), yes)
    SRCREADWRITE_DUMPS= utils_hdf5.f90 utils_dumpfiles_hdf5.f90 readwrite_dumps_common.F90 readwrite_dumps_fortran.F90 readwrite_dumps_hdf5.F90 readwrite_dumps.F90
else
    SRCREADWRITE_DUMPS= readwrite_dumps_common.F90 readwrite_dumps_fortran.F90 readwrite_dumps.F90
endif

ifeq ($(KROME), krome)
    SRCKROME=krome.f90
else
    SRCKROME=
endif

SOURCES= physcon.f90 ${CONFIG} ${SRCKERNEL} io.F90 units.f90 boundary.f90 \
         mpi_utils.F90 dtype_kdtree.F90 utils_omp.F90 utils_cpuinfo.f90 \
         utils_allocate.f90 icosahedron.f90 \
         utils_system.f90 utils_mathfunc.f90 part.F90 ${DOMAIN} utils_timing.f90 mpi_balance.F90 \
         setup_params.f90 timestep.f90 utils_dumpfiles.f90 utils_indtimesteps.F90 utils_infiles.f90 \
         utils_sort.f90 utils_supertimestep.F90 utils_tables.f90 utils_gravwave.f90 \
         utils_sphNG.f90 utils_vectors.f90 utils_datafiles.f90 datafiles.f90 \
         gitinfo.f90 ${SRCFASTMATH} random.f90 ${SRCEOS} cullendehnen.f90 ${SRCNIMHD} ${SRCGR} \
			checkoptions.F90 viscosity.f90 options.f90 cons2primsolver.f90 radiation_utils.f90 cons2prim.f90 \
         centreofmass.f90 ${SRCPOT} damping.f90 checkconserved.f90 \
         partinject.F90 utils_inject.f90 utils_filenames.f90 utils_summary.F90 ${SRCCHEM}\
         ${SRCDUST} dust_formation.F90 ptmass_radiation.F90 ptmass_heating.F90 \
         mpi_dens.F90 mpi_force.F90 mpi_memory.F90 mpi_derivs.F90 kdtree.F90 linklist_kdtree.F90 ${SRCTURB} \
         ${SRCPHOTO} ${SRCINJECT} ${SRCKROME} memory.F90 ${SRCREADWRITE_DUMPS}  \
         quitdump.f90 ptmass.F90 \
         readwrite_infile.F90 dens.F90 force.F90 utils_deriv.f90 deriv.F90 energies.F90 sort_particles.F90 \
         utils_shuffleparticles.F90 evwrite.F90 step_leapfrog.F90 writeheader.F90 ${SRCAN} step_supertimestep.F90 \
         mf_write.f90 evolve.F90 \
         checksetup.F90 initial.F90

OBJECTS1 = $(SOURCES:.f90=.o)
OBJECTS = $(OBJECTS1:.F90=.o)

ifeq ($(KROME), krome)
    .PHONY: all

    all: checksystem krome_setup krome phantom
    include MakeKrome
else
    .PHONY: phantom
endif

phantom: checksystem checkparams $(OBJECTS) phantom.o
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJECTS) phantom.o $(LDFLAGS)
ifeq ($(UNAME), Darwin)
	dsymutil $(BINDIR)/$@
endif


	@echo ""
	@echo "=============== CHEMISTRY ==============="
	@echo ""
ifeq ($(KROME), krome)
	@echo "krome coupling status = enabled"
else
	@echo "krome coupling status = disabled"
endif
	@echo ""
	@echo "========================================="
	@sh ../scripts/phantom_version_gen.sh "$(IDFLAGS)"
	@echo ""
	@echo "The Phantom is here (in $(BINDIR)/phantom)"
	@echo ""

#----------------------------------------------------
# generic target for compiling ALL phantom utilities
# this is used in the nightly build checks
#
utils: phantomsetup phantomanalysis \
       multirun phantom_moddump \
       phantom2divv phantom2divb combinedustdumps \
       diffdumps showheader showarrays ev2mdot phantomevcompare acc2ang \
       phantom2sphNG phantom2gadget testbinary \
       sfutils phantom2pdf-amr splitpart \
       phantom2struct libphantom

cleanutils: cleansetup cleananalysis \
            cleanmultirun cleantestbinary cleanmoddump \
            cleanphantom2divv cleanphantom2divb \
            cleandiffdumps cleanshowheader cleanshowarrays cleanev2mdot cleanacc2ang \
            cleanp2s cleanphantom2gadget \
            cleansfutils cleanphantom2pdf-amr cleansplitpart \
            cleanphantom2struct cleanphantomevcompare cleanlibphantom cleanphantomsinks

#--------------------------------------------------------------
# edit target opens current setup module in the default editor
#
edit: checksetup
	$(EDITOR) ../src/setup/$(SETUPFILE)

#----------------------------------------------------
# these are the sources for anything which uses the readwrite_dumps module
#
SRCDUMP= physcon.f90 ${CONFIG} ${SRCKERNEL} io.F90 units.f90 boundary.f90 mpi_utils.F90 \
         utils_timing.f90 utils_infiles.f90 dtype_kdtree.f90 utils_allocate.f90 part.F90 ${DOMAIN} \
         mpi_dens.F90 mpi_force.F90 \
         mpi_balance.F90 mpi_memory.F90 mpi_derivs.F90 kdtree.F90 linklist_kdtree.F90 \
         utils_dumpfiles.f90 utils_vectors.f90 utils_mathfunc.f90 \
         utils_datafiles.f90 utils_filenames.f90 utils_system.f90 utils_tables.f90 datafiles.f90 gitinfo.f90 \
         centreofmass.f90 \
         timestep.f90 ${SRCEOS} cullendehnen.f90 \
         ${SRCGR} ${SRCPOT} ${SRCPHOTO} \
         memory.F90 \
         utils_sphNG.f90 \
         setup_params.f90 ${SRCFASTMATH} checkoptions.F90 \
         viscosity.f90 options.f90 checkconserved.f90 prompting.f90 ${SRCDUST} \
         ${SRCREADWRITE_DUMPS} \
         utils_sort.f90 sort_particles.F90
OBJDUMP1= $(SRCDUMP:.f90=.o)
OBJDUMP= $(OBJDUMP1:.F90=.o)

#-------------------------------------------------------
# these are the sources for phantom libsetup
# must NOT contain .F90 files or pre-processing options
#
LIBSETUP=$(BINDIR)/libphantomsetup.a
SRCLIBSETUP=physcon.f90 geometry.f90 random.f90 utils_tables.f90 utils_vectors.f90 stretchmap.f90 \
            utils_binary.f90 set_binary.f90 set_flyby.f90 set_unifdis.f90 set_sphere.f90 set_shock.f90 \
            set_dust.f90 libsetup.f90
OBJLIBSETUP=${SRCLIBSETUP:.f90=.o}

libsetup: $(OBJLIBSETUP)
	@echo ""
	@echo "Phantom libsetup built"
	@echo ""
	ar rcs $(LIBSETUP) $^

#----------------------------------------------------
# these are the sources for phantom setup utility
#
.PHONY: phantomsetup
phantomsetup: setup

SRCSETUP= utils_omp.F90 utils_summary.F90 utils_gravwave.f90 \
          set_dust_options.f90 \
          utils_indtimesteps.F90 partinject.F90 \
          ${SRCTURB} ${SRCNIMHD} ${SRCCHEM} \
          ptmass.F90 energies.F90 density_profiles.f90 set_slab.f90 set_disc.F90 \
          set_cubic_core.f90 set_fixedentropycore.f90 set_softened_core.f90 set_star.f90 relax_star.f90 \
          set_vfield.f90 ptmass_radiation.F90 ${SRCINJECT} \
          ${SETUPFILE} checksetup.F90 \
          set_Bfield.f90 damping.f90 readwrite_infile.f90 ${SRCKROME}

OBJSETUP1= $(SRCSETUP:.f90=.o)
ifeq ($(KROME), krome)
    OBJSETUP1 += $(KROME_OBJS)
endif
OBJSETUP= $(OBJDUMP) $(OBJSETUP1:.F90=.o) phantomsetup.o

setup: checksystem checkparams libsetup libphantom $(OBJSETUP)
	$(FC) $(FFLAGS) -o $(BINDIR)/phantomsetup  $(OBJSETUP) $(LIBSETUP) $(LIBPHANTOM) $(LDFLAGS) $(LIBS)
	@echo ""
	@echo "Phantom setup built"
	@echo ""

cleansetup:
	rm -f $(BINDIR)/phantomsetup

#----------------------------------------------------
# phantom test suite
#
SRCTESTS=utils_testsuite.f90 ${TEST_FASTMATH} test_kernel.f90 \
         test_dust.F90 test_growth.F90 test_smol.F90 \
         test_nonidealmhd.F90 directsum.f90 test_gravity.F90 \
         test_derivs.F90 test_cooling.f90 test_eos.f90 test_externf.f90 test_rwdump.f90 \
         test_step.F90 test_indtstep.F90 set_disc.F90 test_setdisc.F90 \
         test_link.F90 test_kdtree.F90 test_part.f90 test_ptmass.f90 test_luminosity.F90\
         test_gnewton.f90 test_corotate.f90 test_geometry.f90 \
         test_sedov.F90 test_radiation.F90 \
         testsuite.F90 phantomtest.f90

ifeq (X$(SRCTEST), X)
SRCTEST=${SRCTESTS}
endif
# replace a few test routines for the GR code
ifeq ($(GR),yes)
    SRCTEST1 = $(SRCTESTS:test_externf.f90=test_externf_gr.f90)
    SRCTEST2 = $(SRCTEST1:test_gnewton.f90=)
    SRCTEST = $(SRCTEST2:test_corotate.f90=test_gr.f90)
endif
OBJTEST1= $(SRCTEST:.f90=.o)
OBJTEST= $(OBJTEST1:.F90=.o)

.PHONY: phantomtest

phantomtest: checksystem checkparams libphantom libsetup $(OBJTEST)
	@echo ""
	@echo "Phantomtest: Getting your life back one test at a time"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/phantomtest $(OBJTEST) $(LDFLAGS) $(LIBS) $(LIBPHANTOM) $(LIBSETUP)

cleantest:
	rm -f $(BINDIR)/phantomtest

#----------------------------------------------------
# update phantom version number
#
config.o: phantom-version.h

phantom-version.h:
	@echo "creating $@"
	@echo "#define PHANTOM_VERSION_MAJOR $(PHANTOM_VERSION_MAJOR)" > $@
	@echo "#define PHANTOM_VERSION_MINOR $(PHANTOM_VERSION_MINOR)" >> $@
	@echo "#define PHANTOM_VERSION_MICRO $(PHANTOM_VERSION_MICRO)" >> $@
	@echo "#define PHANTOM_VERSION_STRING \"$(VERSION)\"" >> $@

#----------------------------------------------------
# Smoluchowsky library
growth_smol.o: checksmol ../src/main/growth_smol.f90
	$(FC) -c $(FFLAGS) ${FPP_PREFIX} $(FPPFLAGS) ../src/main/growth_smol.f90 -I$(SMOL_DIR) -o $@

#----------------------------------------------------
# Probability Distribution Functions from fixed grid
# (produced via splash to grid)
# these are the sources for the grid2pdf utility
#
utils_outputhdf5.o: checkhdf5

write_grid_hdf5.o: checkhdf5

pdfs.o: checksplash $(SPLASH_DIR)/src/pdfs.f90
	$(FC) $(FFLAGS) -o $@ -c $(SPLASH_DIR)/src/pdfs.f90

# In case you need the old pdfs.f90 module located in phantom/src/utils/
# rather than the on located in splash. (e.g. analysis_MWpdf.f90 requires the
# phantom version)
phantom_pdfs.o: ../src/utils/pdfs.f90
	$(FC) $(FFLAGS) -o $@ -c $<

asciiutils.o: checksplash $(SPLASH_DIR)/src/asciiutils.f90
	$(FC) $(FFLAGS) -o $@ -c $(SPLASH_DIR)/src/asciiutils.f90

write_griddata.o: checksplash $(SPLASH_DIR)/src/write_griddata.F90
	$(FC) $(FFLAGS) -o $@ -c $(SPLASH_DIR)/src/write_griddata.F90

OBJG2PDF= io.o utils_filenames.o asciiutils.o write_griddata.o \
          hdf5utils.o read_grid_hdf5.o write_grid_hdf5.o io_grid.o pdfs.o rhomach.o grid2pdf.o

.PHONY: grid2pdf
grid2pdf: checksys checkparams checkhdf5 $(OBJG2PDF)
	@echo "objects are $(OBJG2PDF)"
	$(FC) $(FFLAGS) -o $(BINDIR)/grid2pdf $(OBJG2PDF) $(LDFLAGS) -L$(HDF5ROOT)/lib -lhdf5
	@echo ""
	@echo "Grid2pdf: we are Possibly Dangerously Fanatical"
	@echo ""

cleang2p:
	rm -f $(BINDIR)/grid2pdf

#------------------------------------------------------
# Probability Distribution Functions via adaptive mesh
#
phantom2pdf: phantom2pdf-amr

.PHONY: phantom2pdf-amr
phantom2pdf-amr:
	${MAKE} phantomanalysis ANALYSIS="adaptivemesh.f90 interpolate3D_amr.F90 asciiutils.f90 pdfs.f90 analysis_pdfs.f90"\
        ANALYSISBIN=$@ ANALYSISONLY=yes

cleanphantom2pdf-amr:
	rm -f $(BINDIR)/phantom2struct

analysis_pdfs.o: interpolate3D_amr.o adaptivemesh.o
interpolate3D_amr.o: adaptivemesh.o

#----------------------------------------------------
# these are the sources for the phantom_moddump utility
#
ifndef MODDUMPBIN
MODDUMPBIN=phantommoddump
endif
OBJMOD1 = prompting.o set_Bfield.o density_profiles.o ${MODFILE:.f90=.o} phantom_moddump.o
OBJMOD = ${OBJMOD1:.F90=.o}

phantom_moddump: checksystem checkparams libphantom libsetup $(OBJMOD)
	@echo ""
	@echo "phantom_moddump: we are here to help you"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$(MODDUMPBIN) $(OBJMOD) $(LDFLAGS) $(LIBS) $(LIBPHANTOM) $(LIBSETUP)

moddump: phantom_moddump

cleanmoddump:
	rm -f $(BINDIR)/phantommoddump

# files from MCFOST used by phantommoddump
mess_up_SPH.o: checkmcfost $(MCFOST_DIR)/src/mess_up_SPH.f90
	$(FC) $(FFLAGS) -o $@ -c $(MCFOST_DIR)/src/mess_up_SPH.f90

#----------------------------------------------------
# these are the sources for the phantomanalysis utility
#
OBJAN1= ${ANALYSIS:.f90=.o}
OBJAN= ${OBJAN1:.F90=.o}
OBJA= leastsquares.o solvelinearsystem.o prompting.o \
      utils_disc.o set_dust.o utils_binary.o set_binary.o ${OBJAN} phantomanalysis.o

ifndef ANALYSISBIN
ANALYSISBIN=phantomanalysis
endif

.PHONY: phantomanalysis
phantomanalysis: checksystem checkparams libphantom $(OBJA)
	@echo ""
	@echo "phantomanalysis: we live to serve you"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$(ANALYSISBIN) $(OBJA) $(LDFLAGS) $(LIBS) $(LIBPHANTOM)

analysis: phantomanalysis


cleananalysis:
	rm -f $(BINDIR)/phantomanalysis

.PHONY: libphantom
# extra files used in libphantom that are not in OBJECTS
SRCLIB= # libphantom.F90
OBJLIB1=${SRCLIB:.f90=.o}
OBJLIB=${OBJLIB1:.F90=.o}

LIBPHANTOM=$(BINDIR)/libphantom.a

libphantom.a: checksystem checkparams $(OBJECTS) $(OBJLIB)
	ar rc $(LIBPHANTOM) $(OBJLIB) $(OBJECTS)

libphantom.so: checksystem checkparams phantom ${OBJLIB}
	$(FC) -shared $(FFLAGS) $(FPPFLAGS) $(DBLFLAG) ${OBJLIB} ${OBJECTS} $(LDFLAGS) -o $(BINDIR)/libphantom.so

libphantom: libphantom.a

cleanlibphantom:
	rm -f $(BINDIR)/libphantom.so $(LIBPHANTOM)

.PHONY: pyanalysis
pyanalysis: libphantom.so

#------------------------------------------------------
# Various utilities for computing structure functions
# and manipulating the resulting output
#
.PHONY: phantom2struct
phantom2struct:
	${MAKE} phantomanalysis ANALYSIS="utils_timing.f90 io_structurefn.f90 random.f90 struct_part.f90 analysis_structurefn.f90"\
        ANALYSISBIN=$@ ANALYSISONLY=yes

cleanphantom2struct:
	rm -f $(BINDIR)/phantom2struct

# conversion between structure function file formats
.PHONY: struct2struct
STRUCT2STRUCTOBJ= utils_filenames.o io_structurefn.o struct2struct.o
struct2struct: checksys checkparams ${STRUCT2STRUCTOBJ}
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ ${STRUCT2STRUCTOBJ}

cleanstruct2struct:
	rm -f $(BINDIR)/struct2struct

# structure function slope calculation
.PHONY: get_struct_slope get_struct_slope
GETSLOPESFOBJ=utils_filenames.o io_structurefn.o leastsquares.o get_struct_slope.o
get_slope_sf: get_struct_slope
get_struct_slope: checksys checkparams ${GETSLOPESFOBJ}
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ ${GETSLOPESFOBJ}

cleanget_struct_slope:
	rm -f $(BINDIR)/get_struct_slope

sfutils: structutils
structutils: struct2struct get_slope_sf

cleansfutils: cleanstructutils
cleanstructutils: cleanstruct2struct cleanget_struct_slope

#------------------------------------------------------
# particle splitting utility (this is a moddump compiled as a standalone utility)
#
.PHONY: splitpart
splitpart:
	${MAKE} moddump MODFILE="utils_getneighbours.F90 utils_splitmerge.f90 splitpart.f90 moddump_splitpart.f90"\
        MODDUMPBIN=$@

cleansplitpart:
	rm -f $(BINDIR)/splitpart

#------------------------------------------------------
# particle merging utility (this is a moddump compiled as a standalone utility)
#
.PHONY: mergepart
mergepart:
	${MAKE} moddump MODFILE="utils_getneighbours.F90 utils_splitmerge.f90 splitpart.f90 moddump_mergepart.f90"\
        MODDUMPBIN=$@

cleanmergepart:
	rm -f $(BINDIR)/mergepart

#----------------------------------------------------
# utility to calculate divv from a dump file
# compile using all phantom files
#
phantom2divv: checksys checkparams $(OBJECTS) phantom2divv.o
	@echo ""
	@echo "phantom2divv: divergence is beautiful"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJECTS) phantom2divv.o

cleanphantom2divv:
	rm -f $(BINDIR)/phantom2divv

#----------------------------------------------------
# utility to calculate divB & curlB from a dump file
# compile using all phantom files
#
phantom2divb: checksys checkparams $(OBJECTS) phantom2divb.o
	@echo ""
	@echo "phantom2divb: divergence should be eradicated"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJECTS) phantom2divb.o

cleanphantom2divb:
	rm -f $(BINDIR)/phantom2divb

#----------------------------------------------------
# these are the sources for the diffdumps utility
#
diffdumps: checksys checkparams $(OBJDUMP) utils_testsuite.o diffdumps.o
	@echo ""
	@echo "diffdumps: we welcome you"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJDUMP) utils_testsuite.o diffdumps.o $(LDFLAGS)

cleandiffdumps:
	rm -f $(BINDIR)/phantom2divb

#----------------------------------------------------
# these are the sources for the phantom2sphNG utility
#
phantom2sphNG: checksystem checkparams $(OBJDUMP) phantom2sphNG.o
	@echo ""
	@echo "phantom2sphNG: now why would you want to do that?"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJDUMP) phantom2sphNG.o

p2s: phantom2sphNG

cleanp2s:
	rm -f $(BINDIR)/phantom2sphNG

#----------------------------------------------------
# these are the sources for the phantom2sphNG utility
#
phantom2gadget: checksystem checkparams $(OBJDUMP) phantom2gadget.o
	@echo ""
	@echo "phantom2gadget: now why would you want to do that?"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJDUMP) phantom2gadget.o

p2g: phantom2gadget

cleanphantom2gadget:
	rm -f $(BINDIR)/phantom2gadget

#----------------------------------------------------
# these are the sources for the phantom2mcfost utility
#
.PHONY: phantom2mcfost
phantom2mcfost: checkmcfost
	${MAKE} phantomanalysis ANALYSIS="analysis_mcfost.f90"\
        ANALYSISBIN=$@ ANALYSISONLY=yes LDFLAGS="-L$(MCFOST_DIR)/src -lmcfost $(LIBCXX)"

analysis_mcfost.o: analysis_mcfost.f90
	$(FC) -c $(FFLAGS) -I$(MCFOST_DIR)/src $< -o $@

analysis_mcfost.o: checkmcfost

cleanphantom2mcfost:
	rm -f $(BINDIR)/phantom2mcfost

#----------------------------------------------------
# these are the sources for the phantom2hdf5 utility
#
checkhdf5flags:
   ifneq ($(HDF5), yes)
	@echo "-----------------------------"
	@echo "Need to compile with HDF5=yes"
	@echo "-----------------------------"
	${MAKE} err
   endif

.PHONY: phantom2hdf5
phantom2hdf5: checksystem checkparams checkhdf5flags $(OBJDUMP) readwrite_dumps.o phantom2hdf5.o
	@echo ""
	@echo "phantom2hdf5: welcome to the future"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJDUMP)  phantom2hdf5.o $(LDFLAGS)

p2h: phantom2hdf5

cleanphantom2hdf5:
	rm -f $(BINDIR)/phantom2hdf5

#----------------------------------------------------
# utility to rewrite .ev files using a common header
#
SRCEV=utils_infiles.f90 utils_evfiles.f90 prompting.f90 phantomevcompare.f90
OBJEVC1 = ${SRCEV:.f90=.o}
OBJEVC = ${OBJEVC1:.F90=.o}

.PHONY: phantomevcompare
phantomevcompare: $(OBJEVC)
	@echo ""
	@echo "phantomevcompare: let the graphing begin!"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJEVC)

cleanphantomevcompare:
	rm -f $(BINDIR)/phantomevcompare

#----------------------------------------------------
# utility to extract sink information from a master file
#
.PHONY: phantomsinks
phantomsinks: phantomextractsinks.o
	@echo ""
	@echo "phantomsinks: extracting sink information!"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ phantomextractsinks.o

cleanphantomsinks:
	rm -f $(BINDIR)/phantomextractsinks

#----------------------------------------------------
# these are the sources for the multirun utility
#
SRCMULT = physcon.f90 ${CONFIG} ${SRCKERNEL} io.F90 mpi_utils.F90 utils_timing.f90 dtype_kdtree.F90 ${SRCFASTMATH} \
          units.f90 boundary.f90 utils_allocate.f90 part.F90 timestep.f90 setup_params.f90 \
          utils_filenames.f90 utils_mathfunc.f90 utils_system.f90 utils_vectors.f90 \
          utils_omp.F90 utils_datafiles.f90 datafiles.f90 utils_tables.f90 utils_infiles.f90 \
          ${SRCEOS} viscosity.f90 options.f90 damping.f90 utils_gravwave.f90 \
          utils_dumpfiles.f90 utils_summary.f90 centreofmass.f90 \
          ${SRCCHEM} ${DOMAIN} ${SRCPOT} dust_formation.F90 ptmass_radiation.F90 mpi_balance.F90 mpi_memory.F90 mpi_force.F90 mpi_dens.F90 mpi_derivs.F90 kdtree.F90 linklist_kdtree.F90 ptmass.F90 ${SRCTURB} \
          checkconserved.f90 prompting.f90 ${SRCDUST} ${SRCNIMHD} readwrite_infile.f90 ${MULTIRUNFILE}
OBJM1 = ${SRCMULT:.f90=.o}
OBJMULT = ${OBJM1:.F90=.o}

multirun: checksystem checkparams $(OBJMULT)
	@echo ""
	@echo "multirun: your hope is our desire"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJMULT)

cleanmultirun:
	rm -f $(BINDIR)/multirun

#----------------------------------------------------
# utility to plot orbits based on orbital elements (a,e,i,o,w,f)
#
SRCBIN = prompting.f90 utils_datafiles.f90 datafiles.f90 ${CONFIG} physcon.f90 io.F90 \
         mpi_utils.F90 utils_allocate.f90 utils_binary.f90 set_binary.f90 test_binary.f90 testbinary.f90
OBJBIN1 = ${SRCBIN:.f90=.o}
OBJBIN = ${OBJBIN1:.F90=.o}

.PHONY: testbinary

testbin: testbinary

testbinary: checksys checkparams $(OBJBIN)
	@echo ""
	@echo "test_binary: may your orbits orbit"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/testbinary $(OBJBIN)

cleantestbinary:
	rm -f $(BINDIR)/testbinary

#----------------------------------------------------
# check for anything that depends on HDF5
#
checkhdf5:
   ifeq (X${HDF5ROOT}, X)
	@echo; echo "ERROR: HDF5ROOT should be set before compiling with HDF5 utilities"; echo; ${MAKE} err;
   else
	@if [ -d $$HDF5ROOT ]; then echo; echo "HDF5ROOT=$$HDF5ROOT"; echo; else echo; echo "ERROR: Directory given by HDF5ROOT=$$HDF5ROOT does not exist"; echo; ${MAKE} err; fi;
   endif

#----------------------------------------------------
# these are the sources for the plot_kernel utility
#

OBJPLOTK= physcon.o ${SRCKERNEL:.f90=.o} giza-fortran.o plot_kernel.o

plotkernel: checksys checkparams checksplash $(OBJPLOTK)
	@echo ""
	@echo "plot_kernel: may your kernels be normalised"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJPLOTK) $(LDFLAGS) -L$(SPLASH_DIR)/giza/lib -lgiza

plot_kernel.o: ${SRCKERNEL:.f90=.o}
#giza-fortran.o: ${SPLASH_DIR}/giza/src/$@
#	$(FC) $(FFLAGS) -o $@ -c ${SPLASH_DIR}/giza/interface/giza-fortran.F90

cleanplotkernel:
	rm -f $(BINDIR)/plotkernel

#----------------------------------------------------
# these are the sources for the showheader utility
#
SRCSHOWHEADER= utils_dumpfiles.f90 showheader.f90
OBJSHOWHEADER= $(SRCSHOWHEADER:.f90=.o)
showheader: checksys $(OBJSHOWHEADER)
	@echo ""
	@echo "showheader: show me the header!"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJSHOWHEADER)

cleanshowheader:
	rm -f $(BINDIR)/showheader

#----------------------------------------------------
# these are the sources for the showarrays utility
#
SRCSHOWARRAYS= utils_dumpfiles.f90 showarrays.f90
OBJSHOWARRAYS= $(SRCSHOWARRAYS:.f90=.o)
showarrays: checksys $(OBJSHOWARRAYS)
	@echo ""
	@echo "showarrays: show me the arrays!"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJSHOWARRAYS)

cleanshowarrays:
	rm -f $(BINDIR)/showheader

#----------------------------------------------------
# these are the sources for the evol_dustywaves utility
#
SRCDUSTEVOL= cubicsolve.f90 dustywaves.f90 evol_dustywaves.f90
OBJDUSTEVOL= $(SRCDUSTEVOL:.f90=.o)

evol_dustywaves: checksys $(OBJDUSTEVOL)
	@echo ""
	@echo "dusty wave .ev solutions^TM: All the energy you need."
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJDUSTEVOL)

#----------------------------------------------------
# these are the sources for the ev2mdot utility
#
.PHONY: ev2mdot
ev2mdot: checksys utils_filenames.o utils_infiles.o utils_evfiles.o ev2mdot.o
	@echo ""
	@echo "ev2mdot: Accretion rates R us."
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ ev2mdot.o utils_filenames.o utils_evfiles.o utils_infiles.o

cleanev2mdot:
	rm -f $(BINDIR)/ev2mdot
#----------------------------------------------------
# these are the sources for the ev2kdot utility
#
.PHONY: ev2kdot
ev2kdot: checksys utils_filenames.o utils_infiles.o utils_evfiles.o ev2kdot.o
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ ev2kdot.o utils_filenames.o utils_evfiles.o utils_infiles.o

cleanev2kdot:
	rm -f $(BINDIR)/ev2kdot

#----------------------------------------------------
# these are the sources for the acc2ang utility
#
.PHONY: acc2ang
acc2ang: checksys acc2ang.o
	@echo ""
	@echo "acc2ang: Accreted ang. mom. R us."
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ acc2ang.o

cleanacc2ang:
	rm -f $(BINDIR)/acc2ang

#---------------------------
# sources for the mass_flow utility
#
OBJMF1 = ${ANALYSIS:.f90=.o}
OBJMF2 = ${OBJMF1:.F90=.o}
OBJMF = ${OBJMF2:.f=.o}
OBJM= utils_sort.o leastsquares.o solvelinearsystem.o ${OBJDUMP} ${OBJMF} utils_binary.o set_binary.o mf_write.o

.PHONY: mflow
mflow: checksys $(OBJM)  mflow.o ev2mdot lombperiod
	@echo ""
	@echo "mflow: mass flow R us."
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@  $(OBJM) mflow.o

.PHONY:lombperiod
lombperiod: powerspectrums.o lombperiod.o
	$(FC) $(FFLAGS) -o $(BINDIR)/$@  lombperiod.o powerspectrums.o

#----------------------------------------------------
# these are the sources for the combinedustdumps utility
#
OBJCDD= ${OBJECTS} combinedustdumps.o

combinedustdumps: checksys checkparams $(OBJCDD)
	@echo ""
	@echo "combinedustdumps: many grains make light work"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJCDD) $(LDFLAGS)

cleancombinedustdumps:
	rm -f $(BINDIR)/combinedustdumps


include Makefile_qscripts


#----------------------------------------------------
# unit test for block limits
#
test1: checksystem checkparams $(OBJDUMP) test_blocklimits.o
	$(FC) $(FFLAGS) -o $(BINDIR)/test1 $(OBJDUMP) test_blocklimits.o

#----------------------------------------------------
# run test suite
#
.PHONY: test test2 testcyl testgrav testall
test:
	${MAKE} SETUP=test phantomtest && $(RUNMPI) $(BINDIR)/phantomtest

test2:
	${MAKE} SETUP=test2 phantomtest && $(RUNMPI) $(BINDIR)/phantomtest

testkd:
	${MAKE} SETUP=testkd phantomtest && $(RUNMPI) $(BINDIR)/phantomtest

testcyl:
	${MAKE} SETUP=testcyl phantomtest && $(RUNMPI) $(BINDIR)/phantomtest

testgrav:
	${MAKE} SETUP=testgrav phantomtest && $(RUNMPI) $(BINDIR)/phantomtest gravity

testdust:
	${MAKE} SETUP=testdust phantomtest && $(RUNMPI) $(BINDIR)/phantomtest dust

testgr:
	${MAKE} SETUP=testgr phantomtest && $(RUNMPI) $(BINDIR)/phantomtest gr

testgrowth:
	${MAKE} SETUP=testgrowth phantomtest && $(RUNMPI) $(BINDIR)/phantomtest growth

testnimhd:
	${MAKE} SETUP=testnimhd phantomtest && $(RUNMPI) $(BINDIR)/phantomtest nimhd

testall: test test2 testcyl testgrav

include Makefile_fastmath

LASTSYSTEM = ${shell if [ -e .make_lastsystem ]; then cat .make_lastsystem; fi}
LASTSETUP = ${shell if [ -e .make_lastsetup ]; then cat .make_lastsetup; fi}
LASTFPPFLAGS = ${shell if [ -e .make_lastfppflags ]; then cat .make_lastfppflags; fi}
LASTFFLAGS = ${shell if [ -e .make_lastfflags ]; then cat .make_lastfflags; fi}

include Makefile_checks

giza-fortran.o : $(SPLASH_DIR)/giza/interface/giza-fortran.F90 $(SPLASH_DIR)/giza/lib/libgiza.a
	$(FC) $(FFLAGS) -I$(SPLASH_DIR)/giza/include/ -c $< -o $@

compilers:
	@echo "I suggest one of the following, based on detected Fortran compilers..."; echo;
	@if type -p ifort > /dev/null; then echo "make SYSTEM=ifort"; fi;
	@if type -p pathf90 > /dev/null; then echo "make SYSTEM=pathf90"; fi;
	@if type -p pgf90 > /dev/null; then echo "make SYSTEM=pgf90"; fi;
	@if type -p xlf90_r > /dev/null; then echo "make SYSTEM=ukaff1a [uses xlf90_r]"; fi;
	@if type -p gfortran > /dev/null; then echo "make SYSTEM=gfortran"; fi;
	@if type -p g95 > /dev/null; then echo "make SYSTEM=g95"; fi;
	@echo "(end of possible selections)"; echo;

#----------------------------------------------------
# target to automatically include dependencies in Makefile
# relies on the g95 compiler being present
# (does not have to be used for the main compilation)

depends: clean checksetup
	#@echo '*********************************************************************************'
	#@echo 'First run of Makefile -- creating dependency lines using gfortran, writing to .depends'
	#@echo '*********************************************************************************'
	#@gfortran -M -cpp -c ../src/*/*.*90 > .depends
	#@echo '*************************************************************************'
	#@echo 'If no errors above, then Makefile dependencies were created successfully '
	#@echo ' -- be sure to run "make depends" again if you alter code dependencies'
	#@echo '*************************************************************************'
	#@${MAKE} clean

.depends:
	@if type -p gfortran; then touch .depends; ${MAKE} --quiet SETUP=test depends; else echo "warning: no gfortran so dependencies not calculated"; touch .depends; fi;

include .depends

getdims:
	@echo $(MAXP)

get_setup_opts:
	@echo "${GR:yes=GR} ${METRIC} ${MHD:yes=MHD} ${NONIDEALMHD:yes=non-ideal} ${DUST:yes=dust} ${GRAVITY:yes=self-gravity} ${RADIATION:yes=radiation}  ${H2CHEM:yes=H2_Chemistry} ${DISC_VISCOSITY:yes=disc_viscosity} ${ISOTHERMAL:yes=isothermal} ${PERIODIC:yes=periodic}" | xargs | sed -e 's/ /, /g' -e 's/_/ /g'

get_setup_file:
	@echo "$(SETUPFILE)"

err:
	$(error aborting);

clean: cleanlibphantom
	rm -f *.o *.mod phantom-version.h
	cd $(BINDIR); rm -f phantom phantomsetup
ifeq ($(KROME),krome)
	rm -f ${KROMEPATH}/build/*.o ${KROMEPATH}/build/*.mod
endif

cleanall: clean cleanmathflags
	cd $(BINDIR); rm -f phantom phantomsetup

cleandist: clean cleanall
	rm -f .make_lastsystem .make_lastsetup .make_lastfppflags .depends

cleanmathflags:
	rm -f .make_mathflags bin/getmathflags
