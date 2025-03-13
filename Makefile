# Makefile for RemoTeC

############################################
# Definitions
############################################

# Choose compilation options
############################
# Do not choose parallel=yes for creating synthetic spectra (because of random noise)
# parallel = yes / no
parallel = no
# debug = yes / no
debug = no
# FORT = IUP-gfortran, DKRZ-gfortran, HELIX-gfortran
FORT = HELIX-gfortran
############################

# Directories for source files, object files, modules and dependency files
srcdir = SRC
coredir = $(srcdir)/core
ntdir = $(coredir)/NUM_TOOLS
gsddir = $(coredir)/GAUSDL
lntrndir = $(coredir)/LINTRAN
createdir = $(srcdir)/create
retrievedir = $(srcdir)/retrieve
objdir = ./OBJECTS
moddir = ./MODULES
depdir = ./DEPEND
$(shell [ -d $(objdir) ] || mkdir -p $(objdir))
$(shell [ -d $(moddir) ] || mkdir -p $(moddir))
$(shell [ -d $(depdir) ] || mkdir -p $(depdir))

# Compiler flags
ifeq ($(FORT),IUP-gfortran)
   make = gmake
   fortlib = /home/atmo/software/lib_gnu/
   nc-fortrandir = $(fortlib)/netcdf-fortran-4.4.4/build/
   nc-cdir = $(fortlib)/netcdf-4.6.1/build/
   hdfdir = $(fortlib)/hdf5-1.10.2/build/
   zdir = $(fortlib)/zlib-1.2.11/

   FC = gfortran

   FFLAGS =  -c -m64 -J$(moddir) -I$(moddir)
   FFLAGS += -fimplicit-none -ffree-line-length-none -fdefault-real-8 -fdefault-double-8
   FFLAGS += -I$(hdfdir)/include -I$(nc-fortrandir)/include -I$(nc-cdir)/include -I$(zdir)/include

   LDFLAGS =  -m64 -L$(nc-fortrandir)/lib -lnetcdff -L$(nc-cdir)/lib/ -lnetcdf
   LDFLAGS += -L$(hdfdir)/lib -lhdf5_fortran -lhdf5_hl -lhdf5
   LDFLAGS += -L$(zdir)lib -lz
   LDFLAGS += -ldl -lnsl

   ifeq ($(debug),yes)
      FFLAGS += -g -traceback -debug full -debug-parameters all
      FFLAGS += -fpe0 -check bounds -check pointers -check uninit
      FFLAGS += -check overflow -fp-stack-check -warn unused
   else
      FFLAGS += -O3
   endif
   ifeq ($(parallel),yes)
      FFLAGS +=  -xopenmp
      LDFLAGS +=  -xopenmp
   endif
   FCDEP = $(FC) -m64 -cpp -M -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -J$(depdir) -I$(depdir) -I$(hdfdir)/include -I$(nc-fortrandir)/include -I$(nc-cdir)/include

else ifeq ($(FORT),DKRZ-gfortran)

   make = gmake
   #fortlib = /pf/b/b309168/software/lib_gnu/
   #fortlib = /sw/rhel6-x64/gcc/gcc-7.1.0/
   FC = gfortran
   FFLAGS =  -c -J$(moddir) -I$(moddir)
   FFLAGS += -fimplicit-none -ffree-line-length-none -fdefault-real-8 -fdefault-double-8
   FFLAGS += -I/pf/b/b309168/software/lib/install/include
   LDFLAGS =
   LDFLAGS = -L/pf/b/b309168/software/lib/install/lib -lnetcdff -lnetcdf -lnetcdf -ldl -lm -lhdf5_fortran -lhdf5_hl -lhdf5 -lm -ldl -lz -lcurl


   #FFLAGS += -I$(fortlib)/hdf5-1.8.7/include -I$(fortlib)/netcdf-4.1.3/include
   #LDFLAGS =  -m64 -lcurl -L$(fortlib)/netcdf-4.1.3/lib -lnetcdff -lnetcdf
   #LDFLAGS += -L$(fortlib)/hdf5-1.8.7/lib -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lnsl
   ifeq ($(debug),yes)
#       FFLAGS += -g -traceback -debug full -debug-parameters all
#       FFLAGS += -fpe0 -check bounds -check pointers -check uninit
#       FFLAGS += -check overflow -fp-stack-check -warn unused
      FFLAGS += -g -fbacktrace
      FFLAGS += -fcheck=bounds  -fcheck=pointer -fcheck=do
      FFLAGS += -ffpe-trap=invalid,zero,overflow
      FFLAGS += -Wuninitialized -Wunused -Wunused-dummy-argument -Warray-temporaries
   else
      FFLAGS += -O3
   endif
   ifeq ($(parallel),yes)
      FFLAGS +=  -xopenmp
      LDFLAGS +=  -xopenmp
   endif
   FCDEP = $(FC) -cpp -M -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -J$(depdir) -I$(depdir) -I/pf/b/b309168/software/lib/install/include

else ifeq ($(FORT),HELIX-gfortran)
    make = gmake
    fortlib = /home/hd/hd_hd/hd_gb423/.local/lib64
    nc-fortrandir = $(fortlib)/netcdf-fortran-4.6.1/build/
    nc-cdir = $(fortlib)/netcdf-c-4.9.2/build/
    hdfdir = $(fortlib)/hdf5-1.14.1-2/build/
    zdir = $(fortlib)/zlib-1.2.13/build/

    FC = gfortran

    FFLAGS = -c -J$(moddir) -I$(moddir)
    FFLAGS += -fimplicit-none -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -fallow-invalid-boz -fallow-argument-mismatch
    FFLAGS += -I$(hdfdir)/include -I$(nc-fortrandir)/include -I$(nc-cdir)/include -I$(zdir)/include

    LDFLAGS = -L$(nc-fortrandir)/lib -lnetcdff -L$(nc-cdir)/lib/ -lnetcdf
    LDFLAGS += -L$(hdfdir)/lib -lhdf5_fortran -lhdf5_hl -lhdf5
    LDFLAGS += -L$(zdir)lib -lz -ldl -lcurl

    ifeq ($(debug),yes)
        FFLAGS += -g -traceback -debug full -debug-parameters all
        FFLAGS += -fpe0 -check bounds -check pointers -check uninit
        FFLAGS += -check overflow -fp-stack-check -warn unused
    else
        FFLAGS += -O3
    endif

    ifeq ($(parallel),yes)
        FFLAGS += -xopenmp
        LDFLAGS += -xopenmp
    endif
endif

# Generate dependency files using preprocessor
# add directory of object file
# all prereqs listed will also become command-less, prereq-less targets,
# so that if a prerequisite file is missing/renamed, it is treated as changed
# (so the target will be rebuilt), but not cause an error.
# remove target and continuation backslashes, list one word per line,
# strip leading spaces, add trailing colons
makedep = $(FCDEP) $< > $(depdir)/$*.d;\
mv -f $(depdir)/$*.d $(depdir)/$*.d.tmp; \
sed -e 's|.*:|$@:|' < $(depdir)/$*.d.tmp > $(depdir)/$*.d; \
cp -f $(depdir)/$*.d $(depdir)/$*.d.tmp;\
sed -e 's/.*://' -e 's/\\$$//' < $(depdir)/$*.d.tmp | fmt -1 | \
sed -e 's/^ *//' -e 's/$$/:/' >> $(depdir)/$*.d;\
rm -f $(depdir)/$*.d.tmp;\
rm -f $*.i90

# object files
obj_core =	$(objdir)/header.o\
		$(objdir)/auxiliary_routines.o\
		$(objdir)/read_xsdb.o\
		$(objdir)/optic_molec.o\
		$(objdir)/OpticM_module.o \
		$(objdir)/optic_cirrus.o\
		$(objdir)/aerosol_input_module.o\
		$(objdir)/read_settings.o\
		$(objdir)/depol.o\
		$(objdir)/optic_input_module.o\
		$(objdir)/tau_grid.o\
		$(objdir)/header_gsd.o\
		$(objdir)/perturbation.o \
		$(objdir)/viewing.o \
		$(objdir)/single_scat.o \
		$(objdir)/gautool_acs.o \
		$(objdir)/gaus_mat.o \
		$(objdir)/gaus_pol.o \
		$(objdir)/rad_intf.o \
		$(objdir)/ocean_fresnel_module.o\
		$(objdir)/lintran_interface.o\
		$(objdir)/rad_trans_intf.o \
		$(objdir)/spectral_response.o\
		$(objdir)/forward_model.o\
		$(objdir)/forward_model_noscat.o\
		$(objdir)/pt_inversion.o\
		$(objdir)/adhoc_inversion.o\
		$(objdir)/profile_inversion.o\
		$(objdir)/spectrum_internal.o\
		$(objdir)/atmosphere_internal.o\
		$(objdir)/read_miprep.o\
		$(objdir)/retrieval.o

obj_lntrn = 	$(objdir)/lintran_header.o         	     \
		$(objdir)/lintran_constants_module.o         \
		$(objdir)/lintran_types_module.o             \
		$(objdir)/lintran_grid_module.o              \
		$(objdir)/lintran_geometry_module.o          \
		$(objdir)/lintran_pixel_module.o             \
		$(objdir)/lintran_tlc_module.o               \
		$(objdir)/lintran_fields_module.o            \
		$(objdir)/lintran_matrix_module.o            \
		$(objdir)/lintran_interface_module.o         \
		$(objdir)/lintran_nst1_module.o              \
		$(objdir)/lintran_nst3_module.o              \
		$(objdir)/lintran_nst4_module.o              \
		$(objdir)/lintran_remotec_module.o           \
		$(objdir)/lintran_module.o

obj_create =	$(objdir)/spectrum_interface_create.o\
		$(objdir)/atmosphere_interface_create.o\
		$(objdir)/solar_model_create.o\
		$(objdir)/noise.o\
                $(objdir)/read_synsettings.o\
		$(objdir)/calculate_syn_spectrum.o\
		$(objdir)/main_create.o

obj_retrieve =	$(objdir)/read_errors.o\
		$(objdir)/spectrum_interface_retrieve.o\
		$(objdir)/atmosphere_interface_retrieve.o\
		$(objdir)/solar_model_retrieve.o\
		$(objdir)/synthetic_input.o\
		$(objdir)/diagnostics_retrieve.o\
		$(objdir)/wrapper_retrieve.o\
		$(objdir)/main_retrieve.o

libraries =  $(ntdir)/libNts.a

num_tools =	$(ntdir)/drealft.f90\
		$(ntdir)/dpythag.f90\
		$(ntdir)/dfour1.f90\
		$(ntdir)/dsvdcmp.f90\
		$(ntdir)/spline.f90\
		$(ntdir)/linterp.f90\
		$(ntdir)/ludcmp.f90\
		$(ntdir)/lubksb.f90\
		$(ntdir)/sort.f90\
		$(ntdir)/qpythag.f90\
		$(ntdir)/qsvdcmp.f90\
		$(ntdir)/dtwofft.f90\
		$(ntdir)/indexx.f90


############################################
# Rules
############################################

# Goal
all: create retrieve

# link
create: $(obj_lntrn) $(obj_core) $(obj_create) $(libraries)
	$(FC) $(obj_lntrn) $(obj_core) $(obj_create) $(LDFLAGS) $(libraries) -o RemoTeC_create

# link
retrieve: $(obj_lntrn) $(obj_core) $(obj_retrieve) $(libraries)
	$(FC) $(obj_lntrn) $(obj_core) $(obj_retrieve) $(LDFLAGS) $(libraries) -o RemoTeC_retrieve

# compile and generate dependency files:
$(objdir)/%.o: $(coredir)/%.f90
	$(makedep)
	$(FC) $(FFLAGS) -o $@ $<

$(objdir)/%.o: $(createdir)/%.f90
	$(makedep)
	$(FC) $(FFLAGS) -o $@ $<

$(objdir)/%.o: $(retrievedir)/%.f90
	$(makedep)
	$(FC) $(FFLAGS) -o $@ $<

$(objdir)/%.o: $(lntrndir)/%.f90
	$(makedep)
	$(FC) $(FFLAGS) -o $@ $<

$(objdir)/%.o: $(gsddir)/%.f90
	$(makedep)
	$(FC) $(FFLAGS) -o $@ $<

# Include dependency files
-include $(subst $(objdir), $(depdir), $(obj_core:.o=.d))
-include $(subst $(objdir), $(depdir), $(obj_create:.o=.d))
-include $(subst $(objdir), $(depdir), $(obj_retrieve:.o=.d))
-include $(subst $(objdir), $(depdir), $(obj_lntrn:.o=.d))

# make library
$(ntdir)/libNts.a: $(num_tools)
	cd $(ntdir); $(make) FC="$(FC)" FFLAGS="$(FFLAGS)" -f Make_Lib

.PHONY : clean
clean:
	-rm $(depdir)/*.d
	-rm $(depdir)/*.d.tmp
	-rm $(depdir)/*.mod
	-rm $(objdir)/*.o
	-rm $(moddir)/*.mod
	-rm $(ntdir)/libNts.a
	-rm -r $(depdir)
	-rm -r $(objdir)
	-rm -r $(moddir)
	-rm RemoTeC_create
	-rm RemoTeC_retrieve
