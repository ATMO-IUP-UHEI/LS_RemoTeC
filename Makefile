# ===========================
#  Makefile for RemoTeC
# ===========================

# Check if correct conda environment is used
ifneq ($(notdir $(CONDA_PREFIX)),remotec-env)
    $(error "Run `conda activate remotec-env`.")
endif

# ===========================
#  Compilation Options
# ===========================

# parallel = yes / no
# Do not choose parallel=yes for creating synthetic spectra (because of random noise)
parallel := no

# debug = yes / no
debug    := no

# ===========================
#  Directories
# ===========================

srcdir      := SRC
coredir     := $(srcdir)/core
ntdir       := $(coredir)/NUM_TOOLS
gsddir      := $(coredir)/GAUSDL
lntrndir    := $(coredir)/LINTRAN
createdir   := $(srcdir)/create
retrievedir := $(srcdir)/retrieve
objdir      := ./OBJECTS
moddir      := ./MODULES
depdir      := ./DEPEND

# Ensure necessary directories exist
$(shell mkdir -p $(objdir) $(moddir) $(depdir))

# ===========================
#  Compiler and Flags
# ===========================

FC   := $(CONDA_PREFIX)/bin/gfortran
make := $(CONDA_PREFIX)/bin/make

# Compiler Flags
FFLAGS := -c -J$(moddir) -I$(moddir)
FFLAGS += -fimplicit-none -ffree-line-length-none -fdefault-real-8 -fdefault-double-8
FFLAGS += -fallow-invalid-boz -fallow-argument-mismatch
FFLAGS += -I$(CONDA_PREFIX)/include

# Linker Flags
LDFLAGS := -L$(CONDA_PREFIX)/lib -Wl,-rpath,$(CONDA_PREFIX)/lib
LDFLAGS += -lnetcdff -lnetcdf -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -ldl -lstdc++

# Preprocessor Dependency Generator
FCDEP := $(FC) -cpp -M -ffree-line-length-none -fdefault-real-8 -fdefault-double-8
FCDEP += -fallow-invalid-boz -fallow-argument-mismatch
FCDEP += -J$(depdir) -I$(depdir) -I$(CONDA_PREFIX)/include

# Debugging and Optimization
ifeq ($(debug),yes)
    FFLAGS += -g -traceback -debug full -debug-parameters all
    FFLAGS += -fpe0 -check bounds -check pointers -check uninit
    FFLAGS += -check overflow -fp-stack-check -warn unused
else
    FFLAGS += -O3
endif

# Parallelization
ifeq ($(parallel),yes)
    FFLAGS  += -xopenmp
    LDFLAGS += -xopenmp
endif

# ===========================
#  Depencency Generation
# ===========================

define makedep
    $(FCDEP) $< > $(depdir)/$*.d; \
    mv -f $(depdir)/$*.d $(depdir)/$*.d.tmp; \
    sed -e 's|.*:|$@:|' < $(depdir)/$*.d.tmp > $(depdir)/$*.d; \
    cp -f $(depdir)/$*.d $(depdir)/$*.d.tmp; \
    sed -e 's/.*://' -e 's/\\$$//' < $(depdir)/$*.d.tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(depdir)/$*.d; \
    rm -f $(depdir)/$*.d.tmp; \
    rm -f $*.i90
endef

# ===========================
#  Object Files
# ===========================

# object files
obj_core := $(addprefix $(objdir)/, \
    header.o \
    auxiliary_routines.o \
    read_xsdb.o \
    optic_molec.o \
    OpticM_module.o \
    optic_cirrus.o \
    aerosol_input_module.o \
    read_settings.o \
    depol.o \
    optic_input_module.o \
    tau_grid.o \
    header_gsd.o \
    perturbation.o \
    viewing.o \
    single_scat.o \
    gautool_acs.o \
    gaus_mat.o \
    gaus_pol.o \
    rad_intf.o \
    ocean_fresnel_module.o \
    lintran_interface.o \
    rad_trans_intf.o \
    spectral_response.o \
    forward_model.o \
    forward_model_noscat.o \
    pt_inversion.o \
    adhoc_inversion.o \
    profile_inversion.o \
    spectrum_internal.o \
    atmosphere_internal.o \
    read_miprep.o \
    retrieval.o \
)

obj_create := $(addprefix $(objdir)/, \
    spectrum_interface_create.o \
    atmosphere_interface_create.o \
    solar_model_create.o \
    noise.o \
    read_synsettings.o \
    calculate_syn_spectrum.o \
    main_create.o \
)

obj_retrieve := $(addprefix $(objdir)/, \
    read_errors.o \
    spectrum_interface_retrieve.o \
    atmosphere_interface_retrieve.o \
    solar_model_retrieve.o \
    synthetic_input.o \
    diagnostics_retrieve.o \
    wrapper_retrieve.o \
    main_retrieve.o \
)

obj_lntrn := $(addprefix $(objdir)/, \
    lintran_header.o \
    lintran_constants_module.o \
    lintran_types_module.o \
    lintran_grid_module.o \
    lintran_geometry_module.o \
    lintran_pixel_module.o \
    lintran_tlc_module.o \
    lintran_fields_module.o \
    lintran_matrix_module.o \
    lintran_interface_module.o \
    lintran_nst1_module.o \
    lintran_nst3_module.o \
    lintran_nst4_module.o \
    lintran_remotec_module.o \
    lintran_module.o \
)

num_tools := $(addprefix $(ntdir)/, \
    drealft.f90 \
    dpythag.f90 \
    dfour1.f90 \
    dsvdcmp.f90 \
    spline.f90 \
    linterp.f90 \
    ludcmp.f90 \
    lubksb.f90 \
    sort.f90 \
    qpythag.f90 \
    qsvdcmp.f90 \
    dtwofft.f90 \
    indexx.f90 \
)

libraries := \
    $(ntdir)/libNts.a

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
	-rm -rf $(depdir)
	-rm -rf $(objdir)
	-rm -rf $(moddir)
	-rm -f RemoTeC_create
	-rm -f RemoTeC_retrieve
	-rm $(ntdir)/libNts.a
