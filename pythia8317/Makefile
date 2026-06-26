# Makefile is a part of the PYTHIA event generator.
# Copyright (C) 2026 Torbjorn Sjostrand.
# PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
# Please respect the MCnet Guidelines, see GUIDELINES for details.
# Author: Philip Ilten, October 2014 - November 2017.
#
# This is is the Makefile used to build PYTHIA on POSIX systems.
# Example usage is:
#     make -j2
# For help using the make command please consult the local system documentation,
# i.e. "man make" or "make --help".

################################################################################
# VARIABLES: Definition of the relevant variables from the configuration script
# and the distribution structure.
################################################################################

# Set the shell.
SHELL=/usr/bin/env bash

# Include the configuration and set the local directory structure.
ifeq (,$(findstring clean, $(MAKECMDGOALS)))
  ifneq (Makefile.inc,$(MAKECMDGOALS))
    LOCAL_CFG:=$(shell $(MAKE) Makefile.inc))
    include Makefile.inc
  endif
endif
LOCAL_BIN=bin/
LOCAL_DOCS=AUTHORS COPYING GUIDELINES README ../../examples/Makefile.inc
LOCAL_EXAMPLE=examples
LOCAL_INCLUDE=include
LOCAL_LIB=lib
LOCAL_SHARE=share/Pythia8
LOCAL_SRC=src
LOCAL_TMP=tmp
LOCAL_MKDIRS:=$(shell mkdir -p $(LOCAL_TMP) $(LOCAL_LIB))
CXX_COMMON:=-I$(LOCAL_INCLUDE) $(CXX_COMMON)
OBJ_COMMON:=-MD $(CXX_COMMON) $(OBJ_COMMON)
LIB_COMMON=-pthread -Wl,-rpath,../lib:$(PREFIX_LIB) -ldl $(GZIP_LIB)

# Determine the archiver and archiver flags.
AR      ?= ar
ARFLAGS ?= cr

# PYTHIA.
OBJECTS=$(patsubst $(LOCAL_SRC)/%.cc,$(LOCAL_TMP)/%.o,\
	$(sort $(wildcard $(LOCAL_SRC)/*.cc)))
TARGETS=$(LOCAL_LIB)/libpythia8.a $(LOCAL_LIB)/libpythia8$(LIB_SUFFIX)

# LHAPDF.
ifeq ($(LHAPDF5_USE),true)
  TARGETS+=$(LOCAL_LIB)/libpythia8lhapdf5.so
endif
ifeq ($(LHAPDF6_USE),true)
  TARGETS+=$(LOCAL_LIB)/libpythia8lhapdf6.so
endif

# MG5 matrix element plugins.
ifeq ($(MG5MES_USE),true)
  TARGETS+=mg5mes
endif

# POWHEG (needs directory that contains just POWHEG libraries).
ifeq ($(POWHEG_USE),true)
  TARGETS+=$(LOCAL_LIB)/libpythia8powhegHooks.so
  POWHEG_DIR=$(subst -L,,$(filter -L%,$(POWHEG_LIB)))/
  ifneq ($(POWHEG_DIR),.)
    TARGETS+=$(patsubst $(POWHEG_DIR)lib%.so,\
	     $(LOCAL_LIB)/libpythia8powheg%.so,$(wildcard $(POWHEG_DIR)*))
  endif
endif

# Define RIVET options and fix C++ version, rpath, missing HDF5.
ifeq ($(RIVET_USE),true)
  COMMA=,
  RIVET_VERSION=$(shell $(RIVET_BIN)$(RIVET_CONFIG) --version)
  RIVET_LPATH=$(filter -L%,$(shell $(RIVET_BIN)$(RIVET_CONFIG) --ldflags))
  RIVET_FLAGS=$(subst -L,-Wl$(COMMA)-rpath$(COMMA),$(RIVET_LPATH))
  RIVET_FLAGS+= $(shell $(RIVET_BIN)$(RIVET_CONFIG) --cppflags --libs)
  RIVET_CSTD=c++14
  ifeq ("4.0.0","$(word 1, $(sort 4.0.0 $(RIVET_VERSION)))")
    RIVET_CSTD=c++17
    RIVET_LDIR=$(shell $(RIVET_BIN)$(RIVET_CONFIG) --libdir)
    RIVET_HDF5=$(shell nm $(RIVET_LDIR)/libRivet$(LIB_SUFFIX) | grep H5open)
    ifneq ($(strip $(RIVET_HDF5)),)
      RIVET_FLAGS+= -lhdf5
    endif
    TARGETS+=$(LOCAL_LIB)/libpythia8rivet.so
  endif
  RIVET_OPTS=$(CXX_COMMON:c++11=$(RIVET_CSTD)) $(RIVET_FLAGS) $(CXX_DTAGS)
endif

# Define HepMC3 options.
ifeq ($(HEPMC3_USE),true)
  HEPMC3_OPTS=$(CXX_COMMON) $(HEPMC3_INCLUDE) $(HEPMC3_LIB) -DHEPMC3
  TARGETS+=$(LOCAL_LIB)/libpythia8hepmc3.so
endif

# Python.
ifeq ($(PYTHON_USE),true)
  TARGETS+=python
endif

################################################################################
# RULES: Definition of the rules used to build PYTHIA.
################################################################################

# Rules without physical targets (secondary expansion for documentation).
.SECONDEXPANSION:
.PHONY: all install clean distclean mg5mes python

# All targets.
all: $(TARGETS) $(addprefix $(LOCAL_SHARE)/, $(LOCAL_DOCS))

# The documentation.
$(addprefix $(LOCAL_SHARE)/, $(LOCAL_DOCS)): $$(notdir $$@)
	cp $^ $@

# The Makefile configuration.
Makefile.inc:
	./configure

# Auto-generated (with -MD flag) dependencies.
-include $(LOCAL_TMP)/*.d $(LOCAL_TMP)/mg5mes/*.d

# PYTHIA.
$(LOCAL_TMP)/Pythia.o: $(LOCAL_SRC)/Pythia.cc Makefile.inc
	$(CXX) $< -o $@ -c $(OBJ_COMMON) -DXMLDIR=\"$(PREFIX_SHARE)/xmldoc\"
$(LOCAL_TMP)/FJcore.o: $(LOCAL_SRC)/FJcore.cc
	$(CXX) $< -o $@ -c $(OBJ_COMMON) -DFJCORE_HAVE_LIMITED_THREAD_SAFETY
$(LOCAL_TMP)/Streams.o: $(LOCAL_SRC)/Streams.cc Makefile.inc
	$(CXX) $< -o $@ -c $(OBJ_COMMON)
$(LOCAL_TMP)/%.o: $(LOCAL_SRC)/%.cc
	$(CXX) $< -o $@ -c $(OBJ_COMMON)
$(LOCAL_LIB)/libpythia8.a: $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $^
$(LOCAL_LIB)/libpythia8$(LIB_SUFFIX): $(OBJECTS)
	$(CXX) $^ -o $@ $(CXX_COMMON) $(CXX_SHARED) $(CXX_SONAME)$(notdir $@)\
	  $(LIB_COMMON) $(CXX_DTAGS)

# LHAPDF (turn off all warnings for readability).
$(LOCAL_TMP)/LHAPDF%Plugin.o: $(LOCAL_INCLUDE)/Pythia8Plugins/LHAPDF%.h
	$(CXX) -x c++ $< -o $@ -c -MD -w $(CXX_COMMON) $(LHAPDF$*_INCLUDE)
$(LOCAL_LIB)/libpythia8lhapdf%.so: $(LOCAL_TMP)/LHAPDF%Plugin.o\
	$(LOCAL_LIB)/libpythia8$(LIB_SUFFIX)
	$(CXX) $< -o $@ $(CXX_COMMON) $(CXX_SHARED) $(CXX_SONAME)$(notdir $@)\
	 $(LHAPDF$*_LIB) -lLHAPDF -Llib -lpythia8

# POWHEG.
$(LOCAL_TMP)/LHAPowheg.o: $(LOCAL_INCLUDE)/Pythia8Plugins/LHAPowheg.h
	$(CXX) -x c++ $< -o $@ -c -MD -w $(CXX_COMMON)
$(LOCAL_TMP)/PowhegHooks.o: $(LOCAL_INCLUDE)/Pythia8Plugins/PowhegHooks.h
	$(CXX) -x c++ $< -o $@ -c -MD -w $(CXX_COMMON)
$(LOCAL_LIB)/libpythia8powheg%.so: $(POWHEG_DIR)lib%.so\
	$(LOCAL_TMP)/LHAPowheg.o $(LOCAL_LIB)/libpythia8$(LIB_SUFFIX)
	$(CXX) $(LOCAL_TMP)/LHAPowheg.o -o $@ $(CXX_COMMON) $(CXX_SHARED)\
	 $(CXX_SONAME)$(notdir $@) -Llib -lpythia8\
	 -Wl,-rpath,../lib:$(POWHEG_DIR) -L$(POWHEG_DIR) -l$*
$(LOCAL_LIB)/libpythia8powhegHooks.so: $(LOCAL_TMP)/PowhegHooks.o\
	$(LOCAL_LIB)/libpythia8$(LIB_SUFFIX)
	$(CXX) $< -o $@ $(CXX_COMMON) $(CXX_SHARED) $(CXX_SONAME)$(notdir $@)\
	 -Llib -lpythia8

# RIVET.
$(LOCAL_LIB)/libpythia8rivet.so: $(LOCAL_INCLUDE)/Pythia8Plugins/RivetHooks.h
	$(CXX) -x c++ $< -o $@ -w $(RIVET_OPTS) $(CXX_SHARED)\
	 $(CXX_SONAME)$(notdir $@) -Wl,-undefined,dynamic_lookup

# HepMC3.
$(LOCAL_LIB)/libpythia8hepmc3.so: $(LOCAL_INCLUDE)/Pythia8Plugins/HepMC3Hooks.h
	$(CXX) -x c++ $< -o $@ -w $(HEPMC3_OPTS) $(CXX_SHARED)\
	 $(CXX_SONAME)$(notdir $@) -Wl,-undefined,dynamic_lookup

# MG5 matrix element plugins.
mg5mes:
	cd $(MG5MES_BIN) && $(MAKE)

# Python.
python: $(LOCAL_LIB)/libpythia8$(LIB_SUFFIX)
	cd plugins/python && $(MAKE)

# Install.
install: all
	mkdir -p $(PREFIX_BIN) $(PREFIX_INCLUDE) $(PREFIX_LIB) $(PREFIX_SHARE)
	rsync -a $(LOCAL_BIN)/* $(PREFIX_BIN)
	rsync -a $(LOCAL_INCLUDE)/* $(PREFIX_INCLUDE)
	rsync -a $(LOCAL_LIB)/* $(PREFIX_LIB)
	rsync -a $(LOCAL_SHARE)/* $(PREFIX_SHARE)
	rsync -a $(LOCAL_EXAMPLE) $(PREFIX_SHARE)

# Clean.
clean:
	cd plugins/python && $(MAKE) clean
	cd plugins/mg5mes && $(MAKE) clean
	rm -rf $(LOCAL_TMP) $(LOCAL_LIB)
	cd $(LOCAL_EXAMPLE) && $(MAKE) clean

# Clean all temporary and generated files.
distclean: clean
	cd examples && make clean
	find . -type f -name Makefile.inc -print0 | xargs -0 rm -f
	find . -type f -name "*~" -print0 | xargs -0 rm -f
	find . -type f -name "#*" -print0 | xargs -0 rm -f
	rm -rf $(LOCAL_BIN)
	rm -f $(LOCAL_SHARE)/AUTHORS
	rm -f $(LOCAL_SHARE)/COPYING
	rm -f $(LOCAL_SHARE)/GUIDELINES
	rm -f $(LOCAL_SHARE)/README
