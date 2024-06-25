# $Id: GNUmakefile,v 1.1 1999-01-07 16:05:42 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := hdgeant4
G4TARGET := $(name)
G4EXLIB := true
G4LIB_BUILD_SHARED := true
G4WORKDIR := $(shell pwd)
G4SYSTEM := Linux-g++

ifndef G4ROOT
    $(error G4ROOT is not set, please set it and try again)
endif

ifndef G4SYSTEM
    $(error Geant4 environment not set up, please source $(G4ROOT)/share/Geant4-10.*/geant4make/geant4make.sh and try again)
endif

ifdef DIRACXX_DIR
    DIRACXX_CMAKE := $(shell if [ -f $(DIRACXX_HOME)/CMakeLists.txt ]; then echo true; else echo false; fi)
    ifeq ($(DIRACXX_CMAKE), true)
        CPPFLAGS += -I$(DIRACXX_DIR)/include -DUSING_DIRACXX -L$(DIRACXX_DIR)/lib -lDirac
    else
        CPPFLAGS += -I$(DIRACXX_DIR) -DUSING_DIRACXX -L$(DIRACXX_DIR) -lDirac
    endif
endif

PYTHON_CONFIG = python-config

CPPFLAGS += -I$(HDDS_HOME) -I./src -I./src/G4fixes
CPPFLAGS += -I./src/G4debug
CPPFLAGS += -I$(HALLD_RECON_HOME)/$(BMS_OSNAME)/include
CPPFLAGS += -I$(JANA_HOME)/include
CPPFLAGS += -I$(shell root-config --incdir)
CPPFLAGS += $(shell $(PYTHON_CONFIG) --includes)
CPPFLAGS += -Wno-unused-parameter -Wno-unused-but-set-variable
CPPFLAGS += -DUSE_SSE2
#CPPFLAGS += -I/usr/include/Qt
#CPPFLAGS += -DLINUX_CPUTIME_PROFILING=1
#CPPFLAGS += -DCHECK_OVERLAPS_MM=1e-6
CPPFLAGS += -DBYPASS_DRAWING_CLIPPED_VOLUMES
CPPFLAGS += -DLAYERED_GEOMETRY_PICKING_EXTENSIONS
CPPFLAGS += -DREDUCE_OPTIMIZATION_OF_CDC=1
#CPPFLAGS += -DG4UI_USE_EXECUTIVE
CPPFLAGS += -DG4VIS_BUILD_OPENGL_DRIVER
CPPFLAGS += -DG4VIS_BUILD_OPENGLX_DRIVER
CPPFLAGS += -DG4MULTITHREADED
#CPPFLAGS += -DUSING_BERNARD
#CPPFLAGS += -DMERGE_STEPS_BEFORE_HITS_GENERATION
#CPPFLAGS += -DVERBOSE_RANDOMS=1
#CPPFLAGS += -DFORCE_PARTICLE_TYPE_CHARGED_GEANTINO
#CPPFLAGS += -DBP_DEBUG
#CPPFLAGS += -DMOD_SPONCE
#CPPFLAGS += -DDEBUG_PLACEMENT
#CPPFLAGS += -DDEBUG_SECTIONPLANE
#CPPFLAGS += -DDEBUG_SECTIONPLANE_ZAVE
ifneq (, $(wildcard $(HALLD_RECON_HOME)/src/libraries/DIRC/DDIRCPmtHit.*))
	CPPFLAGS += -DDIRCTRUTHEXTRA
endif

G4LIB_USE_GDML = 1
CPPVERBOSE = 1
#G4DEBUG = 1

hdgeant4_sources := $(filter-out src/CobremsGeneration.cc, $(wildcard src/*.cc))
Cobrems_sources := $(wildcard src/Cobrems*.cc)
G4fixes_sources := $(wildcard src/G4fixes/*.cc)
G4debug_sources := $(wildcard src/G4debug/*.cc)
HDDS_sources := $(HDDS_HOME)/XString.cpp $(HDDS_HOME)/XParsers.cpp $(HDDS_HOME)/hddsCommon.cpp

ROOTLIBS = $(shell root-config --libs) -lGeom -lTMVA -lTreePlayer
ifeq ($(shell test -s $(ROOTSYS)/lib/libtbb.so && echo -n yes),yes)
    ROOTLIBS += $(ROOTSYS)/lib/libtbb.so.2
endif
ifeq ($(shell test -s $(VDTHOME)/lib/libvdt.so && echo -n yes),yes)
    ROOTLIBS += -L$(VDTHOME)/lib -lvdt
endif

ifdef SQLITECPP_VERSION
    SQLITECPP_MAJOR_VERSION := $(shell echo $(SQLITECPP_VERSION) | awk -F. '{print $$1}')
    SQLITECPP_MINOR_VERSION := $(shell echo $(SQLITECPP_VERSION) | awk -F. '{print $$2}')
    SQLITECPP_LIBDIR := $(shell if [[ $(SQLITECPP_MAJOR_VERSION) -ge 3 || $(SQLITECPP_MAJOR_VERSION) -eq 2 && $(SQLITECPP_MINOR_VERSION) -ge 5 ]]; then echo lib64; else echo lib; fi)
else
    SQLITECPP_LIBDIR = lib64
endif

DANALIBS = -L$(HALLD_RECON_HOME)/$(BMS_OSNAME)/lib -lHDGEOMETRY -lDANA \
           -lANALYSIS -lBCAL -lCCAL -lECAL -lCDC -lCERE -lTRD -lDIRC -lFCAL \
           -lFDC -lFMWPC -lHDDM -lPAIR_SPECTROMETER -lPID -lRF \
           -lSTART_COUNTER -lTAGGER -lTOF -lTPOL -lTRACKING \
           -lTRIGGER -lDAQ -lTTAB -lEVENTSTORE -lKINFITTER -lTAC \
           -L$(SQLITECPP_HOME)/$(SQLITECPP_LIBDIR) -lSQLiteCpp -L$(SQLITE_HOME)/lib -Wl,-rpath=$(SQLITE_HOME)/lib -lsqlite3 \
           -lxstream -lbz2 -lz \
           -L/usr/lib64/mysql -lmysqlclient\
           -L$(JANA_HOME)/lib -lJANA \
           -L$(CCDB_HOME)/lib -lccdb \
           -L$(EVIOROOT)/lib -levioxx -levio \
           -lpthread -ldl

ifdef ETROOT
DANALIBS += -L$(ETROOT)/lib -let -let_remote
endif

G4shared_libs := $(wildcard $(G4ROOT)/lib64/*.so)
ifeq ($(PYTHON_GE_3), true)
  BOOST_PYTHON_LIB = -lboost_python$(PYTHON_MAJOR_VERSION)$(PYTHON_MINOR_VERSION) -lpython$(PYTHON_MAJOR_VERSION).$(PYTHON_MINOR_VERSION)
  CPPFLAGS += -DBOOST_NO_AUTO_PTR -DBOOST_BIND_GLOBAL_PLACEHOLDERS
else
  BOOST_PYTHON_LIB = -lboost_python
endif

ifdef HDF5ROOT
CPPFLAGS += -DHDF5_SUPPORT -I ${HDF5ROOT}/include
DANALIBS +=	-L $(HDF5ROOT)/lib -lhdf5_cpp -lhdf5_hl -lhdf5 -lsz -lz -lbz2 -ldl 
endif

INTYLIBS += -Wl,--whole-archive $(DANALIBS) -Wl,--no-whole-archive
INTYLIBS += -fPIC -I$(HDDS_HOME) -I$(XERCESCROOT)/include
INTYLIBS += -L${XERCESCROOT}/lib -lxerces-c
INTYLIBS += -L$(G4TMPDIR) -lhdds
INTYLIBS += $(BOOST_PYTHON_LIB)
INTYLIBS += -L$(shell $(PYTHON_CONFIG) --prefix)/lib
INTYLIBS += $(shell $(PYTHON_CONFIG) --ldflags) $(PYTHON_LIB_OPTION)
INTYLIBS += -L$(G4ROOT)/lib64 $(patsubst $(G4ROOT)/lib64/lib%.so, -l%, $(G4shared_libs))
INTYLIBS += -lgfortran
INTYLIBS += -L/usr/lib64
INTYLIBS += -ltirpc
INTYLIBS += $(ROOTLIBS)

EXTRALIBS += -lG4fixes -lGLU

.PHONY: all
all: hdds cobrems G4fixes_symlink g4fixes sharedlib exe lib bin g4py

include $(G4INSTALL)/config/binmake.gmk

cobrems: $(G4TMPDIR)/libcobrems.so
hdds:  $(G4TMPDIR)/libhdds.so
g4fixes: $(G4TMPDIR)/libG4fixes.so

CXXFLAGS = -g -O4 -fPIC -W -Wall -pedantic -Wno-non-virtual-dtor -Wno-long-long

GCCVERSION = $(shell gcc --version | awk -F'[. ]*' '/gcc/{print $$3}')
ifeq ($(shell test $(GCCVERSION) -ge 8; echo $$?),0)
    CPPFLAGS += -std=c++17
    CXXFLAGS += -std=c++17
else ifeq ($(shell test $(GCCVERSION) -ge 5; echo $$?),0)
    CPPFLAGS += -std=c++14
    CXXFLAGS += -std=c++14
else
    CPPFLAGS += -std=c++11
    CXXFLAGS += -std=c++11
endif

HDDSDIR := $(G4TMPDIR)/hdds
G4FIXESDIR := $(G4TMPDIR)/G4fixes
G4fixes_symlink:
	@if [ ! -L src/G4fixes ]; then\
	    echo "ERROR - symbolic link G4fixes does not exist in directory src,"\
	    "cannot continue!";\
	    false;\
	fi

$(G4TMPDIR)/libcobrems.so: $(Cobrems_sources)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -Wl,--export-dynamic -Wl,-soname,libcobrems.so \
	-shared -o $@ $^ $(G4shared_libs) $(BOOST_PYTHON_LIB)

hdgeant4_objects := $(patsubst src/%.cc, $(G4TMPDIR)/%.o, $(hdgeant4_sources))
G4fixes_objects := $(patsubst src/G4fixes/%.cc, $(G4FIXESDIR)/%.o, $(G4fixes_sources))
G4debug_objects := $(patsubst src/G4debug/%.cc, $(G4FIXESDIR)/%.o, $(G4debug_sources))
sharedlib: $(G4TMPDIR)/libhdgeant4.so

$(G4TMPDIR)/libhdgeant4.so: $(hdgeant4_objects)

$(G4TMPDIR)/libG4fixes.so: $(G4FIXESDIR)/G4fixes.o $(G4fixes_objects) $(G4debug_objects)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -Wl,--export-dynamic -Wl,-soname,libG4fixes.so \
	-shared -o $@ $^ $(G4shared_libs) $(BOOST_PYTHON_LIB)

$(G4FIXESDIR)/G4fixes.o: src/G4fixes.cc
	@mkdir -p $(G4FIXESDIR)
ifdef CPPVERBOSE
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $^
else
	@echo Compiling $*.cc ...
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $^
endif

$(G4fixes_objects): $(G4FIXESDIR)/%.o: src/G4fixes/%.cc
	@mkdir -p $(G4FIXESDIR)
ifdef CPPVERBOSE
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $^
else
	@echo Compiling $*.cc ...
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $^
endif

$(G4debug_objects): $(G4FIXESDIR)/%.o: src/G4debug/%.cc
	@mkdir -p $(G4FIXESDIR)
ifdef CPPVERBOSE
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $^
else
	@echo Compiling $*.cc ...
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $^
endif

$(HDDSDIR)/%.o: $(HDDS_HOME)/%.cpp
	@mkdir -p $(HDDSDIR)
ifdef CPPVERBOSE
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $^
else
	@echo Compiling $*.cc ...
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $^
endif

$(G4TMPDIR)/libhdds.so: $(patsubst $(HDDS_HOME)/%.cpp, $(HDDSDIR)/%.o, $(HDDS_sources))
	@$(CXX) -Wl,-soname,libhdds.so -shared -o $@ $^

exe:
	@mkdir -p $(G4LIBDIR)/$@

g4py: $(G4LIBDIR)/../../../g4py/HDGeant4/libhdgeant4.so \
      $(G4LIBDIR)/../../../g4py/Cobrems/libcobrems.so \
      $(G4LIBDIR)/../../../g4py/G4fixes/libG4fixes.so

$(G4LIBDIR)/../../../g4py/HDGeant4/libhdgeant4.so: $(G4LIBDIR)/libhdgeant4.so
	@rm -f $@
	@cd g4py/HDGeant4 && ln -s ../../tmp/*/hdgeant4/libhdgeant4.so .

$(G4LIBDIR)/../../../g4py/Cobrems/libcobrems.so: $(G4LIBDIR)/libcobrems.so 
	@rm -f $@
	@cd g4py/Cobrems && ln -s ../../tmp/*/hdgeant4/libcobrems.so .

$(G4LIBDIR)/../../../g4py/G4fixes/libG4fixes.so: $(G4LIBDIR)/libG4fixes.so 
	@rm -f $@
	@cd g4py/G4fixes && ln -s ../../tmp/*/hdgeant4/libG4fixes.so .

utils: $(G4BINDIR)/beamtree $(G4BINDIR)/genBH $(G4BINDIR)/adapt $(G4BINDIR)/geneBH $(G4BINDIR)/samplesep

$(G4BINDIR)/beamtree: src/utils/beamtree.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ -L$(G4LIBDIR) -lhdgeant4 $(ROOTLIBS) -Wl,-rpath=$(G4LIBDIR) $(G4shared_libs) $(BOOST_PYTHON_LIB)

$(G4BINDIR)/genBH: src/utils/genBH.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ -L$(G4LIBDIR) -lhdgeant4 $(DANALIBS) $(ROOTLIBS) -Wl,-rpath=$(G4LIBDIR)

$(G4BINDIR)/geneBH: src/utils/geneBH.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ -L$(G4LIBDIR) -lhdgeant4 $(DANALIBS) $(ROOTLIBS) -Wl,-rpath=$(G4LIBDIR)

$(G4BINDIR)/adapt: src/utils/adapt.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ -L$(G4LIBDIR) -lhdgeant4 $(DANALIBS) $(ROOTLIBS) -Wl,-rpath=$(G4LIBDIR)

$(G4BINDIR)/samplesep: src/utils/samplesep.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -O4 -fopenmp -o $@ $^ -L$(G4LIBDIR) -lhdgeant4 $(DANALIBS) $(ROOTLIBS) -Wl,-rpath=$(G4LIBDIR)

show_env:
	@echo PYTHON_VERSION = $(PYTHON_VERSION)
	@echo PYTHON_MAJOR_VERSION = $(PYTHON_MAJOR_VERSION)
	@echo PYTHON_MINOR_VERSION = $(PYTHON_MINOR_VERSION)
	@echo PYTHON_SUBMINOR_VERSION = $(PYTHON_SUBMINOR_VERSION)
	@echo PYTHON_GE_3 = $(PYTHON_GE_3)

diff:
	diff -q -r ../jlab . -x ".[a-z]*" -x tmp -x bin -x "*.pyc" -x "*.so" -x test -x "*-orig"

check_flags:
	@echo CPPFLAGS=$(CPPFLAGS)
	@echo __cpluslplus=$(__cplusplus)