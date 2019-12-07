# $Id: GNUmakefile,v 1.1 1999-01-07 16:05:42 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := hdgeant4
G4TARGET := $(name)
G4EXLIB := true
G4LIB_BUILD_SHARED := true
G4WORKDIR := $(shell pwd)

ifndef G4ROOT
    $(error G4ROOT is not set, please set it and try again)
endif

ifndef G4SYSTEM
    $(error Geant4 environment not set up, please source $(G4ROOT)/share/Geant4-10.2.2/geant4make/geant4make.sh and try again)
endif

ifdef DIRACXX_HOME
    CPPFLAGS += -I$(DIRACXX_HOME) -DUSING_DIRACXX -L$(DIRACXX_HOME) -lDirac
endif

CPPFLAGS += -I$(HDDS_HOME) -I./src -I./src/G4fixes
CPPFLAGS += -I./src/G4debug
CPPFLAGS += -I$(HALLD_RECON_HOME)/$(BMS_OSNAME)/include
CPPFLAGS += -I$(JANA_HOME)/include
CPPFLAGS += -I$(shell root-config --incdir)
CPPFLAGS += -I/usr/include/Qt
CPPFLAGS += $(shell python-config --includes)
CPPFLAGS += -Wno-unused-parameter -Wno-unused-but-set-variable
CPPFLAGS += -DUSE_SSE2 -std=c++11
#CPPFLAGS += -I/usr/include/Qt
#CPPFLAGS += -DLINUX_CPUTIME_PROFILING=1
#CPPFLAGS += -DCHECK_OVERLAPS_MM=1e-6
CPPFLAGS += -DBYPASS_DRAWING_CLIPPED_VOLUMES
CPPFLAGS += -DLAYERED_GEOMETRY_PICKING_EXTENSIONS
#CPPFLAGS += -DG4UI_USE_EXECUTIVE
CPPFLAGS += -DG4VIS_BUILD_OPENGL_DRIVER
CPPFLAGS += -DG4VIS_BUILD_OPENGLX_DRIVER
CPPFLAGS += -DG4MULTITHREADED
#CPPFLAGS += -DVERBOSE_RANDOMS=1
#CPPFLAGS += -DFORCE_PARTICLE_TYPE_CHARGED_GEANTINO
#CPPFLAGS += -DBP_DEBUG
#CPPFLAGS += -DMOD_SPONCE
#CPPFLAGS += -DDEBUG_PLACEMENT
#CPPFLAGS += -DDEBUG_SECTIONPLANE
#CPPFLAGS += -DDEBUG_SECTIONPLANE_ZAVE
#CPPFLAGS += -DG4VERSION_10_04_OR_LATER=1
ifneq (, $(wildcard $(HALLD_RECON_HOME)/src/libraries/DIRC/DDIRCPmtHit.*))
	CPPFLAGS += -DDIRCTRUTHEXTRA
endif

# If you want to build against Geant4.10.03 or greater, you will need this line uncommented
#CPPFLAGS += -DG4VUSERPHYSICSLIST_HAS_GETPARTICLEITERATOR

G4LIB_USE_GDML = 1
CPPVERBOSE = 1
#G4DEBUG = 1

hdgeant4_sources := $(filter-out src/CobremsGeneration.cc, $(wildcard src/*.cc))
Cobrems_sources := $(wildcard src/Cobrems*.cc)
G4fixes_sources := $(wildcard src/G4fixes/*.cc)
G4debug_sources := $(wildcard src/G4debug/*.cc)
HDDS_sources := $(HDDS_HOME)/XString.cpp $(HDDS_HOME)/XParsers.cpp $(HDDS_HOME)/hddsCommon.cpp

ROOTLIBS = $(shell root-config --libs) -lGeom -lTMVA -lTreePlayer

DANALIBS = -L$(HALLD_RECON_HOME)/$(BMS_OSNAME)/lib -lHDGEOMETRY -lDANA \
           -lANALYSIS -lBCAL -lCCAL -lCDC -lCERE -lDIRC -lFCAL \
           -lFDC -lFMWPC -lHDDM -lPAIR_SPECTROMETER -lPID -lRF \
           -lSTART_COUNTER -lTAGGER -lTOF -lTPOL -lTRACKING -lTRD \
           -lTRIGGER -lDAQ -lTTAB -lEVENTSTORE -lKINFITTER -lTAC \
           -L$(SQLITECPP_HOME)/lib -lSQLiteCpp -L$(SQLITE_HOME)/lib -Wl,-rpath=$(SQLITE_HOME)/lib -lsqlite3 \
           -lxstream -lbz2 -lz \
           -L/usr/lib64/mysql -lmysqlclient\
           -L$(JANA_HOME)/lib -lJANA \
           -L$(CCDB_HOME)/lib -lccdb \
           -L$(EVIOROOT)/lib -levioxx -levio \
           $(ROOTLIBS) \
           -lpthread -ldl

ifdef ETROOT
DANALIBS += -L$(ETROOT)/lib -let -let_remote
endif

G4shared_libs := $(wildcard $(G4ROOT)/lib64/*.so)

INTYLIBS += -Wl,--whole-archive $(DANALIBS) -Wl,--no-whole-archive
INTYLIBS += -fPIC -I$(HDDS_HOME) -I$(XERCESCROOT)/include
INTYLIBS += -L${XERCESCROOT}/lib -lxerces-c
INTYLIBS += -L$(G4TMPDIR) -lhdds
INTYLIBS += -lboost_python -L$(shell python-config --prefix)/lib $(shell python-config --ldflags)
INTYLIBS += -L$(G4ROOT)/lib64 $(patsubst $(G4ROOT)/lib64/lib%.so, -l%, $(G4shared_libs))
INTYLIBS += -lgfortran
INTYLIBS += -L/usr/lib64

EXTRALIBS += -lG4fixes

.PHONY: all
all: hdds cobrems G4fixes_symlink g4fixes sharedlib exe lib bin g4py

include $(G4INSTALL)/config/binmake.gmk

cobrems: $(G4TMPDIR)/libcobrems.so
hdds:  $(G4TMPDIR)/libhdds.so
g4fixes: $(G4TMPDIR)/libG4fixes.so

CXXFLAGS = -O4 -fPIC -W -Wall -pedantic -Wno-non-virtual-dtor -Wno-long-long

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
	-shared -o $@ $^ $(G4shared_libs) -lboost_python

hdgeant4_objects := $(patsubst src/%.cc, $(G4TMPDIR)/%.o, $(hdgeant4_sources))
G4fixes_objects := $(patsubst src/G4fixes/%.cc, $(G4FIXESDIR)/%.o, $(G4fixes_sources))
G4debug_objects := $(patsubst src/G4debug/%.cc, $(G4FIXESDIR)/%.o, $(G4debug_sources))
sharedlib: $(G4TMPDIR)/libhdgeant4.so

$(G4TMPDIR)/libhdgeant4.so: $(hdgeant4_objects)

$(G4TMPDIR)/libG4fixes.so: $(G4FIXESDIR)/G4fixes.o $(G4fixes_objects) $(G4debug_objects)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -Wl,--export-dynamic -Wl,-soname,libG4fixes.so \
	-shared -o $@ $^ $(G4shared_libs) -lboost_python

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

utils: $(G4BINDIR)/beamtree $(G4BINDIR)/genBH

$(G4BINDIR)/beamtree: src/utils/beamtree.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ -L$(G4LIBDIR) -lhdgeant4 $(ROOTLIBS) -Wl,-rpath=$(G4LIBDIR)

$(G4BINDIR)/genBH: src/utils/genBH.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ -L$(G4LIBDIR) -lhdgeant4 $(DANALIBS) $(ROOTLIBS) -Wl,-rpath=$(G4LIBDIR)
