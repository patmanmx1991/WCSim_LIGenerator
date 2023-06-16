# $Id: GNUmakefile,v 1.17 2006/09/04 15:43:27 t2k Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

G4DEBUG = 1

name := WCSim
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

ROOTCFLAGS   := $(shell root-config --cflags) -DUSE_ROOT -fPIC
ROOTLIBS     := $(shell root-config --libs)

LIBNAME := WCSim

# NOTE: Geant4.7.0 changes the way Maximum Step size is defined.  
# We need extra code for versions 4.7.0 and above; eventually 
# everyone should upgrade to geant4.7
ifneq (,$(findstring 4.7,$(G4INSTALL)))
GEANT4_7_0 = 0
else
GEANT4_7_0 = 1
endif

ifdef GEANT4_7_0
CPPFLAGS += -DGEANT4_7_0
endif

ifdef GCCVERS296
CPPFLAGS += -DUSE_STRSTREAM
endif

CPPFLAGS  += -I$(ROOTSYS)/include $(ROOTCFLAGS) -std=c++0x            ##for unordered_map 
EXTRALIBS += $(ROOTLIBS)

EXTRA_LINK_DEPENDENCIES := 

DOXYGEN_VERSION := $(shell doxygen --version 2>/dev/null)
ifdef DOXYGEN_VERSION
DOXYGEN_EXISTS = 1
else
DOXYGEN_EXISTS = 0
endif

# Set flag for git version to be used as c++ preprocessor macro
CPPFLAGS += -DGIT_HASH=\"$(shell git describe --always --long --tags --dirty)\"

# Set flag for checking geometry overlap in construction
CPPFLAGS += -DWCSIM_CHECK_GEOMETRY_OVERLAPS=false 

.PHONY: all
all: rootcint shared libWCSim.a lib bin

# Note dependencies not yet set up right yet

ROOTSO    := libWCSimRoot.so

ROOTSRC  := ./src/WCSimRootEvent.cc ./include/WCSimRootEvent.hh ./src/WCSimRootGeom.cc ./include/WCSimRootGeom.hh ./include/WCSimPmtInfo.hh ./src/WCSimEnumerations.cc ./include/WCSimEnumerations.hh ./src/WCSimRootOptions.cc ./include/WCSimRootOptions.hh ./src/TJNuBeamFlux.cc ./include/TJNuBeamFlux.hh ./src/TNRooTrackerVtx.cc ./include/TNRooTrackerVtx.hh ./include/WCSimRootLinkDef.hh

ROOTOBJS  := $(G4WORKDIR)/tmp/$(G4SYSTEM)/WCSim/WCSimRootEvent.o $(G4WORKDIR)/tmp/$(G4SYSTEM)/WCSim/WCSimRootGeom.o $(G4WORKDIR)/tmp/$(G4SYSTEM)/WCSim/WCSimPmtInfo.o $(G4WORKDIR)/tmp/$(G4SYSTEM)/WCSim/WCSimEnumerations.o $(G4WORKDIR)/tmp/$(G4SYSTEM)/WCSim/WCSimRootOptions.o $(G4WORKDIR)/tmp/$(G4SYSTEM)/WCSim/TNRooTrackerVtx.o $(G4WORKDIR)/tmp/$(G4SYSTEM)/WCSim/TJNuBeamFlux.o $(G4WORKDIR)/tmp/$(G4SYSTEM)/WCSim/WCSimRootDict.o 

shared: $(ROOTSRC) $(ROOTOBJS) 
	g++ -shared -O $(ROOTOBJS) -o $(ROOTSO) $(ROOTLIBS)

libWCSim.a : $(ROOTOBJS)
	$(RM) $@
	ar clq $@ $(ROOTOBJS) 

./WCSimRootDict.cxx : $(ROOTSRC)
	rootcint  -f ./WCSimRootDict.cxx -c -I./include -I$(shell root-config --incdir) WCSimRootEvent.hh WCSimRootGeom.hh  WCSimPmtInfo.hh WCSimEnumerations.hh WCSimRootOptions.hh TJNuBeamFlux.hh TNRooTrackerVtx.hh WCSimRootLinkDef.hh

rootcint: ./WCSimRootDict.cxx

doxy:
	@if [ ${DOXYGEN_EXISTS} = 1 ]; \
	then \
		doxygen WCSim_doxygen_config; \
	else\
		echo "Error: doxygen program not found in path. Exiting"; \
	fi

clean_wcsim:
	echo  $(G4WORKDIR); 
	$(RM) -r $(G4WORKDIR); $(RM) *.o *.a *.so *~ */*~ ./WCSimRootDict.*
	@if [ -d "doc/doxygen" ]; \
		then \
		rm -r doc/doxygen; \
	fi	

include $(G4INSTALL)/config/binmake.gmk

$(G4TMPDIR)/%.o: %.cxx
	@echo Compiling $*.cxx ...
	@if [ ! -d $(G4TMPDIR) ] ; then echo mkdir -p $(G4TMPDIR) ; mkdir -p $(G4TMPDIR) ; fi
	@echo $(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $(G4TMPDIR)/$(*F).o $*.cxx
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $(G4TMPDIR)/$(*F).o $*.cxx
