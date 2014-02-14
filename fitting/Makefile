ROOTCFLAGS	=	$(shell root-config --cflags)
ROOTGLIBS	=	$(shell root-config --glibs)

CXX           = g++
CXXFLAGS      = -g -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared

ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)

ifeq ($(PLATFORM),macosx)
# MacOS X with cc (GNU cc 2.95.2 and gcc 3.3)
MACOSX_MINOR := $(shell sw_vers | sed -n 's/ProductVersion://p' | cut -d . -f 2)
CXX           = c++ -lm
CXXFLAGS      = -O2 -pipe -Wall -W -Woverloaded-virtual -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE  
LD           = c++
LDFLAGS       = -O2 -g
# The SOFLAGS will be used to create the .dylib,
# the .so will be created separately
DllSuf        = dylib
ifeq ($(MACOSX_MINOR),9)
UNDEFOPT      = dynamic_lookup
LD            = MACOSX_DEPLOYMENT_TARGET=10.9 c++
else
ifeq ($(MACOSX_MINOR),6)
UNDEFOPT      = dynamic_lookup
LD            = MACOSX_DEPLOYMENT_TARGET=10.6 c++
else
ifeq ($(MACOSX_MINOR),3)
UNDEFOPT      = dynamic_lookup
LD            = MACOSX_DEPLOYMENT_TARGET=10.3 c++
else
UNDEFOPT      = suppress
LD            = c++
endif
endif
endif

endif


CXXFLAGS      += $(ROOTCFLAGS)

NGLIBS		= $(ROOTGLIBS)
NGLIBS		+= -lMathMore -lMinuit -lRooFit -lRooFitCore -lRooStats -lFoam 
GLIBS		= $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR	= ./
CXX		+= -I$(INCLUDEDIR)
OUTLIB		= ./

.SUFFIXES:	.cc,.C,.hh,.h
.PREFIXES:	../

RooHistPdfConv.o: $(INCLUDEDIR)/RooHistPdfConv.cpp
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/libRooHistPdfConv.o $(NGLIBS) $<

Tree2Datasets:	$(INCLUDEDIR)tree2Datasets.cpp
	$(CXX) $(CXXFLAGS) -o Tree2Datasets $< $(GLIBS) 

Tree2DatasetsMC: $(INCLUDEDIR)tree2DatasetsMC.cpp
	$(CXX) $(CXXFLAGS) -o Tree2DatasetsMC $< $(GLIBS) 

Fit2DDataPbPb:	$(INCLUDEDIR)fit2DData_pbpb.cpp
	$(CXX) $(CXXFLAGS) -o Fit2DDataPbPb $< $(OUTLIB)/*.o $(GLIBS)

Fit1DDataPbPb:	$(INCLUDEDIR)fit1DData_pbpb.cpp
	$(CXX) $(CXXFLAGS) -o Fit1DDataPbPb $< $(OUTLIB)/*.o $(GLIBS)

Fit1DDataPbPbSim:	$(INCLUDEDIR)fit1DData_pbpb_simul.cpp
	$(CXX) $(CXXFLAGS) -o Fit1DDataPbPbSim $< $(OUTLIB)/*.o $(GLIBS)

Sideband:	$(INCLUDEDIR)ctau.cpp
	$(CXX) $(CXXFLAGS) -o Sideband $< $(OUTLIB)/*.o $(GLIBS)

clean:
	rm -f $(OUTLIB)*.o $(OUTLIB)*.so

all:	Tree2Datasets, Fit2DData, Sideband
