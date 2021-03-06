# for exception handling need everything compiled via g++
CC=g++
CXX=g++
LIBS=-lm -lstdc++ -lpthread

# the BVECTOR thing is to force g++ to use a normal vector<bool>
# instead of the special one.  this is necessary in order to use
# IndexedVector<*,bool>.  probably not portable.
# you'll need to use it in your project makefile as well.
#  -D_BVECTOR_H=1 -D_STL_BVECTOR_H=1
CFLAGS= -Wall -g $(OTHER_CFLAGS) -IDisplays -I.
# -O3
LDFLAGS = $(LIBS)

# to compile without MPI do 'make fresh MPI=no' ... etc.
MPI=no
ifeq ($(MPI),yes)
LOCALCFLAGS+= -I/usr/local/mpi/include -DMPI
endif

VNL=yes
ifeq ($(VNL),yes)
VXLDIR ?= ../vxl
LOCALCFLAGS+=-DVNL -I$(VXLDIR)/core -I$(VXLDIR)/core/vnl -I$(VXLDIR)/core/vnl/algo -I$(VXLDIR)/vcl -I$(VXLDIR)/lib
#CFLAGS+=-DVNL -I$(VXLDIR)-1.0-beta/vxl -I../vxl-1.0-beta/vxl/config.Linux2-gcc-2.95 -finline-functions
LDFLAGS += -L$(VXLDIR)/lib -lvnl_algo -lvnl -lvcl -lnetlib
#LDFLAGS+=-L../vxl-1.0-beta/vxl/lib/Linux2-gcc-2.95 -L../vxl-1.0-beta/v3p/lib/Linux2-gcc-2.95 -Wl,-rpath,../vxl-1.0-beta/v3p/lib/Linux2-gcc-2.95:../vxl-1.0-beta/vxl/lib/Linux2-gcc-2.95 -lvnl-algo -lnetlib -lvnl -lvcl -lstdc++ -ldl -rdynamic -lm -lc
endif

CVODE=yes
ifeq ($(CVODE),yes)
LOCALCFLAGS+=-DCVODE
LOCALCFLAGS+=-I/usr/include/cvode -I/usr/include/nvector -I/usr/include/sundials
endif

EXECSTREAMDIR ?= ../libexecstream

SOURCES = main.cpp Community.cpp Communicator.cpp Parameters.cpp \
	OutputController.cpp NodeOutputController.cpp Simulation.cpp \
	Integrator.cpp NRIntegrator.cpp Iterator.cpp \
	Indexing.cpp util.cpp genrand2.cpp rand.cpp Algebra.cpp \
	Displays/PopulationDisplay.cpp Displays/DotDisplay.cpp \
	Displays/EvolutionaryTreeDisplay.cpp \
	Displays/ConstrainedPhenoDisplay2d.cpp \
	Displays/TimeSeriesDisplay.cpp
#CVIntegrator.cpp 
OBJS = $(SOURCES:.cpp=.o) $(EXECSTREAMDIR)/exec-stream.o
#CVodeWithStoppingCondition.o 
#l_cvode.o

THIS = GNUmakefile

ARCHIVE = libadap-dyn.a

default: $(ARCHIVE)

lib $(ARCHIVE): $(OBJS) $(THIS)
	$(AR) rs $(ARCHIVE) $(OBJS)

$(EXE): $(OBJS)
	$(CXX) $(CFLAGS) $(LOCALCFLAGS) -o $(EXE) $(OBJS) $(LIBS)

run: $(EXE) clean
	$(EXE)

test-rand: test-rand.cpp $(OBJS)
	$(CXX) $(CFLAGS) $(LOCALCFLAGS) -o test-rand test-rand.cpp $(OBJS)

ALLOBJS = $(OBJS) $(OTHER_OBJS)

FLOTSAM = core test-rand *.idb *.pdb *.exe *.ilk *.obj *~ *\#

clean:
	$(RM) $(ARCHIVE) $(EXE) $(ALLOBJS) $(PFILES) $(FLOTSAM)

tags TAGS:
	etags *.{h,cpp}

fresh over: clean default

# stuff for automatic header file dependencies

%.o : %.cpp
	$(CXX) $(CFLAGS) $(LOCALCFLAGS) -MD -c $< -o $@
	@cp $*.d $*.P; \
	    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	        -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	    rm -f $*.d

%.o : %.c
	$(CC) $(CFLAGS) $(LOCALCFLAGS) -MD -c $< -o $@
	@cp $*.d $*.P; \
	    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	        -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	    rm -f $*.d

%.E : %.cpp
	$(CXX) $(CFLAGS) $(LOCALCFLAGS) -E $< -o $@

PXXS = $(SOURCES:.cpp=.P)
PFILES = $(PXXS:.c=.P)
-include $(PFILES)

