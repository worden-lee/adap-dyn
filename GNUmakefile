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
LOCALCFLAGS+=-DVNL -I../vxl/core -I../vxl/core/vnl -I../vxl/core/vnl/algo -I../vxl/vcl -I../vxl/lib
#CFLAGS+=-DVNL -I../vxl-1.0-beta/vxl -I../vxl-1.0-beta/vxl/config.Linux2-gcc-2.95 -finline-functions
LDFLAGS += -L../vxl/lib -lvnl_algo -lvnl -lvcl -lnetlib
#LDFLAGS+=-L../vxl-1.0-beta/vxl/lib/Linux2-gcc-2.95 -L../vxl-1.0-beta/v3p/lib/Linux2-gcc-2.95 -Wl,-rpath,../vxl-1.0-beta/v3p/lib/Linux2-gcc-2.95:../vxl-1.0-beta/vxl/lib/Linux2-gcc-2.95 -lvnl-algo -lnetlib -lvnl -lvcl -lstdc++ -ldl -rdynamic -lm -lc
endif

CVODE=yes
ifeq ($(CVODE),yes)
LOCALCFLAGS+=-DCVODE
LOCALCFLAGS+=-I/usr/include/cvode -I/usr/include/nvector -I/usr/include/sundials
endif

SOURCES = main.cpp Community.cpp Communicator.cpp Parameters.cpp \
	OutputController.cpp NodeOutputController.cpp Simulation.cpp \
	Integrator.cpp NRIntegrator.cpp CVIntegrator.cpp Iterator.cpp \
	Indexing.cpp util.cpp genrand2.cpp rand.cpp Algebra.cpp \
	Displays/PopulationDisplay.cpp Displays/DotDisplay.cpp \
	Displays/EvolutionaryTreeDisplay.cpp \
	Displays/ConstrainedPhenoDisplay2d.cpp \
	Displays/TimeSeriesDisplay.cpp
OBJS = $(SOURCES:.cpp=.o) CVodeWithStoppingCondition.o \
	../libexecstream/exec-stream.o
#OBJS = $(SOURCES:.cpp=.o) l_cvode.o ../libexecstream/exec-stream.o
#OBJS = $(SOURCES:.cpp=.o) ../libexecstream/exec-stream.o
THIS = GNUmakefile

ARCHIVE = libadap-dyn2.a

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

