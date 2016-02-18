###############################################################
#                                                             #
###############################################################

#-------------------------------------------------------------------------  
# OpenMP compiler switch
 OMP = -fopenmp

#OPT += -DUSE_MPI

 OPT += -DGRID    #gridding visibilities

#------------------------------------------------------------------------- 

SUP_INCL = -I. 
LIB_OPT = 
OPTIMIZE = -O3 -g 

ifeq (MPI,$(findstring MPI,$(OPT)))
  CC  =  mpiCC -g 
 else
  CC  = g++-4.9
 endif

OPTIONS = $(OPTIMIZE) $(OPT)
EXEC = RadioLensfit.x   
 
OBJS  =	 RadioLensfit.o read_coordinates.o distributions.o galaxy_visibilities.o  evaluate_uv_grid.o generate_random_values.o likelihood.o random_gaussian.o marginalise_r.o  
 
INCL   = *.h Makefile
LIB_OPT = -lgsl -lgslcblas -lm

CPPFLAGS = $(OPTIONS) $(SUP_INCL)  $(OMP)

LIBS   = $(LIB_OPT) $(OMP)

.SUFFIXES: .o .cc .cxx .cpp .cu

.cc.o:
	$(CC) -c $(CPPFLAGS) -o "$@" "$<"

.cxx.o:
	$(CC) -c $(CPPFLAGS) -o "$@" "$<"

.cpp.o:
	$(CC) -c $(CPPFLAGS) -o "$@" "$<"


$(EXEC): $(OBJS)
	$(CC)  $(OBJS)  $(OPTIONS) $(LIBS) -o $(EXEC)

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS)

realclean: clean
	rm -f $(EXEC)

