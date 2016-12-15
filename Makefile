F90 = gfortran #CHOOSE YOUR FORTRAN COMPILER HERE

FFLAGS = -O3   #OPTIONAL: OPTIMIZATION FLAG
##########################
# Object Files for build #
##########################

OBJS = \
example.o \
RKHS.o \

all: example.x clean

example.x : $(OBJS)
	 ${F90}  -o $@ $(OBJS)

#######################################
# Object dependencies and compilation #
#######################################
example.o : src/example.f90 \
RKHS.o
	$(F90) -c $(FFLAGS) -o $@ src/example.f90

RKHS.o : src/RKHS.f90
	$(F90) -c $(FFLAGS) -o $@ src/RKHS.f90

.PHONY: clean
clean: 
	rm *.mod *.o


