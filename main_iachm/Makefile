# C++ compiler
C++ = icpc
# SIAM path
SP = ../source
# executable path
RP = ~/RUN_SIAM/run_iachm
FLAGS =#-g # -static
LIBS =# use this if needed 

all : main.o $(SP)/IACHM.o $(SP)/IASIAM.o $(SP)/FFT.o $(SP)/Broyden.o $(SP)/IAGRID.o $(SP)/GRID.o $(SP)/routines.o $(SP)/nrutil.o
	$(C++) $(FLAGS) -o $(RP)/main $(LIBS) main.o $(SP)/IACHM.o $(SP)/IASIAM.o $(SP)/FFT.o $(SP)/Broyden.o $(SP)/IAGRID.o $(SP)/GRID.o $(SP)/routines.o $(SP)/nrutil.o

# main program
main.o : main.cpp $(SP)/IACHM.o $(SP)/IAGRID.o
	$(C++) $(FLAGS) -c -o $@ main.cpp

# Grid utility for initializing and keeping omega grids and all grid dependent routines
$(SP)/IACHM.o : $(SP)/IACHM.cpp $(SP)/IACHM.h
	$(C++) $(FLAGS) -c -o $@ $(SP)/IACHM.cpp

# Grid utility for initializing and keeping omega grids and all grid dependent routines
$(SP)/IASIAM.o : $(SP)/IASIAM.cpp $(SP)/IASIAM.h 
	$(C++) $(FLAGS) -c -o $@ $(SP)/IASIAM.cpp

# Grid utility for initializing and keeping omega grids and all grid dependent routines
$(SP)/IAGRID.o : $(SP)/IAGRID.cpp $(SP)/IAGRID.h 
	$(C++) $(FLAGS) -c -o $@ $(SP)/IAGRID.cpp

# Grid utility for initializing and keeping omega grids and all grid dependent routines
$(SP)/FFT.o : $(SP)/FFT.cpp $(SP)/FFT.h
	$(C++) $(FLAGS) -c -o $@ $(SP)/FFT.cpp

# Broyden solver for systems of equations
$(SP)/Broyden.o : $(SP)/Broyden.h $(SP)/Broyden.cpp $(SP)/routines.o
	$(C++) $(FLAGS) -c -o $@ $(SP)/Broyden.cpp

# Grid utility for initializing and keeping omega grids and all grid dependent routines
$(SP)/GRID.o : $(SP)/GRID.cpp $(SP)/GRID.h 
	$(C++) $(FLAGS) -c -o $@ $(SP)/GRID.cpp

# contains some constants and useful numerical routines
$(SP)/routines.o : $(SP)/routines.cpp $(SP)/routines.h 
	$(C++) $(FLAGS) -c -o $@ $(SP)/routines.cpp

# numerical routines from NumRec
$(SP)/nrutil.o : $(SP)/nrutil.h $(SP)/nrutil.c
	$(C++) $(FLAGS) -c -o $@ $(SP)/nrutil.c

# clean all object and exec files
clean :
	rm -f $(RP)/main main.o $(SP)/*.o
