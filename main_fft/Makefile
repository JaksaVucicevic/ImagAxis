# C++ compiler
C++ = icpc
# SIAM path
SP = ../source
# executable path
RP = ~/RUN_SIAM/run_fft
FLAGS =#-g # -static
LIBS =# use this if needed 

all : main.o $(SP)/FFT.o $(SP)/IAGRID.o $(SP)/GRID.o $(SP)/routines.o
	$(C++) $(FLAGS) -o $(RP)/main $(LIBS) main.o  $(SP)/FFT.o $(SP)/IAGRID.o $(SP)/GRID.o $(SP)/routines.o 

# main program
main.o : main.cpp $(SP)/routines.h $(SP)/routines.o
	$(C++) $(FLAGS) -c -o $@ main.cpp

# Fast fourier transform with corrections specific for green's functions
$(SP)/FFT.o : $(SP)/FFT.cpp $(SP)/FFT.h
	$(C++) $(FLAGS) -c -o $@ $(SP)/FFT.cpp

# Grid utility for initializing and keeping omega grids and all grid dependent routines
$(SP)/IAGRID.o : $(SP)/IAGRID.cpp $(SP)/IAGRID.h 
	$(C++) $(FLAGS) -c -o $@ $(SP)/IAGRID.cpp

# Grid utility for initializing and keeping omega grids and all grid dependent routines
$(SP)/GRID.o : $(SP)/GRID.cpp $(SP)/GRID.h 
	$(C++) $(FLAGS) -c -o $@ $(SP)/GRID.cpp

# contains some constants and useful numerical routines
$(SP)/routines.o : $(SP)/routines.cpp $(SP)/routines.h 
	$(C++) $(FLAGS) -c -o $@ $(SP)/routines.cpp

# clean all object and exec files
clean :
	rm -f $(RP)/main main.o $(SP)/*.o
