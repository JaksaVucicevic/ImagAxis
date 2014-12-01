# C++ compiler
C++ = /opt/intel/bin/icpc
# SIAM path
SP = ../source
# executable path
OP =~/SIAM/obj/
RP =~/SIAM/bin/
main=main_hirschfye

FLAGS = -openmp #-g # -static
LIBS =# use this if needed 

all : $(OP)/$(main).o $(OP)/IAHFCHM.o $(OP)/HirschFye.o $(OP)/FFT.o $(OP)/IAGRID.o $(OP)/GRID.o $(OP)/Resampler.o $(OP)/routines.o $(SP)/arrayInitializers.h
	$(C++) $(FLAGS) -o $(RP)/$(main) $(LIBS) $(OP)/$(main).o $(OP)/HirschFye.o $(OP)/IAHFCHM.o $(OP)/FFT.o $(OP)/IAGRID.o $(OP)/GRID.o $(OP)/Resampler.o $(OP)/routines.o 

# main program
$(OP)/$(main).o : $(main).cpp $(SP)/arrayInitializers.h
	$(C++) $(FLAGS) -c -o $@ $(main).cpp

$(OP)/HirschFye.o : $(SP)/HirschFye.cpp $(SP)/HirschFye.h $(SP)/arrayInitializers.h
	$(C++) $(FLAGS) -c -o $@ $(SP)/HirschFye.cpp

# Grid utility for initializing and keeping omega grids and all grid dependent routines
$(OP)/IAHFCHM.o : $(SP)/IAHFCHM.cpp $(SP)/IAHFCHM.h
	$(C++) $(FLAGS) -c -o $@ $(SP)/IAHFCHM.cpp

# Grid utility for initializing and keeping omega grids and all grid dependent routines
$(OP)/IAGRID.o : $(SP)/IAGRID.cpp $(SP)/IAGRID.h 
	$(C++) $(FLAGS) -c -o $@ $(SP)/IAGRID.cpp

# Grid utility for initializing and keeping omega grids and all grid dependent routines
$(OP)/FFT.o : $(SP)/FFT.cpp $(SP)/FFT.h
	$(C++) $(FLAGS) -c -o $@ $(SP)/FFT.cpp

# Grid utility for initializing and keeping omega grids and all grid dependent routines
$(OP)/GRID.o : $(SP)/GRID.cpp $(SP)/GRID.h 
	$(C++) $(FLAGS) -c -o $@ $(SP)/GRID.cpp

# contains some constants and useful numerical routines
$(OP)/Resampler.o : $(SP)/Resampler.cpp $(SP)/Resampler.h $(SP)/routines.h
	$(C++) $(FLAGS) -c -o $@ $(SP)/Resampler.cpp


# contains some constants and useful numerical routines
$(OP)/routines.o : $(SP)/routines.cpp $(SP)/routines.h 
	$(C++) $(FLAGS) -c -o $@ $(SP)/routines.cpp
# clean all object and exec files
clean :
	rm -f $(RP)/main $(OP)/main.o $(OP)/*.o
