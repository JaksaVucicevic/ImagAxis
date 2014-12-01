# C++ compiler
C++ = icpc
# SIAM path
SP = ../source
# executable path
RP = ../runinsulator
FLAGS =#-g # -O2 -static
LIBS =# use this if needed 

all : main.o $(SP)/CHM.o $(SP)/IDMFT.o $(SP)/SIAM.o $(SP)/GRID.o $(SP)/Broyden.o $(SP)/nrutil.o $(SP)/routines.o 
	$(C++) $(FLAGS) -o $(RP)/main $(LIBS) main.o  $(SP)/CHM.o $(SP)/IDMFT.o $(SP)/SIAM.o $(SP)/GRID.o $(SP)/Broyden.o $(SP)/nrutil.o $(SP)/routines.o 

# main program
main.o : main.cpp $(SP)/Mixer.h $(SP)/CHM.o
	$(C++) $(FLAGS) -c -o $@ main.cpp

# Clean Hubbard model solver. Uses SIAM (below)
$(SP)/CHM.o : $(SP)/CHM.cpp $(SP)/CHM.h $(SP)/SIAM.o $(SP)/GRID.o $(SP)/Mixer.h
	$(C++) $(FLAGS) -c -o $@ $(SP)/CHM.cpp

# Clean Hubbard model solver. Uses SIAM (below)
$(SP)/IDMFT.o : $(SP)/IDMFT.cpp $(SP)/IDMFT.h $(SP)/SIAM.o $(SP)/GRID.o $(SP)/Mixer.h
	$(C++) $(FLAGS) -c -o $@ $(SP)/IDMFT.cpp

# IPT/MPT Single impurity Anderson model solver
$(SP)/SIAM.o : $(SP)/SIAM.cpp $(SP)/SIAM.h $(SP)/GRID.o $(SP)/Broyden.o
	$(C++) $(FLAGS) -c -o $@ $(SP)/SIAM.cpp

# Grid utility for initializing and keeping omega grids and all grid dependent routines
$(SP)/GRID.o : $(SP)/GRID.cpp $(SP)/GRID.h 
	$(C++) $(FLAGS) -c -o $@ $(SP)/GRID.cpp

# Broyden solver for systems of equations
$(SP)/Broyden.o : $(SP)/Broyden.h $(SP)/Broyden.cpp $(SP)/routines.o
	$(C++) $(FLAGS) -c -o $@ $(SP)/Broyden.cpp

# numerical routines from NumRec
$(SP)/nrutil.o : $(SP)/nrutil.h $(SP)/nrutil.c
	$(C++) $(FLAGS) -c -o $@ $(SP)/nrutil.c

# contains some constants and useful numerical routines
$(SP)/routines.o : $(SP)/routines.cpp $(SP)/routines.h 
	$(C++) $(FLAGS) -c -o $@ $(SP)/routines.cpp
# clean all object and exec files
clean :
	rm -f $(RP)/main main.o $(SP)/*.o
