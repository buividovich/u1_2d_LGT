SRC = LGTDescriptor.cpp LGT_QMC.cpp LGT_ExDiag.cpp linalg.cpp
HDR = $(SRC:.cpp=.hpp)
OBJ = $(SRC:.cpp=.o)

HDR += ./ansi_colors.hpp ./observables.hpp

CC = g++ -std=c++20 -O2 -I./ 
CC += -DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP

#This trick is needed for Barkla, because the module system on Barkla only sets LD_LIBRARY_PATH
export LIBRARY_PATH := $(LD_LIBRARY_PATH):$(LIBRARY_PATH)

CC += -fmax-errors=1 -fopenmp

LIB = -lm -lopenblas -lstdc++ -lfftw3 -lboost_program_options

CC += -I /opt/OpenBLAS/include/ -L /opt/OpenBLAS/lib/

%.o: %.cpp
	$(CC) -c $< $(LIB) -o $@


qmc: qmc_main.cpp $(HDR) $(OBJ)
	$(CC) ./$< $(OBJ) $(LIB) -o ./$@

exdiag: exdiag_main.cpp $(HDR) $(OBJ)
	$(CC) ./$< $(OBJ) $(LIB) -o ./$@

clean:
	rm -f -v ./qmc ./exdiag ./*.o