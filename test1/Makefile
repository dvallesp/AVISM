NAME=mock_voids.x

#### This is an example for a local machine ####
ifeq ($(COMP),1)          
 FC=gfortran
 FFLAGS=-O3 -lgcc -mcmodel=medium -fopenmp -mieee-fp -ftree-vectorize -march=native
endif

#### This is an example for a local machine DEBUGGIN ####
ifeq ($(COMP),2)          
 FC=gfortran
 FFLAGS=-O1 -g -fcheck=all -fbounds-check -fbacktrace -fopenmp -mcmodel=medium
endif
##########################################################################

EXEC=$(NAME)

OBJ=mock_voids.o

$(NAME): $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o $(EXEC) $(LIBS)

mock_voids.o: mock_voids.f90
	$(FC) $(FFLAGS) $(INC) -c -o mock_voids.o mock_voids.f90

clean:
	 rm -f $(OBJ) $(EXEC)


