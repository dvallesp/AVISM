#### This is an example for a local machine ####
ifeq ($(COMP),1)          
 FC=gfortran
 FFLAGS=-O2 -mcmodel=medium -fopenmp -mieee-fp -ftree-vectorize -march=native -funroll-loops
 LIBS=
 INC=
endif

#### This is an example for a local machine DEBUGGING ####
ifeq ($(COMP),2)          
 FC=gfortran
 FFLAGS=-O1 -g -mcmodel=medium -fopenmp -mieee-fp -ftree-vectorize -march=native -fcheck=all -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow
 LIBS=
 INC=
endif

##########################################################################

# DIRECTORIES
SRCDIR     := src
BINDIR     := bin

# -J is the directory where the .mod files are stored (BIN)
FFLAGS  += -J $(BINDIR) 

# EXECUTABLE
EXEC=voids.x

# OBJECTS
OBJ=commondata.o kdtree.o particles.o avism.o

# COMPILATION
$(EXEC): $(addprefix $(BINDIR)/, $(OBJ))
	$(FC) $(FFLAGS) $(addprefix $(BINDIR)/, $(OBJ)) -o $(EXEC) $(LIBS)

# Rule for commondata.o
$(BINDIR)/commondata.o: $(SRCDIR)/commondata.f90
	$(FC) $(FFLAGS) $(INC) -c -o $(BINDIR)/commondata.o $(SRCDIR)/commondata.f90

# Rule for kdtree.o
$(BINDIR)/kdtree.o: $(SRCDIR)/kdtree.f90
	$(FC) $(FFLAGS) $(INC) -c -o $(BINDIR)/kdtree.o $(SRCDIR)/kdtree.f90

# Rule for particles.o
$(BINDIR)/particles.o: $(SRCDIR)/particles.f90
	$(FC) $(FFLAGS) $(INC) -c -o $(BINDIR)/particles.o $(SRCDIR)/particles.f90

# Rule for avism.o 
$(BINDIR)/avism.o: $(SRCDIR)/avism.f90 $(BINDIR)/commondata.o $(BINDIR)/kdtree.o $(BINDIR)/particles.o
	$(FC) $(FFLAGS) $(INC) -c -o $(BINDIR)/avism.o $(SRCDIR)/avism.f90

# CLEAN
clean:
	rm -f $(BINDIR)/*.o $(BINDIR)/*.mod $(EXEC) *.mod
