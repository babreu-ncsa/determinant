# Fortran GCC compiler
FC=gfortran

# compilation flags 
FDFLAGS=-llapack

# source file
SOURCE=./det.f90

# binary name
EXEC=a.out

all:
	$(FC) $(FDFLAGS) $(SOURCE) -o $(EXEC)
	@echo -e "----- COMPILATION DONE -----"


clean:
	rm -r $(EXEC)

