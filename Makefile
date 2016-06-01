# Makefile

F90 = gfortran
F90FLAGS = -ffixed-line-length-0 -ffree-form -fopenmp -O2
F90LIBS = -llapack -lblas -lgfortran 
F90INCLUDE = -L/home/bminchew/local/lib 
TARGETS = ./fourDvel
prefix = /home/bminchew/Utilities
PREFIX = $(prefix)

all: $(TARGETS)

./fourDvel: ./fourDvel.f90
	 $(F90) -o $@ $< $(F90FLAGS) $(F90INCLUDE) $(F90LIBS)

.PHONY: install
install:
	 mkdir -p $(PREFIX)
	 install $(TARGETS) $(PREFIX)

.PHONY: uninstall
uninstall:
	 rm -f $(PREFIX)/$(TARGETS)

.PHONY: clean
clean:
	 rm -f $(TARGETS)

again: clean all
 
