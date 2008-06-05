F90=ifort -warn all
F90OPT=

AFNLPATH=/home/alberto/fortran/lib/
AFNLNAME=f90

VPATH=src
SRCDIR=$(VPATH)


all: lib

genetics.o: genetics.f90
	$(F90) -c $^ -o $@ -I$(AFNLPATH) -L$(AFNLPATH) -lf90

evolution.o: evolution.f90
	$(F90) -c $^ -o $@ -I$(AFNLPATH) -L$(AFNLPATH) -lf90

lib: genetics.o evolution.o 
	ar rcs libga.a *.o
#	rm *.o

cleanall:
	rm *.o *.mod src/*~ libga.a
