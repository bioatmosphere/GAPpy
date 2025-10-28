PROG =	UVAFME.exe

SRCS =	Climate.f90 Constants.f90 csv_file.f90                             \
	GenusGroups.f90 Input.f90 IO_Utils.f90 lists.f90                       \
	Model.f90 Output.f90 Parameters.f90 Plot.f90 Random.f90 Site.f90       \
	Sitelist.f90 Soil.f90 Species.f90 Tree.f90                             \
	Utilities.f90 UVAFME.f90 vararray.f90 FileUtils.F90

OBJS =	Climate.o Constants.o csv_file.o                                   \
	GenusGroups.o Input.o IO_Utils.o lists.o Model.o Output.o              \
	Parameters.o Plot.o Random.o Site.o Sitelist.o Soil.o Species.o        \
	Tree.o Utilities.o UVAFME.o vararray.o FileUtils.o

LIBS =	

F90 =  ifort
DBG      =-CB -g -fp-model strict

#F90=pgf90
#DBG =-g -Mbounds

OPT      = -O
F90FLAGS = $(DBG) $(OPT)
LDFLAGS = 
all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

.PHONY: clean
clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90 .F90 .f95

.f90.o .f95.o .F90.o:
	$(F90) $(F90FLAGS) -c $<

Climate.o: Constants.o Random.o
csv_file.o: csv_file_1d.f90 csv_file_2d.f90
GenusGroups.o: Site.o vararray.o
Input.o: Constants.o IO_Utils.o Parameters.o Site.o
IO_Utils.o: Constants.o FileUtils.o Parameters.o csv_file.o
lists.o: Species.o
Model.o: Climate.o Constants.o Input.o Parameters.o Random.o Site.o Soil.o \
	Species.o Tree.o
Output.o: Constants.o IO_Utils.o Species.o Tree.o Utilities.o csv_file.o
Parameters.o: Constants.o
Plot.o: csv_file.o Constants.o Species.o Tree.o
Site.o: Constants.o Plot.o Soil.o Species.o Utilities.o lists.o
Sitelist.o: Constants.o Input.o Site.o
Soil.o: csv_file.o Constants.o FileUtils.o Parameters.o
Species.o: csv_file.o Constants.o FileUtils.o
Tree.o: csv_file.o Constants.o FileUtils.o Random.o Species.o
Utilities.o: Constants.o
UVAFME.o: Constants.o Input.o Model.o Output.o Parameters.o Site.o Sitelist.o \
	Species.o
vararray.o: Constants.o
FileUtils.o: Constants.o Parameters.o
csv_file.o: csv_file_1d.f90 csv_file_1d.f90 csv_file_1d.f90 csv_file_1d.f90 \
	csv_file_2d.f90 csv_file_2d.f90 csv_file_2d.f90 csv_file_2d.f90
