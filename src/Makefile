# src makefile

PROG =	../bin/elbom

SRCS =	ocean_master.f90 

OBJS := $(SRCS:.f90=.o)

MOBJS =	$(MODULES)/mo_nctools.o $(MODULES)/mo_mathtools.o \
	$(MODULES)/mo_typedef.o $(MODULES)/mo_tautsubmodel.o $(MODULES)/mo_prepostproc.o \
	$(MODULES)/mo_cf23model.o $(MODULES)/mo_oceansolver.o
	

all: $(PROG)

$(PROG): $(OBJS) $(MOBJS) 
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(MOBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

#####################################################################################################
ocean_master.o: $(MODULES)/mo_nctools.o $(MODULES)/mo_mathtools.o\
	        $(MODULES)/mo_typedef.o \
	        $(MODULES)/mo_tautsubmodel.o $(MODULES)/mo_prepostproc.o \
                $(MODULES)/mo_cf23model.o $(MODULES)/mo_oceansolver.o


