# modules makefile

SRCS =	mo_nctools.f90 mo_mathtools.f90 mo_typedef.f90 mo_tautsubmodel.f90 mo_prepostproc.f90 \
        mo_cf23model.f90 mo_oceansolver.f90 

OBJS := $(SRCS:.f90=.o)

all: $(OBJS)

clean:
	rm -f $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

################################################################################################

mo_nctools.o:
mo_mathtools.o:
mo_typedef.o: mo_nctools.o mo_mathtools.o
mo_tautsubmodel.o: mo_nctools.o mo_mathtools.o mo_typedef.o
mo_cf23model.o:mo_nctools.o mo_mathtools.o mo_typedef.o
mo_oceansolver.o: mo_nctools.o mo_mathtools.o mo_typedef.o mo_tautsubmodel.o
mo_prepostproc.o: mo_nctools.o mo_mathtools.o mo_typedef.o mo_oceansolver.o

#################################################################################################
