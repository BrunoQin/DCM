export

SHELL = /bin/sh

ARCH  = LINUX

srcdir = .
top_srcdir = .

prefix = .
exec_prefix = ${prefix}

bindir = ${exec_prefix}/bin
sbindir = ${exec_prefix}/sbin
libexecdir = ${exec_prefix}/libexec
datadir = ${prefix}/share
sysconfdir = ${prefix}/etc
libdir = ${exec_prefix}/lib
includedir = ${prefix}/include
oldincludedir = /usr/include
infodir = ${prefix}/info
mandir = ${prefix}/man

sharedstatedir = ${prefix}/com
localstatedir = ${prefix}/var

program_transform_name = s,x,x,

MPIROOT        = /home/jywang/intel/oneapi/mpi/2021.1-beta10
MPI_LIB        = -L$(MPIROOT)/lib -L$(MPIROOT)/lib/debug -lmpi
MPI_INCLUDE    = -I$(MPIROOT)/include

NETCDFROOT     = /home/jywang/netcdf4
NETCDF_LIB     = -L$(NETCDFROOT)/netcdfc/lib -lnetcdf -L$(NETCDFROOT)/netcdf_fortran/lib -lnetcdff -L/home/jywang/FTA/install/lib  -ltorch_wrapper_F
NETCDF_INCLUDE = -I$(NETCDFROOT)/netcdfc/include -I$(NETCDFROOT)/netcdf_fortran/include

LIBS     = $(NETCDF_LIB) $(MPI_LIB)

MODOPT   = -I
MODULES  = ../modules

INCLUDE  = ../include
INCLUDES = $(MODOPT)$(MODULES) -I$(INCLUDE) $(NETCDF_INCLUDE) $(MPI_INCLUDE)

F90      = mpiifort -mkl -qopenmp 
FC       = mpiifort -mkl -qopenmp
CC       = gcc
CPP      = gcc -E
AR       = ar
AS       = as

DEFS     = -DHAVE_CONFIG_H

CFLAGS   = -I../config -O -DNAGf90Fortran
FFLAGS   = -O3 -mp1 -heap-arrays -convert big_endian  -r8
F90FLAGS = $(INCLUDES) -cpp
CPPFLAGS =
ARFLAGS  = crv
LDFLAGS  = -cpp

SRCDIRS = modules src

all:
	@for DIR in $(SRCDIRS) ;\
	  do \
	    back=`pwd`; \
	    cd $$DIR ;\
	    $(MAKE) ; status=$$? ; \
	    if [ $$status != 0 ] ; then \
	      echo "Exit status from make was $$status" ; exit $$status ; \
	    fi ; \
	    cd $$back ; \
	  done

clean:
	@for DIR in $(SRCDIRS) ;\
	  do \
	  (cd $$DIR ;\
	  $(MAKE) clean ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	  done
	-rm -f config.cache
	-rm -f lib/*.a bin/elbom
	-rm -f html/[a-z]*

tar:
	@tarfile=../echam5.f90.`date +%y%m%d`.taz ; gtar zcvf $$tarfile \
	`ls */*.f90 */*.[fhc] */*inc */Makefile Makefile.in Makefile run/hjob*`

index:
	-rm -f html/[a-z]*
	util/f2html.pl -f util/fgenrc -d html \
          modules src
