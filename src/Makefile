# src/Makefile.  Generated from Makefile.in by configure.

# Package-related substitution variables
package	= phase_estimation
version	= 0.1
tarname	= phase_estimation
distdir	= $(tarname)-$(version)

# Prefix-related substitution variables
prefix	 = /usr/local
exec_prefix    = ${prefix}
bindir	 = ${exec_prefix}/bin
libdir	 = ${exec_prefix}/lib

# Tool-related substitution variables
CXX		         = mpic++
CXXFLAGS       =  -Ofast -march=native -Wall
LIBS	         = 
DEFS           = -DHAVE_CONFIG_H
INSTALL	       = /usr/bin/install -c
INSTALL_DATA   = ${INSTALL} -m 644
INSTALL_PROGRAM= ${INSTALL}
INSTALL_SCRIPT = ${INSTALL}
MPI_INC        = -I/global/software/openmpi/openmpi-1.6.5-intel/include
MPI_LIBDIR     = -L/global/software/openmpi/openmpi-1.6.5-intel/lib
MPI_LIBS       = -lmpi
VSL_INC        = 
VSL_LIBDIR     = 
VSL_LIBS       = 
CUDA_INC       = 
CUDA_CFLAGS    = 
CUDA_LIBS      = 
CUDA_LDFLAGS   = 
NVCC       	   = 

# VPATH-related substitution variables
srcdir	 = .


OBJS=main.o candidate.o phase_loss_opt.o io.o problem.o mpi_optalg.o mpi_pso.o mpi_de.o candidate.o rng.o aux_functions.o

ifdef CUDA_LIBS
	OBJS+=rng_gpu.cu.co
else 
  ifdef VSL_LIBS
		OBJS+=rng_vsl.o
	endif
endif

all: phase_estimation
	
phase_estimation: $(OBJS)
	$(CXX) $(DEFS) $(CXXFLAGS) $(MPI_LIBDIR) $(VSL_LIBDIR) $(CUDA_LDFLAGS) -o $@ $^ $(LIBS) $(MPI_LIBS) $(VSL_LIBS) $(CUDA_LIBS)

%.o: %.cpp
		$(CXX) $(DEFS) $(CXXFLAGS) $(MPI_INC) $(CUDA_INC) $(VSL_INC) -I$(srcdir) -I.. -o $@ -c $(srcdir)/$<

%.cu.co: %.cu
		$(NVCC) $(DEFS) $(MPI_INC) $(CUDA_INC) $(CUDA_CFLAGS) -I$(srcdir) -I.. -o $@ -c $(srcdir)/$<

clean:
	-rm -f phase_estimation $(OBJS) 1>/dev/null

install:
	$(INSTALL) -d $(DESTDIR)$(bindir)
	$(INSTALL_PROGRAM) -m 0755 phase_estimation \
	 $(DESTDIR)$(bindir)
	 
uninstall:
	-rm $(DESTDIR)$(bindir)/phase_estimation &>/dev/null

Makefile: Makefile.in ../config.status
	cd .. && ./config.status $@

../config.status: ../configure
	cd .. && ./config.status --recheck

.PHONY: all clean install uninstall
