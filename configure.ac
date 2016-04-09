# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.61])
AC_INIT([phase_estimation], [1.0])
m4_include([m4/ax_mpi_options.m4])
m4_include([m4/ax_mpi_tests.m4])
AC_CONFIG_SRCDIR([src/main.cpp])
AC_CONFIG_HEADERS([config.h])
# Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_ARG_ENABLE([whatever],
[AS_HELP_STRING([--enable-openmp], [Use OpenMP (default: disabled)])])
if test "$enable_openmp" == yes; then
AC_OPENMP
else
enable_openmp=no
fi
CXXFLAGS="${OPENMP_CFLAGS}"
#find out what version we are running
ARCH=`uname -m`
if [[ $ARCH == "x86_64" ]];
then
SUFFIX="64"
else
SUFFIX=""
fi

string=`$CXX --version`
if [[[ $string == *"ICC"* ]]]
then  
  CXXFLAGS="${CXXFLAGS} -fast"
else
  AS_IF([test "x$GXX" = "xyes"],[CXXFLAGS="${CXXFLAGS} -Ofast -march=native -Wall"])
fi

#Setup MPI Paths
# ------------------------------------------------------------------------------
AC_ARG_WITH([mpi],
[ --with-mpi=PATH prefix where MPI is installed])
mpi_enabled=no
if test x"$with_mpi" != x"no" ; then
AX_MPI_OPTIONS
if test x"$MPI_CXX" != x"none"; then
AX_MPI_TESTS
CXX=$MPI_CXX
mpi_enabled=yes
else
AC_MSG_ERROR(MPI is necessary to compile)
fi
else
AC_MSG_ERROR(MPI is necessary to compile)
fi


# Setup CUDA paths
# ------------------------------------------------------------------------------
CUDA_PATH=/usr/local/cuda
AC_ARG_WITH([cuda],
   [  --with-cuda=PATH    prefix where cuda is installed [default=/usr/local/cuda]])

if test x"$with_cuda" != x"no" ; then

  if test -n "$with_cuda"
  then
     CUDA_PATH=$with_cuda
  fi

  AC_CHECK_PROG(NVCC_CHECK,nvcc,yes,no,$CUDA_PATH/bin)
  if test x"$NVCC_CHECK" = x"yes" ; then
     cuda_enabled=yes
     CUDA_INC="-I$CUDA_PATH/include"
     CUDA_CFLAGS="-use_fast_math -Xcompiler \"${CXXFLAGS}\""
     CUDA_LIBS="-lcudart -lcurand"
     CUDA_LDFLAGS="-L$CUDA_PATH/lib$SUFFIX"
     NVCC="$CUDA_PATH/bin/nvcc"
     AC_SUBST(CUDA_INC)
     AC_SUBST(CUDA_CFLAGS)
     AC_SUBST(CUDA_LIBS)
     AC_SUBST(CUDA_LDFLAGS)
     AC_SUBST(NVCC)
     AC_DEFINE([CUDA], 1, [cuda enabled])
  elif test -n "$with_cuda" ; then
     echo "---------------------------------------"
     echo "Unable to find CUDA in $with_cuda."
     echo "Building a version without CUDA.    "
     echo "---------------------------------------"
     cuda_enabled=no
  else 
     cuda_enabled=no
  fi
else
  cuda_enabled=no  
fi
AM_CONDITIONAL([HAVE_CUDA], [test x"$cuda_enabled" = x"yes"])

# Set up VSL
# ------------------------------------------------------------------------------
AC_ARG_WITH([vsl],
   [  --with-vsl=PATH    prefix where Intel MKL/VSL is installed])
if test x"$with_vsl" != x"no" ; then

  if test -n "$with_vsl"
  then
    VSL_PATH=$with_vsl
    vsl_enabled=yes
    VSL_LIBDIR="-L$VSL_PATH/lib/intel$SUFFIX"
    VSL_INC="-I$VSL_PATH/include"
    VSL_LIBS="-lmkl_intel_ilp64 -lmkl_sequential -lmkl_core"
    AC_SUBST(VSL_INC)
    AC_SUBST(VSL_LIBS)
    AC_SUBST(VSL_LIBDIR)
    AC_DEFINE([VSL], 1, [vsl enabled])
  else
    vsl_enabled=no
  fi
fi
AM_CONDITIONAL([HAVE_VSL], [test x"$vsl_enabled" = x"yes"])

AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_OUTPUT
echo \
"-------------------------------------------------
${PACKAGE_NAME} Version ${PACKAGE_VERSION}
Prefix: ${prefix}.
Compiler: ${CXX} ${CXXFLAGS} ${MPI_INC} ${MPI_LIBDIR} ${MPI_LIBS} ${CUDA_INC} ${VSL_INC} ${VSL_LIBDIR} ${VSL_LIBS}
Package features:
OpenMP enabled: ${enable_openmp}
MPI enabled: ${mpi_enabled}
CUDA enabled: ${cuda_enabled}
VSL enabled: ${vsl_enabled}
Now type 'make @<:@<target>@:>@'
where the optional <target> is:
all - build all binaries
install - install everything
test - build tests
--------------------------------------------------"
