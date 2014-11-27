#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.61])
AC_INIT([DE_MPI_loss_100], [0.1])
m4_include([m4/ax_mpi_options.m4])
m4_include([m4/ax_mpi_tests.m4])
AC_CONFIG_SRCDIR([src/main.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
CFLAGS="-O3 -fPIC"
AC_PROG_CC
AC_PROG_INSTALL

AC_ARG_ENABLE([whatever],
  [AS_HELP_STRING([--enable-openmp], [Use OpenMP (default: disabled)])])

if test "$enable_openmp" == yes; then
  AC_OPENMP
else
  enable_openmp=no
fi
CFLAGS="${CFLAGS} ${OPENMP_CFLAGS}"

#find out what version we are running
ARCH=`uname -m`
if [[ $ARCH == "x86_64" ]];
then
  SUFFIX="64"
else
  SUFFIX=""
fi

#Setup MPI Paths
# ------------------------------------------------------------------------------
AC_ARG_WITH([mpi],
   [  --with-mpi=PATH    prefix where MPI is installed])

mpi_enabled=no
if test x"$with_mpi" != x"no" ; then
    AX_MPI_OPTIONS
    if test x"$MPI_CXX" != x"none"; then
      AX_MPI_TESTS
      CC=$MPI_CXX
      mpi_enabled=yes
    else
      AC_MSG_ERROR(MPI is necessary to compile)
    fi
else
  AC_MSG_ERROR(MPI is necessary to compile)
fi

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT

echo \
"-------------------------------------------------

 ${PACKAGE_NAME} Version ${PACKAGE_VERSION}

 Prefix: ${prefix}.
 Compiler: ${CC} ${CFLAGS} ${MPI_INC} ${MPI_LIBDIR} ${MPI_LIBS}

 Package features:
   OpenMP enabled: ${enable_openmp}
   MPI enabled: ${mpi_enabled}
   
 Now type 'make @<:@<target>@:>@'
   where the optional <target> is:
     all                - build all binaries
     install            - install everything

--------------------------------------------------"
