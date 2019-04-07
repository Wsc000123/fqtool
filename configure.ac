# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.69])
AC_INIT([fqtool], [0.0.0], [wljlinksly@gmail.com])
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([src/main.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h inttypes.h limits.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_SUBST(AM_CXXFLAGS, "-std=c++11 -g -Wall -O3 -Wextra -Wno-unknown-pragmas -Wno-unused-parameter -Wno-sign-compare")
AC_SUBST(LDFLAGS, "$LDFLAGS -pthread -lz")
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_MODE_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_MMAP
AC_CHECK_FUNCS([getcwd memset mkdir munmap realpath])

AC_CONFIG_FILES([Makefile
                 src/Makefile])

AC_OUTPUT
