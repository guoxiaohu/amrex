# AMREX_HOME defines the directory in which we will find all the AMReX code
AMREX_HOME = ../../../

DEBUG     = FALSE
USE_MPI   = FALSE
COMP      = gnu
DIM       = 3
PRECISION = DOUBLE
EBASE     = cv_analytic_fp

CVODE_LIB_DIR ?= ${CVODE_LIB}
USE_SUNDIALS_3x4x = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
