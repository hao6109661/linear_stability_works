#Name of executable
noinst_PROGRAMS=hao

#Sources for the executable
hao_SOURCES = hao.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
hao_LDADD = -L@libdir@ -lbeam -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.
AM_CPPFLAGS += -I@includedir@  

