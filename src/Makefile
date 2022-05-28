#Replace with your FORTRAN 95 compiler command.
F77 = gfortran

#Place any local compiler options here.
FFLAGS = -std=f95 -pedantic -Wall -fbounds-check
#Any special compiler flags required for fixed format Fortran source.
F77FLAGS = 

#By default, use the BLAS and/or LAPACK routines provided.
BLAS = -L. blas.o
LAPACK = -L. lapack.o
BLASDP = blas.o
LAPACKDP = lapack.o

#Uncomment the following to use your local system BLAS and/or LAPACK.
#BLAS = -lblas
#LAPACK = -llapack
#LAPACKDP = 
#BLASDP = 

#Most compilers produce lower case file names for the .mod files. 
REALMOD = real_precision.mod
LSMOD = linear_shepard_mod.mod
CSMOD = cshep2d_mod.mod
QS2MOD = qshep2d_mod.mod
QS3MOD = qshep3d_mod.mod
QSDATA = qshepmdata.mod
QSMMOD = qshepmd_mod.mod
SHMOD = sheppack.mod

#The Absoft FORTRAN compiler produces module files with capital letters
#in the file name whereas most other compilers use all lower case.
#Uncomment the following if using the Absoft compiler.
#REALMOD = REAL_PRECISION.mod
#LSMOD = LINEAR_SHEPARD_MOD.mod
#CSMOD = CSHEP2D_MOD.mod
#QS2MOD = QSHEP2D_MOD.mod
#QS3MOD = QSHEP3D_MOD.mod
#QSDATA = QSHEPMDATA.mod
#QSMMOD = QSHEPMD_MOD.mod
#SHMOD = SHEPPACK.mod

all: test_cshep2d test_qshep2d test_qshep3d test_qshepmd test_lshep
	./test_cshep2d
	./test_qshep2d
	./test_qshep3d
	./test_qshepmd
	./test_lshep

blas.o: blas.f
	$(F77) $(FFLAGS) $(F77FLAGS) -c blas.f

$(REALMOD): real_precision.f90
	$(F77) $(FFLAGS) -c real_precision.f90

blas_shep.o: blas_shep.f90
	$(F77) $(FFLAGS) -c blas_shep.f90

lapack.o: lapack.f
	$(F77) $(FFLAGS) $(F77FLAGS) -c lapack.f

$(LSMOD): linear_shepard.f90 $(REALMOD)
	$(F77) $(FFLAGS) -c linear_shepard.f90 

$(CSMOD): cshep2d.f90 $(REALMOD)
	$(F77) $(FFLAGS) -c cshep2d.f90

$(QS2MOD): qshep2d.f90 $(REALMOD)
	$(F77) $(FFLAGS) -c qshep2d.f90

$(QS3MOD): qshep3d.f90 $(REALMOD)
	$(F77) $(FFLAGS) -c qshep3d.f90

$(QSDATA): qshepmdata.f90 $(REALMOD)
	$(F77) $(FFLAGS) -c qshepmdata.f90

$(QSMMOD): qshepmd.f90 $(REALMOD) $(QSDATA)
	$(F77) $(FFLAGS) -c qshepmd.f90 

$(SHMOD): sheppack.f90 $(REALMOD) $(QSMMOD) \
	$(QSDATA) $(QS3MOD) $(QS2MOD) $(CSMOD) \
	$(LSMOD)
	$(F77) $(FFLAGS) -c sheppack.f90

install: $(SHMOD)


test_cshep2d: cshep2d_test.f90 $(SHMOD) blas_shep.o
	$(F77) $(FFLAGS) -o test_cshep2d cshep2d_test.f90 sheppack.o \
	blas_shep.o cshep2d.o

test_qshep2d: qshep2d_test.f90 $(SHMOD) blas_shep.o
	$(F77) $(FFLAGS) -o test_qshep2d qshep2d_test.f90 sheppack.o \
	blas_shep.o qshep2d.o

test_qshep3d: qshep3d_test.f90 $(SHMOD) blas_shep.o
	$(F77) $(FFLAGS) -o test_qshep3d qshep3d_test.f90 sheppack.o \
	blas_shep.o qshep3d.o

test_qshepmd: qshepmd_test.f90 $(SHMOD) blas_shep.o $(QSDATA)
	$(F77) $(FFLAGS) -o test_qshepmd qshepmd_test.f90 sheppack.o \
	blas_shep.o qshepmd.o qshepmdata.f90

test_lshep: linear_shepard_test.f90 $(SHMOD) $(BLASDP) $(LAPACKDP)
	$(F77) $(FFLAGS) -o test_lshep linear_shepard_test.f90 \
	sheppack.o linear_shepard.o $(LAPACK) $(BLAS) 

test_all: test_cshep2d test_qshep2d test_qshep3d test_qshepmd test_lshep
	

clean: 
	rm -f test_lshep test_cshep2d test_qshep2d test_qshep3d \
	test_qshepmd *.o *.mod

