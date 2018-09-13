# programming environment
COMPILER = /usr/lib64/mpi/gcc/mvapich2/bin/mpic++
INCLUDE      := -I$(HOME)/code/LATfield2
LIBFFTW := -L$(HOME)/bin/fftw3/lib -I$(HOME)/bin/fftw3/include
LIBHDF5 := -L$(HOME)/code/hdf5/lib -I$(HOME)/code/hdf5/include
LIB          := -lfftw3 -lm -lhdf5 -lgsl -lgslcblas $(LIBFFTW) $(LIBHDF5)
#-lfftw3 -lm -lhdf5 -lgsl -lgslcblas

# target and source
EXEC         := gevolution
SOURCE       := main.cpp
HEADERS      := $(wildcard *.hpp)

# mandatory compiler settings (LATfield2)
DLATFIELD2   := -DFFT3D -DHDF5

# optional compiler settings (LATfield2)
#DLATFIELD2   += -DH5_HAVE_PARALLEL
#DLATFIELD2   += -DEXTERNAL_IO

# optional compiler settings (gevolution)
DGEVOLUTION  := -DPHINONLINEAR
DGEVOLUTION  += -DBENCHMARK
#DGEVOLUTION  += -DCHECK_B
#DGEVOLUTION  += -DHAVE_CLASS # requires OPT -fopenmp and LIB -lclass

# further compiler options
OPT          := -O3 -std=c++11

$(EXEC): $(SOURCE) $(HEADERS) makefile
	$(COMPILER) $< -o $@ $(OPT) $(DLATFIELD2) $(DGEVOLUTION) $(INCLUDE) $(LIB)

