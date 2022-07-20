[![docs.rs](https://img.shields.io/badge/scientific%20documentation-read-green)](https://gerics.pages.hzdr.de/remo/remo-doc)
[![Documentation Status](https://readthedocs.org/projects/remo-user-guide/badge/?version=latest)](https://remo-user-guide.readthedocs.io/en/latest/?badge=latest)
[![pyremo](https://img.shields.io/badge/pyremo-docs-green)](https://pyremo.readthedocs.io)

# REMO2020 - A modular regional Earth System Model

## Quickstart

A more detailed Documentation will follow. For now try:

```
./setup -auto Remo2020 -objdir=Remo2020
```
This command will browse the `source` directory and link all source code files
required to compile Remo in the 2020 configuration in the `Remo2020` directory. 
It will also create a `Makefile` and automatically choose a header called `Makefile.h`
depending on your machine. If you want to change compiler options, you should do this
in the header file `Makefile.h`.

There are some more helpful options, e.g., to setup the model with debug compiler
options, you can choose:

```
./setup -auto Remo2020 -debug -objdir=Remo2020_debug
```

## Compilation and Running on Levante

REMO is preconfigured if compiled from a login node at [DKRZ HLRE-4 (Levante)](https://docs.dkrz.de/doc/levante/). Different combinations of Fortran and MPI compilers are possible and should be used for compilation and runscripts consistently. Note, that these examples are only snapshots and might change when new compiler versions are available. Note that the `mpif90` compiler command is deprecated, see also [here](https://docs.dkrz.de/doc/levante/code-development/compiling-and-linking.html#mpi-compiler-wrappers).

### OpenMPI4

This configurations uses openmpi and the Intel Fortran compiler as recommended by DKRZ. It is recommended to load those module consistently for compilation and in the runscript like, e.g, 

```bash
module load \
       intel-oneapi-compilers/2022.0.1-gcc-11.2.0 \
       openmpi/4.1.2-intel-2021.5.0 \
       intel-oneapi-mkl/2022.0.1-gcc-11.2.0 \
       netcdf-fortran/4.5.3-openmpi-4.1.2-intel-2021.5.0

export LD_LIBRARY_PATH=/sw/spack-levante/netcdf-fortran-4.5.3-k6xq5g/lib:$LD_LIBRARY_PATH

export OMPI_MCA_pml="ucx"
export OMPI_MCA_btl=self
export OMPI_MCA_osc="pt2pt"
export UCX_IB_ADDR_TYPE=ib_global
# for most runs one may or may not want to disable HCOLL
export OMPI_MCA_coll="^ml,hcoll"
export OMPI_MCA_coll_hcoll_enable="0"
export HCOLL_ENABLE_MCAST_ALL="0"
export HCOLL_MAIN_IB=mlx5_0:1
export UCX_NET_DEVICES=mlx5_0:1
export UCX_TLS=mm,knem,cma,dc_mlx5,dc_x,self
export UCX_UNIFIED_MODE=y
export HDF5_USE_FILE_LOCKING=FALSE
export OMPI_MCA_io="romio321"
export UCX_HANDLE_ERRORS=bt

ulimit -s 409600

```

The mpi compiler for openmpi is called `mpifort`.

### IntelMPI

It is also possible to use the IntelMPI compiler by sourcing the following environment (for compilation and in the runscript):

```bash
source /sw/intel/compilers_and_libraries_2020.2.254/linux/bin/compilervars.sh -ofi_internal=1 intel64
source /sw/intel/compilers_and_libraries_2020.2.254/linux/mkl/bin/mklvars.sh intel64
```
However, there might come updates in the future so you might not want to load specific version but always use the latest.

In your run-script we recommend having the following:
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/spack-levante/netcdf-fortran-4.5.3-r5r3ev/lib/
source /sw/intel/compilers_and_libraries_2020.2.254/linux/bin/compilervars.sh -ofi_internal=1 intel64
source /sw/intel/compilers_and_libraries_2020.2.254/linux/mkl/bin/mklvars.sh intel64
export MKL_DEBUG_CPU_TYPE=5
export MKL_CBWR=COMPATIBLE
export I_MPI_FABRICS=shm:ofi
export I_MPI_PMI=PMI2
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
```
Please update the versions, if needed. To run the model, please use:
```
srun  --ntasks=${NCPU} --mpi=pmi2...
```

The mpi compiler for openmpi is called `mpiifort`.

### IO

The model now contains a stream-like approach for handling internal fields
and meta data for IO for both IEG and NetCDF IO. Right now, the model requires
the NetCDF library to be linked. To choose the IO file formats, you have to
change the `out_filetype` to `NC` in `mo_file.f90` and set the file suffix 
(`out_suffix`) to `NC_SUFFIX` as well.

## Unit Architecture

The code of the model is scattered into several subdirectories in the `source` directory
depending on its purpose. Subdirectories that start with a lower case name, e.g., `physics`
are for organizational matters while units that implement an interface start with an upper case
name, e.g., `Radiation`. This directory should contain an empty `stub` implementation 
that describes the interface to the unit. The `stub` implementation will also make sure
that the rest of the code needs no modification if the unit is not used. 
Each unit typically contains a main implementation called `Main`, e.g., `RadiationMain` that 
can have different implementation of it's own. The original implementation of `REMO2015`
is typically found in the subdirectory `EC4` which contains the legacy (physics)
implementation of `ECHAM4` that `REMO2015` normally used, e.g., `Radiation/RadiationMain/EC4`.
However, if a new radiation physics implementation should be implemented, this would 
be in a new subdirectory of `RadiationMain`, e.g., `RadiationMain/MyNewScheme`.

At setup time, you can then choose the alternative implementation, using, e.g.
```
./setup -auto Remo2020 -objdir=Remo2020 -with-unit=physics/atmosphere/Radiation/RadiationMain/MyNewScheme
```
You will see, that the new radiation scheme is now linked into the object directory.

## Contact

lars.buntemeyer@hzg.de


