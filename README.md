![SOMBRERO](sombrero.png)

A benchmarking utility for high performance computing based on lattice field theory applications.

## Requirements
  
  Requires an MPI library and compiler.

## Set Up

  1. Clone the repository

```
git clone https://github.com/sa2c/sombrero.git
cd sombrero
```

  2. Edit `Make/MkFlags` to set the compiler, `CFLAGS` and linker flags
  3. To compile run
 
```
make
```

## Usage

To run all benchmarks use the `sombrero.sh` script

```
./sombrero.sh -n <num-cores>  [ -w ] [ -s small | medium | large | very_large ]
```

By default this will run a medium scale problem using `<num-cores>` mpi ranks. The problem size independent of the number of processes. This is useful for strong scaling.

All output is printed to stdout and error messages to stderr.


## Parameters

`-n <num-cores>`: </br>
&nbsp;&nbsp;&nbsp; Set the number of mpi ranks. SOMBRERO will attempt to divide the problem accross the ranks automatically. The number of ranks `num-cores` should be a power of 2, `num-cores = 2^n`. This can be multiplied by a factor of 3, `num-cores = 3 * 2^n`.

`-s small | medium | large | very_large:`</br>
&nbsp;&nbsp;&nbsp; Specify the problem size. Compatible with weak or strong scaling.

```
./sombrero.sh -n Np -s large
```

`-w`:</br>
&nbsp;&nbsp;&nbsp; Weak scaling. The problem size will be a multiple of the number of MPI ranks.

```
./sombrero.sh -n Np -w 
```

By default a small local problem size is used, requiring maximal communication. A larger local problem size can be set as in the strong scaling case, for example

```
./sombrero.sh -n Np -w -s medium
```

Note that for the medium and large cases a minimum of 4 MPI ranks is required.

## Example Slurm Script

An example job script for a strong scaling study:

```
#!/bin/bash
#
#SBATCH --ntasks=256
#SBATCH --job-name=sombrero_strong
#SBATCH --time=0-0:20
#SBATCH --ntasks-per-node=32

cd $SLURM_SUBMIT_DIR

./sombrero.sh -n $SLURM_NTASKS -s medium > strong_$n
```

A weak scaling example:

```
#!/bin/bash
#
#SBATCH --ntasks=128
#SBATCH --job-name=sombrero_strong
#SBATCH --time=0-0:20
#SBATCH --ntasks-per-node=32

cd $SLURM_SUBMIT_DIR

./sombrero.sh -n $SLURM_NTASKS -w -s small > weak_$n
```

## Benchmarks

SOMBRERO runs 6 benchmarks, each representing a theory in under active study by the lattice field theory community. For the purposes of benchmarking, the models vary only in the amount of data communicated and the number of floating point operation required. Since the exact numbers depend on the problem size and the number of ranks, they are printed for each case.

The theories are

Case 1: Two color QCD (Quantum Chronodynamics)

Case 2: Two color QCD with adjoint fermions

Case 3: QCD

Case 4: Symplectic four color QCD

Case 5: QCD with sextet fermions

Case 6: Symplectic four color QCD with adjoint fermions

## Miscellaneous

Additional parameters `-l` and `-p` exist for more precise control of the problem size. The benchmark theories are defined on a four dimensional lattice, which is partitioned equally among the MPI ranks.

`-l NxNxNxN`:</br>
&nbsp;&nbsp;&nbsp; Set the lattice size in each direction. Replace each `N` with an integer larger or equal to 4.

`-p NxNxNxN`:</br>
&nbsp;&nbsp;&nbsp; Manually partition the lattice among the mpi ranks. The lattice is divided accross `N` ranks in each direction. A positive integer must be used.

Once compiled, binaries for individual test cases can be found in the sombrero folder. For example the case 3 can be run separately with

```
mpirun -n Np sombrero/sombrero3 -s medium 
```

The verbose mode will produce additional information:

```
mpirun -n Np sombrero/sombrero3 -s medium -v verbose 
```

