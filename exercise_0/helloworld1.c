/*
////////////////////////////////////////////////////////////////////
This is a specific example on how to run a MPI program successfully.
////////////////////////////////////////////////////////////////////
*/

/*
--Log onto DelftBlue

--Create a file (helloworld.c) and type the following contents into helloworld.c 
====================================================================
*/

#include "mpi.h"
#include <stdio.h>

int np, rank;

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  printf("Node %d of %d says: Hello world!\n", rank + 1, np);
  
  MPI_Finalize();
  return 0;
}

/*
====================================================================

--Load the two modules
====================================================================
module load 2022r2 openmpi
====================================================================

--Compile the helloworld program (Be careful with the input order ! ! !)
====================================================================
mpicc -o helloworld helloworld.c
====================================================================

--Run the helloworld program with 2 processes, 4 cores for each process
====================================================================
srun -n 2 -c 4 --mem-per-cpu=1GB  ./helloworld
====================================================================

--Welcome to this fancy MPI world!
*/