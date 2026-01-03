#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main(int argc, char *argv[]){
    int N = 1000;
    if (argc > 1) N = atoi(argv[1]);

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int i, j, k;
    int NRA = N;
    int NCA = N;
    int NCB = N;

    double a[NRA][NCA], b[NCA][NCB], c[NRA][NCB];

    int base_rows = N / size;
    int remainder = N % size;

    // int rows_per_proc = base_rows + (rank < remainder ? 1 : 0);
    int rows_per_proc = base_rows;

    double local_a[rows_per_proc][NCA], local_c[rows_per_proc][NCB];

    /*** Initialize matrices ***/
    
    if (rank==0){
        for (i=0; i<NRA; i++)
            for (j=0; j<NCA; j++)
                a[i][j]= i+j;
    
        for (i=0; i<NCA; i++)
            for (j=0; j<NCB; j++)
                b[i][j]= i*j;
    
        for (i=0; i<NRA; i++)
            for (j=0; j<NCB; j++)
                c[i][j]= 0;
    }

    // Broadcast full B to every single process
    MPI_Bcast(b, NCA*NCB, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Stripe-partition A and distribute over all processes
    MPI_Scatter(a, rows_per_proc * NCA, MPI_DOUBLE, local_a, rows_per_proc * NCA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();

    for (i=0; i<rows_per_proc; i++)
        for(j=0; j<NCB; j++){
            local_c[i][j] = 0;
            for (k=0; k<NCA; k++)
                local_c[i][j] += local_a[i][k] * b[k][j];
        }

    MPI_Barrier(MPI_COMM_WORLD);
    double t1 = MPI_Wtime();
    
    double local_time = t1 - t0;
    double max_time;

    MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    MPI_Gather(local_c, rows_per_proc * NCB, MPI_DOUBLE, c, rows_per_proc * NCB, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // for (i=0; i<NRA; i++)
        //     for (j=0; j<NCB; j++)
        //         printf("Example C[%d][%d] = %f\n", i, j, c[0][0]);
        // printf("Example C[0][0] = %f\n", c[0][0]);
        // printf("Example C[N-1][N-1] = %f\n", c[N-1][N-1]);
        printf("N=%d, p=%d, time=%f seconds\n", N, size, max_time);
    }

    MPI_Finalize();
    return 0;
}