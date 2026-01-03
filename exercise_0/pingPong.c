#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

// Maximum array size 2^20= 1048576 elements
#define MAX_ARRAY_SIZE (1<<20)
#define NUM_TRIALS 10

int main(int argc, char **argv)
{
    // Variables for the process rank and number of processes
    int myRank, numProcs, i;
    MPI_Status status;

    // Initialize MPI, find out MPI communicator size and process rank
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);


    int *myArray = (int *)malloc(sizeof(int)*MAX_ARRAY_SIZE);
    if (myArray == NULL)
    {
        printf("Not enough memory\n");
        exit(1);
    }
    // Initialize myArray
    for (i=0; i<MAX_ARRAY_SIZE; i++)
        myArray[i]=1;

    int numberOfElementsToSend;
    int numberOfElementsReceived;
    double startTime, endTime;

    // PART C
    if (numProcs < 2)
    {
        printf("Error: Run the program with at least 2 MPI tasks!\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /// for (int trial = 0; trial < NUM_TRIALS; trial++){
    // Loop over every message size in powers of 2
    for (int n=1; n<=MAX_ARRAY_SIZE; n*=4) {
        
        numberOfElementsToSend = n;

        if (myRank == 0)
        {
            printf("Rank %2.1i: Sending %i elements\n", myRank, numberOfElementsToSend);

            myArray[0]=myArray[1]+1; // activate in cache (avoids possible delay when sending the 1st element)

            startTime = MPI_Wtime();
            for (i=0; i<5; i++)
            {
                MPI_Send(myArray, numberOfElementsToSend, MPI_INT, 1, 0, MPI_COMM_WORLD);

                // Probe message in order to obtain the amount of data
                MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, MPI_INT, &numberOfElementsReceived);
                
                MPI_Recv(myArray, numberOfElementsReceived, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } // end of for-loop

            endTime = MPI_Wtime();

            printf("Rank %2.1i: Received %i elements\n", myRank, numberOfElementsReceived);

            // average communication time of 1 send-receive (total 5*2 times)
            printf("Ping Pong took %f seconds\n", (endTime - startTime)/10);
        }
        else if (myRank == 1)
        {
            // Probe message in order to obtain the amount of data
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_INT, &numberOfElementsReceived);
           
            for (i=0; i<5; i++)
            {
                MPI_Recv(myArray, numberOfElementsReceived, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                // printf("Rank %2.1i: Received %i elements\n", myRank, numberOfElementsReceived);

                // printf("Rank %2.1i: Sending back %i elements\n", myRank, numberOfElementsToSend);

                MPI_Send(myArray, numberOfElementsToSend, MPI_INT, 0, 0, MPI_COMM_WORLD);
            } // end of for-loop
        }
    } // end of for-loop
//}

    // Finalize MPI
    MPI_Finalize();

    return 0;
}
