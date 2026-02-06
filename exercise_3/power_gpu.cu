// the subroutine for GPU code can be found in several separated text file from the Brightspace. 
// You can add these subroutines to this main code.
////////////////////////////////////////////


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "cuda.h"

#ifndef BLOCK_SIZE
#define BLOCK_SIZE 32
#endif

// Input Array Variables
float* h_MatA = NULL;
float* d_MatA = NULL;

// Output Array
float* h_VecV = NULL;
float* d_VecV = NULL;
float* h_VecW = NULL;
float* d_VecW = NULL;
float* h_NormW = NULL;
float* d_NormW = NULL;
float* d_Lamda = NULL;

// Variables to change
int GlobalSize = 5000;         // this is the dimension of the matrix, GlobalSize*GlobalSize
int BlockSize = BLOCK_SIZE;            // number of threads in each block
const float EPS = 0.000005;    // tolerence of the error
int max_iteration = 100;       // the maximum iteration steps


// Functions
void Cleanup(void);
void InitOne(float*, int);
void UploadArray(float*, int);
float CPUReduce(float*, int);
void  Arguments(int, char**);
void checkCardVersion(void);

// Kernels
__global__ void AvProduct(float* g_MatA, float* g_VecV, float* g_VecW, int N);
__global__ void FindNormW(float* g_VecW, float * g_NormW, int N);
__global__ void NormalizeW(float* g_VecV,float* g_VecW, float* g_NormW, int N);
__global__ void ComputeLamda( float* g_VecV,float* g_VecW, float * g_Lamda, int N);


__global__ void AvProduct(float* g_MatA, float* g_VecV, float* g_VecW, int N)
{
    int b_idx = blockIdx.x;     // block index
    int t_idx = threadIdx.x;    // thread index

    int a_begin = N * b_idx * BLOCK_SIZE;
    int a_end = a_begin + N;

    int step = BLOCK_SIZE;

    int v_begin = 0;      // BLOCK_SIZE * b_idx
    int v_idx = 0;
    int a_idx = 0;
    float w_sub = 0;      // vector w for each block vector subspace

    for (int a = a_begin, v = v_begin; a < a_end; a += step, v += step)
    {
        __shared__ float A_sub[BLOCK_SIZE * BLOCK_SIZE];
        __shared__ float v_sub[BLOCK_SIZE];

        for (int aa = 0; aa < BLOCK_SIZE; aa++)
        {
            a_idx = a + t_idx + aa * N;
            if (a_idx < N * N)
                A_sub[t_idx + aa * BLOCK_SIZE] = g_MatA[a_idx];
            else
                A_sub[t_idx + aa * BLOCK_SIZE] = 0;
        }

        v_idx = t_idx + v;
        if (v_idx < N)
            v_sub[t_idx] = g_VecV[v_idx];
        else
            v_sub[t_idx] = 0;

        __syncthreads();

        for (int k = 0; k < BLOCK_SIZE; k++)
        {
            w_sub += A_sub[k + t_idx * BLOCK_SIZE] * v_sub[k];
        }
        __syncthreads();
    }

    g_VecW[BLOCK_SIZE * b_idx + t_idx] = w_sub;
}

__global__ void FindNormW(float* g_VecW, float * g_NormW, int N)
{
    extern __shared__ float partial_sum[];

    int t_idx = threadIdx.x;
    int g_idx = blockIdx.x * blockDim.x + threadIdx.x;

    // Square every thread's element
    if (g_idx < N)
        partial_sum[t_idx] = g_VecW[g_idx] * g_VecW[g_idx];
    else
        partial_sum[t_idx] = 0;

    __syncthreads();

    // Parallel reduction within each block
    unsigned int n = blockDim.x;
    while (n > 1)
    {
        unsigned int half = n / 2;

        if (t_idx < half)
            partial_sum[t_idx] += partial_sum[t_idx + half];

        __syncthreads();

        if ((n % 2 == 1) && t_idx == 0)
            partial_sum[0] += partial_sum[n - 1];

        __syncthreads();

        n = half;
    }

    if (t_idx == 0)
    {
        atomicAdd(g_NormW, partial_sum[0]);
    }
}

__global__ void NormalizeW(float* g_VecV,float* g_VecW, float* g_NormW, int N)
{
    extern __shared__ float s_data[];

    int t_idx = threadIdx.x;
    int g_idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (t_idx == 0)
        s_data[0] = g_NormW[0];
    __syncthreads();

    if (g_idx < N)
        g_VecV[g_idx] = g_VecW[g_idx] / s_data[0];
}

__global__ void ComputeLamda(float* g_VecV,float* g_VecW, float* g_Lamda, int N)
{
    extern __shared__ float partial_sum[];

    int t_idx = threadIdx.x;
    int g_idx = blockIdx.x * blockDim.x + threadIdx.x;

    // Square every thread's element
    if (g_idx < N)
        partial_sum[t_idx] = g_VecV[g_idx] * g_VecW[g_idx];
    else
        partial_sum[t_idx] = 0;

    __syncthreads();

    // Parallel reduction within each block
    unsigned int n = blockDim.x;
    while (n > 1)
    {
        unsigned int half = n / 2;

        if (t_idx < half)
            partial_sum[t_idx] += partial_sum[t_idx + half];

        __syncthreads();

        if ((n % 2 == 1) && t_idx == 0)
            partial_sum[0] += partial_sum[n - 1];

        __syncthreads();

        n = half;
    }

    if (t_idx == 0)
    {
        atomicAdd(g_Lamda, partial_sum[0]);
    }
}

void CPU_AvProduct()
{
	int N = GlobalSize;
	int matIndex = 0;
    for(int i=0;i<N;i++)
	{
		h_VecW[i] = 0;
		for(int j=0;j<N;j++)
		{
			matIndex = i*N + j;
			h_VecW[i] += h_MatA[matIndex] * h_VecV[j];
		}
	}
}

void CPU_NormalizeW()
{
	int N = GlobalSize;
	float normW=0;
	for(int i=0;i<N;i++)
		normW += h_VecW[i] * h_VecW[i];
	
	normW = sqrt(normW);
	for(int i=0;i<N;i++)
		h_VecV[i] = h_VecW[i]/normW;
}

float CPU_ComputeLamda()
{
	int N = GlobalSize;
	float lamda =0;
	for(int i=0;i<N;i++)
		lamda += h_VecV[i] * h_VecW[i];
	
	return lamda;
}

void RunCPUPowerMethod()
{
	printf("*************************************\n");
	float oldLamda =0;
	float lamda=0;
	
	//AvProduct
	CPU_AvProduct();
	
	//power loop
	for (int i=0;i<max_iteration;i++)
	{
		CPU_NormalizeW();
		CPU_AvProduct();
		lamda= CPU_ComputeLamda();
		printf("CPU lamda at %d: %f \n", i, lamda);
		// If residual is lass than epsilon break
		if(abs(oldLamda - lamda) < EPS)
			break;
		oldLamda = lamda;	
	
	}
	printf("*************************************\n");
	
}

// Host code
int main(int argc, char** argv)
{
    double memcpy_time = 0;
    double memset_time = 0;
    double memalloc_time = 0;
    struct timespec t1, t2;

    struct timespec t_start, t_end;
    double runtime;
    double cpu_runtime;
    double cpu_malloc_time = 0;
    Arguments(argc, argv);
		
    int N = GlobalSize;
    printf("Matrix size %d X %d \n", N, N);
    size_t vec_size = N * sizeof(float);
    size_t mat_size = N * N * sizeof(float);
    size_t norm_size = sizeof(float);
  
    clock_gettime(CLOCK_REALTIME,&t_start);
    // Allocate normalized value in host memory
    h_NormW = (float*)malloc(norm_size);
    // Allocate input matrix in host memory
    h_MatA = (float*)malloc(mat_size);
    // Allocate initial vector V in host memory
    h_VecV = (float*)malloc(vec_size);
    // Allocate W vector for computations
    h_VecW = (float*)malloc(vec_size);
    clock_gettime(CLOCK_REALTIME,&t_end);

    cpu_malloc_time += (t_end.tv_sec - t_start.tv_sec) + 1e-9*(t_end.tv_nsec - t_start.tv_nsec);


    // Initialize input matrix
    UploadArray(h_MatA, N);
    InitOne(h_VecV,N);

    printf("Power method in CPU starts\n");	   
    clock_gettime(CLOCK_REALTIME,&t_start);
    RunCPUPowerMethod();   // the lamda is already solved here
    clock_gettime(CLOCK_REALTIME,&t_end);
    cpu_runtime = (t_end.tv_sec - t_start.tv_sec) + 1e-9*(t_end.tv_nsec - t_start.tv_nsec);
    printf("CPU: run time = %f secs.\n", cpu_runtime);
    // printf("CPU: malloc time = %f secs.\n", cpu_malloc_time);
    printf("Power method in CPU is finished\n");
    
    
    /////////////////////////////////////////////////
    // This is the starting points of GPU
    printf("\nPower method in GPU starts\n");
    checkCardVersion();

    // Initialize input matrix
    InitOne(h_VecV,N);
    
    clock_gettime(CLOCK_REALTIME, &t_start);  // Here I start to count

    // Set the kernel arguments
    int threadsPerBlock = BlockSize;   
    int sharedMemSize = threadsPerBlock * threadsPerBlock * sizeof(float); // in per block, the memory is shared   
    int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;

    // Allocate matrix and vectors in device memory
    clock_gettime(CLOCK_REALTIME, &t1);
    cudaMalloc((void**)&d_MatA, mat_size); 
    cudaMalloc((void**)&d_VecV, vec_size); 
    cudaMalloc((void**)&d_VecW, vec_size); // This vector is only used by the device
    cudaMalloc((void**)&d_NormW, norm_size); 
    cudaMalloc((void**)&d_Lamda, norm_size);
    clock_gettime(CLOCK_REALTIME, &t2);

    memalloc_time += (t2.tv_sec - t1.tv_sec) + 1e-9 * (t2.tv_nsec - t1.tv_nsec);

    //Copy from host memory to device memory
    clock_gettime(CLOCK_REALTIME, &t1);
    cudaMemcpy(d_MatA, h_MatA, mat_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_VecV, h_VecV, vec_size, cudaMemcpyHostToDevice);
    clock_gettime(CLOCK_REALTIME, &t2);

    memcpy_time += (t2.tv_sec - t1.tv_sec) + 1e-9 * (t2.tv_nsec - t1.tv_nsec);

	// cutilCheckError(cutStopTimer(timer_mem));

    //Power method loops
    printf("*************************************\n");


    float oldLamda = 0;
    float lamda = 0;

    AvProduct<<<blocksPerGrid, threadsPerBlock>>>(d_MatA, d_VecV, d_VecW, N);
    cudaDeviceSynchronize();

    for (int i=0; i < max_iteration; i++)
    {
        clock_gettime(CLOCK_REALTIME, &t1);
        cudaMemset(d_NormW, 0, norm_size);
        cudaMemset(d_Lamda, 0, norm_size);
        clock_gettime(CLOCK_REALTIME, &t2);

        memset_time += (t2.tv_sec - t1.tv_sec) + 1e-9 * (t2.tv_nsec - t1.tv_nsec);

        FindNormW<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(d_VecW, d_NormW, N);
        cudaDeviceSynchronize();

        clock_gettime(CLOCK_REALTIME, &t1);
        cudaMemcpy(h_NormW, d_NormW, norm_size, cudaMemcpyDeviceToHost);
        clock_gettime(CLOCK_REALTIME, &t2);

        memcpy_time += (t2.tv_sec - t1.tv_sec) + 1e-9 * (t2.tv_nsec - t1.tv_nsec);

        *h_NormW = sqrt(*h_NormW);

        clock_gettime(CLOCK_REALTIME, &t1);
        cudaMemcpy(d_NormW, h_NormW, norm_size, cudaMemcpyHostToDevice);
        clock_gettime(CLOCK_REALTIME, &t2);

        memcpy_time += (t2.tv_sec - t1.tv_sec) + 1e-9 * (t2.tv_nsec - t1.tv_nsec);

        NormalizeW<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(d_VecV, d_VecW, d_NormW, N);
        cudaDeviceSynchronize();

        AvProduct<<<blocksPerGrid, threadsPerBlock>>>(d_MatA, d_VecV, d_VecW, N);
        cudaDeviceSynchronize();
    
        ComputeLamda<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(d_VecV, d_VecW, d_Lamda, N);
        cudaDeviceSynchronize();

        clock_gettime(CLOCK_REALTIME, &t1);
        cudaMemcpy(&lamda, d_Lamda, norm_size, cudaMemcpyDeviceToHost);
        clock_gettime(CLOCK_REALTIME, &t2);

        memcpy_time += (t2.tv_sec - t1.tv_sec) + 1e-9 * (t2.tv_nsec - t1.tv_nsec);

        printf("GPU lamda at %d: %f\n", i, lamda);

		if(abs(oldLamda - lamda) < EPS)
			break;
		oldLamda = lamda;	
    }
    printf("*************************************\n");
    
    clock_gettime(CLOCK_REALTIME, &t_end);
    runtime = (t_end.tv_sec - t_start.tv_sec) + 1e-9*(t_end.tv_nsec - t_start.tv_nsec);
    double kernel_runtime = runtime - memcpy_time;
    printf("GPU: Total runtime (incl. Memory copies) = %f secs.\n", runtime);
    // printf("GPU: Kernel runtime = %f secs.\n", kernel_runtime);
    // printf("GPU: Memcpy runtime = %f secs.\n", memcpy_time);
    // printf("GPU: Memset runtime = %f secs.\n", memset_time);
    // printf("GPU: Memalloc runtime = %f secs.\n", memalloc_time);
    printf("GPU: Speedup (incl. Memory copies) = %f.\n", cpu_runtime / runtime);
    printf("GPU: Speedup (excl. Memory copies) = %f.\n", cpu_runtime / kernel_runtime);
    // printf("Overall CPU Execution Time: %f (ms) \n", cutGetTimerValue(timer_CPU));

    Cleanup();
}

void Cleanup(void)
{
    // Free device memory
    if (d_MatA)
        cudaFree(d_MatA);
    if (d_VecV)
        cudaFree(d_VecV);
    if (d_VecW)
        cudaFree(d_VecW);
	if (d_NormW)
		cudaFree(d_NormW);
    if (d_Lamda)
        cudaFree(d_Lamda);
		
    // Free host memory
    if (h_MatA)
        free(h_MatA);
    if (h_VecV)
        free(h_VecV);
    if (h_VecW)
        free(h_VecW);
    if (h_NormW)
        free(h_NormW);
    
    exit(0);
}

// Allocates an array with zero value.
void InitOne(float* data, int n)
{
    for (int i = 0; i < n; i++)
        data[i] = 0;
	data[0]=1;
}

void UploadArray(float* data, int n)
{
   int total = n*n;
   int value=1;
    for (int i = 0; i < total; i++)
    {
    	data[i] = (int) (rand() % (int)(101));//1;//value;
	    value ++; if(value>n) value =1;
      // data[i] = 1;
    }
}

// Obtain program arguments
void Arguments(int argc, char** argv)
{
    for (int i = 0; i < argc; ++i) 
    {
        if (strcmp(argv[i], "--size") == 0 || strcmp(argv[i], "-size") == 0)
        {
            GlobalSize = atoi(argv[i+1]);
		    i = i + 1;
        }
        if (strcmp(argv[i], "--max_iteration") == 0 || strcmp(argv[i], "-max_iteration") == 0)
        {
            max_iteration = atoi(argv[i+1]);
		    i = i + 1;
        }
    }
}


void checkCardVersion()
{
   cudaDeviceProp prop;
   
   cudaGetDeviceProperties(&prop, 0);
   
   printf("This GPU has major architecture %d, minor %d \n",prop.major,prop.minor);
   if(prop.major < 2)
   {
      fprintf(stderr,"Need compute capability 2 or higher.\n");
      exit(1);
   }
}
