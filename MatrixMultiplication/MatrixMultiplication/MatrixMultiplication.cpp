#define _CRT_SECURE_NO_WARNINGS
#include <locale>
#include <stdio.h>
#include <omp.h>
#include "mpi.h"

#define TYPE_FLOAT

#ifdef TYPE_FLOAT
#define TYPE float
#define MPI_TYPE MPI_FLOAT
#elif TYPE_INT
#define TYPE int
#define MPI_TYPE MPI_INT
#elif TYPE_DOUBLE
#define TYPE double
#define MPI_TYPE MPI_DOUBLE
#endif

typedef TYPE** Matrix;

int myid, count;
int outN;
int TEST_N = 3; // Кол-во доп прогонок
//double *resPosl, *resParal;
double *runtimePosl, *runtimeMpi;
const int MIN = 100, MAX = 1000, STEP = 100;
const int SIZE = (MAX - MIN) / STEP + 1;

// Рандомная генерация матрицы
static Matrix GenerateMatrix(int N, int M)
{
	TYPE* data = (TYPE*)malloc(N * M * sizeof(TYPE));
	Matrix m = (Matrix)malloc(N * sizeof(TYPE*));
	for (int i = 0; i < N; i++)
	{
		m[i] = &(data[M * i]);
		for (int j = 0; j < M; j++)
			m[i][j] = ((TYPE)rand()) / RAND_MAX;
	}
	return m;
}

// Матрица из массива N*M
static Matrix FromRawMatrix(TYPE* data, int N, int M)
{
	Matrix m = (Matrix)malloc(N * sizeof(TYPE*));
	for (int i = 0; i < N; i++)
		m[i] = &(data[M * i]);
	return m;
}

// Пустая матрица
static Matrix EmptyMatrix(int N, int M)
{
	TYPE* data = (TYPE*)malloc(N * M * sizeof(TYPE));
	Matrix m = (Matrix)malloc(N * sizeof(TYPE*));
	for (int i = 0; i < N; i++)
		m[i] = &(data[M * i]);
	return m;
}

static void FreeMatrix(Matrix m, int N, int M)
{
	if (m == NULL)
		return;

	free(m[0]);
	free(m);
}

static void PrintMatrix(Matrix m, int N, int M)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
			printf("%f ", m[i][j]);
		printf("\n");
	}
}

static Matrix BroadcastMatrix(Matrix orig, int N, int M)
{
	Matrix m;
	if (myid != 0)
		m = EmptyMatrix(N, M);
	else
		m = orig;

	MPI_Bcast(m[0], M*N, MPI_TYPE, 0, MPI_COMM_WORLD);
	return m;
}

static Matrix Multiply(Matrix m1, int N1, int M1, Matrix m2, int N2, int M2)
{
	if (M1 != N2)
		return NULL;

	Matrix ret = EmptyMatrix(N1, M2);
	int i, j, k;

	for (i = 0; i < N1; i++)
		for (j = 0; j < M2; j++)
		{
			ret[i][j] = 0;
			for (k = 0; k < M1; k++) {
				ret[i][j] += m1[i][k] * m2[k][j];
			}
		}
	return ret;
}

static Matrix MultiplyOmp(Matrix m1, int N1, int M1, Matrix m2, int N2, int M2){
	if (M1 != N2)
		return NULL;

	Matrix ret = EmptyMatrix(N1, M2);
	int i, j, k;
	#pragma omp parallel for shared(m1,m2,ret,N1,M2) private (i,j,k)
	for (i = 0; i < N1; i++)
		for (j = 0; j < M2; j++)
		{
			ret[i][j] = 0;
			for (k = 0; k < M1; k++) {
				ret[i][j] += m1[i][k] * m2[k][j];
			}
		}
	return ret;
}

static Matrix MultiplyMPI(Matrix m1, int N1, int M1, Matrix m2, int N2, int M2)
{
	if (M1 != N2)
		return NULL;

	int chunk = N1 / count;
	if (myid + 1 < count)
		outN = chunk;
	else
		outN = chunk + N1 % count;
	int x0 = chunk * myid;

	Matrix ret = EmptyMatrix(outN, M2);
	for (int i = x0; i < x0 + outN; i++)
		for (int j = 0; j < M2; j++)
		{
			ret[i - x0][j] = 0;
			for (int k = 0; k < M1; k++)
				ret[i - x0][j] += m1[i][k] * m2[k][j];
		}
	return ret;
}

// Генерирует две квадратные матрицы и перемножает их
static void CalcMatrix(int N, int M, bool output = false)
{
	Matrix m1 = EmptyMatrix(N,M);
	Matrix m2 = EmptyMatrix(N,M);
	if (myid == 0) {
		m1 = GenerateMatrix(N, M);
		m2 = GenerateMatrix(M, N);
		if (output)
			PrintMatrix(m1, N, M);
	}

	m1 = BroadcastMatrix(m1, N, M);
	m2 = BroadcastMatrix(m2, M, N);
	Matrix myAns = MultiplyMPI(m1, N, M, m2, M, N);

	TYPE *buf = (TYPE*)malloc(N * M * sizeof(TYPE));
	int*counted = (int*)malloc(sizeof(int));
	int *disp = (int*)malloc(sizeof(int));
	if (myid == 0) {
		counted = (int*)malloc(count * sizeof(int));
		disp = (int*)malloc(count * sizeof(int));
		for (int i = 0; i + 1 < count; i++) {
			counted[i] = N / count * M;
			disp[i] = counted[i] * i;
		}
		counted[count - 1] = (N / count + N % count) *M;
		disp[count - 1] = counted[count - 2] * (count - 1);
	}

	MPI_Gatherv(myAns[0], outN * M, MPI_TYPE, buf, counted, disp, MPI_TYPE, 0, MPI_COMM_WORLD);
	if (myid == 0) {
		Matrix totalAns = FromRawMatrix(buf, N, M);
		if (output)
		{
			printf("==========Answer=========\n");
			PrintMatrix(totalAns, N, M);
		}
		FreeMatrix(totalAns, N, M);
	}

	FreeMatrix(m1, N, M);
	FreeMatrix(m2, M, N);
	FreeMatrix(myAns, outN, M);
}

static void MPITask()
{
	int dim[7] = {10,100,200,400,800,1000,2000};
	for (int i = 0; i <7; i ++)
	{
		if (myid == 0)
			printf("Matrix size: %d ", dim[i]);
			double start = 0, end = 0;
			if (myid == 0)
				start = omp_get_wtime();
			CalcMatrix(dim[i], dim[i]);
			if (myid == 0)
			{
				end = omp_get_wtime();
				runtimeMpi[i] = (end - start);
			}
		if (myid == 0)
			printf("Runtime: %f \n", runtimeMpi[i]);
	}
}

static double PoslTask(int N1, int M1, int N2, int M2)
{
	//Генерируем матрицы
	Matrix m1 = GenerateMatrix(N1, M1);
	Matrix m2 = GenerateMatrix(N2, M2);

	//Подсчитаем результат
	double time;
		double start, end;
		start = omp_get_wtime();
		Matrix res = Multiply(m1, N1, M1, m2, N2, M2);
		end = omp_get_wtime();
		time = (end - start);
		FreeMatrix(res, M1, N2);
	// Освобождаем матрицы
	FreeMatrix(m1, N1, M1);
	FreeMatrix(m2, N2, M2);
	return time;
}

static double OmpTask(int N1, int M1, int N2, int M2){
	//Генерируем матрицы
	Matrix m1 = GenerateMatrix(N1, M1);
	Matrix m2 = GenerateMatrix(N2, M2);

	//Подсчитаем результат
	double time;
		double start, end;
		start = omp_get_wtime();
		Matrix res = MultiplyOmp(m1, N1, M1, m2, N2, M2);
		end = omp_get_wtime();
		time = (end - start);
		FreeMatrix(res, M1, N2);
	// Освобождаем матрицы
	FreeMatrix(m1, N1, M1);
	FreeMatrix(m2, N2, M2);
	return time;
}

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &count);
	int dim[7] = {10,100,200,400,800,1000,2000};
	
	double runtimeOmp[7]={0,0,0,0,0,0,0};
	
	if (myid == 0)
	{
		//resPosl = new double[SIZE];
		//resParal = new double[SIZE];
		runtimePosl=new double[7];
		runtimeMpi = new double[7];
		printf("==========Posl=========\n");
		
		for (int i = 0; i < 7; i++ )
		{
			printf("Size: %d ", dim[i]);
			runtimePosl[i]=PoslTask(dim[i],dim[i],dim[i],dim[i]);
			printf("Runtime: %f \n", runtimePosl[i]);
		}
		printf("==========OMP=========\n");
		for (int i = 0; i < 7; i++ )
		{
			printf("Size: %d ", dim[i]);
			runtimeOmp[i]=PoslTask(dim[i],dim[i],dim[i],dim[i]);
			printf("Runtime: %f \n", runtimeOmp[i]);
		}
		printf("==========MPI=========\n");
	}
	MPITask();

	// Подсчитываем ускорение
	if (myid == 0)
	{
		double acmRes = 0;
		for (int i = 0; i < 7; i++){
			printf("Size: %d ",dim[i]);
			printf("Uscor Posl-Omp: %f\n", runtimePosl[i] / runtimeOmp[i]);
			printf("Uscor Posl-Mpi: %f\n", runtimePosl[i] / runtimeMpi[i]);
		}
	}

	MPI_Finalize();
	return 0;
}
