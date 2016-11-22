// MatrixMultiplication.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <iostream>
#include <ctime>
#include <omp.h>
using namespace std;

int** matrixInt(int n){
	int **m = new int* [n];
	for (int i = 0; i < n; i++){
		m[i] = new int[n];
	}
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			m[i][j] = rand() % 10;
		}
	}
	return m;
}

double** matrixDouble(int n){
	double** m = new double* [n];
	for (int i = 0; i < n; i++){
		m[i] = new double[n];
	}
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			m[i][j] = rand() % 10;
		}
	}
	return m;

}

float** matrixFloat(int n){
	float** m = new float* [n];
	for (int i = 0; i < n; i++){
		m[i] = new float[n];
	}
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			m[i][j] = rand() % 10;
		}
	}
	return m;
}


void printMatrix(int** m, int n){
	cout << "Matrix";
	cout << endl;
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			cout <<m[i][j]<< " ";
		}
		cout << endl;
	}
	cout << endl;
}

void freeMemoryInt(int** m, int n){
	for (int i = 0; i < n; i++)
    {
        delete[] m[i]; 
    }
    delete [] m; 
}

void freeMemoryDouble(double** m, int n){
	for (int i = 0; i < n; i++)
    {
        delete[] m[i]; 
    }
    delete [] m;
}

void freeMemoryFloat(float** m, int n){
	for (int i = 0; i < n; i++)
    {
        delete[] m[i]; 
    }
    delete [] m;
}

int _tmain(int argc, _TCHAR* argv[])
{
	int n = 0;
	int i,j,k;
	cout<<"Size of matrix n: ";
	cin>>n;
	//int
	int **a = matrixInt(n);
	int **b = matrixInt(n);
	int **c = matrixInt(n);
	double** a2 = matrixDouble(n);
	double** b2 = matrixDouble(n);
	double** c2 = matrixDouble(n);
	float** a3 = matrixFloat(n);
	float** b3 = matrixFloat(n);
	float** c3 = matrixFloat(n);
	//последовательный алгоритм
	double start = time(NULL);
	for (i=0; i<n; i++){
		for (j=0; j<n; j++){
			for (k=0; k<n; k++){
				c[i][j]+=a[i][k]*b[k][j];
			}
		}
	}
	time_t end = time(NULL);
	double runtime1 = end-start; 
	cout<<"Simple program int "<<runtime1<<endl;
	freeMemoryInt(c,n);
	//параллельный алгоритм
	c = matrixInt(n);
	start = time(NULL);
		#pragma omp parallel for shared(a,b,c,n) private (i,j,k)
			for (i=0; i<n; i++){
				for (j=0; j<n; j++){
					for (k=0; k<n; k++){
						c[i][j]+=a[i][k]*b[k][j];
					}
				}
			}
			//cout<<"threads "<<omp_get_num_threads()<<endl;
	end = time(NULL);
	double runtime2 = end-start;
		cout<<"Parallel program int "<<runtime2<<"  "<<runtime1/runtime2<<endl;
	freeMemoryInt(a,n);
	freeMemoryInt(b,n);
	freeMemoryInt(c,n);
	
	//double
	//последовательный алгоритм
	start = time(NULL);
	for (i=0; i<n; i++){
		for (j=0; j<n; j++){
			for (k=0; k<n; k++){
				c2[i][j]+=a2[i][k]*b2[k][j];
			}
		}
	}
	end = time(NULL);
	runtime1 = end-start;
	cout<<"Simple program double "<<runtime1<<endl;
	freeMemoryDouble(c2,n);
	//параллельный алгоритм
	c2 = matrixDouble(n);
	start = time(NULL);
		#pragma omp parallel for shared(a2,b2,c2,n) private (i,j,k)
			for (i=0; i<n; i++){
				for (j=0; j<n; j++){
					for (k=0; k<n; k++){
						c2[i][j]+=a2[i][k]*b2[k][j];
					}
				}
			}

	end = time(NULL);
	runtime2 = end-start;
	cout<<"Parallel program double "<<runtime2<<"  "<<runtime1/runtime2<<endl;
	freeMemoryDouble(a2,n);
	freeMemoryDouble(b2,n);
	freeMemoryDouble(c2,n);

	//float
	//последовательный алгоритм
	start = time(NULL);
	for (i=0; i<n; i++){
		for (j=0; j<n; j++){
			for (k=0; k<n; k++){
				c3[i][j]+=a3[i][k]*b3[k][j];
			}
		}
	}
	end = time(NULL);
	runtime1 = end - start;
	cout<<"Simple program float "<<runtime1<<endl;
	freeMemoryFloat(c3,n);
	//параллельный алгоритм
	c3 = matrixFloat(n);
	start = time(NULL);
		#pragma omp parallel for shared(a3,b3,c3,n) private (i,j,k)
			for (i=0; i<n; i++){
				for (j=0; j<n; j++){
					for (k=0; k<n; k++){
						c3[i][j]+=a3[i][k]*b3[k][j];
					}
				}
			}
	end = time(NULL);
	runtime2 = end - start;
	cout<<"Parallel program float "<<end-start<<"  "<<runtime1/runtime2<<endl;
	freeMemoryFloat(a3,n);
	freeMemoryFloat(b3,n);
	freeMemoryFloat(c3,n);
	
	return 0;
}