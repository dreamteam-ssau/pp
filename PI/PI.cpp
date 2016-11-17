#include "stdafx.h"



//using namespace std;

void buffon(int n);
void buffonOMP(int n);

int main()
{
	
	int n;
	std::cout << "number of trials ";
	std::cin >> n;
	buffon(n);
	buffonOMP(n);
	system("pause");

	return 0;
}

void buffon(int n) {
	srand(time(NULL));
	double x; // x coordinate of needle's center 
	double k; // angle between vertical position and needle
	int p = 0; // positive trials
	double y; // sin(angle) * l
	double pi;
	time_t start = time(NULL);
	double l = 1.0;
	for (int i = 0; i < n; i++)
	{
		k = (double)rand() / (RAND_MAX) * 2 * M_PI;       // random angle

		x = (double)rand() / (RAND_MAX * 2);         // random x (0 do 1)

		y = (l / 2) * sin(k);


		if (x <= y)
		{
			p++;
		}

	}
	time_t end = time(NULL);
	pi = (l*n) / (p);

	std::cout << "n = ";
	std::cout << n << std::endl;
	std::cout << "p = ";
	std::cout << p << std::endl;
	std::cout << "net OMP: " << end - start << std::endl;
	std::cout  << pi << std::endl << std::endl;
}

void buffonOMP(int n) {//29 cек на 100000000
	srand(time(NULL));
	double x; // x coordinate of needle's center 
	double k; // angle between vertical position and needle
	int p = 0; // positive trials
	double y; // sin(angle) * l
	double pi;
	time_t start = time(NULL);
	double l = 1.0;
	#pragma omp parallel for reduction(+:p)
	for (int i = 0; i < n; i++)
	{
		k = (double)rand() / (RAND_MAX) * 2 * M_PI;       // random angle
		x = (double)rand() / (RAND_MAX * 2);         // random x (0 do 1)

		y = (l / 2) * sin(k);
		//std::cout << omp_get_num_threads();
		if (x <= y)
		{
			p++;
		}

	}
	time_t end = std::time(NULL);
	pi = (l*n) / (p);

	std::cout << "n = ";
	std::cout << n << std::endl;
	std::cout << "p = ";
	std::cout << p << std::endl;
	std::cout << "OMP: " << end - start << std::endl;
	std::cout << pi << std::endl;
}