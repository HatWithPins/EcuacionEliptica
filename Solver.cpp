#include "Solver.h"

double pi = 3.14159265359;
//Definición de b
//Define b
double* defineB(bool isPolar, double boundaryBottom = 0.0, double boundaryArc = 1.0, int N = 3, int M = 6)
{
	int sizeB = (N - 1) * (M - 1);
	int col = 0;
	int row = 0;
	double h = 1.0 / N;
	double k = 2.0 / M;
	double r_i = 0.0;
	double x_i = 0.0;
	double y_i = 0.0;
	bool outOfBounds = false;

	double* b = new double[sizeB];

	if (isPolar)
	{
		for (int i = 0; i < sizeB; i++)
		{
			row = i / M;
			col = i % M;
			r_i = (N - 1.0 - row)/N;
			k = r_i * pi * (M - 1.0 - col) / M;
			b[i] = (col == 0 || col == M - 2)*(1 / (h * h) - 1 / (2 * h * r_i)) * boundaryBottom + 1 / (r_i * r_i * k * k) * boundaryArc*(row == 0) + (row == N - 2)* (1 / (h * h) - 1 / (2 * h * r_i))* boundaryBottom;
		}
	}
	else
	{
		for (int i = 0; i < sizeB; i++)
		{
			row = i / (M-1);
			col = i % (M-1);
			x_i = k*(col+1)-1.0;
			y_i = h * (N-row-1);
			outOfBounds = (x_i * x_i + y_i * y_i) > 1;
			b[i] = outOfBounds + !outOfBounds*((col == 0 || col == M-2)*(1/(h*h))*boundaryArc + (row == 0)*(1/(k*k))*boundaryArc + (row == N - 2)*(1/(k*k))*boundaryBottom);
		}
	}


	return b;
}
//Matriz A
//A matrix
double* defineA(bool isPolar, int N, int M)
{
	double* A;

	return A;
}
//Algoritmo SOR
//SOR algotithm
double* SOR(double* A, double* b, double* x, int N, int M)
{
	double* x_i;

	return x_i;
}
//Función para escribir los resultados.
//Funtion to write results.
void writeResults(bool isPolar, double* x, int N, int M)
{

}