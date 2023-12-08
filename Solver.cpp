#include "Solver.h"

double pi = 3.14159265359;
//Definición de b
//Define b
double* defineB(bool isPolar, double boundaryBottom, double boundaryArc, int N, int M)
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
			row = i / (M-1);
			col = i % (M-1);
			r_i = (N - 1.0 - row)/N;
			k = r_i * pi / M;
			b[i] = (col == 0)*(1 / (h * h) - 1 / (2 * h * r_i)) * boundaryBottom + (col == M - 2)* (1 / (h * h) + 1 / (2 * h * r_i)) * boundaryBottom + 1 / (r_i * r_i * k * k) * boundaryArc * (row == 0) + (row == N - 2) * (1 / (h * h) - 1 / (2 * h * r_i)) * boundaryBottom;
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
	int cols = (N - 1) * (M - 1);
	int sizeA = cols * cols;
	int col = 0;
	int row = 0;
	double h = 1.0 / N;
	double k = 2.0 / M;
	double r_i = 0.0;
	double x_i = 0.0;
	double y_i = 0.0;
	bool outOfBounds = false;

	double* A = new double[sizeA];

	if (isPolar)
	{
		for (int i = 0; i < sizeA; i++)
		{
			row = i / cols;
			col = i % cols;
			r_i = (N - 1.0 - row) / N;
			k = r_i * pi / M;
			A[i] = (col == row) * (2 * (1 / (h * h) + 1 / (r_i * r_i * k * k))) //Diagonal element.
				+ (row - col == 1 && row % (M - 1) != 0) * (1 / (2 * h * r_i) - 1 / (h * h)) //Left point.
				- (col - row == 1 && row % (M - 2) != 0) * (1 / (2 * h * r_i) + 1 / (h * h)) //Right point.
				- (row > M - 2 && row - col == M - 1) * (1 / (r_i * r_i * k * k)) //Upper point
				- (row < M - 1 && col - row == M - 1) * (1 / (r_i * r_i * k * k));//Lower point
		}
	}
	else
	{
		for (int i = 0; i < sizeA; i++)
		{
			row = i / cols;
			col = i % cols;
			x_i = k * (col + 1) - 1.0;
			y_i = h * (N - row - 1);
			outOfBounds = (x_i * x_i + y_i * y_i) > 1;
			A[i] = !outOfBounds * (
				(col == row) * (2 * (1 / (h * h) + 1 / (k * k))) //Diagonal element.
				+ (row - col == 1 && row % (M - 1) != 0) * (1 / (h * h)) //Left point.
				- (col - row == 1 && row % (M - 2) != 0) * (1 / (h * h)) //Right point.
				- (row > M - 2 && row - col == M - 1) * (1 / (k * k)) //Upper point
				- (row < M - 1 && col - row == M - 1) * (1 / (k * k)));//Lower point
		}
	}

	return A;
}
//Algoritmo SOR
//SOR algotithm
double* SOR(double* A, double* b, double* x, int N, int M)
{
	int sizeX = (N - 1) * (M - 1);
	double* x_i;

	return x_i;
}
//Función para escribir los resultados.
//Funtion to write results.
void writeResults(bool isPolar, double* x, int N, int M)
{

}