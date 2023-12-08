#include "Solver.h"
#include <cmath>
#include <fstream>
#include <string>

#define PI  3.14159265359

using namespace std;
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
			k = r_i * PI / M;
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
			k = r_i * PI / M;
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
double* SOR(double* A, double* b, double* x0, int N, int M)
{
	int sizeX = (N - 1) * (M - 1);
	double* x1 = new double[sizeX];

	double omega = 4 / (2 + sqrt(4 - (cos(PI / M) + cos(PI / N) * (cos(PI / M) + cos(PI / N)))));
	double sigma;
	double difference;
	double TOL = 0.000001;
	int maxIterations = 1000;

	for (int k = 0; k < maxIterations; k++)
	{
		for (int i = 0; i < sizeX; i++)
		{
			sigma = 0.0;
			for (int j = 0; j < sizeX; j++)
			{
				if (j != i)
				{
					sigma += A[j+i*sizeX] * x0[j];
				}
			}
			x1[i] = (1 - omega) * x0[i] + omega * (b[i] - sigma) / A[i + i * sizeX];
		}
		difference = 0.0;
		for (int i = 0; i < sizeX; i++)
		{
			difference += (x1[i] - x0[i]) * (x1[i] - x0[i]);
		}
		difference = sqrt(difference);
		if (difference < TOL) break;
	}

	return x1;
}
//Función para calcular la serie de Fourier.
//Function for calculating Fourier series.
double* fourier(bool isPolar, int N, int M)
{
	int sizeX = (N - 1) * (M - 1);
	double* x = new double[sizeX];
	int maxN = 100;
	int col = 0;
	int row = 0;
	double r_i = 0.0;
	double theta_i = 0.0;
	double x_i = 0.0;
	double y_i = 0.0;
	double h = 1.0 / N;
	double k = 2.0 / M;
	bool outOfBounds = false;

	if (isPolar)
	{
		for (int i = 0; i < sizeX; i++)
		{
			x[i] = 0.0;
			row = i / (M - 1);
			col = i % (M - 1);
			r_i = (N - 1.0 - row) / N;
			theta_i = PI*(M - 1.0 - col) / M;
			for (int j = 1; j <= maxN; j++)
			{
				x[i] += pow(r_i, j) * 2 / (j * PI) * (1 - cos(j * theta_i)) * sin(j * theta_i);
			}
		}
	}
	else {
		for (int i = 0; i < sizeX; i++)
		{
			x[i] = 0.0;
			row = i / (M - 1);
			col = i % (M - 1);
			x_i = k * (col + 1) - 1.0;
			y_i = h * (N - row - 1);
			outOfBounds = (x_i * x_i + y_i * y_i) > 1;
			theta_i = atan(y_i / x_i);
			for (int j = 1; j <= maxN; j++)
			{
				x[i] += !outOfBounds*pow(x_i*x_i+y_i*y_i, j/2.0) * 2 / (j * PI) * (1 - cos(j * theta_i)) * sin(j * theta_i);
			}
		}
	}

	return x;
}
//Función para escribir los resultados.
//Funtion to write results.
void writeResults(bool isPolar, double* x, int N, int M)
{
	int sizeX = (N - 1) * (M - 1);
	int col = 0;
	int row = 0;
	double r_i = 0.0;
	double theta_i = 0.0;
	double x_i = 0.0;
	double y_i = 0.0;
	double h = 1.0 / N;
	double k = 2.0 / M;

	if (isPolar)
	{
		double* x_fourier = fourier(isPolar, N, M);

		ofstream file{ "results/polar-N-" + to_string(N) + "-M-" + to_string(M) + ".csv"};

		file << "r_i,theta_i,u,fourier,diff\n";

		for (int i = 0; i < sizeX + 1; i++)
		{
			row = i / (M - 1);
			col = i % (M - 1);
			r_i = (N - 1.0 - row) / N;
			theta_i = PI * (M - 1.0 - col) / M;
			file << to_string(r_i) + "," + to_string(theta_i) + "," + to_string(x[i]) + "," + to_string(x_fourier[i]) + "," + to_string(abs(x[i]-x_fourier[i])) + "\n";
		}

		file.close();
	}
	else
	{
		double* x_fourier = fourier(isPolar, N, M);

		ofstream file{ "results/cartesians-N-" + to_string(N) + "-M-" + to_string(M) + ".csv" };

		file << "x_i,y_i,u,fourier,diff\n";

		for (int i = 0; i < sizeX + 1; i++)
		{
			row = i / (M - 1);
			col = i % (M - 1);
			x_i = k * (col + 1) - 1.0;
			y_i = h * (N - row - 1);
			file << to_string(x_i) + "," + to_string(y_i) + "," + to_string(x[i]) + "," + to_string(x_fourier[i]) + "," + to_string(abs(x[i] - x_fourier[i])) + "\n";
		}

		file.close();
	}
}