#include "Solver.h"
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>

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
	double k;
	double r_i = 0.0;
	double x_i = 0.0;
	double y_i = 0.0;
	double x_2;
	double y_2;
	double theta_1;
	double theta_2;
	bool closeToBorder = false;

	double* b = new double[sizeB];

	if (isPolar)
	{
		k = PI / M;
		for (int i = 0; i < sizeB; i++)
		{
			row = i / (M-1);
			col = i % (M-1);
			r_i = (N-row-1) * h;
			b[sizeB - i - 1] = (col == 0) * (2 * h * h / (k * k)) * boundaryBottom //"Left" boundary.
				+ (col == M - 2) * (2 * h * h / (k * k)) * boundaryBottom //"Right" boundary.
				+ (row == 0) * (2 * r_i * r_i + h * r_i) * boundaryArc  //Top boundary.
				+ (row == N-2) * (2 * r_i * r_i - h * r_i) * boundaryBottom; //Boundary bottom.
		}
	}
	else
	{
		k = 2.0 / M;
		for (int i = 0; i < sizeB; i++)
		{
			row = i / (M-1);
			col = i % (M-1);
			x_i = k * (col % (M - 1) + 1) - 1.0;
			y_i = h * (col % (M - 1) + 1) - 1.0;
			x_2 = sqrt(1 - y_i * y_i);
			y_2 = sqrt(1 - x_i * x_i);
			theta_1 = (k - x_2) / k;
			theta_2 = (h - y_2) / h;
			closeToBorder = (1 - x_2 < k) || (1 - y_2 < h);
			b[i] = (col == 0 || col == M - 2) * boundaryArc //Left or right point.
				+ (row == 0) * (h * h / (k * k)) * boundaryArc //Top point.
				+ (row == N - 2) * (1 / (k * k)) * boundaryBottom; //Bottom point.
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
	int pointsSet = 0;
	double h = 1.0 / N;
	double k = 0.0;
	double r_i = 0.0;
	double x_i = 0.0;
	double y_i = 0.0;
	double x_2;
	double y_2;
	double theta_1;
	double theta_2;
	bool closeToBorder = false;

	double* A = new double[sizeA];

	if (isPolar)
	{
		k = PI / M;

		for (int i = 0; i < sizeA; i++)
		{
			row = i / cols;
			col = i % cols;
			pointsSet = row / (M - 1) + 1;
			r_i = pointsSet * h;
			A[i] = (col == row) * (4 * (r_i * r_i + h * h / (k * k))) //Diagonal element.
				- (row - col == 1 && row % (M - 1) != 0) * (2 * h * h / (k * k)) //Left point.
				- (col - row == 1 && row % (M - 1) != M - 2) * (2 * h * h / (k * k)) //Right point.
				- (row - col == M - 1) * (2 * r_i * r_i - h * r_i) //Upper point
				- (col - row == M - 1) * (2 * r_i * r_i + h * r_i); //Lower point
		}
	}
	else
	{
		k = 2.0 / M;

		for (int i = 0; i < sizeA; i++)
		{
			row = i / cols;
			col = i % cols;
			pointsSet = col % (M - 1) + 1;
			x_i = k * pointsSet - 1.0;
			y_i = h * pointsSet - 1.0;
			x_2 = sqrt(1 - y_i * y_i);
			y_2 = sqrt(1 - x_i * x_i);
			theta_1 = (k - x_2) / k;
			theta_2 = (h - y_2) / h;
			closeToBorder = (1 - x_2 < k) || (1 - y_2 < h);
			A[i] = (col == row) * (2 * (h * h / (k * k) + 1)) //Diagonal element.
				- (row - col == 1 && row % (M - 1) != 0) //Left point.
				- (col - row == 1 && row % (M - 1) != M - 2) //Right point.
				- (row - col == M - 1) * (h * h / (k * k)) //Upper point.
				- (col - row == M - 1) * (h * h / (k * k)); //Bottom point.
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

	double omega = 4 / (2 + sqrt(4 - (cos(PI / M) + cos(PI / N)) * (cos(PI / M) + cos(PI / N))));
	omega = 0.7;
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
		if (difference < TOL)
		{
			cout << "Solved after " << k << " iterations" << endl;
			break;
		}
		for (int i = 0; i < sizeX; i++)
		{
			x0[i]=x1[i];
		}
	}

	return x1;
}
//Función para calcular la serie de Fourier.
//Function for calculating Fourier series.
double* fourier(bool isPolar, int N, int M)
{
	int sizeX = (N - 1) * (M - 1);
	double* x = new double[sizeX];
	int maxN = 100000;
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
				x[i] += pow(r_i, j) * 2.0 / (j * PI) * (1 - cos(j * PI)) * sin(j * theta_i);
			}
		}
	}
	else {
		for (int i = 0; i < sizeX; i++)
		{
			x[i] = 0.0;
			row = i / (M - 1);
			col = i % (M - 1);
			x_i = k * (col % (M - 1) + 1) - 1.0;
			y_i = h * (N - row - 1);
			outOfBounds = (x_i * x_i + y_i * y_i) > 1;
			theta_i = abs(atan(y_i / x_i));
			for (int j = 1; j <= maxN; j++)
			{
				x[i] += pow(x_i*x_i+y_i*y_i, j/2.0) * 2.0 / (j * PI) * (1 - cos(j * PI)) * sin(j * theta_i);
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

		for (int i = 0; i < sizeX; i++)
		{
			row = i / (M - 1);
			col = i % (M - 1);
			r_i = (N - 1.0 - row) / N;
			theta_i = 180 * (M - 1.0 - col) / M;
			file << to_string(r_i) + "," + to_string(theta_i) + "," + to_string(x[sizeX - i - 1]) + "," + to_string(x_fourier[i]) + "," + to_string(abs(x[sizeX - i - 1]-x_fourier[i])) + "\n";
		}

		file.close();
	}
	else
	{
		double* x_fourier = fourier(isPolar, N, M);

		ofstream file{ "results/cartesians-N-" + to_string(N) + "-M-" + to_string(M) + ".csv" };

		file << "x_i,y_i,u,fourier,diff\n";

		for (int i = 0; i < sizeX; i++)
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