#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

//Definición de b
//Define b
double* defineB(bool isPolar, double boundaryBottom = 0.0, double boundaryArc = 1.0, int N = 3, int M = 6);
//Matriz A
//A matrix
double* defineA(bool isPolar, int N, int M);
//Algoritmo SOR
//SOR algotithm
double* SOR(double* A, double* b, double* x, int N, int M);
//Función para escribir los resultados.
//Funtion to write results.
void writeResults(bool isPolar, double* x, int N, int M);
#endif