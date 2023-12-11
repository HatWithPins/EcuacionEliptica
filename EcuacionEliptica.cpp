#include "EcuacionEliptica.h"

using namespace std;

int main(int argc, char** argv)
{
	string argument;
	size_t pos;
	size_t check;
	double boundaryBottom = 0.0;
	double boundaryArc = 1.0;
	int N = 3;
	int M = 6;
	int expectedArguments = 5;
	vector<string> expectedArgumentsList = { "boundary_bottom=", "boundary_arc=", "N=", "M"};

	//Se comprueba el número de parámetros recibidos. 
	//Checking the number of received parameters.
	if (argc > expectedArguments)
	{
		cout << "Error, expected " << expectedArguments - 1 << ", but received " << argc - 1 << endl;
		return -1;
	}

	//Se comprueban los nombres de los parámetros.
	//Checking if parameter names are correct.
	try
	{
		for (int i = 1; i < argc; i++)
		{
			argument = argv[i];
			check = argument.find(expectedArgumentsList[i - 1]);
			if (check < 0 || check > argument.size())
			{
				vector<string> exceptionVector = { argument, expectedArgumentsList[i - 1] };
				throw(exceptionVector);
			}
			pos = argument.find("=");

			try
			{
				if (argument.substr(0,pos) == "boundary_bottom")
				{
					boundaryBottom = stod(argument.substr(pos + 1));
				}
				else if (argument.substr(0,pos) == "boundary_arc")
				{
					boundaryArc = stod(argument.substr(pos + 1));
				}
				else if (argument.substr(0,pos) == "N")
				{
					N = stoi(argument.substr(pos + 1));
				}
				else if (argument.substr(0,pos) == "M")
				{
					M = stoi(argument.substr(pos + 1));
				}
			}
			catch (const std::exception& e)
			{
				std::cerr << e.what() << '\n';
				return -1;
			}
		}
	}
	catch (vector<string> errorVector)
	{
		cout << "Error, expected " << errorVector[1] << "something, but received " << errorVector[0] << endl;
		return -1;
	}

	int sizeX = (N - 1) * (M - 1);
	double* x = new double[sizeX];
	for (int i = 0; i < sizeX; i++)
	{
		x[i] = 1.0;
	}
	double* polarB = defineB(true, boundaryBottom, boundaryArc, N, M);
	double* polarA = defineA(true, N, M);
	double* polarSolution = SOR(polarA, polarB, x, N, M);
	double* cartesianB = defineB(false, boundaryBottom, boundaryArc, N, M);
	double* cartesianA = defineA(false, N, M);
	double* cartesianSolution = SOR(cartesianA, cartesianB, x, N, M);

	writeResults(true, polarSolution, N, M);
	writeResults(false, cartesianSolution, N, M);

	return 0;
}
