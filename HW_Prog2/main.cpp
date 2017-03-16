#include <iostream>
#include <sstream>
#include <time.h>
#include <limits.h>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace std;

//Standard Matrix Multiplication
void RunStandardMatrixMultiplication(const int& dimension, const vector<int>& matA, const vector<int>& matB, vector<int>& matC) {
	int i=0, j=0, k=0, temp=0;
	for (k = 0; k<dimension; k++) {
		for (i = 0; i<dimension; i++) {
			temp = matA[i*dimension + k];
			for (j = 0; j<dimension; j++)
				matC[i*dimension + j] += temp * matB[k*dimension + j];
		}
	}
}

void ProduceOutput(const vector<int>& matC, const int& dimension) {
	for (int i = 0; i < dimension; i++) {
		cout << matC[i*dimension + i] <<endl;
	}
}

int main(int argc, char** argv) {
	if (argc != 4) {
		cout << "Incorrect input args list" << endl;
	}

	//Read the input file for matrix elements.
	ifstream file;
	file.open(argv[3]);
	if (!file.is_open()) return -1;
	vector<int> matrixElements;
	int elementInput;
	while (file >> elementInput)
	{
		matrixElements.push_back(elementInput);
	}
	file.close();

	//Read the matrix dimension from arguments.
	istringstream ss1(argv[2]);
	int dimension;
	if (!(ss1 >> dimension))
		cerr << "Invalid number for dimension" << argv[2] << '\n';
	
	vector<int> matrixA(dimension*dimension);
	vector<int> matrixB(dimension*dimension);
	vector<int> matrixC(dimension*dimension);

	//Initialise the matrices from the input file
	vector<int>::iterator rangeIterator;
	rangeIterator = matrixElements.begin() + (dimension*dimension);
	int  i = 0;
	for (vector<int>::iterator it = matrixElements.begin(); it != rangeIterator; ++it) {
		matrixA[i++] = *it;
	}
	i = 0;
	for (vector<int>::iterator it = rangeIterator; it != matrixElements.end(); ++it) {
		matrixB[i++] = *it;
	}

	//Standard Matrix Multilication
	RunStandardMatrixMultiplication(dimension, matrixA, matrixB, matrixC);

	//Output
	ProduceOutput(matrixC, dimension);
	return 0;
}