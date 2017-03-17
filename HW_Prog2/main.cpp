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
#include "func.h"

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

vector<int> RunStrassenMatrixMultiplication(const int& originalDimension, int& dimension, int level, vector<int>& matA, int rowAS, int colAS, 
										  vector<int>& matB, int rowBS, int colBS) {
	vector<int> matC(dimension*dimension);
	if (dimension == 1) {
		matC[0] = matA[rowAS*originalDimension + colAS] * matB[rowBS*originalDimension + colBS];
	}
	else {
		int newDimension = dimension / 2;
		level++;
		
		int nMatrixCoords = originalDimension / ((int)pow(2, level));
		vector<int> tempC1(newDimension*newDimension);
		tempC1 = RunStrassenMatrixMultiplication(originalDimension, newDimension, level, matA, rowAS, colAS, matB, rowBS, colBS) +
			RunStrassenMatrixMultiplication(originalDimension, newDimension, level, matA, rowAS, colAS+nMatrixCoords, matB, rowBS+nMatrixCoords, colBS);

		vector<int> tempC2(newDimension*newDimension);
		tempC2 = RunStrassenMatrixMultiplication(originalDimension, newDimension, level, matA, rowAS, colAS, matB, rowBS, colBS+nMatrixCoords) +
			RunStrassenMatrixMultiplication(originalDimension, newDimension, level, matA, rowAS, colAS+nMatrixCoords, matB, rowBS+nMatrixCoords, colBS+nMatrixCoords);

		vector<int> tempC3(newDimension*newDimension);
		tempC3 = RunStrassenMatrixMultiplication(originalDimension, newDimension, level, matA, rowAS+nMatrixCoords, colAS, matB, rowBS, colBS) +
			RunStrassenMatrixMultiplication(originalDimension, newDimension, level, matA, rowAS+nMatrixCoords, colAS+nMatrixCoords, matB, rowBS+nMatrixCoords, colBS);

		vector<int> tempC4(newDimension*newDimension);
		tempC4 = RunStrassenMatrixMultiplication(originalDimension, newDimension, level, matA, rowAS+nMatrixCoords, colAS, matB, rowBS, colBS+nMatrixCoords) +
			RunStrassenMatrixMultiplication(originalDimension, newDimension, level, matA, rowAS+nMatrixCoords, colAS+nMatrixCoords, matB, rowBS+nMatrixCoords, colBS+nMatrixCoords);
		
		for (int i = 0; i < newDimension; i++) {
			for (int j = 0; j < newDimension; j++) {
				matC[i*dimension + j] = tempC1[i*newDimension + j];
			}
		}
		for (int i = 0; i < newDimension; i++) {
			for (int j = newDimension; j < dimension; j++) {
				matC[i*dimension + j] = tempC2[i*newDimension + (j-newDimension)];
			}
		}
		for (int i = newDimension; i < dimension; i++) {
			for (int j = 0; j < newDimension; j++) {
				matC[i*dimension + j] = tempC3[(i-newDimension)*newDimension + j];
			}
		}
		for (int i = newDimension; i < dimension; i++) {
			for (int j = newDimension; j < dimension; j++) {
				matC[i*dimension + j] = tempC4[(i - newDimension)*newDimension + (j-newDimension)];
			}
		}
	}
	return matC;
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
	RunStandardMatrixMultiplication(dimension,matrixA,matrixB,matrixC);

	//Strassen Matrix Multiplication
	vector<int> returnMatrixC = RunStrassenMatrixMultiplication(dimension,dimension,0,matrixA,0,0,matrixB,0,0);

	//Output
	ProduceOutput(returnMatrixC, dimension);
	return 0;
}