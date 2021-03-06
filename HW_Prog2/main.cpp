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

//Crossoversize for transition from Strassen to Standard multiplication
int crossoverSize = 32;

//Standard Matrix Multiplication
void RunStandardMatrixMultiplication(const int& dimension, const vector<int>& matA, const vector<int>& matB, vector<int>& matC) {
	int i = 0, j = 0, k = 0, temp = 0;
	for (k = 0; k<dimension; k++) {
		for (i = 0; i<dimension; i++) {
			temp = matA[i*dimension + k];
			for (j = 0; j<dimension; j++)
				matC[i*dimension + j] += temp * matB[k*dimension + j];
		}
	}
}

//Strassen Matrix Multiplication
vector<int> RunStrassenMatrixMultiplication(const int& originalDimension, int& dimension, int level, vector<int>& matA, vector<int>& matB) {
	vector<int> matC(dimension*dimension);
	if (dimension <= crossoverSize) {
		RunStandardMatrixMultiplication(dimension, matA, matB, matC);
	}
	else {
		int newDimension = dimension / 2;
		level++;

		//Divide the matrices A and B
		vector<int> A11(newDimension*newDimension); vector<int> A12(newDimension*newDimension);
		vector<int> A21(newDimension*newDimension); vector<int> A22(newDimension*newDimension);
		for (int i = 0; i < newDimension; i++) {
			for (int j = 0; j < newDimension; j++) {
				A11[i*newDimension + j] = matA[i*dimension + j];
			}
		}
		for (int i = 0; i < newDimension; i++) {
			for (int j = 0; j < newDimension; j++) {
				A12[i*newDimension + j] = matA[((i*dimension) + newDimension) + j];
			}
		}
		for (int i = 0; i < newDimension; i++) {
			for (int j = 0; j < newDimension; j++) {
				A21[i*newDimension + j] = matA[((i + newDimension)*dimension) + j];
			}
		}
		for (int i = 0; i < newDimension; i++) {
			for (int j = 0; j < newDimension; j++) {
				A22[i*newDimension + j] = matA[(((i + newDimension)*dimension) + newDimension) + j];
			}
		}
		vector<int> B11(newDimension*newDimension); vector<int> B12(newDimension*newDimension);
		vector<int> B21(newDimension*newDimension); vector<int> B22(newDimension*newDimension);
		for (int i = 0; i < newDimension; i++) {
			for (int j = 0; j < newDimension; j++) {
				B11[i*newDimension + j] = matB[i*dimension + j];
			}
		}
		for (int i = 0; i < newDimension; i++) {
			for (int j = 0; j < newDimension; j++) {
				B12[i*newDimension + j] = matB[((i*dimension) + newDimension) + j];
			}
		}
		for (int i = 0; i < newDimension; i++) {
			for (int j = 0; j < newDimension; j++) {
				B21[i*newDimension + j] = matB[((i + newDimension)*dimension) + j];
			}
		}
		for (int i = 0; i < newDimension; i++) {
			for (int j = 0; j < newDimension; j++) {
				B22[i*newDimension + j] = matB[(((i + newDimension)*dimension) + newDimension) + j];
			}
		}
		vector<int> S1 = B12 - B22;
		vector<int> S2 = A11 + A12;
		vector<int> S3 = A21 + A22;
		vector<int> S4 = B21 - B11;
		vector<int> S5 = A11 + A22;
		vector<int> S6 = B11 + B22;
		vector<int> S7 = A12 - A22;
		vector<int> S8 = B21 + B22;
		vector<int> S9 = A11 - A21;
		vector<int> S10 = B11 + B12;

		vector<int> P1 = RunStrassenMatrixMultiplication(originalDimension, newDimension, level, A11, S1);

		vector<int> P2 = RunStrassenMatrixMultiplication(originalDimension, newDimension, level, S2, B22);

		vector<int> P3 = RunStrassenMatrixMultiplication(originalDimension, newDimension, level, S3, B11);

		vector<int> P4 = RunStrassenMatrixMultiplication(originalDimension, newDimension, level, A22, S4);

		vector<int> P5 = RunStrassenMatrixMultiplication(originalDimension, newDimension, level, S5, S6);

		vector<int> P6 = RunStrassenMatrixMultiplication(originalDimension, newDimension, level, S7, S8);

		vector<int> P7 = RunStrassenMatrixMultiplication(originalDimension, newDimension, level, S9, S10);

		//C11=P5+P4-P2+P6
		vector<int> tempC1(newDimension*newDimension);
		tempC1 = P5 + P4 - P2 + P6;
		
		//C12=P1+P2
		vector<int> tempC2(newDimension*newDimension);
		tempC2 = P1 + P2;
	
		//C21=P3+P4
		vector<int> tempC3(newDimension*newDimension);
		tempC3 = P3 + P4;

		//C22=P5+P1-P3-P7
		vector<int> tempC4(newDimension*newDimension);
		tempC4 = P5 + P1 - P3 - P7;
		
		for (int i = 0; i < newDimension; i++) {
			for (int j = 0; j < newDimension; j++) {
				matC[i*dimension + j] = tempC1[i*newDimension + j];
			}
		}
		for (int i = 0; i < newDimension; i++) {
			for (int j = newDimension; j < dimension; j++) {
				matC[i*dimension + j] = tempC2[i*newDimension + (j - newDimension)];
			}
		}
		for (int i = newDimension; i < dimension; i++) {
			for (int j = 0; j < newDimension; j++) {
				matC[i*dimension + j] = tempC3[(i - newDimension)*newDimension + j];
			}
		}
		for (int i = newDimension; i < dimension; i++) {
			for (int j = newDimension; j < dimension; j++) {
				matC[i*dimension + j] = tempC4[(i - newDimension)*newDimension + (j - newDimension)];
			}
		}
	}
	return matC;
}

void ProduceOutput(const vector<int>& matC, const int& dimension) {
	for (int i = 0; i < dimension; i++) {
		cout << matC[i*dimension + i] << endl;
	}
}

void checkIfDimensionPowerofTwo(int dim, int checkNum, bool* retVal, int* modifiedDim) {
	if (dim == checkNum * 2)
	{
		*retVal = true;
		return;
	}
	else if (dim < checkNum * 2)
	{
		*modifiedDim = checkNum * 2;
		*retVal = false;
		return;
	}
	else
		checkIfDimensionPowerofTwo(dim, checkNum * 2, retVal, modifiedDim);

	return;
}

int main(int argc, char** argv) {
	if (argc != 4) {
		cout << "Incorrect input args list" << endl;
	}
	//using namespace chrono;
	//chrono::steady_clock::time_point tStart;
	//chrono::steady_clock::time_point tEnd;

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

	int modifiedDimension = dimension;
	bool isPowerofTwo = false;
	checkIfDimensionPowerofTwo(dimension, 1, &isPowerofTwo, &modifiedDimension);

	vector<int> matrixA(modifiedDimension*modifiedDimension);
	vector<int> matrixB(modifiedDimension*modifiedDimension);

	//Initialise the matrices from the input file
	//Pad with zeros when dimension is not power of 2.
	if (modifiedDimension != dimension) {
		vector<int>::iterator rangeIterator;
		rangeIterator = matrixElements.begin() + (dimension*dimension);
		int  i = 0;
		int counter = 0;
		for (vector<int>::iterator it = matrixElements.begin(); it != rangeIterator; ++it) {
			matrixA[i++] = *it;
			counter++;
			if (counter%dimension == 0) {
				for (int j = dimension; j < modifiedDimension; j++) {
					matrixA[i++] = 0;
				}
			}
		}
		i = 0;
		counter = 0;
		for (vector<int>::iterator it = rangeIterator; it != matrixElements.end(); ++it) {
			matrixB[i++] = *it;
			counter++;
			if (counter%dimension == 0) {
				for (int j = dimension; j < modifiedDimension; j++) {
					matrixB[i++] = 0;
				}
			}
		}
		for (i = dimension; i < modifiedDimension; i++) {
			for (int j = 0; j < modifiedDimension; j++) {
				matrixA[i*modifiedDimension + j] = 0;
				matrixB[i*modifiedDimension + j] = 0;
			}
		}
	}
	else {
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
	}

	//Strassen Matrix Multiplication
	//tStart = steady_clock::now();
	vector<int> returnMatrixC = RunStrassenMatrixMultiplication(modifiedDimension, modifiedDimension, 0, matrixA, matrixB);
	//tEnd = steady_clock::now();
	/*cout << "Time for STRASSEN multiplication with dimension: " << dimension << " crossover: " << crossoverSize
		<< " : " << chrono::duration_cast<chrono::milliseconds>(tEnd - tStart).count() << " seconds" << endl;*/

	//Unpadding the zeros in case when dimension is not power of 2
	if (modifiedDimension != dimension) {
		vector<int> matrixC(dimension*dimension);
		for (int i = 0; i < dimension; i++) {
			for (int j = 0; j < modifiedDimension; j++) {
				if( j<dimension)
					matrixC[i*dimension + j] = returnMatrixC[i*modifiedDimension + j];
			}
		}
		ProduceOutput(matrixC, dimension);
	}
	else
		ProduceOutput(returnMatrixC, dimension);
	return 0;
}