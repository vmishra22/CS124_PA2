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

int main(int argc, char** argv) {
	if (argc != 4) {
		cout << "Incorrect input args list" << endl;
	}

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

	istringstream ss1(argv[2]);
	int dimension;
	if (!(ss1 >> dimension))
		cerr << "Invalid number for dimension" << argv[2] << '\n';
	


	return 0;
}