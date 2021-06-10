#include <string>
#include <vector>
#include <fstream> 
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "utils.h"

using namespace std;

void readZ(string file, arma::mat& z, string * snpNames) {
	int index = 0;
	double score = 0;
	string name = "";
	string line = "";
	ifstream fin(file.c_str(), ifstream::in);
	
	while( getline(fin, line) ){
		istringstream iss(line);
		iss >> name;
		iss >> score;
		snpNames[index] = name; 
		z(index,0) = score;
		index++;
	}
	fin.close();
}

void writeOutput(string file, string * snpNames, arma::vec& post1, arma::vec& post2, double p1, double p2) {
	double p0 = 1-p1-p2;
	int size = post1.n_elem;
	arma::vec post = post1 % post2;
	double norm = accu(post);

	ofstream fout(file);

	fout << "postGWAS_raw\tposteQTL_raw\tCLPP_raw\tCLPP_norm\tgene0col\tgene1col\tgene2col\n";
	for (int i=0; i<size; i++) {
		fout << post1(i) << "\t" << post2(i) << "\t" << post(i) << "\t" << post(i)/norm << "\t" << 1-p1-p2 << "\t" << p1 << "\t" << p2 << "\n";
	}

	fout.close();
}