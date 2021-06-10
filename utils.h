#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <vector>
#include <armadillo>

using namespace std;

void readZ(string file, arma::mat& z, string * snpNames);
void writeOutput(string file, string * snpNames, arma::vec& post1, arma::vec&post2, double p1, double p2);

#endif