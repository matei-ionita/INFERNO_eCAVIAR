#include "Model.h"
#include "utils.h"
#include<iostream>
#include<vector>
#include<armadillo>
#include<math.h>

using namespace arma;
using namespace std;

Model::Model(string ldFile, string zFile1, string zFile2, string outFile, int maxCausal, bool variance){
	// constructor for the class: just initializing a bunch of attributes
	this->sigma.load(ldFile);
	this->nSNP = this->sigma.n_rows;

	this->maxCausal = maxCausal;
	this->outFile = outFile;

	this->snpNames = new string [this->nSNP];

	mat z1(nSNP,1), z2(this->nSNP,1);
	readZ(zFile1, z1, snpNames);
	readZ(zFile2, z2, snpNames);
	this->z1 = z1;
	this->z2 = z2;

	// log likelihood initialized with -inf
	this->lik1 = zeros<mat>(this->nSNP, this->nSNP) - datum::inf;
	this->lik2 = zeros<mat>(this->nSNP, this->nSNP) - datum::inf;

	this->post1 = zeros<vec>(this->nSNP);
	this->post2 = zeros<vec>(this->nSNP);
	this->postSNP = zeros<vec>(this->nSNP);
}

void Model::getPosteriors(double NCP) {
	int size = this->nSNP;
	this->NCP = NCP; // sigma^2 from the eCAV paper
	this->sigmaReg = pow(NCP,-1) * eye(size,size)  + this->sigma; // matrix in exponent of Bayes factor

	// GWAS
	double logTotal1 = iterateConfig(size, this->z1, this->post1, this->lik1);
	this->post1 = exp(this->post1 - logTotal1); // SNP posteriors in GWAS
	this->lik1 = exp(this->lik1 - logTotal1); // lik1[i,j] = posterior of configuration with i,j causal in GWAS

	// eQTL
	double logTotal2 = iterateConfig(size, this->z2, this->post2, this->lik2);
	this->post2 = exp(this->post2 - logTotal2); // SNP posteriors in GWAS
	this->lik2 = exp(this->lik2 - logTotal2); // lik1[i,j] = posterior of configuration with i,j causal in eQTL

	// compute posterior prob that either 1 or 2 variants are colocalized
	mat prodLik = this->lik1 * this->lik2;
	mat elemwiseLik = this->lik1 % this->lik2;

	double p2 = ( trace(prodLik) - trace(elemwiseLik) )/2;
	double p1 = accu(prodLik) - 2*p2;

	writeOutput(this->outFile, this->snpNames, this->post1, this->post2, p1, p2);
}

double Model::iterateConfig(int size, mat& z, vec& post, mat& likMat) {
	// iterate through all configs with 1 or 2 causal variants
	int nCausal;
	int last = size;
	uvec ids;
	double logTotal = size * log(1-pr);

	for (int i=0; i < size; i++) {
		if (this->maxCausal==1)
			last = i+1;

		for (int j=i; j < last; j++) {
			if (j==i) {
				nCausal = 1;
				ids = {i};
			}
			else {
				nCausal = 2;
				ids = {i,j};
			}

			mat subZ = z.rows(ids);
			mat subSigma = this->sigmaReg.submat(ids,ids);

			double logPost = getLogLik(subSigma,subZ); 
			logPost += nCausal * log(pr) + (nSNP - nCausal) * log(1-pr); // priors as in eCAV paper

			// save the (log) posterior prob of this configuration
			likMat(i,j) = logPost;
			if (j!=i) {
				likMat(j,i) = logPost;
			}

			logTotal = logAdd(logTotal, logPost);
			for (int k=0; k < nCausal; k++) {
				post(ids(k)) = logAdd(post(ids(k)), logPost);
			}
		}
	}

	return(logTotal);
}

double Model::getLogLik(mat& subSigma, mat& subZ) {
	// computes (log of) a Bayes factor: logLik of this config - logLik of null config
	mat res = trans(subZ) * inv_sympd(subSigma) * subZ;
	double exponent =  res(0,0) / 2;
	return(exponent - log(det(this->NCP * subSigma))/2 );
}

double Model::logAdd(double x, double y) {
	// computes log(e^x + e^y), careful with numeric precision
	if (x == 0)
		return y;
	if (y == 0)
		return x;

	double m = max(x,y);
	if (m - min(x,y) > 500)
		return m;

	return (m + log(1 + exp(min(x,y)-m)) );
}
