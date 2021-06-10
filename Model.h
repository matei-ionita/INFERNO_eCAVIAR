#include <armadillo>

using namespace std;
using namespace arma;

class Model{
public:
	Model(string, string, string, string, int, bool);
	void getPosteriors(double NCP = 5.2);
	~Model() {}

private:
	int nSNP, maxCausal;
	double pr = 0.01, NCP;
	string * snpNames;
	mat sigma, sigmaReg, z1, z2, lik1, lik2;
	vec post1, post2, postSNP;
	string outFile;

	double iterateConfig(int size, mat& z, vec& post, mat& likMat);
	double getLogLik(mat& subSigma, mat& subZ);
	double logAdd(double x, double y);
};
