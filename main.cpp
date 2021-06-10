#include <string>
#include <unistd.h> 
#include "Model.h"

using namespace std;

int main(int argc, char *argv[]) {

	int opt = 0;
	int maxCausal = 2;
	string ldFile = "", zFile1 = "", zFile2 = "", outFile = "";
	bool variance = false;

	while ((opt = getopt(argc, argv, "l:z:c:o:")) != -1) {
		switch(opt) {
			case 'l':
				ldFile = string(optarg);
				break;
			case 'z':
				if (zFile1 == "")
					zFile1 = string(optarg);
				else
					zFile2 = string(optarg);
				break;
			case 'c':
				maxCausal = atoi(optarg);
				break;
			case 'o':
				outFile = string(optarg);
		}
	}

	Model model(ldFile, zFile1, zFile2, outFile, maxCausal, variance);
	model.getPosteriors();

	return 0;
}