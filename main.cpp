#include <iostream>
#include "qif"
#include <vector>
#include <list>
#include <math.h>

#include "point.h"
#include "cluster.h"
#include "clustering.h"

using namespace qif;

using namespace std;

typedef qif::MatrixEntry<double> ME;



//static double Pi[numberOfPoints];

static list<ME> entries; 
static LinearProgram<double> lp;

/*template <int N>
class templateclass {
public:
	void blabla(){
		cout << N << endl;
	}
};
*/

int main(int argc, char const *argv[])
{

/*templateclass<5788> a;
a.blabla();	
*/
	
	initSetOfLocations();

	//displaySetOfLocations();	
	initDistanceMatrix();
	//displayDistanceMatrix();

	vector<Cluster> tempClus = naiveClustering();
	initDistanceMatrixForClusters(tempClus);
	//displayDistanceMatrixForClusters();

	initNeighbourMatrixOfClusters(1.0);
	//displayNeighbourMatrixOfClusters();
	initAmountOfNeighbour();
	//displayAmountOfNeighbourVector();
	//displayMapPointsBelongsToClusters();

	
	initT1(tempClus);
	initT2(tempClus);
	//cout << endl << numberOfVariables(tempClus)<<endl;

	//displayT1();
	//displayT2();

	
	//cout << endl << l1(4, 3, tempClus) << endl;
	
	//cout << endl << l2(4, 3, tempClus) << endl;

	//cout << endl << l3(4, 3, tempClus) << endl;
	
	//cout << endl << l4(4, 3, tempClus) << endl;
	
	initPriors();
	int numTrivConst = initTrivialConstraints(tempClus);
	initFirstConstraints(tempClus);
	initSecondConstraints(tempClus);
	initThirdConstraints(tempClus);
	initFourthConstraints(tempClus);
	
	int c_size, numVariabl, numConst, b_size, sense_size; 
	c_size = numVariabl = numberOfVariables(tempClus);
	numConst = b_size = sense_size = initPrivacyConstraints(tempClus);

	initParametersOfLP(c_size, b_size, sense_size);

	fillVector_B_and_Sence(numConst, numTrivConst);

	fillVector_C(numVariabl);

	solve();




	//fillMatrixA();

	

 	//initParametersOfLP(5, 9, 10);
 	//fillVector_B_and_Sence(5, 2);

	/*
	vector<Cluster> tempClus = naiveClustering(16, 16);
	for (int i = 0; i < tempClus.size(); ++i)
	{
		cout<<"Number § "<<i<<" Amount of Points "<< tempClus[i].size()<<endl;
		tempClus[i].displayAllPoints();
	}
	*/
	

	return 0;
}