#include <iostream>
#include "qif"
#include <vector>
#include <list>
#include <math.h>

#include "point.h"
#include "cluster.h"
#include "clustering.h"

#include "UF.h"


using namespace qif;

using namespace std;


typedef qif::lp::MatrixEntry<double> ME;

static list<ME> entries; 
static qif::lp::LinearProgram<double> lp;


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


/*
templateclass<5788> a;
a.blabla();	
*/
	UF UnionFindSet(100);
	UnionFindSet.merge(0, 1);
	UnionFindSet.merge(0, 2);
	UnionFindSet.merge(0, 3);
	UnionFindSet.merge(0, 4);
	UnionFindSet.merge(5, 6);
	UnionFindSet.merge(5, 7);
	UnionFindSet.merge(5, 8);
//cout<<UnionFindSet.size(5)<<endl;

	


	cout << "---------------+++++++++++++-----------++++++++++++++----------------" << endl;
	double neighbourThreshold = 1.5;
	vector<Cluster> vClusters = InitBasic(neighbourThreshold);

	vector<Cluster> V = PreClustering();

	int numConst = InitConstraints(vClusters);
	int numVariabl = numberOfVariables(vClusters);
	initParametersOfLP(numVariabl, numConst);
	fill_Objects_Of_LP(numConst);
	double executionTime = solve("min", "simplex_dual");
	cout << " Time: "<< executionTime << endl;
	showInfo();

	return 0;
}

