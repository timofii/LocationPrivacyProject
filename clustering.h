#ifndef ADD_CLUSTERINGH
#define ADD_CLUSTERINGH

void initSetOfLocations();
void displaySetOfLocations();
std::vector<Cluster> naiveClustering();//(int lengthOfCluster, int widthOfCluster);


//precomputations.cpp
double distance(MyPoint i, MyPoint j);
double distanceBetweenClusters(Cluster first, Cluster second);
void initDistanceMatrix();
void displayDistanceMatrix();


void initDistanceMatrixForClusters( std::vector<Cluster> vClusters);
void displayDistanceMatrixForClusters();

void initNeighbourMatrixOfClusters(double distanceThreshold);
void displayNeighbourMatrixOfClusters();

void initAmountOfNeighbour();
void displayAmountOfNeighbourVector();

void displayMapPointsBelongsToClusters();

void initT1(std::vector<Cluster> vClusters);
void displayT1();

int numberOfVariables(std::vector<Cluster> vClusters);

void initT2(std::vector<Cluster> vClusters);
void displayT2();

int indexOfPointInCluster(int index, std::vector<Cluster> vClusters);

int l0(int i, int j);
int l1(int i, int k, std::vector<Cluster> vClusters);
int l2(int i, int k, std::vector<Cluster> vClusters);
int l3(int i, int k, std::vector<Cluster> vClusters);
int l4(int i, int k, std::vector<Cluster> vClusters);

//Check on Monday
// All return number of the last constraint that was introduced.
int initTrivialConstraints(std::vector<Cluster> vClusters);
int initFirstConstraints(std::vector<Cluster> vClusters);
int initSecondConstraints(std::vector<Cluster> vClusters);
int initThirdConstraints(std::vector<Cluster> vClusters);
int initFourthConstraints(std::vector<Cluster> vClusters);

int initPrivacyConstraints(std::vector<Cluster> vClusters);

void initParametersOfLP(int c_size, int b_size, int sense_size);

void fillVector_B_and_Sence(int numberOfConstraints, int numberOfTrivialConstraints);

void initPriors();

void fillVector_C(int numberOfVariables);
void fillMatrixA();
void solve();
#endif