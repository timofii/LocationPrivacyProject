#ifndef ADD_CLUSTERINGH
#define ADD_CLUSTERINGH

void initSetOfLocations();

int getNumberOfPointInVectorOfPoints(MyPoint p);
void displaySetOfLocations();
std::vector<Cluster> naiveClustering();//(int lengthOfCluster, int widthOfCluster);
std::vector<Cluster> PreClustering();



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


// All return number of the last constraint that was introduced.
int initTrivialConstraints(std::vector<Cluster> vClusters);
int initFirstConstraints(std::vector<Cluster> vClusters);
int initSecondConstraints(std::vector<Cluster> vClusters);
int initThirdConstraints(std::vector<Cluster> vClusters);
int initFourthConstraints(std::vector<Cluster> vClusters);

int initPrivacyConstraints(std::vector<Cluster> vClusters);

int InitConstraints(std::vector<Cluster> vClusters);

void initParametersOfLP(int numberOfVariables, int numberOfConstraints);


void fillVector_B_and_Sence(int numberOfConstraints);

void initPriors();
void displayPriors();

void fillVector_C();
vector<Cluster> InitBasic(double distanceThresholdForNeighbourClusters);
void fill_Objects_Of_LP(int numberOfConstraints);
double solve(string extremum, string method);
void showInfo();

#endif