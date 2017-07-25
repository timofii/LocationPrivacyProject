// Clustering naive

#include <iostream>
#include <vector>
#include <math.h>
#include <unordered_map>
#include <string>

#include "qif"
#include "CONSTANTS.h"
#include "point.h"
#include "cluster.h"
#include "clustering.h"
#include "UF.h"

using namespace std;

//remove it from here
//start

static std::vector<MyPoint> vPoints(NUMBER_OF_POINTS);

static std::vector<int> vAmountOfNeighbour(NUMBER_OF_CLUSTERS);

static double distanceMatrix[NUMBER_OF_POINTS][NUMBER_OF_POINTS];

static double distanceMatrixForClusters[NUMBER_OF_CLUSTERS][NUMBER_OF_CLUSTERS];

static bool neighbourMatrixOfClusters[NUMBER_OF_CLUSTERS][NUMBER_OF_CLUSTERS]; 

static std::unordered_map <int, int> mPointsBelongsToClusters;


//T1(k-1) = t1(k)
static std::vector<int> vT1(NUMBER_OF_CLUSTERS + 1);

static int aT2[NUMBER_OF_CLUSTERS][NUMBER_OF_CLUSTERS];
//end

typedef qif::lp::MatrixEntry<double> ME;

static list<ME> entries; 

static double Pi[NUMBER_OF_POINTS];

static qif::lp::LinearProgram<double> lp;


void initSetOfLocations()
{	/*
	 Initialize Points of the big Grid of locations 
	 */
	//cout << "Start Initialization of Lotations." << endl;
	int sizeOfBigGrid = sqrt(NUMBER_OF_POINTS);

	for(int i = 0; i < sizeOfBigGrid; ++i)
	{
		for(int j = 0; j < sizeOfBigGrid; ++j)
		{
			vPoints[ (sizeOfBigGrid * i) + j ].setX(i);
			vPoints[ (sizeOfBigGrid * i) + j ].setY(j);
			
		}
	}
	//cout << "Finish Initialization of Lotations!" << endl;
}

int getNumberOfPointInVectorOfPoints(MyPoint p)
{
	/*
	 First Initialize Points of the big Grid of locations // use  initSetOfLocations();
	 */
	int sizeOfBigGrid = sqrt(NUMBER_OF_POINTS);
	return (sizeOfBigGrid * p.getX()) + p.getY();
}

void displaySetOfLocations()
{
	for (int i = 0; i < vPoints.size(); ++i)
	{
		vPoints[ i ].display();
	}
}

//issue with width and lenght
std::vector<Cluster> naiveClustering()//(int lengthOfCluster, int widthOfCluster)
{	/*
	 Clustered the Grid by pices which are paralelepipeds.
	 Does not exist check for a size of clusters, therefore carefully choosing parameters of a function to avare an error!!!
	 */
	//cout << "Start Naive Clustering." << endl;
	//int numberOfClusters = NUMBER_OF_POINTS / (lengthOfCluster * widthOfCluster);
	std::vector<Cluster> vClusters(NUMBER_OF_CLUSTERS);

	int x = 0; 
	int y = 0;

	int sizeOfBigGrid = sqrt(NUMBER_OF_POINTS);
	int coeficientOfLinearization = ((sizeOfBigGrid - 1) / LENGTH_OF_CLUSTERS) + 1;
	for (int i = 0; i < NUMBER_OF_POINTS; ++i)
	{
		x = vPoints[ i ].getX() / WIDTH_OF_CLUSTERS;
		y = vPoints[ i ].getY() / LENGTH_OF_CLUSTERS;

		vClusters[ (coeficientOfLinearization * x) + y ].push_back(vPoints[ i ]);
		//fill unordered map:  
		//   key-> #point;
		//   	value-> #cluster
		mPointsBelongsToClusters[ i ] = (coeficientOfLinearization * x) + y;
	}
	//cout << "End Naive Clustering!" << endl;
	return vClusters;
}

std::vector<Cluster> PreClustering()
{	
	//Some parameters
	bool flag = false;
	double epsilon = 1.0;
	double A = 6;
	//-------------

	UF *pClusters;
	while(not flag)
	{
		flag = true;
		UF Clusters(NUMBER_OF_POINTS);
		for (int i = 0; (i < NUMBER_OF_POINTS - 1) and (flag); ++i){
			//cout << "position " << i << endl;
		 	for (int j = i + 1; (j < NUMBER_OF_POINTS) and (flag) ; ++j){
		 		//cout << "position " << j << endl;
		 		if(distanceMatrix[ i ][ j ] <= epsilon)
		 		{
		 			Clusters.merge(i, j);
		 			//cout << Clusters.size(i)<<  " " << Clusters.size(j) << "A = " << A << endl; 
		 			if(Clusters.size(i) > A)
		 			{
		 				epsilon = epsilon / 2.0;
		 				flag = false;

		 			}
		 		}
		
		 	}
		}

		if(flag)
			{pClusters = new UF(Clusters);}

	}
//	for(int i = 0; i < NUMBER_OF_POINTS; ++i)
//		cout << "pClusters->size " << pClusters->size(i) << endl;
//	cout << "pClusters->count " << pClusters->count() << endl;

	std::vector<Cluster> vClusters(NUMBER_OF_CLUSTERS);
	return vClusters;
}

///////////////////////////////////////
//precomputations.cpp
//////////////////////////////////////

double distance(MyPoint i, MyPoint j)
{	/*
	Calculate euclidian distance between points.
	Now without scaling factor.
	 */

	return sqrt(((i.getX() - j.getX())*(i.getX() - j.getX())) + ((i.getY() - j.getY())*(i.getY() - j.getY())));
}


void initDistanceMatrix()
{	/*
	Calculate distance between all points and fill the Distance Matrix
	 */

	//cout << "Start Initialization of Distance Matrix Between Points." << endl;
	for (int i = 0; i < NUMBER_OF_POINTS; ++i)
	{
		for (int j = 0; j < NUMBER_OF_POINTS; ++j)
		{
			distanceMatrix[ i ][ j ] = distance(vPoints[ i ], vPoints[ j ]);
		}
	}
	//cout << "Finish Initialization of Distance Matrix Between Points!" << endl;
}

void displayDistanceMatrix()
{
	cout<<"----------------Matrix of Distance For All Points ---------------"<<endl;
	for (int i = 0; i < NUMBER_OF_POINTS; ++i)
	{
		for (int j = 0; j < NUMBER_OF_POINTS; ++j)
		{
			cout<<"Point_"<<i<<" Point_"<<j<<", dist = "<<distanceMatrix[ i ][ j ]<<endl;
		}
	}
}


void initDistanceMatrixForClusters( std::vector<Cluster> vClusters)
{	/*
	Calculate distance between all clusters and fill the Distance Matrix for clusters
	 */
	double min = DBL_MAX;
	
	for (int k = 0; k < vClusters.size(); ++k)
	{
		for (int p = 0; p < vClusters.size(); ++p)
		{
			min = DBL_MAX; 
			for (int i = 0; i < vClusters[ k ].size() ; ++i)
			{
				for (int j = 0; j < vClusters[ p ].size(); ++j)
				{

					if ( distance(vClusters[ k ].getPoint(i), vClusters[ p ].getPoint(j)) <= min)
					{	
						min = distance(vClusters[ k ].getPoint(i), vClusters[ p ].getPoint(j));
					}
					//cout << "i " << i <<" j "<< j << min;
				}				
			}
			
			distanceMatrixForClusters[ k ][ p ] = min;	
		}
	}
	
}

void displayDistanceMatrixForClusters()
{
	cout<<"----------------Matrix of Distance Between Clusters ---------------"<<endl;
	for (int i = 0; i < NUMBER_OF_CLUSTERS; ++i)
	{
		for (int j = 0; j < NUMBER_OF_CLUSTERS; ++j)
		{
			cout<<"Cluster_"<<i<<" Cluster_"<<j<<", d = "<<distanceMatrixForClusters[ i ][ j ]<<endl;
		}
	}
}


void initNeighbourMatrixOfClusters(double distanceThreshold)
{	/*
	Fill the Neighbour Matrix using as the parameter the distance Threshhold.
	 */
	// initDistanceMatrixForClusters( std::vector<Cluster> vClusters);
	for (int i = 0; i < NUMBER_OF_CLUSTERS; ++i)
	{
		for (int j = 0; j < NUMBER_OF_CLUSTERS; ++j)
		{
			if( distanceMatrixForClusters[ i ][ j ] <= distanceThreshold)
				{
					neighbourMatrixOfClusters[ i ][ j ] = true;
				}
			else
				{
					neighbourMatrixOfClusters[ i ][ j ] = false;
				}
		}
	}
}

void displayNeighbourMatrixOfClusters()
{
	cout<<"----------------Adjacent Matrix ---------------"<<endl;
	for (int i = 0; i < NUMBER_OF_CLUSTERS; ++i)
	{
		for (int j = 0; j < NUMBER_OF_CLUSTERS; ++j)
		{
			cout<<"Cluster_"<<i<<" Cluster_"<<j<<", d = "<<neighbourMatrixOfClusters[ i ][ j ]<<endl;
		}
	}
}


void initAmountOfNeighbour()
{	/*
	Calculate the amount of neighbours for each cluster and fill the vector which corresponds to this amount
	 */
	// initDistanceMatrixForClusters( std::vector<Cluster> vClusters);
	// initNeighbourMatrixOfClusters(double distanceThreshold);
	int amount = 0;
	for (int i = 0; i < NUMBER_OF_CLUSTERS; ++i)
	{
		for (int j = 0; j < NUMBER_OF_CLUSTERS; ++j)
		{
			amount += neighbourMatrixOfClusters[ i ][ j ];
		}
		vAmountOfNeighbour[ i ] = amount;
		amount = 0;
	}
}

void displayAmountOfNeighbourVector()
{
	cout<<"----------------Vector of Neighbours---------------"<<endl;	
	for (int i = 0; i < NUMBER_OF_CLUSTERS; ++i)
	{
		cout<<"Cluster_"<<i<<" "<<" amount = "<<vAmountOfNeighbour[ i ]<<endl;
	}
}

void displayMapPointsBelongsToClusters()
{
	cout<<"----------------Map of Points which Belongs To Clusters---------------"<<endl;	
	for (int i = 0; i < NUMBER_OF_POINTS; ++i)
	{
		cout<<"Point_"<<i<<" #cluster: "<<mPointsBelongsToClusters[ i ]<<endl;
	}
}



//////////////////////////////////////////////////////////////
////////////////////Calculation of Variables and Linearization
/////////////////////////////////////////////////////////////
////////// Ostorojno poshel govno kod!!!!

// new function from 0 to K
void initT1(std::vector<Cluster> vClusters)
{
	vT1[ 0 ] = 0;
	for (int k = 0; k < NUMBER_OF_CLUSTERS; ++k)
	{
		vT1[ k + 1 ] = vT1[ k ] + (vClusters[ k ].size() * (NUMBER_OF_CLUSTERS - vAmountOfNeighbour[ k ]));
	}	
}


void displayT1()
{
	cout<<"----------------T1---------------"<<endl;	
	for (int i = 0; i < NUMBER_OF_CLUSTERS + 1; ++i)
	{
		cout << vT1[i] << endl;
	}
}


int numberOfVariables(std::vector<Cluster> vClusters)
{
	/*
	Have to be initialized before throw runing  initT1(std::vector<Cluster> vClusters)!!!!!!!!
	*/
	return NUMBER_OF_POINTS_2 + (2 * vT1[ NUMBER_OF_CLUSTERS ]) + (2 * NUMBER_OF_CLUSTERS * NUMBER_OF_POINTS); 

}

void initT2(std::vector<Cluster> vClusters)
{
	/*
	Be careful aT2 is NOT symetric!!!! 
	Right order is k == k0, k' == k1
	Have to be initialized before throw runing initNeighbourMatrixOfClusters()!!!!!!!!
	*/
	int counter = 0;
	for (int k0 = 0; k0 < NUMBER_OF_CLUSTERS; ++k0)
	{
		for (int k1 = 0; k1 < NUMBER_OF_CLUSTERS; ++k1)
		{
			aT2[k0][k1] = 0;
			for (int i = 0; i < k1; ++i)
			{
				if(neighbourMatrixOfClusters[ k0 ][ i ] == false)
					counter++;
			}
			aT2[k0][k1] = vClusters[ k0 ].size() * counter  ;
			counter = 0;
		}
	}
}

void displayT2()
{
	/*
	Be careful aT2 is NOT symetric!!!!
	Right order is k == k0, k' == k1
	*/
	cout<<"----------------T2---------------"<<endl;	
	for (int i = 0; i < NUMBER_OF_CLUSTERS; ++i)
	{
		for (int j = 0; j < NUMBER_OF_CLUSTERS; ++j)
		{
			cout << aT2[i][j] << " ";
		}
		cout << endl;
	}
}

int indexOfPointInCluster(int index, std::vector<Cluster> vClusters)
{
	/*
	REWRITE using unordered map with custom class
	*/
	if(index == NUMBER_OF_POINTS)
		cout<< "ERROR: parameter index is not valid in function indexOfPointInCluster()" << endl;
	int x0 = vPoints[index].getX();
	int y0 = vPoints[index].getY();

	//cout << "index " << index << "index "<<vPoints.size()<< endl;	
	int indexOfCluster = mPointsBelongsToClusters[ index ];
	//cout << "--"<<indexOfCluster << "--" <<endl;
	int x1 = 0;
	int y1 = 0;

	for (int i = 0; i < vClusters[ indexOfCluster ].size(); ++i)
	{
		//cout << x0 << " " << y0 << " == " << x1 << " "<< y1 << endl;
		x1 = vClusters[ indexOfCluster ].getPoint(i).getX();
		y1 = vClusters[ indexOfCluster ].getPoint(i).getY();

		if( (x0 == x1) && (y0 == y1) )
			return i;
	}
	cout << " Problem with indexOfPointInCluster !!! Do not find Point in the for loop !!!" << endl;
	return -1;
}

int l0(int i, int j)
{ 
	return ((NUMBER_OF_POINTS * i) + j);
}

int l1(int i, int k1, std::vector<Cluster> vClusters)
{
	if(i == NUMBER_OF_POINTS)
		cout<< "ERROR: problem with parameter i in function l1()" << endl;
	/*
	Be careful aT2 is NOT symetric!!!!
	Right order is k == k0, k' == k1
	*/
	int k0 = mPointsBelongsToClusters[i];

	return NUMBER_OF_POINTS_2 + vT1[ k0 ] + aT2[ k0 ][ k1 ] + indexOfPointInCluster(i, vClusters);
}

int l2(int i, int k1, std::vector<Cluster> vClusters)
{
	int k0 = mPointsBelongsToClusters[i];
	return NUMBER_OF_POINTS_2 + vT1[ NUMBER_OF_CLUSTERS ] 
							  + vT1[ k0 ] + aT2[ k0 ][ k1 ] 
							  + indexOfPointInCluster(i, vClusters);
}

int l3(int i, int k, std::vector<Cluster> vClusters)
{

	return  NUMBER_OF_POINTS_2 
							  + (2 * vT1[ NUMBER_OF_CLUSTERS ]) 
							  + (k * NUMBER_OF_POINTS) 
							  + i;
}

int l4(int i, int k, std::vector<Cluster> vClusters)
{

	return NUMBER_OF_POINTS_2 
							  + (2 * vT1[ NUMBER_OF_CLUSTERS ]) 
							  + (NUMBER_OF_CLUSTERS * NUMBER_OF_POINTS)
							  + (k * NUMBER_OF_POINTS) 
							  + i;
}

int initTrivialConstraints(std::vector<Cluster> vClusters)
{
	/*
	Check this function !!!
	*/
	//cout << "Init Trivial Constraints! " << endl;
	int amountOfVariables = numberOfVariables(vClusters);

	for (int x = 0; x < NUMBER_OF_POINTS; ++x)
	{
		for (int y = 0; y < amountOfVariables; ++y)
		{
			if( ( (x * NUMBER_OF_POINTS) <= y) && (((x + 1) * NUMBER_OF_POINTS) > y) )
			{
				entries.push_back( ME(x, y, 1) );
			}
			else
			{
				entries.push_back( ME(x, y, 0) );
			}
		}
	}
	
	return NUMBER_OF_POINTS;
}


 int initFirstConstraints(std::vector<Cluster> vClusters)
 {
 	//cout << "Init 1th Constraints! " << endl;
 	int i = 0, j = 0, k0 = 0, k1 = 0;
 	for (int x = NUMBER_OF_POINTS; x < NUMBER_OF_POINTS + NUMBER_OF_POINTS_2; ++x)
 	{
 		i = ((x - NUMBER_OF_POINTS)/ NUMBER_OF_POINTS) + 1;
 		j = x - (i * NUMBER_OF_POINTS);
 		
 		k0 = mPointsBelongsToClusters[ i - 1 ];
 		k1 = mPointsBelongsToClusters[ j ];

 		if(!neighbourMatrixOfClusters[ k0 ][ k1 ])
 		{
 			entries.push_back( ME(x, x - NUMBER_OF_POINTS, -1) );
 			entries.push_back( ME(x, l1(i - 1, k1, vClusters), 1) );
 		}
 		else
 		{
 			entries.push_back( ME(x, x - NUMBER_OF_POINTS, -1) );
 			entries.push_back( ME(x, l1(i - 1, k1, vClusters), 1) );
 		}
 		//if( (x - NUMBER_OF_POINTS) == l1(i - 1, k1, vClusters))
 		//	{cout<< "----------HHHH--------------" << endl;}
 	}
 	return NUMBER_OF_POINTS + NUMBER_OF_POINTS_2;
 }

 int initSecondConstraints(std::vector<Cluster> vClusters)
 {
 	//cout << "Init 2nd Constraints! " << endl;
 	int i = 0, j = 0, k0 = 0, k1 = 0;
 	for (int x = NUMBER_OF_POINTS + NUMBER_OF_POINTS_2; x < NUMBER_OF_POINTS + (2 * NUMBER_OF_POINTS_2) ; ++x)
 	{

 		i = ((x - NUMBER_OF_POINTS_2 - NUMBER_OF_POINTS) / NUMBER_OF_POINTS) + 1;
 		j = x - (i * NUMBER_OF_POINTS) - NUMBER_OF_POINTS_2;

 		k0 = mPointsBelongsToClusters[ i - 1 ];
 		k1 = mPointsBelongsToClusters[ j ];

 		if(!neighbourMatrixOfClusters[ k0 ][ k1 ])
 		{
 			entries.push_back( ME(x, x - NUMBER_OF_POINTS_2 - NUMBER_OF_POINTS, 1) );
 			entries.push_back( ME(x, l2(i - 1, k1, vClusters), -1) );
 		}
 		else
 		{
 			entries.push_back( ME(x, x - NUMBER_OF_POINTS_2 - NUMBER_OF_POINTS, 0) );
 			entries.push_back( ME(x, l2(i - 1, k1, vClusters), 0) );
 		}
 	//	 if( (x - NUMBER_OF_POINTS_2 - NUMBER_OF_POINTS) == l2(i - 1, k1, vClusters))
 	//		{cout<< "----------AAA--------------" << endl;}
 	
 	}
 	return NUMBER_OF_POINTS + (2 * NUMBER_OF_POINTS_2);
 }

int initThirdConstraints(std::vector<Cluster> vClusters)
{
	//cout << "Init 3rd Constraints! " << endl;
	int i = 0, j = 0, k1 = 0;
	for (int x = NUMBER_OF_POINTS + (2 * NUMBER_OF_POINTS_2); x <  NUMBER_OF_POINTS + (3 * NUMBER_OF_POINTS_2); ++x)
 	{
 		i = ((x - (2 * NUMBER_OF_POINTS_2) - NUMBER_OF_POINTS) / NUMBER_OF_POINTS) + 1;
 		j = x - (i * NUMBER_OF_POINTS) - (2 * NUMBER_OF_POINTS_2);

 		entries.push_back( ME(x, x - (2 * NUMBER_OF_POINTS_2) - NUMBER_OF_POINTS, -1) );
 		
 		k1 = mPointsBelongsToClusters[ i - 1 ];
 		entries.push_back( ME(x, l3(j, k1, vClusters), 1) );

 		//if( (x - (2 * NUMBER_OF_POINTS_2) - NUMBER_OF_POINTS) == l3(j, k1, vClusters))
 		//	{cout<< "----------ZZZ--------------" << endl;}
 	}
 	return NUMBER_OF_POINTS + (3 * NUMBER_OF_POINTS_2);
}

int initFourthConstraints(std::vector<Cluster> vClusters)
{
	//cout << "Init 4th Constraints! " << endl;
	int i = 0, j = 0, k1 = 0;

	for (int x = NUMBER_OF_POINTS + (3 * NUMBER_OF_POINTS_2); x <  NUMBER_OF_POINTS + (4 * NUMBER_OF_POINTS_2); ++x)
 	{
 		i = ((x - (3 * NUMBER_OF_POINTS_2) - NUMBER_OF_POINTS) / NUMBER_OF_POINTS) + 1;
 		j = x - (3 * NUMBER_OF_POINTS_2) - (i * NUMBER_OF_POINTS);
 		

 		entries.push_back( ME(x, x - (3 * NUMBER_OF_POINTS_2) - NUMBER_OF_POINTS, 1) );

 		k1 = mPointsBelongsToClusters[ i - 1 ];
 		entries.push_back( ME(x, l4(j, k1, vClusters), -1) );

 		if( (x - (3 * NUMBER_OF_POINTS_2) - NUMBER_OF_POINTS) == l4(j, k1, vClusters))
 			{cout<< "----------XXX--------------" << endl;}
 	}
 	return NUMBER_OF_POINTS + (4 * NUMBER_OF_POINTS_2);
}

// Check the correctnes according to formulas. Finish writing the code with out of distanse
// Write down the distance for points and clusters for exp(i,i') etc.
//Run and Debug this shit!!!!!!!!!  
int initPrivacyConstraints(std::vector<Cluster> vClusters)
{
	//cout << "Init Privacy Constraints! " << endl;
	int x0 = 0,
		x1 = 0,
		y = 0;
	int amount_of_constraints = NUMBER_OF_POINTS + (4 * NUMBER_OF_POINTS_2);
	for (int k0 = 0; k0 < NUMBER_OF_CLUSTERS; ++k0)
	{
		for (int k1 = 0; k1 < NUMBER_OF_CLUSTERS; ++k1)
		{
			for (int k2 = 0; k2 < NUMBER_OF_CLUSTERS; ++k2)
			{

				if (neighbourMatrixOfClusters[ k0 ][ k1 ] && 
					(neighbourMatrixOfClusters[ k0 ][ k2 ] || neighbourMatrixOfClusters[ k1 ][ k2 ]))
				{

					for (int i0 = 0; i0 < vClusters[ k0 ].size(); ++i0)
					{
						for (int i1 = 0; i1 < vClusters[ k1 ].size(); ++i1)
						{
							for (int j = 0; j < vClusters[ k2 ].size(); ++j)
							{
								amount_of_constraints++;

									x0 = getNumberOfPointInVectorOfPoints(vClusters[ k0 ].getPoint(i0));
									x1 = getNumberOfPointInVectorOfPoints(vClusters[ k1 ].getPoint(i1));
									y = getNumberOfPointInVectorOfPoints(vClusters[ k2 ].getPoint(j));

								if( l0(x0, y) != l0(x1, y))
								{
									entries.push_back( ME(amount_of_constraints, l0(x0, y), 1) );
									entries.push_back( ME(amount_of_constraints, l0(x1, y), -exp(distanceMatrix[x0][x1]) ) );
								}
								//Problems HERE!!!!!!!
								//if( l0(i0, j) == l0(i1, j))
 								//	{cout<< "------------$$$$$$$------------" << endl;}
								//cout << "ttt:" << x0 << " " << y << " = " << l0(x1, y) << endl;
								//cout << "ppp:" << l0(x1, y) << endl;

							}
						}
					}
				}
				else if(neighbourMatrixOfClusters[ k0 ][ k1 ] && 
				   (! neighbourMatrixOfClusters[k0][k2]) && (! neighbourMatrixOfClusters[k1][k2]))
				{

					for (int i0 = 0; i0 < vClusters[ k0 ].size(); ++i0)
					{
						for (int i1 = 0; i1 < vClusters[ k1 ].size(); ++i1)
						{
							amount_of_constraints++;
							x0 = getNumberOfPointInVectorOfPoints(vClusters[ k0 ].getPoint(i0));
							x1 = getNumberOfPointInVectorOfPoints(vClusters[ k1 ].getPoint(i1));

							entries.push_back( ME(amount_of_constraints, l2(x0, k2, vClusters), 1) );
							entries.push_back( ME(amount_of_constraints, l1(x1, k2, vClusters), -exp(distanceMatrix[x0][x1]) ) );

							if( l2(i0, k2, vClusters) == l1(i1, k2, vClusters))
 							{cout<< "++++----++++++++---++++++" << endl;}
							
						}
					}
				}
				else
				{

					for (int j = 0; j < vClusters[ k2 ].size(); ++j)
					{
						amount_of_constraints++;

						y = getNumberOfPointInVectorOfPoints(vClusters[ k2 ].getPoint(j));
						entries.push_back( ME(amount_of_constraints, l4(y, k0, vClusters), 1) );
						entries.push_back( ME(amount_of_constraints, l3(y, k1, vClusters), -exp(distanceMatrixForClusters[k0][k1]) ) );
						
						if( l4(j, k0, vClusters) == l3(j, k1, vClusters))
 							{cout<< "/////////////*****////////////" << endl;}

					}
				}
			}
		}
	}
	cout << " CONSTRAINTS: " << amount_of_constraints << endl;
	// Return the amount of constraints == b_size == sense_size
	return amount_of_constraints;
}

int InitConstraints(std::vector<Cluster> vClusters)
{
		initTrivialConstraints(vClusters);
		initFirstConstraints(vClusters);
		initSecondConstraints(vClusters);
		initThirdConstraints(vClusters);
		initFourthConstraints(vClusters);
		int numberOfConstraints =  initPrivacyConstraints(vClusters);
		return numberOfConstraints;
}

void initParametersOfLP(int numberOfVariables, int numberOfConstraints)
{
	//cout << "Init Parameters of LP! " << endl;
	lp.c.set_size(numberOfVariables);	
	lp.b.set_size(numberOfConstraints);
	lp.sense.set_size(numberOfConstraints);
}

void fillVector_B_and_Sence(int numberOfConstraints)
{
	/* 
	Before run initParametersOfLP(int c_size, int b_size, int sense_size)!!!
	*/
	//cout << "Fill vector B and SENCE in LP! " << endl;
	for (int i = 0; i < numberOfConstraints; ++i)
	{
		if (i < NUMBER_OF_POINTS)
		{
			lp.b(i) = 1;
			lp.sense(i) = '=';

		}
		else
		{
			lp.b(i) = 0;
			lp.sense(i) = '<';
		}
	}
}

void initPriors()
{
	//cout << "Init Priors!" << endl;
    for(int i = 0; i< NUMBER_OF_POINTS; i++)
    {
    	Pi[i] = 1.0 / NUMBER_OF_POINTS;
    }
    cout << endl;
}

void displayPriors()
{
	//cout << "Priors:" << endl;
	for (int i = 0; i < NUMBER_OF_POINTS; ++i)
	{
		cout << Pi[i] << endl;
	}
	
}

void fillVector_C()
{
	//cout << "Fill vector C in LP! " << endl;

	for (int i = 0; i < NUMBER_OF_POINTS; ++i)
	{
		for (int j = 0; j < NUMBER_OF_POINTS; ++j)
		{
			lp.c(l0(i, j)) = (Pi[ i ] * distanceMatrix[ i ][ j ]);
		}
	}	
	/*
	for (int i = NUMBER_OF_POINTS_2; i < numberOfVariables; ++i)
	{
		lp.c(i) = 0;
	}
	*/
}

void fill_Objects_Of_LP(int numberOfConstraints)
{
	fillVector_C();
	fillVector_B_and_Sence(numberOfConstraints);
	lp.fill_A(entries);
}


vector<Cluster> InitBasic(double distanceThresholdForNeighbourClusters)
{
	/*
		For all functions F() you have void displayF()
	*/
	initSetOfLocations();
	initPriors();
	initDistanceMatrix();
	vector<Cluster> vClusters = naiveClustering();
	initDistanceMatrixForClusters(vClusters);
	initNeighbourMatrixOfClusters(distanceThresholdForNeighbourClusters);
	initAmountOfNeighbour();
	initT1(vClusters);
	initT2(vClusters);
	return vClusters;
}

double solve(string extremum, string method)
{
	cout << " Method: " << method << endl;
	cout << " ::::::::: " << extremum << " ::::::::" << endl; 
	if(extremum == "max")	
	{
		lp.maximize = true;
	}
	else
	{
		lp.maximize = false;
	}
	if(method == "simplex_dual")
	{
		lp.method = qif::lp::method_t::simplex_dual;
	}
	else if(method == "simplex_primal")
	{
		lp.method = qif::lp::method_t::simplex_primal;
	}
	else
	{
		lp.method = qif::lp::method_t::interior;
	}
	clock_t tStart = clock();
	bool solved = lp.solve();
	double executionTime = (double)(clock() - tStart)/CLOCKS_PER_SEC;
	cout <<" Solved: " << solved << endl;

    return executionTime;
}
void showInfo()
{
 	cout <<" A.n_rows: "<< lp.A.n_rows << endl;
    cout <<" b.n_elem: "<< lp.b.n_elem << endl;
    cout <<" sense.n_elem: "<< lp.sense.n_elem << endl;
    cout <<" A.n_cols: "<< lp.A.n_cols << endl;
    cout <<" c.n_elem: "<< lp.c.n_elem << endl;
    cout <<" Status: "<< lp.status << endl;
    cout <<" lp.x.n_elem: "<< lp.x.n_elem <<endl;
    cout << " Optimum: "<< lp.optimum() << endl;
	
}

/*
// In case that vector B and Sence will be calculated separately!
void fillVectorB(int b_size, int numberOfTrivialConstraints)
{
	//Before run initParametersOfLP(int c_size, int b_size, int sense_size)!!!
	for (int i = 0; i < b_size; ++i)
	{
		if(i < numberOfTrivialConstraints)
			lp.b(i) = 1;
		else
			lp.b(i) = 0;
	}
}
void fillSence(int sence_size, int numberOfTrivialConstraints)
{
	//Before run initParametersOfLP(int c_size, int b_size, int sense_size)!!!
	for (int i = 0; i < sence_size; ++i)
	{
		if(i < numberOfTrivialConstraints)
			lp.sense(i) = '=';
		else
			lp.sense(i) = '<';
	}	
	
}
*/
//brut0303 Brut0303 brut_0303 Brut_0303
//R,eE)g>M



