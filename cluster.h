#ifndef ADD_CLH
#define ADD_CLH


#include <vector>
#include "point.h"
//#include <unordered_map>

class Cluster 
{
private:
	std::vector<MyPoint> vPoints;
//	std::unordered_map <MyPoint, int> mPointsToIndex;

public:
	void push_back(MyPoint p) 
	{ 
		vPoints.push_back(p);
//		mPointsToIndex[ p ] = vPoints.size() - 1;

	}

	void pop_back() 
	{ 
//		MyPoint p =  vPoints[ vPoints.size() - 1 ];
//		mPointsToIndex.erase(p);
		vPoints.pop_back(); 
	}

	int size() { return vPoints.size(); }

	MyPoint getPoint(int index) 
	{ 
		return vPoints.at(index); 
	}

	void setPoint(MyPoint p, int index) 
	{ 
		vPoints.at(index) = p; 
		//mPointsToIndex[ p ] = index;
	}
/*
	int getIndex(MyPoint p)
	{
	 	return mPointsToIndex[ p ]; 
	}
*/
	void displayAllPoints() 
	{
	 	for (int i = 0; i < vPoints.size(); ++i) 
	 		{ 
	 			vPoints[i].display(); 
	 		} 
	 }

};

#endif