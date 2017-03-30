#ifndef ADD_PH
#define ADD_PH

#include <iostream>
#include "CONSTANTS.h"

using namespace std;
class MyPoint 
{ 
private:
	int x;
	int y;

public:
	void setX(int x){ this->x = x; }
	int getX(){ return this->x; }

	void setY(int y){ this->y = y; }
	int getY(){ return this->y; }

	void setXY(int x, int y){ this->x = x; this->y = y; }

	void display(){ cout<<"x = "<< this->x <<" y = "<< this->y<<endl; }

};


#endif
