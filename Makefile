
all:

	g++ main.cpp clustering.cpp -std=c++11 -lqif -larmadillo -lglpk -o location_privacy
	