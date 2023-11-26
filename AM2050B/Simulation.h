#pragma once
#include "VectorMath.h"
#include <vector>
#define G 1
//0.000000000066743

class PointMass {
public:
	PointMass(double mass);
	PointMass(const PointMass& obj);
	PointMass(Vector x, Vector vel, double mass);

	Vector r;	
	Vector v;	//Note for leapfrog model, this is v_{i-1/2} 
	double m;

	friend bool operator==(const PointMass& left, const PointMass& right) { return (left.m == right.m && left.r == right.r && left.v == right.v); }
};

class Simulation {
public:
		Simulation (double dt);
		Simulation (const Simulation& sim);
		Simulation (double dt, std::vector<PointMass> objects);

		std::vector<PointMass> objects;
		void update();
private:
	Vector calcAccel(PointMass object);
	double dt;
};
