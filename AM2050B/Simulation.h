#pragma once
#include "VectorMath.h"
#include <vector>
#define G 0.000000000066743 //Gravitational constant, 10^-11
//0.000000000066743

class PointMass {
public:
	PointMass(double mass, bool canMove = true);
	PointMass(const PointMass& obj);
	PointMass(Vector x, Vector vel, double mass, bool canMove = true);

	Vector r;	
	Vector v;	//Note for leapfrog model, this is v_{i-1/2} 
	Vector a = Vector(3);	//Simple workaround to avoid needless initalizations and destructions
	double m;

	friend bool operator==(const PointMass& left, const PointMass& right) 
		{ return (left.mobile == right.mobile && left.m == right.m && left.r == right.r && left.v == right.v); }
	bool isMobile() const { return mobile; }
protected:
	bool mobile;
};

class Ellipsoid : PointMass {
public:
	Ellipsoid(double mass, long double a, long double b, long double c, bool canMove = true);
	Ellipsoid(double mass, Vector ellipticparams, bool canMove = true);
	Ellipsoid(const Ellipsoid& obj);
	Ellipsoid(Vector x, Vector vel, Vector ellipparams, double mass, bool canMove = true);

	bool isInside(Vector pos);

	//friend bool isInside(Ellipsoid obj, Vector r);
	//friend bool isInside(Ellipsoid obj, PointMass otherobj);

	long double a; long double b; long double c;
};

class Simulation {
public:
	Simulation (long double dt, long double T0 = 0.0);
	Simulation (const Simulation& sim);
	Simulation (long double dt, std::vector<PointMass> objects, long double T0 = 0.0);

	std::vector<PointMass> objects;
	void update();
	void rewind();
	Vector calcAccel(PointMass object);
	void ShiftInitVelsByHalfStep();
	void simpleSave(std::ostream &output);

	long double getTime();
private:
	long double dt;
	long double T;
};
