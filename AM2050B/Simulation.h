#pragma once
#include "VectorMath.h"
#include <Eigen/Dense>
#include <vector>
#define G 0.000000000066743 //Gravitational constant, 10^-11
typedef Eigen::Vector3<long double> Vector3ld;

//0.000000000066743

class PointMass {
public:
	PointMass(double mass, bool canMove = true);
	PointMass(const PointMass& obj);
	PointMass(Vector3ld x, Vector3ld vel, double mass, bool canMove = true);

	Vector3ld r;
	Vector3ld v;	//Note for leapfrog model, this is v_{i-1/2} 
	Vector3ld a = Vector3ld::Zero();	//Simple workaround to avoid needless initalizations and destructions
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
	Ellipsoid(Vector3ld x, Vector3ld vel, Vector ellipparams, double mass, bool canMove = true);

	bool isInside(Vector3ld pos);

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
	Vector3ld calcAccel(PointMass object);
	void ShiftInitVelsByHalfStep();
	void simpleSave(std::ostream &output);

	long double getTime();
private:
	long double dt;
	long double T;
};
