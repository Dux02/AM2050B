#pragma once
#define _USE_MATH_DEFINES
#include "VectorMath.h"
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <queue>
#include <string>
#include <cmath>
#include <zlib.h>

#define G 0.000000000066743 //Gravitational constant, 10^-11
constexpr long double dtcoeff = 12 * 3600;
const long double AU = 149597870700;		//1 AU in meters
//0.000000000066743

struct rawkeplerdata {
	std::string name;
	long double da;		//semi major axis in AU
	long double de;
	long double di;		//in degrees
	long double dOmega;	//in degrees
	long double domega;	//in degrees
	long double dL;		//in degrees
};
struct keplerinfo {
	long double e;					// Eccentricity
	long double semmajaxs;			// Semi-major axis
	Matrix3x3 RotMatrix;
	long double T;					// Orbital period
	//long double L;				// Mean Longitude
	long double M;					// Mean Anomaly
	long double E;
	long double M0;					// Mass of center of orbit (e.g. Sun)
	long double iM;					// Initial Mean Anomaly
	int o;							// Number of complete orbits around the sun
};

keplerinfo convertData(rawkeplerdata rkd, long double T, long double M0);

class PointMass {
public:
	PointMass(long double mass, bool canMove = true);
	PointMass(const PointMass& obj);
	PointMass(Vector x, Vector vel, long double mass, bool canMove = true);
	PointMass(keplerinfo* ki, long double mass);

	Vector r;	
	Vector v;	//Note for leapfrog model, this is v_{i-1/2} 
	Vector a = Vector(3);	//Simple workaround to avoid needless initalizations and destructions
	long double m;
	bool isKepler() const { return keplerian; }
	//long double dt_i;		//Time spent since start of an era (!!)

	friend bool operator==(const PointMass& left, const PointMass& right) 
		{ return (left.mobile == right.mobile && left.m == right.m && left.r == right.r && left.v == right.v); }
	bool isMobile() const { return mobile; }
	void update(long double dt);

	keplerinfo* kepinfo;
protected:
	bool mobile;
	bool keplerian;
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


class ThreadPool {
public:
	ThreadPool();
	void Start();
	void QueueJob(const std::function<void()>& job);
	void Stop();
	void waitUntilCompleted();
	bool busy();

private:
	void ThreadLoop();

	bool should_terminate = false;           // Tells threads to stop looking for jobs
	std::mutex queue_mutex;                  // Prevents data races to the job queue
	std::condition_variable mutex_condition; // Allows threads to wait on new jobs or termination 
	std::vector<std::thread> threads;
	std::queue<std::function<void()>> jobs;
	std::mutex main_mutex;
	std::condition_variable main_condition;
	std::atomic<int> njobs_pending;
};

class Simulation {
public:
	Simulation (long double dt, long double T0 = 0.0);
	Simulation (const Simulation& sim);
	Simulation (long double dt, std::vector<PointMass> objects, long double T0 = 0.0);
	~Simulation();

	std::vector<PointMass> objects;
	void update();
	void advancedUpdate(); //Uses the Block step era method
	void rewind();
	Vector calcAccel(PointMass object);
	long double calcTimeStep(PointMass object);
	void ShiftInitVelsByHalfStep();
	void simpleSave(std::ostream &output);
	bool assertValidKepler();
	void forceKeplerUpdate();

	long double getTime();
private:
	//void init();
	long double dt;
	long double T;
	Vector temp;
	
};
