#include "Simulation.h"
#include "VectorMath.h"
#include <algorithm>
#include <iostream>
#include <cmath>
using std::cerr;
using std::endl;

PointMass::PointMass(double mass, bool canMove)
{
	m = mass;
	r = Vector3ld::Zero();
	v = Vector3ld::Zero();
	mobile = canMove;
}
PointMass::PointMass(const PointMass& obj)
{
	m = obj.m; r = obj.r; v = obj.v; mobile = obj.mobile;
}
PointMass::PointMass(Vector3ld x, Vector3ld vel, double mass, bool canMove)
{
	r = x; v = vel; 
	m = mass; mobile = canMove;
}


Simulation::Simulation(long double dt, long double T0) {
	objects = std::vector<PointMass>();
	this->dt = dt; T = T0;
}
Simulation::Simulation(const Simulation& sim) {
	objects = sim.objects;
	dt = sim.dt; T = sim.T;
}
Simulation::Simulation(long double dt, std::vector<PointMass> objects, long double T0) {
	this->dt = dt; this->objects = objects; T = T0;
}


void Simulation::update()
{
	std::vector<Vector3ld> accelerations;
	for (auto& obj : objects) {
		if (!obj.isMobile()) { accelerations.push_back(Vector3ld::Zero()); continue; }
		accelerations.push_back(calcAccel(obj));
	}
	for (int i = 0; i < accelerations.size(); i++) {
		objects[i].v += accelerations[i] * dt;
		objects[i].r += objects[i].v * dt;
	}
	T += dt;
}

Vector3ld Simulation::calcAccel(PointMass object) { 
	bool found = false; object.a.setZero(); auto acceleration = Vector3ld(0,0,0);
	for (auto& obj : objects) {
		if (!found && obj == object) { found = true; continue; }
		//Vector r_ji = obj.r - object.r; //j-> other object, i-> the object we calculate acceleration for	
		acceleration = (G * obj.m * std::powl ((obj.r - object.r).squaredNorm(), -1.5)) * (obj.r - object.r);
	}
	return acceleration;
}

void Simulation::ShiftInitVelsByHalfStep() {
	for (auto& obj : objects) {
		obj.v -= (dt / 2) * calcAccel(obj);
	}
}

void Simulation::rewind()
{
	std::vector<Vector3ld> accelerations;
	for (auto& obj : objects) {
		obj.r -= dt * obj.v; //x_{i} = x_{i+1} - v_{i + 1/2} * dt
		if (!obj.isMobile()) { accelerations.push_back(Vector3ld::Zero()); continue; }
		accelerations.push_back(calcAccel(obj));
	}
	for (int i = 0; i < accelerations.size(); i++) {
		objects[i].v -= accelerations[i] * dt;
	}
	T -= dt;
}

long double Simulation::getTime()
{
	return T;
}

void Simulation::simpleSave(std::ostream &output) {
	for (int i = 0; i < objects.size(); i++) {
		output << objects[i].r;
		if (i < objects.size() - 1) {
			output << ", ";
		}
		else {
			output << std::endl;
		}
	}
}

Ellipsoid::Ellipsoid(double mass, long double a, long double b, long double c, bool canMove) : PointMass(mass, canMove) {
	this->a = a; this->b = b; this->c = c; 
}

Ellipsoid::Ellipsoid(double mass, Vector ellipticparams, bool canMove) : PointMass(mass, canMove) {
	if (ellipticparams.size() != 3) {
		cerr << "Wrong sized vector" << endl; //TODO
		exit(1);
	}
	a = ellipticparams[0]; b = ellipticparams[1]; c = ellipticparams[2];
}

Ellipsoid::Ellipsoid(const Ellipsoid& obj) : PointMass(obj) {
	a = obj.a; b = obj.b; c = obj.c;
}

Ellipsoid::Ellipsoid(Vector3ld x, Vector3ld vel, Vector ellipparams, double mass, bool canMove) : PointMass(x,vel,mass, canMove) {
	if (ellipparams.size() != 3) {
		cerr << "Wrong sized vector" << endl; //TODO
		exit(1);
	}
	a = ellipparams[0]; b = ellipparams[1]; c = ellipparams[2];
}

bool Ellipsoid::isInside(Vector3ld pos)
{
	return Vector3ld((r - pos)(0) / a, (r - pos)(1) / b, (r - pos)(2) / c).norm() <= 1;
}
