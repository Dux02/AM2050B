#include "Simulation.h"
#include "VectorMath.h"
#include <algorithm>
#include <iostream>
#include <cmath>
using std::cerr;
using std::endl;

PointMass::PointMass(double mass)
{
	m = mass;
	r = Vector(3);
	v = Vector(3);
}
PointMass::PointMass(const PointMass& obj)
{
	m = obj.m; r = obj.r; v = obj.v;
}
PointMass::PointMass(Vector x, Vector vel, double mass)
{
	r = x; v = vel; m = mass;
}


Simulation::Simulation(double dt) {
	objects = std::vector<PointMass>();
	this->dt = dt;
}
Simulation::Simulation(const Simulation& sim) {
	objects = sim.objects;
	dt = sim.dt;
}
Simulation::Simulation(double dt, std::vector<PointMass> objects) {
	this->dt = dt; this->objects = objects;
}


void Simulation::update()
{
	std::vector<Vector> accelerations;
	for (auto& obj : objects) {
		accelerations.push_back(calcAccel(obj));
	}
	for (int i = 0; i < accelerations.size(); i++) {
		objects[i].v += accelerations[i] * dt;	
		objects[i].r += objects[i].v * dt;
	}
}

Vector Simulation::calcAccel(PointMass object) {
	Vector accel(3); bool found = false;
	for (auto& obj : objects) {
		if (!found && obj == object) { found = true; continue; }
		Vector r_ij = obj.r - object.r;
		accel += (G * obj.m * std::pow(Norm(r_ij), -3)) * r_ij;
	}
	return accel;
}


