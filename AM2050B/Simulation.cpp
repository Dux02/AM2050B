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
	r = Vector(3);
	v = Vector(3);
	mobile = canMove; keplerian = false;
}
PointMass::PointMass(const PointMass& obj)
{
	m = obj.m; r = obj.r; v = obj.v; mobile = obj.mobile; keplerian = false;
}
PointMass::PointMass(Vector x, Vector vel, double mass, bool canMove)
{
	r = x; v = vel; 
	m = mass; mobile = canMove; keplerian = false;
}

ThreadPool::ThreadPool() { Start(); }

void ThreadPool::Start() {
	njobs_pending = 0;
	const auto number_of_threads = std::thread::hardware_concurrency();
	for (int i = 0; i < number_of_threads; i++) {
		threads.emplace_back(std::thread(&ThreadPool::ThreadLoop, this));
	}
}
void ThreadPool::ThreadLoop() {
	while (true) {
		std::function<void()> job;
		{
			std::unique_lock<std::mutex> lock(queue_mutex);
			mutex_condition.wait(lock, [this] {
				return !jobs.empty() || should_terminate;
			});
			if (should_terminate) {
				return;
			}

			job = jobs.front();
			jobs.pop();
			//std::cout << "Threadin ova here" << std::endl;
		}
		job();
		if (--njobs_pending == 0) {
			std::unique_lock<std::mutex> lock(main_mutex);
			main_condition.notify_one();
		}
	}
}
void ThreadPool::QueueJob(const std::function<void()>& job) {
	{
		njobs_pending++;
		std::unique_lock<std::mutex> lock(queue_mutex);
		jobs.push(job);
	}
	mutex_condition.notify_one();
}
bool ThreadPool::busy() {
	bool poolbusy;
	{
		std::unique_lock<std::mutex> lock(queue_mutex);
		poolbusy = !jobs.empty();
	}
	return poolbusy;
}
void ThreadPool::Stop() {
	{
		std::unique_lock<std::mutex> lock(queue_mutex);
		should_terminate = true;
	}
	mutex_condition.notify_all();
	for (std::thread& active_thread : threads) {
		active_thread.join();
	}
	threads.clear();
}
	
void ThreadPool::waitUntilCompleted() {
	std::unique_lock<std::mutex> lock(main_mutex);
	main_condition.wait(lock, [this] {return njobs_pending == 0; });
}


Simulation::Simulation(long double dt, long double T0) : dt(dt), objects(std::vector<std::shared_ptr<PointMass>>()), T(T0) {}
Simulation::Simulation(const Simulation& sim) : dt(sim.dt), objects(sim.objects), T(sim.T) {}
Simulation::Simulation(long double dt, std::vector<std::shared_ptr<PointMass>> objects, long double T0) : dt(dt), objects(objects), T(T0) {}

Simulation::~Simulation() { }


void Simulation::update()
{
	std::vector<Vector> accelerations;
	for (auto& obj : objects) {
		if (!obj->isMobile() || obj->isKepler()) { accelerations.push_back(Vector(3)); continue; }
		accelerations.push_back(calcAccel(*obj));
	} 
	for (int i = 0; i < accelerations.size(); i++) {
		if (objects[i]->isKepler()) {
			std::cout << "Hello I am Mr. Kepler" << std::endl;
			objects[i]->update(dt);
		}
		objects[i]->v += accelerations[i] * dt;
		objects[i]->r += objects[i]->v * dt;
	}
	T += dt;
}



void Simulation::advancedUpdate() {
	return;
}

Vector Simulation::calcAccel(PointMass object) { 
	bool found = false; object.a.SetZero();
	for (auto& obj : objects) {
		if (!found && *obj == object) { found = true; continue; }
		object.a += (G * obj->m * std::powl (Norm(obj->r - object.r), -3.0)) * (obj->r - object.r);
	}
	return object.a;
}

long double Simulation::calcTimeStep(PointMass object)
{
	bool found = false; long double min = dt / dtcoeff;
	for (auto& obj : objects) {
		if (!found && *obj == object) { found = true; continue; }
		long double tim = Norm(obj->r - object.r) / Norm(obj->v - object.v);
		if (tim < min) { min = tim; }
	}
	return min*dtcoeff;
}

void Simulation::ShiftInitVelsByHalfStep() {
	for (auto& obj : objects) {
		obj->v -= (dt / 2) * calcAccel(*obj);
	}
}

void Simulation::rewind()
{
	std::vector<Vector> accelerations;
	for (auto& obj : objects) {
		obj->r -= dt * obj->v; //x_{i} = x_{i+1} - v_{i + 1/2} * dt
		if (!obj->isMobile()) { accelerations.push_back(Vector(3)); continue; }
		accelerations.push_back(calcAccel(*obj));
	}
	for (int i = 0; i < accelerations.size(); i++) {
		objects[i]->v -= accelerations[i] * dt;
	}
	T -= dt;
}

long double Simulation::getTime()
{
	return T;
}

void Simulation::simpleSave(std::ostream &output) {
	for (int i = 0; i < objects.size(); i++) {
		output << objects[i]->r;
		if (i < objects.size() - 1) {
			output << ", ";
		}
		else {
			output << "\n";
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

Ellipsoid::Ellipsoid(Vector x, Vector vel, Vector ellipparams, double mass, bool canMove) : PointMass(x,vel,mass, canMove) {
	if (ellipparams.size() != 3) {
		cerr << "Wrong sized vector" << endl; //TODO
		exit(1);
	}
	a = ellipparams[0]; b = ellipparams[1]; c = ellipparams[2];
}

bool Ellipsoid::isInside(Vector pos)
{
	Vector lazy = r - pos; lazy[0] /= a; lazy[1] /= b; lazy[2] /= c;
	return Norm(lazy) <= 1;
}

long double KeplersEq(long double E, long double e, long double M) {
	return E - e * std::sinl(E) - M;
}
long double KeplersEqDfn(long double E, long double e) {
	return 1 - e * std::cosl(E);
}

long double KeplerSolver(long double e, long double M) {
	long double accuracy = 1e-8; //Beyond e-16 trig terms may actually cause problems
	int maxIterations = 100;
	long double E = e > 0.8 ? M_PI : M;
	for (int i = 1; i < maxIterations; i++) {
		//Newton Raphson: E_(N+1) = E_N - ( f(E_N) / f'(E_N) )
		long double nextVal = E - (KeplersEq(E, e, M) / KeplersEqDfn(E, e));
		long double difference = std::abs(E - nextVal);
		E = nextVal;
		if (difference < accuracy) {
			//cout << "Completed KeplerSolver in " << i << " iterations\n";
			break;
		}
	}
	return E;
}


KeplerObject::KeplerObject(keplerinfo data, double m, long double M0, long double T) : PointMass(m) {
	keplerian = true;
	this->m = m; this->M0 = M0; this->mobile = true; this->T = T;
	//Step 0. Set everything to radians
	auto L = data.dL * M_PI / 180.0; auto i = data.di * M_PI / 180.0;
	auto Omega = data.dOmega * M_PI / 180.0; auto omega = data.domega * M_PI / 180.0;
	//And get the other data types too:
	auto e = data.de; auto a = data.da;

	//Step 1. From Mean Longitude (L) to true anomaly (nu)
	//1.1 Find mean anomaly (M)
	auto M = L - Omega - omega;
	//1.2 Find eccentric anomaly (E)
	auto E = KeplerSolver(e, M);
	//1.3 Hence find true anomaly
	auto nu = 2 * std::atanl(std::sqrtl((1 + e) / (1 - e)) * std::tanl(E / 2.0));

	//Step 2. From Semi-major axis (a) to Orbital angular momentum (h)
	auto mu = G * (m + M0);		//Standard gravitational parameter
	auto h = std::sqrtl(a * mu * (1 - e * e));

	//Step 3. From these parameters determine initial positions and velocities:
	auto r = a * (1 - e * std::cosl(E));

	Vector Position(3);		//In Orbital reference frame
	//Position[0] = (h * h * std::cosl(nu)) / (mu * (1 + e * std::cosl(nu)));
	//Position[1] = Position[0] * std::tanl(nu);
	Position[0] = r * std::cosl(nu); Position[1] = r * std::sinl(nu);
	Vector Velocity(3);
	//Velocity[0] = -mu * std::sinl(nu) / h;
	//Velocity[1] = mu * (e + std::cosl(mu)) / h;
	Velocity[0] = -std::sqrtl(mu * a) * std::sinl(E) / r;
	Velocity[1] = std::sqrtl(mu * a * (1 - e * e)) * std::cosl(E) / r;


	//Step 4. Rotate from orbital reference frame to Earth-orbital reference frame (standard)
	Matrix3x3 RotMatrix = XRotation(-Omega) * YRotation(-i) * ZRotation(-omega);
	this->RotMatrix = RotMatrix;
	this->r = RotMatrix * Position; this->v = RotMatrix * Velocity;
	this->e = e; this->M = M;
	this->a = Vector(3); this->semmajaxs = a;
};

void KeplerObject::update(long double dt)
{
	M += dt * M_2_PI / T;
	auto E = KeplerSolver(e, M);
	auto nu = 2 * std::atanl(std::sqrtl((1 + e) / (1 - e)) * std::tanl(E / 2.0));
	auto radius = semmajaxs * (1 - e * std::cosl(E));
	
	r[0] = radius * std::cosl(nu); r[1] = radius * std::sinl(nu);
	r = RotMatrix * r;
}
