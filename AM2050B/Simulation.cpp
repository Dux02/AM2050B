#include "Simulation.h"
#include "VectorMath.h"
#include <algorithm>
#include <iostream>
#include <cmath>

using std::cerr;
using std::endl;


long double KeplersEq(long double E, long double e, long double M) {
	return E - e * std::sinl(E) - M;
}
long double KeplersEqDfn(long double E, long double e) {
	return 1 - e * std::cosl(E);
}

long double KeplerSolver(long double e, long double M, long double E_init = -1) {
	long double accuracy = 1e-12; //Beyond e-16 trig terms may actually cause problems
	int maxIterations = 100;
	long double E = E_init == -1 ? (e > 0.8 ? M_PI : M) : E_init;
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

PointMass::PointMass(long double mass, bool canMove) : m(mass), r(3), v(3), mobile(canMove), keplerian(false), kepinfo(nullptr) {}
PointMass::PointMass(const PointMass& obj) : m(obj.m), r(obj.r), v(obj.v), mobile(obj.mobile), keplerian(obj.keplerian), kepinfo(obj.kepinfo) {}
PointMass::PointMass(Vector x, Vector vel, long double mass, bool canMove) : m(mass), r(x), v(vel), mobile(canMove), keplerian(false), kepinfo(nullptr) {}

PointMass::PointMass(keplerinfo * ki, long double mass) : m(mass), keplerian(true), mobile(true), kepinfo(ki) {
	r = Vector(3); v = Vector(3);
	auto E = ki->E;
	auto nu = 2 * std::atanl(std::sqrtl((1 + ki->e) / (1 - ki->e)) * std::tanl(E / 2));
	auto rad = ki->semmajaxs * (1 - ki->e * std::cosl(E));
	//rad = ki->semmajaxs;
	r[0] = rad * std::cosl(nu); r[1] = rad * std::sinl(nu);
	r = ki->RotMatrix * r;
	auto mu = G * (ki->M0 + m);
	v[0] = -std::sqrtl(mu * ki->semmajaxs) * std::sinl(E) / rad;
	v[1] = std::sqrtl(mu * ki->semmajaxs * (1 - ki->e * ki->e)) * std::cosl(E) / rad;
	v = ki->RotMatrix * v;
}


//Note: dt can be negative here!
void PointMass::update(long double dt) {
	if (!keplerian) { return; }
	if (!kepinfo) { cerr << "Something went wrong! Expected Keplerian PointMass object but didn't find keplerinfo!" << endl; exit(1); }

	kepinfo->M += dt * 2.0 * M_PI / kepinfo->T;
	// Keeps M as a parameter between 0 and 2pi to avoid potential numerical issues if |M| grows too large.
	if (kepinfo->M > 2.0 * M_PI) {
		kepinfo->M -= 2.0 * M_PI; kepinfo->o++;
	} else if (kepinfo->M < 0.0) {
		kepinfo->M += 2.0 * M_PI; kepinfo->o--;
	}
	auto e = kepinfo->e;
	auto E = KeplerSolver(e, kepinfo->M, kepinfo->E); kepinfo->E = E;
	//std::cout << E << endl;
	auto nu = 2 * std::atanl(std::sqrtl((1 + e) / (1 - e)) * std::tanl(E / 2.0));
	auto radius = kepinfo->semmajaxs * (1.0 - e * std::cosl(E));
	//radius = kepinfo->semmajaxs;

	auto x =radius * std::cosl(nu); auto y = radius * std::sinl(nu); r[2] = 0.0;
	Matrix3x3* m = &kepinfo->RotMatrix;
	r[0] = m->m[0][0] * x + m->m[0][1] * y;
	r[1] = m->m[1][0] * x + m->m[1][1] * y;
	r[2] = m->m[2][0] * x + m->m[2][1] * y;
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


Simulation::Simulation(long double dt, long double T0) : dt(dt), objects(std::vector<PointMass>()), T(T0) { temp = Vector(3); }
Simulation::Simulation(const Simulation& sim) : dt(sim.dt), objects(sim.objects), T(sim.T) { temp = Vector(3); }
Simulation::Simulation(long double dt, std::vector<PointMass> objects, long double T0) : dt(dt), objects(objects), T(T0) { 
	temp = Vector(3); for (int i = 0; i < objects.size(); i++) {
		if (objects[i].isKepler()) { 
			auto thingy = new keplerinfo;
			//thingy->RotMatrix = Matrix3x3();
			*thingy = *objects[i].kepinfo;
			this->objects[i].kepinfo = thingy;
			//std::memcpy(this->objects[i].kepinfo, objects[i].kepinfo, sizeof(keplerinfo)); 
		}	// Ensures we don't copy memory addresses!
	}
}

Simulation::~Simulation() { }


void Simulation::update()
{
	for (auto& obj : objects) {
		if (!obj.isMobile() || obj.isKepler()) {  continue; }
		obj.v += calcAccel(obj) * dt;
	} 
	for (auto& obj : objects) {
		if (obj.isKepler()) {
			obj.update(dt); continue;
		}
		obj.r += obj.v * dt;
	}
	T += dt;
}



void Simulation::advancedUpdate() {
	return;
}

Vector Simulation::calcAccel(PointMass object) { 
	bool found = false; object.a.SetZero();
	for (auto& obj : objects) {
		if (!found && obj == object) { found = true; continue; }
		temp = obj.r; temp -= object.r; auto d = Norm(temp);
		object.a[0] += (G * obj.m /d/d/d) * temp[0];
		object.a[1] += (G * obj.m /d/d/d) * temp[1];
		object.a[2] += (G * obj.m /d/d/d) * temp[2];
		//*/
		//object.a += (G * obj.m * std::powl(Norm(obj.r - object.r), -3)) * (obj.r - object.r);
	}
	return object.a;
}

long double Simulation::calcTimeStep(PointMass object)
{
	bool found = false; long double min = dt / dtcoeff;
	for (auto& obj : objects) {
		if (!found && obj == object) { found = true; continue; }
		long double tim = Norm(obj.r - object.r) / Norm(obj.v - object.v);
		if (tim < min) { min = tim; }
	}
	return min*dtcoeff;
}

void Simulation::ShiftInitVelsByHalfStep() {
	for (auto& obj : objects) {
		if (obj.isKepler() || !obj.isMobile()) { continue; }
		obj.v -= (dt / 2) * calcAccel(obj);
	}
}
	
void Simulation::rewind()
{
	for (auto& obj : objects) {
		if (obj.isKepler()) { obj.update(-dt); continue; }
		obj.r -= dt * obj.v; //x_{i} = x_{i+1} - v_{i + 1/2} * dt
	}
	for (auto & obj : objects) {
		if (obj.isKepler() || !obj.isMobile()) { continue; }
		obj.v -= calcAccel(obj) * dt;
	}
	T -= dt;
}

long double Simulation::getTime()
{
	return T;
}

void Simulation::simpleSave(std::ostream &output) {
	auto N = objects.size();
	for (int i = 0; i < N; i++) {
		output << objects[i].r;
		if (i < N - 1) {
			output << ", ";
		}
		else {
			output << "\n";
		}
	}
}

bool Simulation::assertValidKepler()
{
	bool valid = true;
	for (auto& obj : objects) {
		if (!obj.isKepler()) { continue; }
		auto dM = obj.kepinfo->M + 2 * M_PI * obj.kepinfo->o - obj.kepinfo->iM;
		if (abs(dM - 2.0 * M_PI * T / obj.kepinfo->T) > 1e-8) {
			std::cout << "Mr. eccentricity: " << obj.kepinfo->e << endl;
			std::cout << "delta_M: " << dM << ", Expected: " << 2.0 * M_PI * T / obj.kepinfo->T << endl;
			obj.kepinfo->M = obj.kepinfo->iM + 2.0 * M_PI * T / obj.kepinfo->T;
			valid = false; // break;
		}
	}
	return valid;
}

void Simulation::forceKeplerUpdate()
{
	for (auto& obj : objects) {
		if (obj.isKepler()) { obj.update(0.0); }
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

keplerinfo convertData(rawkeplerdata rkd, long double T, long double M0) {
	auto x = keplerinfo();
	x.e = rkd.de;
	x.M = (rkd.dL - rkd.dOmega - rkd.domega) * M_PI / 180;
	x.M0 = M0; x.T = T;
	x.semmajaxs = rkd.da * AU;
	x.RotMatrix = ZRotation(-rkd.domega * M_PI / 180) * XRotation(-rkd.di * M_PI / 180) * ZRotation(-rkd.dOmega * M_PI / 180);
	x.E = KeplerSolver(x.e, x.M);
	x.iM = x.M;
	if (x.iM < 0.0) { x.M += 2 * M_PI; x.o = -1; }
	else { x.o = 0; }
	return x;
};