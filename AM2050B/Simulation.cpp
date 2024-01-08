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
	mobile = canMove;
}
PointMass::PointMass(const PointMass& obj)
{
	m = obj.m; r = obj.r; v = obj.v; mobile = obj.mobile;
}
PointMass::PointMass(Vector x, Vector vel, double mass, bool canMove)
{
	r = x; v = vel; 
	m = mass; mobile = canMove;
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


Simulation::Simulation(long double dt, long double T0) : dt(dt), objects(std::vector<PointMass>()), T(T0) {}
Simulation::Simulation(const Simulation& sim) : dt(sim.dt), objects(sim.objects), T(sim.T) {}
Simulation::Simulation(long double dt, std::vector<PointMass> objects, long double T0) : dt(dt), objects(objects), T(T0) {}

Simulation::~Simulation() { workers.Stop(); }

/*
void Simulation::init() {
	numThreadsFinished = 0;
	for (int i = 0; i < numThreads; ++i) {
		threads.emplace_back([this, i]() {
			while (true) {
				// Wait for new task
				//std::unique_lock<std::mutex> lock(accelMutex);
				accelReady.wait(lock, [this]() { return numThreadsFinished < numThreads || terminate; });

				if (terminate) {
					break;  // Exit the thread if termination is signaled
				}

				// Calculate accelerations for a subset of objects
				for (size_t j = i; j < this->objects.size(); j += numThreads) {
					Vector accel = calcAccel(this->objects[j]);

					// Make sure to update acceleration atomically to avoid race conditions
					//std::lock_guard<std::mutex> accelLock(accelMutex);
					accelerations[j] = accel;
				}

				// Notify that this thread has finished
				++numThreadsFinished;
				//lock.unlock();
				accelReady.notify_one();
			}
		});
	}

	accelerations.resize(objects.size());
}
*/


void Simulation::update()
{
	std::vector<Vector> accelerations;
	/* Multi threading LMAO - This is where we use the treads
	std::mutex reality_land;
	for (auto& obj : objects) {
		if (!obj.isMobile()) { std::unique_lock<std::mutex> fake_land(reality_land); 
		accelerations.push_back(Vector(3)); continue; }
		workers.QueueJob([this, obj, &accelerations, &reality_land] {
			std::unique_lock<std::mutex> fake_land(reality_land);
			accelerations.push_back(calcAccel(obj)); //Make sure there are no data races to acceleration
			//std::cout << accelerations.size() << std::endl;
		});
		
	}
	workers.waitUntilCompleted(); //And then wait until all the workers are completed	
	//std::cout << "We waited until completed" << std::endl;
	//As long as threads complete before here we are good to go! */
	///* Single threaded
	for (auto& obj : objects) {
		if (!obj.isMobile()) { accelerations.push_back(Vector(3)); continue; }
		accelerations.push_back(calcAccel(obj));
	} //*/
	for (int i = 0; i < accelerations.size(); i++) {
		objects[i].v += accelerations[i] * dt;
		objects[i].r += objects[i].v * dt;
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
		object.a += (G * obj.m * std::powl (Norm(obj.r - object.r), -3.0)) * (obj.r - object.r);
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
		obj.v -= (dt / 2) * calcAccel(obj);
	}
}

void Simulation::rewind()
{
	std::vector<Vector> accelerations;
	for (auto& obj : objects) {
		obj.r -= dt * obj.v; //x_{i} = x_{i+1} - v_{i + 1/2} * dt
		if (!obj.isMobile()) { accelerations.push_back(Vector(3)); continue; }
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
