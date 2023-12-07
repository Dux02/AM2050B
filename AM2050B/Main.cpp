#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <chrono>
#include "date.h" //Standard library for chrono formatting, not-self-made	

#include "VectorMath.h"
#include "Simulation.h"
using std::cout; using std::endl;

const auto M0 = 3.955 * std::powl(10, 30); // Mass of the sun (kg)
const auto ME = 5.972 * std::powl(10, 24); // Mass of the earth (kg)

//Seed random generator
std::random_device rd;
std::mt19937 mt(rd());
std::uniform_real_distribution<double> unitdistr(0.0, 1.0);

double ran() {
	return unitdistr(mt);
}
Vector ran2Dv(double normlim, double angmin = 0.0, double angmax = 2.0*M_PI) {
	double norm = ran() * normlim;
	double angle = ran() * (angmax - angmin) + angmin;
	return Vector(norm * std::cos(angle), norm * std::sin(angle), 0.0);
}
std::ofstream openDataFile(std::string basename = "",std::string ending = ".txt") {
	auto now = std::chrono::system_clock::now();
	auto filename = "data/" + basename + date::format("%m-%d %H-%M", now) + ending;
	std::ofstream file(filename);
	//cout << filename << endl;
	return file;
}
static auto stopwatch(std::chrono::steady_clock::time_point &timer) {
	auto str = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - timer).count();
	timer = std::chrono::steady_clock::now();
	return str;
}

static int StableOrbitTest() {
	auto file = openDataFile("sot ",".txt");
	
	const long double v = 10;
	const long double R = 1000;
	const long double M = R * v * v / G;
	const auto DT = 1; //Note DT = 10 is still stable.
	const auto N = 10000;

	PointMass Sun(M, false);
	PointMass Planet(Vector3ld(R, 0.0, 0.0), Vector3ld(0.0, v, 0.0), 1);
	
	Simulation sim(DT, { Sun, Planet });
	sim.ShiftInitVelsByHalfStep();
	for (int i = 0; i < N; i++) {
		sim.simpleSave(file);
		sim.update();
	}
	std::cout << "Finished forward Stable Orbit Test simulation run" << std::endl;
	for (int i = 0; i < N; i++) {
		sim.simpleSave(file);
		sim.rewind();
	}
	std::cout << "Finished backward run. Initial position difference:" << (Vector3ld(R, 0, 0) - sim.objects[1].r).norm() << std::endl;
	sim.simpleSave(file);
	file.close();
	return 0;
}

static int SunEarthMoonSys() {
	auto file = openDataFile("sems ",".txt");
	const long double DT = 3600;
	const long double R = 149.60 * std::powl(10.0, 9);
	const long double Rem = 384400000.0;
	const long double M = 7.347 * std::powl(10, 22);
	const long double velE = 29780.0;
	const long double velM = 1023.055556;

	//Need to use calculated stable velocities to ensure proper orbit
	const long double velcE = std::sqrt(G * M0 / R);
	const long double velcM = std::sqrt(G * ME / Rem);

	cout << "Comparison of real and stable velocities: " << velE << " " << velcE << endl;
	cout << "And for moon: " << velM << " " << velcM << endl;
 
	PointMass Moon(Vector3ld(R + Rem,0, 0), Vector3ld(0, velcE + velcM, 0), M);
	PointMass Earth(Vector3ld(R, 0, 0), Vector3ld(0, velcE, 0), ME);
	PointMass Sun(Vector3ld(3), Vector3ld(3), M0, false);
	Simulation sim(DT, { Earth, Moon, Sun });

	auto N = 20000;
	for (int i = 0; i < N; i++) {
		sim.simpleSave(file);
		sim.update();
	}
	file.close();
	cout << "Finished Sun-Earth-Moon System simulation run" << endl;
	return 0;
}

static int EfficiencyTest() {
	auto file = openDataFile("efft", ".txt");
	auto N = 100;
	auto percent = 100;
	auto T = percent*100;
	const long double M = 1 * std::powl(10,25);
	const long double R = 1000000000; //1 mil km (10^9)
	auto DT = 1;
	auto timer = std::chrono::steady_clock::now();
	auto prev = std::chrono::steady_clock::now();

	std::vector<PointMass> planets;

	for (int i = 1; i <= N; i++) {
		planets.push_back(PointMass(Vector3ld(R * i, 0, 0), Vector3ld(0, sqrt(G * M / R), 0), M * (N+1-i)));
	}
	planets.push_back(PointMass(Vector3ld::Zero(), Vector3ld::Zero(), M * N * 10000, false));
	Simulation sim(DT, planets);
	sim.ShiftInitVelsByHalfStep();
	cout << "Finished setting up Efficiency test in " << stopwatch(timer) << " ms" << endl;

	for (int i = 0; i < T; i++) {
		sim.simpleSave(file);
		sim.update();
		if (i % percent == (percent - 1)) {
			cout << "Looped through " << ceil(100 * i / T) << "% in " << stopwatch(timer) << " ms" << endl;
		}
	}

	cout << "Finished effieciency test in " << std::chrono::duration_cast<std::chrono::seconds>
		(prev - std::chrono::steady_clock::now()).count() << " seconds" << endl;
	file.close();
	return 0;
}

static int CollisionRewind() {
	return 0;
}

int main() {
	
	//*/
	/*
	PointMass planetA(Vector(-R/2, 0, 0), ran2Dv(5), M * (1.1 - 0.2 * ran()));
	PointMass planetB(Vector(R/2, 0, 0), ran2Dv(5), M * (1.1 - 0.2 * ran()));
	PointMass planetC(Vector(0, R*std::sqrt(3)/2, 0), ran2Dv(5), 1000000*M * (1.1 - 0.2 * ran()));
	
	std::cout << planetA.m << ", " << planetB.m << ", " << planetC.m << std::endl;
	*/
	
	return EfficiencyTest();
}