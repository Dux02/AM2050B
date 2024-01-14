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

const auto M0 = 1.989e30;	// Mass of the sun (kg)
const auto ME = 5.972e24;	// Mass of the earth (kg)
const long double AU = 149597870700;		//1 AU in meters

struct planetstats {
	long double mass;			//kg
	long double period;			//days
	long double distfromsun;	//m
	long double orbitalvel;		//m/s
	std::string name;
}; 

const keplerinfo PLANETS_KEP[9] = {
{"Mercury", 0.38709893, 0.20563069, 7.00487, 48.33167, 77.45645, 252.25084},
{"Venus", 0.72333199,    0.00677323,    3.39471    ,76.68069,    131.53298,    181.97973},
{"Earth",    1.00000011,    0.01671022,    0.00005,    -11.26064,    102.94719,    100.46435},
{"Mars",    1.52366231,    0.09341233,    1.85061,    49.57854,    336.04084,    355.45332},
{"Jupiter",    5.20336301,    0.04839266,    1.30530,    100.55615,    14.75385,    34.40438},
{"Saturn",    9.53707032,    0.05415060,    2.48446,    113.71504,    92.43194,    49.94432},
{"Uranus",    19.19126393,    0.04716771,    0.76986,    74.22988,    170.96424,    313.23218},
{"Neptune",    30.06896348,    0.00858587,    1.76917,    131.72169,    44.97135,    304.88003},
{"Pluto",    39.48168677,    0.24880766,    17.14175,    110.30347,    224.06676,    238.92881}
};
/*const planetstats_cy PLANETS[9] = {{"Mercury",    0.00000066,    0.00002527,    -23.51,    -446.30,    573.57,    538101628.29},{
"Venus",    0.00000092,    -0.00004938,    -2.86,    -996.89,    -108.80,    210664136.06},{
"Earth",    -0.00000005,    -0.00003804,    -46.94,    -18228.25,    1198.28,    129597740.63},{
"Mars",    -0.00007221,    0.00011902,    -25.47,    -1020.19,    1560.78,    68905103.78},{
"Jupiter",    0.00060737,    -0.00012880,    -4.15,    1217.17,    839.93,    10925078.35},{
"Saturn",    -0.00301530,    -0.00036762,    6.11 - 1591.05,    -1948.89,    4401052.95},{
"Uranus",    0.00152025,    -0.00019150,    -2.09,    -1681.40,    1312.56,    1542547.79},{
"Neptune",    -0.00125196,    0.0000251,    -3.64, -151.25,    -844.43,    786449.21},{
"Pluto",    -0.00076912,    0.00006465,    11.07,    -37.33,    -132.25,    522747.90}, };*/

const planetstats PLANETS[11] = {
{3.955e30, 0, 0, 0, "Sun"},
{0.33e24, 88.0, 57.9e9, 47400, "Mercury"},
{4.87e24, 224.7, 108.2e9, 35000, "Venus"},
{5.97e24, 365.2, 149.6e9, 29800, "Earth"},
{0.073e24, 27.3, 0.384e9, 1000, "Moon"},
{0.642e24, 687, 228e9, 24100, "Mars"},
{1898e24, 4331, 778.5e9, 13100, "Jupiter"},
{568e24, 10747, 1432e9, 9700, "Saturn"},
{86.8e24, 30589, 2867e9, 6800, "Uranus"},
{102e24, 59800, 4515e9, 5400, "Neptune"},
{0.0130e24, 90560, 5906.4e9, 4700, "Pluto"},
};
const int SECS_IN_DAY = 24 * 60 * 60;
/*
PointMass KeplerInitializer(keplerinfo data, planetstats planet, long double M_sun) {
	//Step 0. Set everything to radians
	auto L = data.dL * M_PI / 180.0; auto i = data.di * M_PI / 180.0;
	auto Omega = data.dOmega * M_PI / 180.0; auto omega = data.domega * M_PI / 180.0;
	//And get the other data types too:
	auto e = data.de; auto a = data.da * AU;

	//Step 1. From Mean Longitude (L) to true anomaly (nu)
	//1.1 Find mean anomaly (M)
	auto M = L - Omega - omega;
	//1.2 Find eccentric anomaly (E)
	auto E = KeplerSolver(e, M);
	//1.3 Hence find true anomaly
	auto nu = 2 * std::atanl(std::sqrtl((1 + e) / (1 - e)) * std::tanl(E / 2.0));
	
	//Step 2. From Semi-major axis (a) to Orbital angular momentum (h)
	auto mu = G * (planet.mass + M_sun);		//Standard gravitational parameter
	auto h = std::sqrtl( a * mu * (1 - e * e));

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

	cout << "FOR THE PLANET " << planet.name << "\n";
	//cout << "Calculated velocity: " << Norm(Velocity) << " v. Real vel: " << planet.orbitalvel << endl;
	cout << "MANLET ALERT? " << std::round(100 * (Norm(Velocity) / planet.orbitalvel - 1)) << "% \n";

	//Step 4. Rotate from orbital reference frame to Earth-orbital reference frame (standard)
	Matrix3x3 RotMatrix = XRotation(-Omega) * YRotation(-i) * ZRotation(-omega);
	return PointMass(RotMatrix * Position, RotMatrix * Velocity, planet.mass);
}
*/
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



/// <summary>
/// This checks if the integrator is indeed time symmetric by going forward and backward in time and comparing
/// the positions
/// </summary>
static int StableOrbitTest() { 
	auto file = openDataFile("sot ",".txt");
	
	const long double v = 10;
	const long double R = 1000;
	const long double M = R * v * v / G;
	const auto DT = 1; //Note DT = 10 is still stable.
	const auto N = 10000;

	auto Sun = std::make_shared<PointMass>(M, false);
	auto Planet = std::make_shared<PointMass>(Vector(R, 0, 0), Vector(0, v, 0), 1);
	
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
	std::cout << "Finished backward run. Initial position difference:" << Norm(Vector(R, 0, 0) - sim.objects[1]->r) << std::endl;
	sim.simpleSave(file);
	file.close();
	return 0;
}

/*
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
 
	PointMass Moon(Vector(R + Rem,0, 0), Vector(0, velcE + velcM, 0), M);
	PointMass Earth(Vector(R, 0, 0), Vector(0, velcE, 0), ME);
	PointMass Sun(Vector(3), Vector(3), M0, false);
	Simulation sim(DT, { &Sun, &Earth, &Moon});

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

	std::vector<PointMass *> planets;

	for (int i = 1; i <= N; i++) {
		planets.push_back(&PointMass(Vector(R * i, 0, 0), Vector(0, sqrt(G * M / R), 0), M * (N+1-i)));
	}
	planets.push_back(&PointMass(Vector(3), Vector(3), M * N * 10000,false));
	Simulation sim(DT, planets);
	sim.ShiftInitVelsByHalfStep();
	cout << "Finished setting up Efficiency test in " << stopwatch(timer) << " ms" << endl;

	for (int i = 0; i < T; i++) {
		//sim.simpleSave(file);
		sim.update();
		if (i % percent == (percent - 1)) {
			cout << "Looped through " << ceil(100 * i / T) << "% in " << stopwatch(timer) << " ms" << endl;
		}
	}

	cout << "Finished efficiency test in " << std::chrono::duration_cast<std::chrono::seconds>
		(prev - std::chrono::steady_clock::now()).count() << " seconds" << endl;
	file.close();
	return 0;
}

static int CollisionRewind() {
	return 0;
}

static int SolarSys() {
	auto file = openDataFile("ssys ", ".txt");
	const long double DT = 10;

	//Need to use calculated stable velocities to ensure proper orbit
	//const long double velcE = std::sqrt(G * M0 / R);
	//const long double velcM = std::sqrt(G * ME / Rem);

	PointMass Sun(Vector(3), Vector(3), PLANETS[0].mass, false);
	std::vector<PointMass*> solarsystem = { &Sun };
	for (int i = 1; i < 11; i++) {
		if (PLANETS[i].name == "Moon") {
			solarsystem.push_back(&PointMass(Vector(PLANETS[i - 1].distfromsun, PLANETS[i].distfromsun, 0),
				Vector(PLANETS[i].orbitalvel, PLANETS[i - 1].orbitalvel, 0), PLANETS[i].mass, true));
		}
		solarsystem.push_back(&PointMass(Vector(PLANETS[i].distfromsun, 0, 0), 
			Vector(0, PLANETS[i].orbitalvel, 0), PLANETS[i].mass, true));
	}
	cout << solarsystem.size() << endl;
	Simulation sim(DT, solarsystem);
	cout << "Starting Solar System simulation run" << endl;

	auto N = 20000;
	for (int i = 0; i < N*100; i++) {
		sim.simpleSave(file);
		sim.update();
		if (i % N == N - 1) {
			cout << "Finished " << ceil(i / N) << "% of simulation" << endl;
		}
	}
	file.close();
	cout << "Finished Solar System simulation run" << endl;
	return 0;
}
*/
static int KeplerSunEarth() {
	auto file = openDataFile("kse ", ".txt");
	const long double DT = 3600;
	auto Sun = std::make_shared<PointMass>(Vector(3), Vector(3), M0, false);
	std::vector<std::shared_ptr<PointMass>> solarsys = { Sun }; 
	solarsys.push_back(std::make_shared<KeplerObject>(PLANETS_KEP[2], PLANETS[3].mass, M0, PLANETS[3].period));
	Simulation sim(DT, solarsys);
	auto N = 20000;
	for (int i = 0; i < N; i++) {
		sim.simpleSave(file);
		sim.update();
	}
	file.close();
	cout << "Finished Sun-Earth simulation with Kepler initializer" << endl;
	return 0;
}
/*
static int KeplerSolarSys() {
	//auto file = openDataFile("kss ", ".txt");
	const long double DT = 3600;
	PointMass Sun(Vector(3), Vector(3), M0, false);
	std::vector<PointMass> solarsystem = { Sun };
	auto diff = 1;
	auto timer = std::chrono::steady_clock::now();

	for (int i = 1; i < 11; i++) {
		if (PLANETS[i].name == "Moon") {
			diff++; continue;
		}
		solarsystem.push_back(KeplerObject(PLANETS_KEP[i - diff], PLANETS[i].mass, M0, PLANETS[i].period));
	}
	cout << "Initialized Kepler system in " << stopwatch(timer) << " ms, for " << solarsystem.size() - 1 << " planets\n";
	Simulation sim(DT, solarsystem);

	auto N = 200;
	for (int i = 0; i < 100 * N; i++) {
		//sim.simpleSave(file);
		sim.update();
		if (i % N == N - 1) {
			cout << "Finished " << ceil(i / N) << "% of simulation" << endl;
		}
	}
	//file.close();
	cout << "Finished Solar System simulation run" << endl;
	return 0;
}
*/
int main() {
	
	//*/
	/*
	PointMass planetA(Vector(-R/2, 0, 0), ran2Dv(5), M * (1.1 - 0.2 * ran()));
	PointMass planetB(Vector(R/2, 0, 0), ran2Dv(5), M * (1.1 - 0.2 * ran()));
	PointMass planetC(Vector(0, R*std::sqrt(3)/2, 0), ran2Dv(5), 1000000*M * (1.1 - 0.2 * ran()));
	
	std::cout << planetA.m << ", " << planetB.m << ", " << planetC.m << std::endl;
	*/
	return StableOrbitTest();	
}