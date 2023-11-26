#include <iostream>
#include <fstream>
#include <cmath>
#include "VectorMath.h"
#include "Simulation.h"

int main() {
	std::ofstream file;
	file.open("data.txt");
	const float R = 1000.0;
	const float M = 4000.0;
	const float vel = std::sqrt(G * M / 4 * R);
	PointMass planetA(Vector(0, R, 0), Vector(vel, 0, 0), M);
	PointMass planetB(Vector(0, -R, 0), Vector(-vel, 0, 0), M);
	std::vector<PointMass> planets = { planetA, planetB };
	Simulation sim(1, planets);

	auto N = 10000;
	for (int i = 0; i < N; i++) {
		for (int it = 0; it < sim.objects.size(); it++) {
			file << sim.objects[it].r;
			if (it < sim.objects.size() - 1) {
				file << ", ";
			}
			else {
				file << std::endl;
			}
		}
		sim.update();
		//std::cout << sim.objects[0].r << "," << sim.objects[1].r << std::endl;
	}
	file.close();
	return 0;
}