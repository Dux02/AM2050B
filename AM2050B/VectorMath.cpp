#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <iomanip>

#include "VectorMath.h"

using std::cerr;
using std::endl;


// Constructors for Vector class
Vector::Vector() : n(0) {
	v = 0;
}

Vector::Vector(int size) : n(size) {
	v = new double[size];
	for (int i = 0; i < size; i++) v[i] = 0.0;
}

Vector::Vector(const Vector& V) : n(V.n) {
	v = new double[V.n];
	for (int i = 0; i < V.n; i++) v[i] = V.v[i];
}

Vector::Vector(double x, double y, double z) : n(3) {
	v = new double[3];
	v[0] = x; v[1] = y; v[2] = z;
}

// Destructor for Vector class
Vector::~Vector() {
	delete[] v;	//Free up unused memory
}

Vector& Vector::operator=(const double value)
{
	for (int i = 0; i < n; i++) v[i] = value;
	return (*this);
}
Vector& Vector::operator=(const Vector& V)
{
	if (this == &V) return (*this); //Don't bother copying
	//Copy vector if this vector is empty
	if (n == 0) {
		n = V.n;
		v = new double[V.n];
	}
	if (n != V.n) {
		cerr << "ERROR: Vectors not the same size in Vector operator=(Vector), left: " << n << "right: " << V.n << endl;
		exit(1);
	}
	for (int i = 0; i < n; i++) v[i] = V.v[i];
	return (*this);
}

void Vector::operator+=(const Vector& V)
{
	if (n != V.n) {
		cerr << "ERROR: Vectors not the same size in Vector operator+=(Vector), left: " << n << "right: " << V.n << endl;
		exit(1);
	}
	for (int i = 0; i < n; i++) v[i] += V.v[i];
}
void Vector::operator-=(const Vector& V)
{
	if (n != V.n) {
		cerr << "ERROR: Vectors not the same size in Vector operator+=(Vector), left: " << n << "right: " << V.n << endl;
		exit(1);
	}
	for (int i = 0; i < n; i++) v[i] -= V.v[i];
}

bool operator==(const Vector& left, const Vector& right)
{
	if (left.n != right.n) return false;
	for (int i = 0; i < left.n; i++) if (left.v[i] != right.v[i]) return false;
	return true;
}

// Vector addition and subtraction
Vector operator+(const Vector& left, const Vector& right)
{
	if (left.n != right.n) {
		cerr << "ERROR: Vectors not the same size in +(Vector,Vector), left: " << left.n << "right: " << right.n << endl;
		exit(1);
	}
	Vector Aux(left.n);
	for (int i = 0; i < left.n; i++) Aux.v[i] = left.v[i] + right.v[i];
	return Aux;
}
Vector operator-(const Vector& left, const Vector& right)
{
	if (left.n != right.n) {
		cerr << "ERROR: Vectors not the same size in -(Vector,Vector), left: " << left.n << "right: " << right.n << endl;
		exit(1);
	}
	Vector Aux(left.n);
	for (int i = 0; i < left.n; i++) Aux.v[i] = left.v[i] - right.v[i];
	return Aux;
}

//Vector operations
double Dot(const Vector& left, const Vector& right)
{
	if (left.n != right.n) {
		cerr << "ERROR: Vectors not the same size in Dot(Vector,Vector), left: " << left.n << "right: " << right.n << endl;
		exit(1);
	}
	double Sum = 0.0;
	for (int i = 0; i < left.n; i++) Sum += left.v[i] * right.v[i];
	return Sum;
}
double Norm(const Vector& V)
{
	return std::sqrt(Dot(V, V));
}
Vector Cross(const Vector& left, const Vector& right)
{
	if (left.n != 3 || right.n != 3) {
		cerr << "ERROR: Vectors not 3 dimensional in Cross(Vector,Vector), left: " << left.n << "right: " << right.n << endl;
		exit(1);
	}
	Vector Result(3);
	Result.v[0] = left.v[1] * right.v[2] - left.v[2] * right.v[1];
	Result.v[1] = left.v[2] * right.v[0] - left.v[0] * right.v[2];
	Result.v[2] = left.v[0] * right.v[1] - left.v[1] * right.v[0];
	return Result;
}

// Scalar multiplication and division
Vector operator*(double value, const Vector& V)
{
	Vector Aux(V.n);
	for (int i = 0; i < V.n; i++) Aux.v[i] = value * V.v[i];
	return Aux;
}
Vector operator*(const Vector& V, double value)
{
	return value * V;
}
Vector operator/(const Vector& V, double value)
{
	Vector Aux(V.n);
	for (int i = 0; i < V.n; i++) Aux.v[i] = V.v[i] / value;
	return Aux;
}



// Create a Vector by using polar coordinates
Vector VecPolar(double r, double theta, double phi) {
	return Vector(r * cos(theta) * cos(phi), r * cos(theta) * sin(phi), r * sin(theta));
}

void CalcPolarCoords(const Vector& Vec, double& r, double& theta, double& phi)
{
	const double rhoSqr = Vec[0] * Vec[0] + Vec[1] * Vec[1];
	r = std::sqrt(rhoSqr + Vec[2] * Vec[2]);

	if (Vec[0] == 0.0 && Vec[1] == 0.0) {
		phi = 0.0;
	} else {
		phi = std::atan2(Vec[1], Vec[0]);
	}
	if (phi < 0.0) phi += 2.0 * M_PI;

	const double rho = std::sqrt(rhoSqr);
	if (Vec[2] == 0.0 && rho == 0.0)
		theta = 0.0;
	else
		theta = std::atan2(Vec[2], rho);
}

std::ostream& operator<<(std::ostream& os, const Vector& Vec)
{
	//int w = os.width();
	for (int i = 0; i < Vec.n; i++)
		os << ' ' << Vec[i];
	return os;
}
