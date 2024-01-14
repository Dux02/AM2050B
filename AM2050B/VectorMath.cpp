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
	v = new long double[size];
	for (int i = 0; i < size; i++) v[i] = 0.0;
}
Vector::Vector(const Vector& V) : n(V.n) {
	v = new long double[V.n];
	for (int i = 0; i < V.n; i++) v[i] = V.v[i];
}
Vector::Vector(long double x, long double y, long double z) : n(3) {
	v = new long double[3];
	v[0] = x; v[1] = y; v[2] = z;
}

// Destructor for Vector class
Vector::~Vector() {
	delete[] v;	//Free up unused memory
}

Vector& Vector::operator=(const long double value)
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
		v = new long double[V.n];
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
long double Dot(const Vector& left, const Vector& right)
{
	if (left.n != right.n) {
		cerr << "ERROR: Vectors not the same size in Dot(Vector,Vector), left: " << left.n << "right: " << right.n << endl;
		exit(1);
	}
	long double Sum = 0.0;
	for (int i = 0; i < left.n; i++) Sum += left.v[i] * right.v[i];
	return Sum;
}
long double Norm(const Vector& V)
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
Vector operator*(long double value, const Vector& V)
{
	Vector Aux(V.n);
	for (int i = 0; i < V.n; i++) Aux.v[i] = value * V.v[i];
	return Aux;
}
Vector operator*(const Vector& V, long double value)
{
	return value * V;
}
Vector operator/(const Vector& V, long double value)
{
	Vector Aux(V.n);
	for (int i = 0; i < V.n; i++) Aux.v[i] = V.v[i] / value;
	return Aux;
}

// To make print statements for vectors
std::ostream& operator<<(std::ostream& os, const Vector& Vec)
{
	//int w = os.width();
	for (int i = 0; i < Vec.n; i++)
		if (i == 0) {
			os << " " << Vec[i];
		}
		else {
			os << ", " << Vec[i];
		}
	return os;
}

// Create a Vector by using polar coordinates
Vector VecPolar(long double r, long double theta, long double phi) {
	return Vector(r * cos(theta) * cos(phi), r * cos(theta) * sin(phi), r * sin(theta));
}

// Find Polar coordinates of a vector
void CalcPolarCoords(const Vector& Vec, long double& r, long double& theta, long double& phi)
{
	const long double rhoSqr = Vec[0] * Vec[0] + Vec[1] * Vec[1];
	r = std::sqrt(rhoSqr + Vec[2] * Vec[2]);

	if (Vec[0] == 0.0 && Vec[1] == 0.0) {
		phi = 0.0;
	} else {
		phi = std::atan2(Vec[1], Vec[0]);
	}
	if (phi < 0.0) phi += 2.0 * M_PI;

	const long double rho = std::sqrt(rhoSqr);
	if (Vec[2] == 0.0 && rho == 0.0)
		theta = 0.0;
	else
		theta = std::atan2(Vec[2], rho);
}


void Vector::SetZero() {
	for (int i = 0; i < n; i++) { v[i] = 0; }
}
// Generate a matrix based off a rotation about the axis `axis` with an angle `angle`
Matrix3x3 RotationMatrix(long double angle, Vector axis) {
	if (axis.size() != 3) {
		cerr << "ERROR: Expected a size 3 vector for axis in RotationMatrix, instead got a size " << axis.size() << "vector!" << endl;
		exit(1);
	}
	auto L = Norm(axis);
	if (L < 1e-9) {
		cerr << "ERROR: Norm of axis vector was too small (or 0) in RotationMatrix\n";
		exit(1);
	}
	auto mat = Matrix3x3();
	auto c = std::cosl(angle);
	auto s = std::sinl(angle);
	auto x = axis[0] / L; auto y = axis[1] / L; auto z = axis[2] / L;

	mat[0][0] = c + x * x * (1 - c); mat[0][1] = x * y * (1 - c) - z * s; mat[0][2] = x * z * (1 - c) + y * s;
	mat[1][0] = y * x * (1 - c) + z * s; mat[1][1] = c + y * y * (1 - c); mat[1][2] = y * z * (1 - c) - x * s;
	mat[2][0] = z * x * (1 - c) - y * s; mat[2][1] = z * y * (1 - c) + x * s; mat[2][2] = c + z * z * (1 - c);
	return mat;
}

// Simpler X-Y-Z Rotation matrices coded below

Matrix3x3 XRotation(long double angle) {
	Matrix3x3 mat = Matrix3x3();
	mat[0][0] = 1.0; mat[1][1] = std::cosl(angle); mat[2][2] = mat[1][1];
	mat[2][1] = std::sinl(angle); mat[1][2] = -mat[2][1];
	return mat;
}
Matrix3x3 YRotation(long double angle) {
	Matrix3x3 mat = Matrix3x3();
	mat[1][1] = 1.0; mat[0][0] = std::cosl(angle); mat[2][2] = mat[0][0];
	mat[0][2] = std::sinl(angle); mat[2][0] = -mat[0][2];
	return mat;
}
Matrix3x3 ZRotation(long double angle) {
	Matrix3x3 mat = Matrix3x3();
	mat[2][2] = 1.0; mat[0][0] = std::cosl(angle); mat[1][1] = mat[0][0];
	mat[1][0] = std::sinl(angle); mat[0][1] = -mat[1][0];
	return mat;
}

// Matrix constructors
Matrix3x3::Matrix3x3()
{
	m = new long double* [3];
	for (int i = 0; i < 3; i++) {
		m[i] = new long double[3];
		for (int j = 0; j < 3; j++) {
			m[i][j] = 0.0;
		}
	}
}
Matrix3x3::Matrix3x3(long double entries[3][3])
{
	m = new long double* [3];
	for (int i = 0; i < 3; i++) {
		m[i] = new long double[3];
		for (int j = 0; j < 3; j++) {
			m[i][j] = entries[i][j];
		}
	}
}
Matrix3x3::Matrix3x3(long double entries[9])
{
	m = new long double* [3];
	for (int i = 0; i < 3; i++) {
		m[i] = new long double[3];
		for (int j = 0; j < 3; j++) {
			m[i][j] = entries[3*i + j];
		}
	}
}

// Matrix destructor
Matrix3x3::~Matrix3x3() { // HIGHLY NECESSARY TO AVOID DATA LEAKS!
	for (int i = 0; i < 3; i++) {
		delete[] m[i];
	}
	delete[] m;
}

// Scalar multiplication
Matrix3x3 operator *(long double scalar, const Matrix3x3& m) {
	Matrix3x3 Aux = Matrix3x3();
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			Aux[i][j] = scalar * m[i][j];
		}
	}
	return Aux;
}
// Matrix multiplication
Matrix3x3 operator*(const Matrix3x3& left, const Matrix3x3& right) {
	Matrix3x3 Aux;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				Aux[i][j] += left[i][k] * right[k][j];
			}
		}
	}
	return Aux;
}
// Matrix-vector multiplication
Vector operator*(const Matrix3x3& m, const Vector& v) {
	return Vector(m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
		m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
		m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2]);
}

// To make print statements with matrices semi-readable
std::ostream& operator<<(std::ostream& os, const Matrix3x3& mat) {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			os << mat[i][j] << " ";	//TODO: Maybe better padding?
		}
		os << endl;
	}
	return os;
}
