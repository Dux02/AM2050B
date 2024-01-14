#pragma once
#include <iostream>
class Vector {
	
	public:
		// Constructors
		Vector();								// Empty vector
		Vector(int Size);						// Null vector with given size
		Vector(const Vector& V);				// Copy from another vector
		Vector(long double x, long double y, long double z);	// Make 3D array

		// Destructor
		~Vector();

		// Assignment
		Vector& operator=(const long double value);
		Vector& operator=(const Vector& V);

		int size() const { return n; };
		// Vector equality
		friend bool operator == (const Vector& left, const Vector& right);

		// Coordinate access
		long double operator [] (int i) const { return v[i]; }
		long double& operator [] (int i)		 { return v[i]; }

		// Vector addition and subtraction
		friend Vector operator + (const Vector& left, const Vector& right);
		friend Vector operator - (const Vector& left, const Vector& right);

		// Vector addition/subtraction with assignment
		void operator += (const Vector& V);
		void operator -= (const Vector& V);

		// Vector operations
		friend long double Dot(const Vector& left, const Vector& right);
		friend long double Norm(const Vector& V);
		friend Vector Cross(const Vector& left, const Vector& right);

		//Scalar multiplication and division
		friend Vector operator * (long double value, const Vector& V);
		friend Vector operator * (const Vector& V, long double value);
		friend Vector operator / (const Vector& V, long double value);

		friend Vector VecPolar(long double r, long double theta, long double phi);	//Make a vector from polar coordinates. Theta is angle w/ z-axis
		friend void CalcPolarCoords(const Vector& Vec, long double& r, long double& theta, long double& phi); //Calculate polar coordinates. Uses pointers
		void SetZero();
		

		// Output
		friend std::ostream& operator << (std::ostream& os, const Vector& Vec);

		friend class Matrix3x3;

private:
	int n;			// Dimension
	long double* v;		// Vector v[n]
};

class Matrix3x3 {
public:
	long double** m;
	Matrix3x3();		//Zero-matrix
	Matrix3x3(long double entries[3][3]);
	Matrix3x3(long double entries[9]);

	~Matrix3x3();

	// Coordinate access
	long double* operator [] (int i) const { return m[i]; }

	friend Matrix3x3 operator * (long double scalar, const Matrix3x3& m);
	friend Matrix3x3 operator * (const Matrix3x3& left, const Matrix3x3& right);
	friend Vector operator * (const Matrix3x3& m, const Vector& v);

	friend std::ostream& operator << (std::ostream& os, const Matrix3x3& mat);
};
Matrix3x3 RotationMatrix(long double angle, Vector axis);
Matrix3x3 XRotation(long double angle);
Matrix3x3 YRotation(long double angle);
Matrix3x3 ZRotation(long double angle);