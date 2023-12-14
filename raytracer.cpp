#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

// vector class for 3D point, normal vector and color.
struct Vec {
    // coordinates, rgb values for color.
    double x, y, z;
    // constructor
    Vec(double x = 0, double y = 0, double z = 0) {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    // overload "+"
    Vec operator+(const Vec &v) const {
        return Vec(x + v.x, y + v.y, z + v.z);
    }
    // overload "-"
    Vec operator-(const Vec &v) const {
        return Vec(x - v.x, y - v.y, z - v.z);
    }
    // overload "*", scalar multiplication
    Vec operator*(double n) const {
        return Vec(x * n, y * n, z * n);
    }
    // Hadamard product ()
    Vec mult(const Vec &v) const {
        return Vec(x * v.x, y * v.y, z * v.z);
    }
    // normalize a vector
    Vec& norm() {
        return *this = *this * (1 / sqrt(x * x + y * y + z * z));
    }
    // dot product
    double dot(const Vec &v) const {
        return x * v.x + y * v.y + z * v.z;
    }
    // cross product
    Vec cross(const Vec &v) const {
        return Vec(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }
};

int main(int argc, char* argv[]) {
    // Create two vectors
    Vec vec1(1.0, 2.0, 3.0);
    Vec vec2(4.0, 5.0, 6.0);

    // Perform vector operations
    Vec sum = vec1 + vec2;
    Vec difference = vec1 - vec2;
    Vec scalarProduct = vec1 * 2.0;
    Vec hadamardProduct = vec1.mult(vec2);
    double dotProduct = vec1.dot(vec2);
    Vec crossProduct = vec1.cross(vec2);

    // Print the results
    std::cout << "Sum: (" << sum.x << ", " << sum.y << ", " << sum.z << ")\n";
    std::cout << "Difference: (" << difference.x << ", " << difference.y << ", " << difference.z << ")\n";
    std::cout << "Scalar Product: (" << scalarProduct.x << ", " << scalarProduct.y << ", " << scalarProduct.z << ")\n";
    std::cout << "Hadamard Product: (" << hadamardProduct.x << ", " << hadamardProduct.y << ", " << hadamardProduct.z << ")\n";
    std::cout << "Dot Product: " << dotProduct << "\n";
    std::cout << "Cross Product: (" << crossProduct.x << ", " << crossProduct.y << ", " << crossProduct.z << ")\n";

    return 0;
}