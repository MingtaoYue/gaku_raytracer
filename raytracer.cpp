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

// ray class
struct Ray {
    // origin
    Vec o;
    // direction
    Vec d;
    // constructor
    Ray(Vec o, Vec d) {
        this->o = o;
        this->d = d;
    }
};

// reflection type
enum Refl_t {
    // diffuse
    DIFF,
    // specular
    SPEC,
    // refractive
    REFR
};

// sphere class
struct Sphere {
    // radius
    double rad;
    // center
    Vec p;
    // emission
    Vec e;
    // color
    Vec c;
    // reflection type
    Refl_t refl;
    // constructor
    Sphere(double rad, Vec p, Vec e, Vec c, Refl_t refl) {
        this->rad = rad;
        this->p = p;
        this->e = e;
        this->c = c;
        this->refl = refl;
    }
    // calculate intersect of sphere and ray
    double intersect(const Ray &r) const {
        // ray origin -> sphere center
        Vec op = p - r.o;
        // fudge factor
        double eps = 1e-4;
        // actually -1/2 b
        double b = op.dot(r.d);
        // (b^2 - 4ac) / 4, note that a = 1 as ray is normalized
        double det = b * b - (op.dot(op) - rad * rad);
        // ray misses sphere
        if (det < 0)
            return 0;
        det = sqrt(det);
        // return smaller positive t
        double t1 = b - det;
        if (t1 > eps)
            return t1;
        double t2 = b + det;
        if (t2 > eps)
            return t2;
        return 0;
    }
};

int main(int argc, char* argv[]) {
    // Create a sphere
    Sphere sphere(1.0, Vec(0.0, 0.0, 0.0), Vec(0.0, 0.0, 0.0), Vec(1.0, 0.0, 0.0), DIFF);

    // Create a ray (adjust the values as needed)
    Ray ray(Vec(0.0, 0.0, -5.0), Vec(0.0, 0.0, 1.0));

    // Test intersection
    double intersectionDistance = sphere.intersect(ray);

    std::cout << intersectionDistance << std::endl;

    return 0;
}