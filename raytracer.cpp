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

// Scene definition, mainly made of spheres.
// Note that when the radius is big enough, the sphere becomes a plane, which forms the background.
Sphere spheres[] = {
    // left wall
    Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), DIFF),
    // right wall
    Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), DIFF),
    // back wall
    Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFF), 
    // front wall
    Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(), DIFF), 
    // bottom wall
    Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFF), 
    // top wall
    Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75), DIFF),
    // mirror sphere
    Sphere(16.5, Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1) * .999, SPEC),
    // glass sphere
    Sphere(16.5, Vec(73, 16.5, 78), Vec(), Vec(1, 1, 1) * .999, REFR),
    // light source
    Sphere(600, Vec(50, 681.6 - .27, 81.6),Vec(12, 12, 12), Vec(), DIFF)
};

// clamp the number to be in range [0, 1]
double clamp(double x) {
    if (x < 0)
        return 0;
    if (x > 1)
        return 1;
    return x;
}

// convert colors to displayable range [0, 255]
int toInt(double x) {
    // clamp the number to be in range [0, 1]
    double clamped_x = clamp(x);
    // gamma correction of 2.2
    double gamma_x = pow(clamped_x, 1 / 2.2);
    // convert to range [0, 255]
    return int(gamma_x * 255 + .5);
}

// check intersection with the scene
bool intersect(const Ray &r, double &t, int &id) {
    // # of objects/spheres in the scene
    double n = sizeof(spheres) / sizeof(Sphere);
    // init t to a very big number
    double t = 1e20;
    // current minimum t
    double t_min;
    // check each of the sphere, choose the minimum distance
    for (int i = 0; i < n; i++) {
        t_min = spheres[i].intersect(r);
        if (t_min > 0 && t_min < t) {
            t = t_min;
            // first intersected sphere id
            id = i;
        }
    }
    return t < 1e20;
}

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

