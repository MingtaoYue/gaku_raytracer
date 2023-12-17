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
    t = 1e20;
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

Vec radiance(const Ray &r, int depth, unsigned short *Xi) {
    // distance to the intersection
    double t;
    // id of the intersected object
    int id = 0;
    // if no intersection, return black
    if (!intersect(r, t, id)) 
        return Vec();
    // first intersected object
    const Sphere &obj = spheres[id];
    // intersection point
    Vec x = r.o + r.d * t;
    // normal vector at the intersection point
    Vec n = (x - obj.p).norm();
    // surface normal that is oppose to the ray direction
    Vec nl = n.dot(r.d) < 0 ? n : n * -1;
    // object color
    Vec f = obj.c;
    // take the max reflectivity
    double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z;
    // increase level of recursion
    depth++;
    if (depth > 5) {
        // russian roulette
        if (erand48(Xi) < p) 
            f = f * (1 / p);
        else return obj.e;
    }
    // ideal diffuse reflection
    if (obj.refl == DIFF){
        // random generation
        double r1 = 2 * M_PI * erand48(Xi), r2 = erand48(Xi), r2s = sqrt(r2);
        // reflection direction
        Vec w = nl;
        // perpendicular to the surface normal and lies in the plane of reflection
        Vec u = ((fabs(w.x)>.1?Vec(0,1):Vec(1)).cross(w)).norm();
        // perpendicular to w and u
        Vec v = w.cross(u);
        // direction of the reflected ray
        Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm(); 
        return obj.e + f.mult(radiance(Ray(x, d), depth, Xi));
    }
    // ideal specular reflection 
    else if (obj.refl == SPEC)
        return obj.e + f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth, Xi));
    return Vec(0, 0, 0);
}

int main(int argc, char* argv[]) {
    // image size
    int w = 1024, h = 768;
    // # of samples, default is 1
    int samps = (argc == 2) ? atoi(argv[1]) / 4 : 1;
    // colors of samples
    Vec r;
    // init image
    Vec *img = new Vec[w * h];

    // camera position, gaze direction
    Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm());
    // x, y direction of camera, 0.5135 defines filed of view angle
    Vec cx = Vec(w * .5135 / h, 0, 0), cy = (cx.cross(cam.d)).norm() * .5135;
    
    // OpenMP
#pragma omp parallel for schedule(dynamic, 1) private(r)
    // loop over rows
    for (int y = 0; y < h; y++) {
        // random seed
        unsigned short Xi[3] = {0, 0, y * y * y};
        // loop over columns
        for (int x = 0; x < w; x++) {
            // current pixel
            int i = (h - y - 1) * w + x;
            // loop over 2x2 sub-pixel rows
            for (int sy = 0; sy < 2; sy++) {
                for (int sx = 0; sx < 2; sx++) {
                    // init radiance
                    r = Vec();
                    // loop over samples
                    for (int s = 0; s < samps; s++) {
                        // generate a random position within the pixel
                        double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1: 1 - sqrt(2 - r1); 
                        double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1: 1 - sqrt(2 - r2);
                        // compute direction
                        Vec d = cx *(((sx + .5 + dx) / 2 + x) / w - .5) + cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d; 
                        // average the radiance over one sub-pixel
                        r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0, Xi) * (1. / samps);
                    }
                    // average the radiance over one pixel
                    img[i] = img[i] + Vec(clamp(r.x),clamp(r.y),clamp(r.z))*.25; 
                }
            }
        }   
    }
    // write image
    FILE *f = fopen("image.ppm", "w");
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255); 
    for (int i = 0; i < w * h; i++) 
        fprintf(f,"%d %d %d ", toInt(img[i].x), toInt(img[i].y), toInt(img[i].z));
    return 0;
}

