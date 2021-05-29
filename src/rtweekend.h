#ifndef RTWEEKEND_H
#define RTWEEKEND_H

#include <cmath>
#include <limits>
#include <random>

// Common Headers

#include "ray.h"
#include "vec3.h"

// Usings

using std::sqrt;

// Constants

extern const double infinity;
extern const double PI;
// extern const double sin_table[360];
// extern const double cos_table[360];

// Utility Functions

/*inline double sintb(double x) {
	return sin_table[static_cast<int>(x * 180.0 / PI)];
}

inline double costb(double x) {
	return cos_table[static_cast<int>(x * 180.0 / PI)];
}*/

inline double degrees_to_radians(double degrees) {
    return degrees * PI / 180.0;
}

inline double clamp(double x, double min, double max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

double clamp(double x);

double gamma(double x);

double random_double();

double random_double(double min, double max);

int random_int(int min, int max);

Vec3 random_in_unit_disk();

Vec3 random_in_unit_sphere();

Vec3 random_in_unit_semisphere(Vec3 normal);

Vec3 reflect(const Vec3& v, const Vec3& n);

Vec3 refract(const Vec3 in_dir, const Vec3& n, double rior); // rior == reflective index of refraction

#endif
