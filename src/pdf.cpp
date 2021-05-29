#include "pdf.h"

Vec3 random_cosine_direction() {
    double phi = 2 * PI * random_double();
    double r2 = random_double();
    double z = sqrt(1 - r2); // z = cos(arcsin(x))
    double x = cos(phi) * sqrt(r2);
    double y = sin(phi) * sqrt(r2);
    return Vec3(x, y, z);
}

Vec3 random_to_sphere(double radius, double distance_squared) {
    double phi = 2 * PI * random_double();
    double r2 = random_double();
    double z = 1 + r2 * (sqrt(1 - radius * radius / distance_squared) - 1);
    double x = cos(phi) * sqrt(1 - z * z);
    double y = sin(phi) * sqrt(1 - z * z);
    return Vec3(x, y, z);
}