#ifndef RAY_H
#define RAY_H

#include "vec3.h"

class Ray {
public:
    Ray() {}
    Ray(Point3 origin, Vec3 direction)
        : orig(origin), dir(direction)
    {
        dir.normalize();
    }

    Point3 origin() const { return orig; }
    Vec3 direction() const { return dir; }

    Point3 at(double t) const {
        return orig + t * dir;
    }

private:
    Point3 orig;
    Vec3 dir;
};

#endif
