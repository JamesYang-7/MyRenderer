#ifndef AABB_H
#define AABB_H

#include "vec3.h"
#include "ray.h"
#include <cmath>

class AABB { // axis-aligned-bounding-box
public:
	AABB() {}
	AABB(const Point3& a, const Point3& b) { minimum = a; maximum = b; }

	Point3 min() const { return minimum; }
	Point3 max() const { return maximum; }

	bool hit(const Ray& ray, double t_min, double t_max) const;
	static AABB surrounding_box(AABB box0, AABB box1);
public:
	Point3 minimum;
	Point3 maximum;
};



#endif
