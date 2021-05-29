#include "aabb.h"

bool AABB::hit(const Ray& ray, double t_min, double t_max) const { // determine whether ray hits aabb
	for (int i = 0; i < 3; ++i) {
		double invD = 1.0 / ray.direction()[i]; // inverse of dir[i]
		/* // primary version
		double t0 = fmin((minimum[i] - ray.origin()[i]) / ray.direction()[i],
			(maximum[i] - ray.origin()[i]) / ray.direction()[i]);
		double t1 = fmax((minimum[i] - ray.origin()[i]) / ray.direction()[i],
			(maximum[i] - ray.origin()[i]) / ray.direction()[i]);
		*/
		double t0 = (minimum[i] - ray.origin()[i]) * invD;
		double t1 = (maximum[i] - ray.origin()[i]) * invD;
		if (invD < 0.0f) { std::swap(t0, t1); }
		t_min = fmax(t0, t_min);
		t_max = fmin(t1, t_max);
		if (t_max <= t_min) { return false; }
	}
	return true;
}

AABB AABB::surrounding_box(AABB box0, AABB box1) {
	Point3 point_min(fmin(box0.min().x(), box1.min().x()),
		fmin(box0.min().y(), box1.min().y()),
		fmin(box0.min().z(), box1.min().z()));
	Point3 point_max(fmax(box0.max().x(), box1.max().x()),
		fmax(box0.max().x(), box1.max().x()),
		fmax(box0.max().x(), box1.max().x()));
	return AABB(point_min, point_max);
}