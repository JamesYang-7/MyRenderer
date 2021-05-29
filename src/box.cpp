#include "box.h"

Box::Box(const Point3& p0, const Point3& p1, Material* ptr) {
	box_min = p0;
	box_max = p1;
	sides.add(new XY_Rect(p0.x(), p1.x(), p0.y(), p1.y(), p1.z(), ptr));
	sides.add(new XY_Rect(p0.x(), p1.x(), p0.y(), p1.y(), p0.z(), ptr));

	sides.add(new XZ_Rect(p0.x(), p1.x(), p0.z(), p1.z(), p1.y(), ptr));
	sides.add(new XZ_Rect(p0.x(), p1.x(), p0.z(), p1.z(), p0.y(), ptr));

	sides.add(new YZ_Rect(p0.y(), p1.y(), p0.z(), p1.z(), p1.x(), ptr));
	sides.add(new YZ_Rect(p0.y(), p1.y(), p0.z(), p1.z(), p0.x(), ptr));
}

bool Box::hit(const Ray& ray, double t_min, double t_max, HitRecord& rec) const {
	return sides.hit(ray, t_min, t_max, rec);
}