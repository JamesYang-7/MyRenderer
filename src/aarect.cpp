#include "aarect.h"

bool XY_Rect::hit(const Ray& ray, double t_min, double t_max, HitRecord& rec) const {
	double t = (k - ray.origin().z()) / ray.direction().z();
	if (t < t_min || t > t_max) { return false; }
	double x = ray.origin().x() + t * ray.direction().x();
	double y = ray.origin().y() + t * ray.direction().y();

	if (x < x0 || x > x1 || y < y0 || y > y1) { return false; }
	rec.u = (x - x0) / (x1 - x0);
	rec.v = (y - y0) / (y1 - y0);
	rec.t = t;
	Vec3 outward_normal = Vec3(0, 0, 1);
	rec.set_face_normal(ray, outward_normal);
	rec.mat_ptr = mp;
	rec.p = ray.at(t);
	return true;
}

bool XZ_Rect::hit(const Ray& ray, double t_min, double t_max, HitRecord& rec) const {
	double t = (k - ray.origin().y()) / ray.direction().y();
	if (t < t_min || t > t_max) { return false; }
	double x = ray.origin().x() + t * ray.direction().x();
	double z = ray.origin().z() + t * ray.direction().z();

	if (x < x0 || x > x1 || z < z0 || z > z1) { return false; }
	rec.u = (x - x0) / (x1 - x0);
	rec.v = (z - z0) / (z1 - z0);
	rec.t = t;
	Vec3 outward_normal = Vec3(0, -1, 0);
	rec.set_face_normal(ray, outward_normal);
	rec.mat_ptr = mp;
	rec.p = ray.at(t);
	return true;
}

bool YZ_Rect::hit(const Ray& ray, double t_min, double t_max, HitRecord& rec) const {
	double t = (k - ray.origin().x()) / ray.direction().x();
	if (t < t_min || t > t_max) { return false; }
	double y = ray.origin().y() + t * ray.direction().y();
	double z = ray.origin().z() + t * ray.direction().z();

	if (y < y0 || y > y1 || z < z0 || z > z1) { return false; }
	rec.u = (y - y0) / (y1 - y0);
	rec.v = (z - z0) / (z1 - z0);
	rec.t = t;
	Vec3 outward_normal = Vec3(1, 0, 0);
	rec.set_face_normal(ray, outward_normal);
	rec.mat_ptr = mp;
	rec.p = ray.at(t);
	return true;
}