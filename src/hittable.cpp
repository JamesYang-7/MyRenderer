#include "hittable.h"

bool Translate::hit(const Ray& r, double t_min, double t_max, HitRecord& rec) const {
	Ray moved_r(r.origin() - offset, r.direction());
	if (!ptr->hit(moved_r, t_min, t_max, rec)) { return false; }
	rec.p += offset;
	rec.set_face_normal(moved_r, rec.normal);

	return true;
}

bool Translate::bounding_box(double time0, double time1, AABB& output_box) const {
	if (!ptr->bounding_box(time0, time1, output_box)) { return false; }
	output_box = AABB(output_box.min() + offset, output_box.max() + offset);
	return true;
}

Rotate_Y::Rotate_Y(Hittable* p, double angle) : ptr(p) {
	double radians = degrees_to_radians(angle);
	sin_theta = sin(radians);
	cos_theta = cos(radians);
	hasbox = ptr->bounding_box(0, 1, bbox);

	Point3 min(infinity, infinity, infinity);
	Point3 max(-infinity, -infinity, -infinity);

	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			for (int k = 0; k < 2; ++k) {
				double x = i * bbox.max().x() + (1 - i) * bbox.min().x();
				double y = i * bbox.max().y() + (1 - i) * bbox.min().y();
				double z = i * bbox.max().z() + (1 - i) * bbox.min().z();

				double newx = cos_theta * x + sin_theta * z;
				double newz = -sin_theta * x + cos_theta * z;
				Vec3 tester(newx, y, newz);
				for (int c = 0; c < 3; ++c) {
					min[c] = fmin(min[c], tester[c]);
					max[c] = fmax(max[c], tester[c]);
				}
			}
		}
	}
	bbox = AABB(min, max);
}

bool Rotate_Y::hit(const Ray& ray, double t_min, double t_max, HitRecord& rec) const {
	Vec3 origin = ray.origin();
	Vec3 direction = ray.direction();

	origin[0] = cos_theta * ray.origin()[0] - sin_theta * ray.origin()[2];
	origin[2] = sin_theta * ray.origin()[0] + cos_theta * ray.origin()[2];

	direction[0] = cos_theta * ray.direction()[0] - sin_theta * ray.direction()[2];
	direction[2] = sin_theta * ray.direction()[0] + cos_theta * ray.direction()[2];

	Ray rotated_r(origin, direction);
	if (!ptr->hit(rotated_r, t_min, t_max, rec)) { return false; }

	Vec3 p = rec.p;
	Vec3 normal = rec.normal;

	p[0] = cos_theta * rec.p[0] + sin_theta * rec.p[2];
	p[2] = -sin_theta * rec.p[0] + cos_theta * rec.p[2];

	normal[0] = cos_theta * rec.normal[0] + sin_theta * rec.normal[2];
	normal[2] = -sin_theta * rec.normal[0] + cos_theta * rec.normal[2];

	rec.p = p;
	rec.set_face_normal(rotated_r, normal);
	return true;
}