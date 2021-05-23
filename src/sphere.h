#ifndef SPHERE_H
#define SPHERE_H

#include "hittable.h"
#include "vec3.h"
#include "ray.h"
#include "onb.h"
#include "pdf.h"

class Sphere : public Hittable {
public:
	Sphere(Vec3 cen, double r, Material* m) :
		center(cen), radius(r), mat_ptr(m)
	{}

	virtual bool hit(const Ray& ray, double t_min, double t_max, HitRecord& rec) const override;
	virtual bool bounding_box(double time0, double time1, AABB& output_box) const override;

	double intersect(const Ray& ray) const { // sphere intersect with ray, return distance from ray.orig to intersection
		double t = 0;
		double eps = 1e-4;
		// Solve t^2*d.*d + 2*t*(o-c).*d + (o-c).*(o-c)-R^2 = 0
		// using |d|==1 and simplified some coefficients
		Vec3 op = ray.origin() - center;
		double b = dot(op, ray.direction());
		double delta = b * b - dot(op, op) + radius * radius;
		if (delta < 0) return 0;
		else { // notice that the ray is unidirectional, so only pick up the smaller positive solution
			double delta_sqrt = sqrt(delta);
			if ((t = -b - delta_sqrt) > eps) {
				return t;
			}
			else if ((t = -b + delta_sqrt) > eps) {
				return t;
			}
			else {
				return 0;
			}
		}
	}

	virtual double pdf_value(const Point3& origin, const Vec3& direction) const override;
	virtual Vec3 generate_random(const Vec3& origin) const override;
	
public:
	double radius;
	Vec3 center;
	Material* mat_ptr;
private:
	static void get_sphere_uv(const Point3& p, double& u, double& v) {
		double theta = acos(-p.y());
		double phi = atan2(-p.z(), p.x()) + PI;
		u = phi / (2 * PI);
		v = theta / PI;
	}
};

bool Sphere::hit(const Ray& ray, double t_min, double t_max, HitRecord& rec) const {
	double root = 0;
	// Solve t^2*d.*d + 2*t*(o-c).*d + (o-c).*(o-c)-R^2 = 0
	double eps = 0;
	Vec3 oc = ray.origin() - center;
	double b = dot(oc, ray.direction());
	double delta = b * b - dot(oc, oc) + radius * radius;
	if (delta < 0) return 0;
	double delta_sqrt = sqrt(delta);
	root = -b - delta_sqrt;
	if (root < t_min || root > t_max) {
		root = -b + delta_sqrt;
		if (root < t_min || root > t_max) {
			return false;
		}
	}
	rec.t = root;
	rec.p = ray.at(root);
	rec.normal = (rec.p - center) / radius;
	rec.mat_ptr = mat_ptr;
	Vec3 outward_normal = (rec.p - center) / radius;
	rec.set_face_normal(ray, outward_normal);
	get_sphere_uv(outward_normal, rec.u, rec.v);
	return true;
}

bool Sphere::bounding_box(double time0, double time1, AABB& output_box) const {
	double absr = fabs(radius);
	output_box = AABB(center - Vec3(absr, absr, absr), center + Vec3(absr, absr, absr));
	return true;
}

double Sphere::pdf_value(const Point3& origin, const Vec3& direction) const {
	HitRecord hrec;
	if (!this->hit(Ray(origin, direction), 1e-4, infinity, hrec)) {
		return 0;
	}
	double cos_theta_max = sqrt(1 - radius * radius / (center - origin).length_squared());
	double solid_angle = 2 * PI * (1 - cos_theta_max);
	return 1 / solid_angle;
}

Vec3 Sphere::generate_random(const Vec3& origin) const {
	Vec3 direction = center - origin;
	auto distance_squared = direction.length_squared();
	ONB uvw;
	uvw.build_from_w(direction);
	return uvw.local(random_to_sphere(radius, distance_squared));
}

#endif
