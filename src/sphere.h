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

#endif
