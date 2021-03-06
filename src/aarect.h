#ifndef AARECT_H
#define AARECT_H

#include "rtweekend.h"
#include "hittable.h"
#include "pdf.h"

class XY_Rect : public Hittable {
public:
	XY_Rect() : mp(nullptr), x0(0), x1(0), y0(0), y1(0), k(0) {}
	XY_Rect(double x0_, double x1_, double y0_, double y1_, double k_, Material* mat)
		: x0(x0_), x1(x1_), y0(y0_), y1(y1_), k(k_), mp(mat)
	{};
	virtual bool hit(const Ray& ray, double t_min, double t_max, HitRecord& rec) const override;
	virtual bool bounding_box(double time0, double time1, AABB& output_box) const override {
		output_box = AABB(Point3(x0, y0, k - 1e-4), Point3(x1, y1, k + 1e-4));
		return true;
	}
public:
	Material* mp;
	double x0, x1, y0, y1;
	double k;
};

class XZ_Rect : public Hittable {
public:
	XZ_Rect() : mp(nullptr), x0(0), x1(0), z0(0), z1(0), k(0) {}
	XZ_Rect(double x0_, double x1_, double z0_, double z1_, double k_, Material* mat)
		: x0(x0_), x1(x1_), z0(z0_), z1(z1_), k(k_), mp(mat)
	{};
	virtual bool hit(const Ray& ray, double t_min, double t_max, HitRecord& rec) const override;

	virtual bool bounding_box(double time0, double time1, AABB& output_box) const override {
		output_box = AABB(Point3(x0, k - 1e-4, z0), Point3(x1, k + 1e-4, z1));
		return true;
	}

	virtual double pdf_value(const Point3& org, const Vec3& dir) const override {
		HitRecord rec;
		Ray r = Ray(org, dir);
		if (!hit(r, 1e-4, infinity, rec)) { return 0; }
		double area = (x1 - x0) * (z1 - z0);
		double distance_squared = rec.t * rec.t;
		double cosine = fabs(dot(rec.normal, r.direction()));
		return distance_squared / (cosine * area);
	}

	virtual Vec3 generate_random(const Vec3& origin) const override {
		Point3 random_point = Point3(random_double(x0, x1), k, random_double(z0, z1));
		return random_point - origin;
	}

	Ray generate_emitted_ray() const {
		Point3 random_point = Point3(random_double(x0, x1), k, random_double(z0, z1));
		CosinePDF cos_pdf(Vec3(0, -1, 0));
		Vec3 random_direction = cos_pdf.generate();
		return Ray(random_point, random_direction);
	}

	Point3 random_point() const {
		return Point3(random_double(x0, x1), k, random_double(z0, z1));
	}

public:
	Material* mp;
	double x0, x1, z0, z1;
	double k;
};

class YZ_Rect : public Hittable {
public:
	YZ_Rect() : mp(nullptr), y0(0), y1(0), z0(0), z1(0), k(0) {}
	YZ_Rect(double y0_, double y1_, double z0_, double z1_, double k_, Material* mat)
		: y0(y0_), y1(y1_), z0(z0_), z1(z1_), k(k_), mp(mat)
	{};
	virtual bool hit(const Ray& ray, double t_min, double t_max, HitRecord& rec) const override;
	virtual bool bounding_box(double time0, double time1, AABB& output_box) const override {
		output_box = AABB(Point3(k - 1e-4, y0, z0), Point3(k + 1e-4, y1, z1));
		return true;
	}
public:
	Material* mp;
	double y0, y1, z0, z1;
	double k;
};

#endif
