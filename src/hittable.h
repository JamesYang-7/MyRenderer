#ifndef HITTABLE_H
#define HITTABLE_H

#include "aabb.h"
#include "rtweekend.h"

class Material;
class PDF;

struct HitRecord {
	Point3 p;
	Vec3 normal;
	Material* mat_ptr = nullptr;
	double t = 0; // ray distance
	double u = 0; // texture coordinate
	double v = 0; // texture coordinate
	bool front_face = true;

	inline void set_face_normal(const Ray& ray, const Vec3& outward_normal) {
		front_face = dot(ray.direction(), outward_normal) < 0;
		normal = front_face ? outward_normal : -outward_normal;
	}
};

struct ScatterRecord {
	Ray scattered;
	Color albedo;
	bool is_specular;
	PDF* p_pdf;

	ScatterRecord() : is_specular(false), p_pdf(nullptr) {}
	~ScatterRecord() { if (p_pdf != nullptr) delete p_pdf; }
};

class Hittable {
public:
	virtual bool hit(const Ray& ray, double t_min, double t_max, HitRecord& rec) const = 0;
	virtual bool bounding_box(double time0, double time1, AABB& output_box) const = 0;
	virtual double pdf_value(const Point3& origin, const Vec3& direction) const { 
		return 0.0;
	}
	virtual Vec3 generate_random(const Vec3& origin) const {
		return Vec3(1, 0, 0);
	}
};

class Translate : public Hittable {
public:
	Translate(Hittable* p, const Vec3& displacement)
		: ptr(p), offset(displacement) {}

	virtual bool hit(const Ray& r, double t_min, double t_max, HitRecord& rec) const override;
	virtual bool bounding_box(double time0, double time1, AABB& output_box) const override;

public:
	Hittable* ptr;
	Vec3 offset;
};

class Rotate_Y : public Hittable {
public:
	Rotate_Y(Hittable* p, double angle);
	virtual bool hit(const Ray& ray, double t_min, double t_max, HitRecord& rec) const override;
	virtual bool bounding_box(double time0, double time1, AABB& output_box) const override {
		output_box = bbox;
		return hasbox;
	}
public:
	Hittable* ptr;
	double sin_theta;
	double cos_theta;
	bool hasbox;
	AABB bbox;
};

#endif
