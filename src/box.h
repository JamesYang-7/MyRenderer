#ifndef BOX_H
#define BOX_H

#include "rtweekend.h"
#include "aarect.h"
#include "hittable_list.h"

class Box : public Hittable {
public:
	Box() {}
	Box(const Point3& p0, const Point3& p1, Material* ptr);

	virtual bool hit(const Ray& ray, double t_min, double t_max, HitRecord& rec) const override;
	virtual bool bounding_box(double time0, double time1, AABB& output_box) const override {
		output_box = AABB(box_min, box_max);
		return true;
	}
public:
	Point3 box_min;
	Point3 box_max;
	HittableList sides;
};

#endif