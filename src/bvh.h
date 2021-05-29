#ifndef BVH_H
#define BVH_H

#include <algorithm>
#include "rtweekend.h"
#include "hittable.h"
#include "hittable_list.h"

inline bool box_compare(const Hittable* a, const Hittable* b, int axis) {
	AABB box_a;
	AABB box_b;

	if (!a->bounding_box(0, 0, box_a) || !b->bounding_box(0, 0, box_b)) {
		std::cerr << "No bounding box in bvh_node constructor.\n";
	}
	return box_a.min().e[axis] < box_b.min().e[axis];
}

bool box_x_compare(const Hittable* a, const Hittable* b);

bool box_y_compare(const Hittable* a, const Hittable* b);

bool box_z_compare(const Hittable* a, const Hittable* b);

class BVH_Node : public Hittable {
public:
	BVH_Node() : left(nullptr), right(nullptr) {};
	BVH_Node(const HittableList& list, double time0, double time1) :
		BVH_Node(list.objects, 0, list.objects.size(), time0, time1) {}
	BVH_Node(const std::vector<Hittable*>& src_objects, size_t start, size_t end, double time0, double time1);

	virtual bool hit(const Ray& r, double t_min, double t_max, HitRecord& rec) const override;
	virtual bool bounding_box(double time0, double time1, AABB& output_box) const override;
public:
	Hittable* left;
	Hittable* right;
	AABB box;
};

#endif