#include "bvh.h"

bool box_x_compare(const Hittable* a, const Hittable* b) {
	return box_compare(a, b, 0);
}

bool box_y_compare(const Hittable* a, const Hittable* b) {
	return box_compare(a, b, 1);
}

bool box_z_compare(const Hittable* a, const Hittable* b) {
	return box_compare(a, b, 2);
}

bool BVH_Node::bounding_box(double time0, double time1, AABB& output_box) const {
	output_box = box;
	return true;
}

bool BVH_Node::hit(const Ray& r, double t_min, double t_max, HitRecord& rec) const {
	if (!box.hit(r, t_min, t_max)) { return false; }
	bool hit_left = left->hit(r, t_min, t_max, rec);
	bool hit_right = right->hit(r, t_min, (hit_left ? rec.t - 1e-4 : t_max), rec);
	return hit_left || hit_right;
}

BVH_Node::BVH_Node(const std::vector<Hittable*>& src_objects, size_t start, size_t end, double time0, double time1) {
	std::vector<Hittable*> objects = src_objects;

	int axis = random_int(0, 2);
	auto comparator = (axis == 0) ? box_x_compare
		: (axis == 1) ? box_y_compare
		: box_z_compare;

	size_t object_span = end - start;
	if (object_span == 1) {
		left = right = objects[start];
	}
	else if (object_span == 2) {
		if (comparator(objects[start], objects[start + 1])) {
			left = objects[start];
			right = objects[start + 1];
		}
		else {
			left = objects[start + 1];
			right = objects[start];
		}
	}
	else {
		std::sort(objects.begin() + start, objects.begin() + end, comparator);
		size_t mid = start + object_span / 2;
		left = new BVH_Node(objects, start, mid, time0, time1);
		right = new BVH_Node(objects, mid, end, time0, time1);
	}

	AABB box_left, box_right;
	if (!left->bounding_box(time0, time1, box_left) || !right->bounding_box(time0, time1, box_right)) {
		std::cerr << "Error: No bounding box in bvh_node constructor.\n";
	}
	box = AABB::surrounding_box(box_left, box_right);
}