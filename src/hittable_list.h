#ifndef Hittable_LIST_H
#define Hittable_LIST_H

#include "hittable.h"

#include <vector>

class HittableList : public Hittable {
public:
    HittableList() {}
    HittableList(Hittable* object) { add(object); }

    void clear() { objects.clear(); }
    void add(Hittable* object) { objects.push_back(object); }

    virtual bool hit(const Ray& r, double t_min, double t_max, HitRecord& rec) const override;
    virtual bool bounding_box(double time0, double time1, AABB& output_box) const override;
    virtual double pdf_value(const Point3& origin, const Vec3& direction) const override;
    virtual Vec3 generate_random(const Vec3& origin) const override;

public:
    std::vector<Hittable*> objects;
};

bool HittableList::hit(const Ray& r, double t_min, double t_max, HitRecord& rec) const {
    HitRecord temp_rec;
    bool hit_anything = false;
    double closest_so_far = t_max;

    for (const auto& object : objects) {
        if (object->hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t - 1e-4;
            rec = temp_rec;
        }
    }

    return hit_anything;
}

bool HittableList::bounding_box(double time0, double time1, AABB& output_box) const {
    if (objects.empty()) return false;
    AABB temp_box;
    bool first_box = true;

    for (const Hittable* object : objects) {
        if (!object->bounding_box(time0, time1, temp_box)) { return false; }
        output_box = first_box ? temp_box : surrounding_box(output_box, temp_box);
        first_box = false;
    }
    return true;
}

double HittableList::pdf_value(const Point3& origin, const Vec3& direction) const
{
    double weight = 1.0 / objects.size();
    double sum = 0;
    for (const Hittable* object : objects) {
        sum += weight * object->pdf_value(origin, direction);
    }
    return sum;
}

Vec3 HittableList::generate_random(const Vec3& origin) const
{
    int n = objects.size();
    return objects[random_int(0, n - 1)]->generate_random(origin);
}

#endif
