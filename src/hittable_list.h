#ifndef Hittable_LIST_H
#define Hittable_LIST_H

#include "hittable.h"
#include <vector>

class HittableList : public Hittable {
public:
    /*~HittableList() {
        for (Hittable* obj_ptr : objects) {
            if (obj_ptr != nullptr) {
                delete obj_ptr;
            }
        }
    }*/
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

#endif
