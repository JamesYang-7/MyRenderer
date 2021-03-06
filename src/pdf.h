#ifndef PDF_H
#define PDF_H

#include "rtweekend.h"
#include "hittable.h"
#include "onb.h"

class PDF {
public:
    virtual ~PDF() {}
    virtual double value(const Vec3& direction) const = 0;
    virtual Vec3 generate() const = 0;
};

Vec3 random_cosine_direction();
Vec3 random_to_sphere(double radius, double distance_squared);

class CosinePDF : public PDF {
public:
    ~CosinePDF() {}
    CosinePDF(const Vec3& w) { uvw.build_from_w(w); }

    virtual double value(const Vec3& direction) const override {
        auto cosine = dot(unit_vector(direction), uvw.w());
        return (cosine <= 0) ? 0 : cosine / PI;
    }

    virtual Vec3 generate() const override {
        return uvw.local(random_cosine_direction());
    }

public:
    ONB uvw;
};

class HittableCustomPDF : public PDF {
public:
    ~HittableCustomPDF() {}
    HittableCustomPDF(Hittable* ptr_, const Point3& origin) : ptr(ptr_), o(origin) {}
    virtual double value(const Vec3& direction) const override {
        return ptr->pdf_value(o, direction);
    }
    virtual Vec3 generate() const override {
        return ptr->generate_random(o);
    }
public:
    Point3 o;
    Hittable* ptr;
};

class MixturePDF : public PDF {
public:
    ~MixturePDF() {} // be aware of memory leak, must delete p[i] somewhere else

    MixturePDF(PDF* pdf_important, PDF* pdf_normal, double imp_ratio) : importance_ratio(imp_ratio) {
        if (importance_ratio < 0 || importance_ratio > 1) importance_ratio = 0;
        p[0] = pdf_important;
        p[1] = pdf_normal;
    }

    virtual double value(const Vec3& direction) const override {
        return importance_ratio * p[0]->value(direction) + (1 - importance_ratio) * p[1]->value(direction);
    }

    virtual Vec3 generate() const override {
        if (random_double() < importance_ratio) {
            return p[0]->generate();
        }
        else {
            return p[1]->generate();
        }
    }
public:
    PDF* p[2];
    double importance_ratio = 0; // probability of choosing p[0]
};

#endif