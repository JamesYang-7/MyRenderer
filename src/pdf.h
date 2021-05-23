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

Vec3 random_cosine_direction() {
    double phi = 2 * PI * random_double();
    double r2 = random_double();
    double z = sqrt(1 - r2); // z = cos(arcsin(x))
    double x = cos(phi) * sqrt(r2);
    double y = sin(phi) * sqrt(r2);
    return Vec3(x, y, z);
}

Vec3 random_to_sphere(double radius, double distance_squared) {
    double phi = 2 * PI * random_double();
    double r2 = random_double();
    double z = 1 + r2 * (sqrt(1 - radius * radius / distance_squared) - 1);
    double x = cos(phi) * sqrt(1 - z * z);
    double y = sin(phi) * sqrt(1 - z * z);
    return Vec3(x, y, z);
}

class CosinePDF : public PDF {
public:
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
    ~MixturePDF() {
        delete p[0];
        delete p[1];
    }

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