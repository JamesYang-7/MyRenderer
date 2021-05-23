#ifndef MATERIAL_H
#define MATERIAL_H

#include "rtweekend.h"
#include "texture.h"
#include "pdf.h"

class Material {
public:
    virtual Color emitted(const Ray& r_in, const HitRecord& hrec, double u, double v, const Point3& p) const {
        return Color();
    }
    virtual bool scatter(const Ray& r_in, const HitRecord& hrec, ScatterRecord& srec) const = 0;   
    virtual double scattering_pdf(const Ray& ray_in, const HitRecord& hrec, const Ray& scattered) const {
        return 0;
    }
};

class Diffuse : public Material {
public:
    Diffuse(const Color& a) : albedo(new SolidColor(a)) {}
    Diffuse(Texture* a) : albedo(a){}

    virtual bool scatter(const Ray& r_in, const HitRecord& hrec, ScatterRecord& srec) const override {
        srec.is_specular = false;
        srec.p_pdf = new CosinePDF(hrec.normal);
        Vec3 scatter_direction = srec.p_pdf->generate();
        srec.scattered = Ray(hrec.p, scatter_direction);
        srec.albedo = albedo->value(hrec.u, hrec.v, hrec.p);
        return true;
    }

    virtual double scattering_pdf(const Ray& ray_in, const HitRecord& hrec, const Ray& scattered) const {
        double cosine = dot(hrec.normal, scattered.direction());
        return cosine < 0 ? 0 : cosine / PI;
    }
public:
    Texture* albedo;
};

class Specular : public Material {
public:
    Specular(const Color& a) : albedo(a) {}

    virtual bool scatter(const Ray& r_in, const HitRecord& hrec, ScatterRecord& srec) const override {
        srec.is_specular = true;
        Vec3 reflected = reflect(r_in.direction(), hrec.normal);
        srec.scattered = Ray(hrec.p, reflected);
        srec.albedo = albedo;
        return dot(srec.scattered.direction(), hrec.normal) > 0;
    }

public:
    Color albedo;
};

class DiffuseLight : public Material {
public:
    DiffuseLight(Texture* a) : emit(a) {}
    DiffuseLight(Color c) : emit(new SolidColor(c)) {}

    virtual bool scatter(const Ray& r_in, const HitRecord& hrec, ScatterRecord& srec) const override {
        return false;
    }

    virtual Color emitted(const Ray& r_in, const HitRecord& hrec, double u, double v, const Point3& p) const override { // *unfinished
        // should be relative to rec.frontface
        if (hrec.front_face) {
            return emit->value(u, v, p);
        }
        else {
            return Vec3();
        }
    }

public:
    Texture* emit;
};

class Dielectric : public Material {
public:
    Dielectric(Texture* a, double idx_rfr) : albedo(a), index_of_refr(idx_rfr) {};
    virtual bool scatter(const Ray& r_in, const HitRecord& hrec, ScatterRecord& srec) const override {
        srec.is_specular = true;
        double refr_ratio = hrec.front_face ? (1.0 / index_of_refr) : index_of_refr;
        double cos_theta = fmin(dot(-r_in.direction(), hrec.normal), 1.0);
        double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
        Vec3 scatter_direction;
        if (refr_ratio * sin_theta >= 1 || reflectance(cos_theta, refr_ratio) > random_double()) { // if total internal reflection or
            scatter_direction = reflect(r_in.direction(), hrec.normal);
        }
        else {
            scatter_direction = refract(r_in.direction(), hrec.normal, refr_ratio);
        }
        srec.scattered = Ray(hrec.p, scatter_direction);
        srec.albedo = albedo->value(hrec.u, hrec.v, hrec.p);
        return true;
    }

public:
    double index_of_refr;
    Texture* albedo;
private:
    static double reflectance(double cosine, double refl_idx) {
        // calculate Fresnel term using Schlick's approximation
        double r0 = (1 - refl_idx) / (1 + refl_idx);
        r0 = r0 * r0;
        return r0 + (1 - r0) * pow((1 - cosine), 5);
    }
};

#endif
