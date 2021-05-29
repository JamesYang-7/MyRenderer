#include "material.h"

bool Diffuse::scatter(const Ray& r_in, const HitRecord& hrec, ScatterRecord& srec) const {
    srec.is_specular = false;
    srec.p_pdf = new CosinePDF(hrec.normal);
    Vec3 scatter_direction = srec.p_pdf->generate();
    srec.scattered = Ray(hrec.p, scatter_direction);
    srec.albedo = albedo->value(hrec.u, hrec.v, hrec.p);
    return true;
}

bool Dielectric::scatter(const Ray& r_in, const HitRecord& hrec, ScatterRecord& srec) const {
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