#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <random>
#include <ctime>
#include <string>
#include "rtweekend.h"
#include "color.h"
#include "hittable_list.h"
#include "sphere.h"
#include "camera.h"
#include "material.h"
#include "aarect.h"
#include "box.h"
#include "pdf.h"
#include "bvh.h"
#include "photon_map.h"

using namespace std;

enum class ParamConfig {DEFALUT = 1, CUSTOM};
enum class RenderingMethod {PATH_TRACING = 1, PHOTON_MAPPING};
string filename = "default.ppm";

Color path_tracing_radiance(const Ray& r, const Color& background, const Hittable& world, Hittable* lights, double imp_ratio, int depth) {
    Color color = Color(0, 0, 0);
    if (depth < 0) {
        return Color();
    }
    HitRecord hrec;
    if (!world.hit(r, 1e-4, infinity, hrec)) { // if not hit
        return background;
    }
    ScatterRecord srec;
    Color emitted = hrec.mat_ptr->emitted(r, hrec, hrec.u, hrec.v, hrec.p);

    if (!hrec.mat_ptr->scatter(r, hrec, srec)) {
        return emitted;
    }
    else if (srec.is_specular) {
        return emitted + srec.albedo * path_tracing_radiance(srec.scattered, background, world, lights, imp_ratio, depth - 1);
    }
    else {
        PDF* pdf_imp = new HittableCustomPDF(lights, hrec.p);
        MixturePDF* mix_pdf = new MixturePDF(pdf_imp, srec.p_pdf, imp_ratio);
        srec.scattered = Ray(hrec.p, mix_pdf->generate());
        double pdf_val = mix_pdf->value(srec.scattered.direction());
        color = emitted + srec.albedo * path_tracing_radiance(srec.scattered, background, world, lights, imp_ratio, depth - 1)
            * hrec.mat_ptr->scattering_pdf(r, hrec, srec.scattered) / pdf_val;

        delete mix_pdf; // must release memory
        delete pdf_imp;

        return color;
        // For lambertian surface, BRDF = albedo/PI. Value of 'albedo * scattering_pdf' here is equal to
        // the 'BRDF * cos(omega)' part in the rendering equation
    }
    return color;
}

Color direct_lighting_radiance(const Ray& r, const Color& background, const Hittable& world, Hittable* lights, int depth) {
    if (depth < 0) {
        return Color();
    }
    HitRecord hrec;
    if (!world.hit(r, 1e-4, infinity, hrec)) { // if not hit
        return background;
    }
    ScatterRecord srec;
    Color emitted = hrec.mat_ptr->emitted(r, hrec, hrec.u, hrec.v, hrec.p);

    if (!hrec.mat_ptr->scatter(r, hrec, srec)) {
        return emitted;
    }
    else if (srec.is_specular) {
        return Color();
    }
    else {
        PDF* pdf = new HittableCustomPDF(lights, hrec.p);
        srec.scattered = Ray(hrec.p, pdf->generate());
        double pdf_val = pdf->value(srec.scattered.direction());

        delete pdf; // release memory

        return emitted + srec.albedo * direct_lighting_radiance(srec.scattered, background, world, lights, depth - 1)
            * hrec.mat_ptr->scattering_pdf(r, hrec, srec.scattered) / pdf_val;
    }
}

Color direct_lighting_radiance(const Ray& r, const Color& background, const Hittable& world, Hittable* lights) {
    return direct_lighting_radiance(r, background, world, lights, 1);
}

Color global_photon_radiance(const Ray& r, const Color& background, const Hittable& world, int depth, PhotonMap* global_map) { // used for global photon map
    if (depth <= 0) {
        return Color();
    }
    HitRecord hrec;
    if (!world.hit(r, 1e-4, infinity, hrec)) { // if not hit
        return background;
    }
    ScatterRecord srec;
    Color emitted = hrec.mat_ptr->emitted(r, hrec, hrec.u, hrec.v, hrec.p);
    if (!hrec.mat_ptr->scatter(r, hrec, srec)) {
        return emitted;
    }
    else if (srec.is_specular) {
        return emitted + srec.albedo * global_photon_radiance(srec.scattered, background, world, depth - 1, global_map);
    }
    else {
        Vec3 ellipsoid_param = Vec3(400, 400, 400);
        return emitted + global_map->radiance_estimate(hrec.p, srec.albedo / PI, ellipsoid_param, hrec.normal);
    }
}

Color indirect_radiance(const Ray& r, const Color& background, const Hittable& world, PhotonMap* map) {
    HitRecord hrec;
    if (!world.hit(r, 1e-4, infinity, hrec)) { // if not hit
        return background;
    }
    ScatterRecord srec;
    if (!hrec.mat_ptr->scatter(r, hrec, srec)) {
        return Color();
    }
    else if (srec.is_specular) {
        return Color();
    }
    else {
        Vec3 ellipsoid_param = Vec3(400, 400, 400);
        return map->radiance_estimate(hrec.p, srec.albedo / PI, ellipsoid_param, hrec.normal);
    }
}

Color hybrid_complete_radiance(const Ray& r, const Color& background, const Hittable& world, Hittable* lights,
    PhotonMap* caustics_map, PhotonMap* indirect_map, int depth) {
    if (depth < 0) {
        return Color();
    }
    HitRecord hrec;
    if (!world.hit(r, 1e-4, infinity, hrec)) { // if not hit
        return background;
    }
    ScatterRecord srec;
    Color emitted = hrec.mat_ptr->emitted(r, hrec, hrec.u, hrec.v, hrec.p);

    if (!hrec.mat_ptr->scatter(r, hrec, srec)) {
        return emitted;
    }
    else if (srec.is_specular) {
        return emitted + srec.albedo * hybrid_complete_radiance(srec.scattered, background, world, lights, caustics_map, indirect_map, depth - 1);
    }
    else {
        // return direct_lighting_radiance + indirect_radiance + caustics_radiance
        Color direct;
        Color indirect;
        Color caustics;
        Vec3 ellipsoid_param = Vec3(400, 400, 400);
        PDF* pdf = new HittableCustomPDF(lights, hrec.p);
        srec.scattered = Ray(hrec.p, pdf->generate());
        double pdf_val = pdf->value(srec.scattered.direction());
        direct = emitted + srec.albedo * direct_lighting_radiance(srec.scattered, background, world, lights, 0)
            * hrec.mat_ptr->scattering_pdf(r, hrec, srec.scattered) / pdf_val;
        indirect = indirect_map->radiance_estimate(hrec.p, srec.albedo / PI, ellipsoid_param, hrec.normal);
        caustics = caustics_map->radiance_estimate(hrec.p, srec.albedo / PI, ellipsoid_param, hrec.normal);

        delete pdf; // release memory

        return direct + indirect + caustics;
    }
}

Color specular_radiance(const Ray& r, const Color& background, const Hittable& world, Hittable* lights, Hittable* spec_list,
    PhotonMap* caustics_map, PhotonMap* indirect_map, int depth) {
    if (depth < 0) {
        return Color();
    }
    HitRecord hrec;
    if (!spec_list->hit(r, 1e-4, infinity, hrec)) { // if not hit specular objects
        return background;
    }
    ScatterRecord srec;

    if (!hrec.mat_ptr->scatter(r, hrec, srec)) {
        return Color();
    }
    else if (srec.is_specular) {
        return srec.albedo * hybrid_complete_radiance(srec.scattered, background, world, lights, caustics_map, indirect_map, depth - 1);
    }
    else {
        return Color();
    }
}

Color photon_brightness(const Ray& r, const Hittable& world, PhotonMap* map) {
    HitRecord hrec;
    if (!world.hit(r, 1e-4, infinity, hrec)) { // if not hit
        return Color();
    }
    ScatterRecord srec;
    if (!hrec.mat_ptr->scatter(r, hrec, srec)) {
        return Color();
    }
    else if (srec.is_specular) {
        return Color();
    }
    Vec3 ellipsoid_param = Vec3(0.5, 0.5, 0.5);
    return map->show_photon(hrec.p, ellipsoid_param, hrec.normal);
}

HittableList cornell_box() {
    HittableList objects;
    Material* white = new Diffuse(Color(0.75, 0.75, 0.75));
    Material* red = new Diffuse(Color(0.65, 0.05, 0.05));
    Material* green = new Diffuse(Color(0.12, 0.45, 0.15));
    Material* light = new DiffuseLight(Color(16, 16, 16));
    Material* white_mirror = new Specular(Color(0.75, 0.75, 0.75));
    Material* glass = new Dielectric(new SolidColor(1, 1, 1), 1.6);

    // add walls
    objects.add(new YZ_Rect(0, 555, 0, 555, 555, red));
    objects.add(new YZ_Rect(0, 555, 0, 555, 0, green));
    objects.add(new XZ_Rect(213, 343, 227, 332, 554, light)); // light
    objects.add(new XZ_Rect(0, 555, 0, 555, 0, white)); // floor
    objects.add(new XZ_Rect(0, 555, 0, 555, 555, white)); // ceil
    objects.add(new XY_Rect(0, 555, 0, 555, 555, white)); // front

    // add boxes
    Hittable* box1 = new Box(Point3(0, 0, 0), Point3(165, 330, 165), white); // high box
    Hittable* box2 = new Box(Point3(0, 0, 0), Point3(165, 165, 165), white); // low box
    box1 = new Rotate_Y(box1, 15);
    box1 = new Translate(box1, Vec3(265, 0, 295));
    box2 = new Rotate_Y(box2, -18);
    box2 = new Translate(box2, Vec3(130, 0, 65));
    objects.add(box1);
    objects.add(box2);
    return objects;
}

HittableList cornell_box_spheres(HittableList* imp_list, HittableList* spec_list, HittableList* light_list) {
    HittableList objects;
    Material* white = new Diffuse(Color(0.75, 0.75, 0.75));
    Material* red = new Diffuse(Color(0.65, 0.05, 0.05));
    Material* green = new Diffuse(Color(0.12, 0.45, 0.15));
    Material* light = new DiffuseLight(Color(20, 20, 20));
    Material* white_mirror = new Specular(Color(0.8, 0.8, 0.8));
    Material* glass = new Dielectric(new SolidColor(1, 1, 1), 1.6);
    Hittable* light_o = new XZ_Rect(213, 343, 227, 332, 554, light);

    // add walls
    objects.add(new YZ_Rect(0, 555, 0, 555, 555, red));
    objects.add(new YZ_Rect(0, 555, 0, 555, 0, green));
    objects.add(light_o); // light
    objects.add(new XZ_Rect(0, 555, 0, 555, 0, white)); // floor
    objects.add(new XZ_Rect(0, 555, 0, 555, 555, white)); // ceil
    objects.add(new XY_Rect(0, 555, 0, 555, 555, white)); // front

    Hittable* sph1 = new Sphere(Point3(415, 100, 360), 100, white_mirror); // mirror sphere
    Hittable* sph2 = new Sphere(Point3(140, 100, 190), 100, glass); // glass sphere

    objects.add(sph1);
    objects.add(sph2);

    light_list->add(light_o);

    imp_list->add(light_o);
    imp_list->add(sph2);

    spec_list->add(sph1);
    spec_list->add(sph2);

    return objects;
}

HittableList earth() {
    // ImageTexture* earth_texture = new ImageTexture("earthmap.jpg"); // error in texture.h
    SolidColor* earth_texture = new SolidColor(Color(1, 1, 1));
    Diffuse* earth_surface = new Diffuse(earth_texture);
    Sphere* globe = new Sphere(Point3(0, 0, 0), 2, earth_surface);
    return HittableList(globe);
}

HittableList simple_scene() {
    HittableList objects;
    Material* glass = new Diffuse(new SolidColor(0.75, 0.75, 0.75));
    Material* blue = new Diffuse(Color(0.2, 0.2, 0.9));
    Material* green = new Diffuse(Color(0.12, 0.65, 0.15));
    objects.add(new Sphere(Vec3(1, 1, 0), 1, glass));
    objects.add(new Sphere(Vec3(-1, 1, 0), 1, blue));
    objects.add(new Sphere(Vec3(0, -100, 0), 100, green)); // ground
    return objects;
}

HittableList test_scene() {
    HittableList objects;
    Material* light = new DiffuseLight(Color(20, 20, 20));
    Material* white = new Diffuse(Color(0.75, 0.75, 0.75));
    objects.add(new XZ_Rect(-1, 1, -1, 1, 3, light));
    objects.add(new XZ_Rect(-2, 2, -2, 2, 0, white));
    objects.add(new Sphere(Vec3(1, 1, 1), 1, white));
    return objects;
}

HittableList empty_cornell_box() {
    HittableList objects;
    Material* white = new Diffuse(Color(0.75, 0.75, 0.75));
    Material* red = new Diffuse(Color(0.65, 0.05, 0.05));
    Material* green = new Diffuse(Color(0.12, 0.45, 0.15));
    Material* light = new DiffuseLight(Color(20, 20, 20));
    Material* white_mirror = new Specular(Color(0.8, 0.8, 0.8));
    Material* glass = new Dielectric(new SolidColor(1, 1, 1), 1.5);

    //// add walls
    objects.add(new YZ_Rect(0, 555, 0, 555, 555, red));
    objects.add(new YZ_Rect(0, 555, 0, 555, 0, green));
    objects.add(new XZ_Rect(213, 343, 227, 332, 554, light)); // light
    objects.add(new XZ_Rect(0, 555, 0, 555, 0, white)); // floor
    objects.add(new XZ_Rect(0, 555, 0, 555, 555, white)); // ceil
    objects.add(new XY_Rect(0, 555, 0, 555, 555, white)); // front
    return objects;
}

void change_file_name(string& s, RenderingMethod method) {
    size_t length = s.length();
    string s1 = s.substr(0, length - 4);
    string s2 = s.substr(length - 4);
    string mid = method == RenderingMethod::PATH_TRACING ? "_pt" : "_pm";
    s = s1 + mid + s2;
}

int main() {
    // image
    // customizable configuration
    int scene_case = 0;
    int max_depth = 16;
    int samples_per_pixel = 32;
    double importance_ratio = 0;
    // internal configuration
    int subpixel_level = 2;
    int subpixel_level2 = subpixel_level * subpixel_level;
    double subpixel_length = 1.0 / subpixel_level;
    double aspect_ratio = 16.0/9;
    int image_width = 400;
    int param_config = static_cast<std::underlying_type<ParamConfig>::type>(ParamConfig::CUSTOM);
    int method_config_in = static_cast<std::underlying_type<RenderingMethod>::type>(RenderingMethod::PATH_TRACING);
    RenderingMethod method_config = RenderingMethod::PATH_TRACING;
    Color global_photon_energy = Color(21, 21, 21);
    Color caustics_photon_energy = Color(4.5, 4.5, 4.5);

    // world
    HittableList world;
    HittableList* imp_list = new HittableList(); // do not delete this, otherwise it will cause redundantly release
    HittableList* spec_list = new HittableList(); // do not delete this, otherwise it will cause redundantly release
    HittableList* light_list = new HittableList(); // do not delete this, otherwise it will cause redundantly release
    Color background(0, 0, 0);

    // camera
    Point3 lookfrom;
    Point3 lookat;
    double vfov = 40.0;
    double aperture = 0.0;
    Vec3 vup = Vec3(0, 1, 0);
    double dist_to_focus = 10.0;

    // others
    clock_t start_time;
    clock_t end_time;

    // get input
    std::cout << "please input scene number: "
        << endl << "1. simple_scene"
        << endl << "2. earth"
        << endl << "3. cornell_box"
        << endl << "4. cornell_box_spheres"
        << endl;
    cin >> scene_case;
    std::cout << "please select rendering method\n1.path tracing\n2.photon mapping\n";
    cin >> method_config_in;
    std::cout << "please select configuration\n1.default\n2.custom\n";
    cin >> param_config;
    if (param_config == static_cast<std::underlying_type<ParamConfig>::type>(ParamConfig::CUSTOM)) {
        std::cout << "please input max reflection depth: ";
        cin >> max_depth;
        std::cout << "please input samples per pixel: ";
        cin >> samples_per_pixel;
        std::cout << "please input importance sampling ratio (between[0, 1]) : ";
        cin >> importance_ratio;
    }
    method_config = static_cast<RenderingMethod>(method_config_in);
    max_depth = static_cast<int>(clamp(max_depth, 1, 32));
    samples_per_pixel = static_cast<int>(clamp(samples_per_pixel, 1, 1000));
    importance_ratio = clamp(importance_ratio, 0, 1);
    // scene
    switch (scene_case)
    {
    case 1:
        world = simple_scene();
        background = Color(0.5, 0.7, 1.0);
        lookfrom = Point3(3, 3, 1);
        lookat = Point3(-1, 1, 0);
        vfov = 60.0;
        break;
    case 2:
        world = earth();
        background = Color(.75, .75, .75);
        lookfrom = Point3(-10 / sqrt(2), 2, -16 / sqrt(2));
        lookat = Point3(0, 0, 0);
        vfov = 20.0;
        break;
    case 3:
        filename = "cornell_box.ppm";
        world = cornell_box();
        aspect_ratio = 1.0;
        image_width = 600;
        background = Color();
        lookfrom = Point3(278, 278, -800);
        lookat = Point3(278, 278, 0);
        vfov = 40.0;
        imp_list->add(new XZ_Rect(213, 343, 227, 332, 554, nullptr)); // light
        break;
    case 4:
        filename = "cornell_spheres.ppm";
        world = cornell_box_spheres(imp_list, spec_list, light_list);
        aspect_ratio = 1.0;
        image_width = 600;
        background = Color();
        lookfrom = Point3(278, 278, -800);
        lookat = Point3(278, 278, 0);
        vfov = 40.0;
        break;
    default:
    case 0:
        filename = "test_scene.ppm";
        world = test_scene();
        background = Color(.25, .25, .25);
        lookfrom = Point3(0, 1, 2);
        lookat = Point3(0, 1, 0);
        vfov = 40.0;
        imp_list->add(new XZ_Rect(-1, 1, -1, 1, 3, nullptr));
        break;
    };
    int image_height = static_cast<int>(image_width / aspect_ratio);

    // camera
    Camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus);
    HittableList empty_box = empty_cornell_box(); // for rendering photon map

    // render
    change_file_name(filename, method_config);
    if (method_config_in == static_cast<std::underlying_type<RenderingMethod>::type>(RenderingMethod::PATH_TRACING)) {
        ofstream outfile;
        outfile.open(filename);
        outfile << "P3\n" << image_width << ' ' << image_height << "\n255\n";
        start_time = clock();
        for (int j = image_height - 1; j >= 0; --j) {
            std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
            for (int i = 0; i < image_width; ++i) {
                Color pixel_color(0, 0, 0);
                for (int s = 0; s < samples_per_pixel; ++s) { // stratified sampling
                    int subpixel_u = s % subpixel_level; // column
                    int subpixel_v = (s % subpixel_level2) / subpixel_level; // row
                    double u = (i + subpixel_u * subpixel_length + random_double() / subpixel_level) / (image_width - 1.0);
                    double v = (j + subpixel_v * subpixel_length + random_double() / subpixel_level) / (image_height - 1.0);
                    Ray r = cam.get_ray(u, v);
                    pixel_color += path_tracing_radiance(r, background, world, light_list, importance_ratio, max_depth); // imp_list ignored
                }
                write_color(outfile, pixel_color, samples_per_pixel);
            }
        }
        end_time = clock();
        outfile.close();
        std::cout << "rendering time is: " << end_time - start_time << endl;
    }
    else if (method_config_in == static_cast<std::underlying_type<RenderingMethod>::type>(RenderingMethod::PHOTON_MAPPING)) {
        XZ_Rect* photon_map_light = new XZ_Rect(213, 343, 227, 332, 554, nullptr);
        /*PhotonRatioSampler sampler = PhotonRatioSampler(100000, world, photon_map_light);
        double photon_ratio = sampler.get_ratio();*/
        // PhotonMap* map = new PhotonMap(20000, world, photon_map_light);
        IndirectPhotonMap* indirect_map = new IndirectPhotonMap(40000, 20.0, global_photon_energy, world, photon_map_light);
        CausticsPhotonMap* caustics_map = new CausticsPhotonMap(20000, 20.0, caustics_photon_energy , world, photon_map_light, spec_list);
        samples_per_pixel = 1; //+

        // render
        ofstream outfile;

        /*samples_per_pixel = 16;
        outfile.open("total_scene.ppm");
        start_time = clock();
        outfile << "P3\n" << image_width << ' ' << image_height << "\n255\n";
        for (int j = image_height - 1; j >= 0; --j) {
            std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
            for (int i = 0; i < image_width; ++i) {
                Color pixel_color(0, 0, 0);
                for (int s = 0; s < samples_per_pixel; ++s) { // stratified sampling
                    int subpixel_u = s % subpixel_level; // column
                    int subpixel_v = (s % subpixel_level2) / subpixel_level; // row
                    double u = (i + subpixel_u * subpixel_length + random_double() / subpixel_level) / (image_width - 1.0);
                    double v = (j + subpixel_v * subpixel_length + random_double() / subpixel_level) / (image_height - 1.0);
                    Ray r = cam.get_ray(u, v);
                    pixel_color += hybrid_complete_radiance(r, background, world, light_list, caustics_map, indirect_map, max_depth);
                }
                write_color(outfile, pixel_color, samples_per_pixel);
            }
        }
        end_time = clock();
        outfile.close();
        std::cout << "\nrendering time is: " << end_time - start_time << endl;*/

        // direct lighting
        /*samples_per_pixel = 64;
        outfile.open("direct_lighting.ppm");
        start_time = clock();
        outfile << "P3\n" << image_width << ' ' << image_height << "\n255\n";
        for (int j = image_height - 1; j >= 0; --j) {
            std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
            for (int i = 0; i < image_width; ++i) {
                Color pixel_color(0, 0, 0);
                for (int s = 0; s < samples_per_pixel; ++s) { // stratified sampling
                    int subpixel_u = s % subpixel_level; // column
                    int subpixel_v = (s % subpixel_level2) / subpixel_level; // row
                    double u = (i + subpixel_u * subpixel_length + random_double() / subpixel_level) / (image_width - 1.0);
                    double v = (j + subpixel_v * subpixel_length + random_double() / subpixel_level) / (image_height - 1.0);
                    Ray r = cam.get_ray(u, v);
                    pixel_color += direct_lighting_radiance(r, background, world, light_list, 1);
                }
                write_color(outfile, pixel_color, samples_per_pixel);
            }
        }
        end_time = clock();
        outfile.close();
        std::cout << "\nrendering time is: " << end_time - start_time << endl;*/

        // caustics
        samples_per_pixel = 4;
        outfile.open("caustics.ppm");
        start_time = clock();
        outfile << "P3\n" << image_width << ' ' << image_height << "\n255\n";
        for (int j = image_height - 1; j >= 0; --j) {
            std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
            for (int i = 0; i < image_width; ++i) {
                Color pixel_color(0, 0, 0);
                for (int s = 0; s < samples_per_pixel; ++s) { // stratified sampling
                    int subpixel_u = s % subpixel_level; // column
                    int subpixel_v = (s % subpixel_level2) / subpixel_level; // row
                    double u = (i + subpixel_u * subpixel_length + random_double() / subpixel_level) / (image_width - 1.0);
                    double v = (j + subpixel_v * subpixel_length + random_double() / subpixel_level) / (image_height - 1.0);
                    Ray r = cam.get_ray(u, v);
                    pixel_color += indirect_radiance(r, background, world, caustics_map);
                }
                write_color(outfile, pixel_color, samples_per_pixel);
            }
        }
        end_time = clock();
        outfile.close();
        std::cout << "\nrendering time is: " << end_time - start_time << endl;

        // indirect lighting
        samples_per_pixel = 4;
        outfile.open("indirect_lighting.ppm");
        start_time = clock();
        outfile << "P3\n" << image_width << ' ' << image_height << "\n255\n";
        for (int j = image_height - 1; j >= 0; --j) {
            std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
            for (int i = 0; i < image_width; ++i) {
                Color pixel_color(0, 0, 0);
                for (int s = 0; s < samples_per_pixel; ++s) { // stratified sampling
                    int subpixel_u = s % subpixel_level; // column
                    int subpixel_v = (s % subpixel_level2) / subpixel_level; // row
                    double u = (i + subpixel_u * subpixel_length + random_double() / subpixel_level) / (image_width - 1.0);
                    double v = (j + subpixel_v * subpixel_length + random_double() / subpixel_level) / (image_height - 1.0);
                    Ray r = cam.get_ray(u, v);
                    pixel_color += indirect_radiance(r, background, world, indirect_map);
                }
                write_color(outfile, pixel_color, samples_per_pixel);
            }
        }
        end_time = clock();
        outfile.close();
        std::cout << "\nrendering time is: " << end_time - start_time << endl;

        // specular
        /*samples_per_pixel = 32;
        outfile.open("specular.ppm");
        start_time = clock();
        outfile << "P3\n" << image_width << ' ' << image_height << "\n255\n";
        for (int j = image_height - 1; j >= 0; --j) {
        // for (int j = 150; j >= 0; --j) {
            //- std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
            for (int i = 0; i < image_width; ++i) {
                std::cerr << "\rScanlines remaining: " << j << ' ' << i << ' ' << std::flush; //+
                Color pixel_color(0, 0, 0);
                for (int s = 0; s < samples_per_pixel; ++s) { // stratified sampling
                    int subpixel_u = s % subpixel_level; // column
                    int subpixel_v = (s % subpixel_level2) / subpixel_level; // row
                    double u = (i + subpixel_u * subpixel_length + random_double() / subpixel_level) / (image_width - 1.0);
                    double v = (j + subpixel_v * subpixel_length + random_double() / subpixel_level) / (image_height - 1.0);
                    Ray r = cam.get_ray(u, v);
                    pixel_color += specular_radiance(r, background, world, light_list, spec_list, caustics_map, indirect_map, max_depth);
                }
                write_color(outfile, pixel_color, samples_per_pixel);
            }
        }
        end_time = clock();
        outfile.close();
        std::cout << "\nrendering time is: " << end_time - start_time << endl;*/

        // render photon map
        /*outfile.open("photon_map.ppm");
        image_width = image_height = 1200;
        start_time = clock();
        outfile << "P3\n" << image_width << ' ' << image_height << "\n255\n";
        for (int j = image_height - 1; j >= 0; --j) {
            std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
            for (int i = 0; i < image_width; ++i) {
                Color pixel_color(0, 0, 0);
                for (int s = 0; s < 1; ++s) { // stratified sampling
                    double u = (i + 0.5) / (image_width - 1.0);
                    double v = (j + 0.5) / (image_height - 1.0);
                    Ray r = cam.get_ray(u, v);
                    pixel_color += photon_brightness(r, empty_box, indirect_map);
                }
                write_color(outfile, pixel_color, 1);
            }
        }
        end_time = clock();
        outfile.close();
        std::cout << "\nrendering time is: " << end_time - start_time << endl;
        std::cout << "max photons is: " << indirect_map->maxphotons << endl; // debug code*/
    }
}