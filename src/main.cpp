#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <random>
#include <ctime>
#include "rtweekend.h"
#include "color.h"
#include "hittable_list.h"
#include "sphere.h"
#include "camera.h"
#include "material.h"
#include "aarect.h"
#include "box.h"
#include "pdf.h"
#include "photon_map.h"
#include "bvh.h"

using namespace std;
int max_photons = 0;
enum class ParamConfig {DEFALUT = 1, CUSTOM};
enum class RenderingMethodConfig {PATH_TRACING = 1, PHOTON_MAPPING};
string filename = "default.ppm";

Color path_tracing_radiance(const Ray& r, const Color& background, const Hittable& world, Hittable* lights, double imp_ratio, int depth) {
    if (depth < 0) {
        return Vec3();
    }
    HitRecord hrec;
    if (!world.hit(r, 1e-4, infinity, hrec)) { // if not hit
        return background;
    }
    ScatterRecord srec;
    double pdf_val = 0;
    Color emitted = hrec.mat_ptr->emitted(r, hrec, hrec.u, hrec.v, hrec.p);

    if (!hrec.mat_ptr->scatter(r, hrec, srec)) {
        return emitted;
    }
    else if (srec.is_specular) {
        return emitted + srec.albedo * path_tracing_radiance(srec.scattered, background, world, lights, imp_ratio, depth - 1);
    }
    else {
        MixturePDF* mix_pdf = new MixturePDF(new HittableCustomPDF(lights, hrec.p), srec.p_pdf, imp_ratio);
        srec.scattered = Ray(hrec.p, mix_pdf->generate());
        pdf_val = mix_pdf->value(srec.scattered.direction());
        delete mix_pdf;
        return emitted + srec.albedo * path_tracing_radiance(srec.scattered, background, world, lights, imp_ratio, depth - 1)
            * hrec.mat_ptr->scattering_pdf(r, hrec, srec.scattered) / pdf_val;
        // For lambertian surface, BRDF = albedo/PI. Value of 'albedo * scattering_pdf' here is equal to
        // the 'BRDF * cos(omega)' part in the rendering equation
    }
}

Color direct_lighting_radiance(const Ray& r, const Color& background, const Hittable& world, Hittable* lights) {
    return path_tracing_radiance(r, background, world, lights, 1, 1);
}

Color photon_radiance(const Ray& r, const Color& background, const Hittable& world, Hittable* lights,
    double imp_ratio, int depth, PhotonMap* global_map) {
    if (depth <= 0) {
        return Vec3();
    }
    HitRecord hrec;
    if (!world.hit(r, 1e-4, infinity, hrec)) { // if not hit
        return Color();
    }
    ScatterRecord srec;
    Color emitted = hrec.mat_ptr->emitted(r, hrec, hrec.u, hrec.v, hrec.p);
    if (!hrec.mat_ptr->scatter(r, hrec, srec)) {
        return emitted;
    }
    else if (srec.is_specular) {
        return emitted + srec.albedo * photon_radiance(srec.scattered, background, world, lights, imp_ratio, depth - 1, global_map);
    }
    else {
        Vec3 ellipsoid_param = Vec3(400, 400, 400);
        return emitted + global_map->radiance_estimate(hrec.p, srec.albedo / PI, ellipsoid_param, hrec.normal);
    }
}

Color photon_brightness(const Ray& r, const Hittable& world, PhotonMap* map) {
    HitRecord hrec;
    if (!world.hit(r, 1e-4, infinity, hrec)) { // if not hit
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
    Material* glass = new Dielectric(new SolidColor(1, 1, 1), 1.5);

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

HittableList cornell_box_spheres() {
    HittableList objects;
    Material* white = new Diffuse(Color(0.75, 0.75, 0.75));
    Material* red = new Diffuse(Color(0.65, 0.05, 0.05));
    Material* green = new Diffuse(Color(0.12, 0.45, 0.15));
    Material* light = new DiffuseLight(Color(20, 20, 20));
    Material* white_mirror = new Specular(Color(0.8, 0.8, 0.8));
    Material* glass = new Dielectric(new SolidColor(1, 1, 1), 1.5);

    // add walls
    objects.add(new YZ_Rect(0, 555, 0, 555, 555, red));
    objects.add(new YZ_Rect(0, 555, 0, 555, 0, green));
    objects.add(new XZ_Rect(213, 343, 227, 332, 554, light)); // light
    objects.add(new XZ_Rect(0, 555, 0, 555, 0, white)); // floor
    objects.add(new XZ_Rect(0, 555, 0, 555, 555, white)); // ceil
    objects.add(new XY_Rect(0, 555, 0, 555, 555, white)); // front

    Hittable* sph1 = new Sphere(Point3(415, 100, 360), 100, white_mirror);
    Hittable* sph2 = new Sphere(Point3(140, 100, 190), 100, glass); // glass sphere

    objects.add(sph1);
    objects.add(sph2);
    return objects;
}

HittableList earth() {
    ImageTexture* earth_texture = new ImageTexture("earthmap.jpg");
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
    int method_config = static_cast<std::underlying_type<RenderingMethodConfig>::type>(RenderingMethodConfig::PATH_TRACING);

    // world
    HittableList world;
    HittableList* lights = new HittableList();
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
    cout << "please input scene number: "
        << endl << "1. simple_scene"
        << endl << "2. earth"
        << endl << "3. cornell_box"
        << endl << "4. cornell_box_spheres"
        << endl;
    cin >> scene_case;
    cout << "please select rendering method\n1.path tracing\n2.photon mapping\n";
    cin >> method_config;
    cout << "please select configuration\n1.default\n2.custom\n";
    cin >> param_config;
    if (param_config == static_cast<std::underlying_type<ParamConfig>::type>(ParamConfig::CUSTOM)) {
        cout << "please input max reflection depth: ";
        cin >> max_depth;
        cout << "please input samples per pixel: ";
        cin >> samples_per_pixel;
        cout << "please input importance sampling ratio (between[0, 1]) : ";
        cin >> importance_ratio;
    }
   
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
        lights->add(new XZ_Rect(213, 343, 227, 332, 554, nullptr)); // light
        break;
    case 4:
        filename = "cornell_box_spheres.ppm";
        world = cornell_box_spheres();
        aspect_ratio = 1.0;
        image_width = 600;
        background = Color();
        lookfrom = Point3(278, 278, -800);
        lookat = Point3(278, 278, 0);
        vfov = 40.0;
        lights->add(new XZ_Rect(213, 343, 227, 332, 554, nullptr)); // light
        lights->add(new Sphere(Point3(190, 90, 190), 90, nullptr)); // glass ball
        break;
    default:
    case 0:
        filename = "test_scene.ppm";
        world = test_scene();
        background = Color(.25, .25, .25);
        lookfrom = Point3(0, 1, 2);
        lookat = Point3(0, 1, 0);
        vfov = 40.0;
        lights->add(new XZ_Rect(-1, 1, -1, 1, 3, nullptr));
        break;
    };
    int image_height = static_cast<int>(image_width / aspect_ratio);

    // camera
    Camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus);
    HittableList empty_box = empty_cornell_box();

    // render
    if (method_config == static_cast<std::underlying_type<RenderingMethodConfig>::type>(RenderingMethodConfig::PATH_TRACING)) {
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
                    pixel_color += path_tracing_radiance(r, background, world, lights, importance_ratio, max_depth);
                }
                write_color(outfile, pixel_color, samples_per_pixel);
            }
        }
        end_time = clock();
        outfile.close();
        cout << "rendering time is: " << end_time - start_time << endl;
    }
    else if (method_config == static_cast<std::underlying_type<RenderingMethodConfig>::type>(RenderingMethodConfig::PHOTON_MAPPING)) {
        XZ_Rect* photon_map_light = new XZ_Rect(213, 343, 227, 332, 554, nullptr);
        PhotonMap* map = new PhotonMap(200000, world, photon_map_light);
        std::cout << "Photon map constructed\n";
        samples_per_pixel = 4; //+

        // render
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
                    pixel_color += photon_radiance(r, background, world, lights, importance_ratio, max_depth, map);
                }
                write_color(outfile, pixel_color, samples_per_pixel);
            }
        }
        end_time = clock();
        outfile.close();
        cout << "\nrendering time is: " << end_time - start_time << endl;

        // render photon map
        outfile.open("photon_map.ppm");
        image_width = image_height = 1200;
        outfile << "P3\n" << image_width << ' ' << image_height << "\n255\n";
        for (int j = image_height - 1; j >= 0; --j) {
            std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
            for (int i = 0; i < image_width; ++i) {
                Color pixel_color(0, 0, 0);
                for (int s = 0; s < 1; ++s) { // stratified sampling
                    double u = (i + 0.5) / (image_width - 1.0);
                    double v = (j + 0.5) / (image_height - 1.0);
                    Ray r = cam.get_ray(u, v);
                    pixel_color += photon_brightness(r, world, map);
                }
                write_color(outfile, pixel_color, 1);
            }
        }
        outfile.close();
        // cout << "max photons is: " << map->maxphotons << endl; // debug code
    }
}