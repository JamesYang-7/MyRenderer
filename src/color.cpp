#include "color.h"

void write_color(std::ostream& out, Color pixel_color, int samples_per_pixel) {
    double r = pixel_color.x();
    double g = pixel_color.y();
    double b = pixel_color.z();

    double scale = 1.0 / samples_per_pixel;
    r = gamma(scale * r);
    g = gamma(scale * g);
    b = gamma(scale * b);

    // Write the translated [0,255] value of each color component.
    out << static_cast<int>(255.9 * clamp(r, 0, 1)) << ' '
        << static_cast<int>(255.9 * clamp(g, 0, 1)) << ' '
        << static_cast<int>(255.9 * clamp(b, 0, 1)) << '\n';
}