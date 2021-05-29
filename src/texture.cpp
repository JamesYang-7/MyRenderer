#include "texture.h"

//Color ImageTexture::value(double u, double v, const Vec3& p) const {
//	if (data == nullptr) { return Color(0, 1, 1); }
//
//	u = clamp(u);
//	v = 1.0 - clamp(v);
//
//	int i = static_cast<int>(u * width);
//	int j = static_cast<int>(v * height);
//
//	if (i >= width) i = width - 1;
//	if (j >= height) j = height - 1;
//
//	const double color_scale = 1.0 / 255.0;
//	unsigned char* pixel = data + j * bytes_per_scanline + i * bytes_per_pixel;
//	return Color(color_scale * pixel[0], color_scale * pixel[1], color_scale * pixel[2]);
//}