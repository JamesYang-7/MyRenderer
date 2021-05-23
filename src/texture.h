#ifndef TEXTURE_H
#define TEXTURE_H

#include "rtweekend.h"
#include "rtw_stb_image.h"

class Texture {
public:
	virtual Color value(double u, double v, const Point3& p) const = 0;
};

class SolidColor : public Texture {
public:
	SolidColor() {}
	SolidColor(Color c) : color_value(c) {}
	SolidColor(double red, double green, double blue) : SolidColor(Color(red, green, blue)) {}
	virtual Color value(double u, double v, const Point3& p) const override {
		return color_value;
	}
private:
	Color color_value;
};

class ImageTexture : public Texture {
public:
	const static int bytes_per_pixel = 3;

	ImageTexture() : data(nullptr), width(0), height(0), bytes_per_scanline(0) {}

	ImageTexture(const char* filename) {
		int components_per_pixel = bytes_per_pixel;
		data = stbi_load(filename, &width, &height, &components_per_pixel, components_per_pixel);
		if (!data) {
			std::cerr << "ERROR: Could not load texture image file '" << filename << "'.\n";
			width = height = 0;
		}
		bytes_per_scanline = bytes_per_pixel * width;
	}

	~ImageTexture() {
		delete data;
	}

	virtual Color value(double u, double v, const Vec3& p) const override {
		if (data == nullptr) { return Color(0, 1, 1); }

		u = clamp(u);
		v = 1.0 - clamp(v);

		int i = static_cast<int>(u * width);
		int j = static_cast<int>(v * height);

		if (i >= width) i = width - 1;
		if (j >= height) j = height - 1;

		const double color_scale = 1.0 / 255.0;
		unsigned char* pixel = data + j * bytes_per_scanline + i * bytes_per_pixel;
		return Color(color_scale * pixel[0], color_scale * pixel[1], color_scale * pixel[2]);
	}
private:
	unsigned char* data;
	int width;
	int height;
	int bytes_per_scanline;
};
#endif // !TEXTURE_H

