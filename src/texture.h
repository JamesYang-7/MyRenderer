#ifndef TEXTURE_H
#define TEXTURE_H

#include "rtweekend.h"
//#include "rtw_stb_image.h" // multi-defination error because there exists function defination in hearder file external/rtw_stb_image.h
// need to put function definations in .cpp file to eliminate errors

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

//class ImageTexture : public Texture { 
//public:
//	const static int bytes_per_pixel = 3;
//
//	ImageTexture() : data(nullptr), width(0), height(0), bytes_per_scanline(0) {}
//
//	ImageTexture(const char* filename) {
//		int components_per_pixel = bytes_per_pixel;
//		data = stbi_load(filename, &width, &height, &components_per_pixel, components_per_pixel);
//		if (!data) {
//			std::cerr << "ERROR: Could not load texture image file '" << filename << "'.\n";
//			width = height = 0;
//		}
//		bytes_per_scanline = bytes_per_pixel * width;
//	}
//
//	~ImageTexture() {
//		delete data;
//	}
//
//	virtual Color value(double u, double v, const Vec3& p) const override;
//private:
//	unsigned char* data;
//	int width;
//	int height;
//	int bytes_per_scanline;
//};
#endif // !TEXTURE_H

