#include "stdafx.h"
#include "Image.h"
#include "tinyexr.h"


Image::Image(int _width, int _height)
	:width(_width), height(_height)
{
	for (int i = 0; i < 3; i++)
	{
		images[i].resize(width * height);
		std::fill_n(images[i].begin(), images[i].size(), 0.f);
	}
}


Image::~Image()
{
}

void Image::write(const char* fileName)
{
	EXRImage image;
	InitEXRImage(&image);

	image.num_channels = 3;

	// Must be BGR(A) order, since most of EXR viewers expect this channel order.
	const char* channel_names[] = { "B", "G", "R" }; // "B", "G", "R", "A" for RGBA image

	float* image_ptr[3];
	image_ptr[0] = &(images[2].at(0)); // B
	image_ptr[1] = &(images[1].at(0)); // G
	image_ptr[2] = &(images[0].at(0)); // R

	image.channel_names = channel_names;
	image.images = (unsigned char**)image_ptr;
	image.width = width;
	image.height = height;

	image.pixel_types = (int *)malloc(sizeof(int) * image.num_channels);
	image.requested_pixel_types = (int *)malloc(sizeof(int) * image.num_channels);
	for (int i = 0; i < image.num_channels; i++) {
		image.pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT; // pixel type of input image
		image.requested_pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT; // pixel type of output image to be stored in .EXR
	}

	const char* err;
	int ret = SaveMultiChannelEXRToFile(&image, fileName, &err);
	if (ret != 0) {
		throw std::string("Save EXR err: ") + err;
	}
	printf("Saved EXR file. [ %s ] \n", fileName);
}

void Image::write(int row, int col, const Array3f& val)
{
	images[0][row * width + col] = val.x();
	images[1][row * width + col] = val.y();
	images[2][row * width + col] = val.z();
}
