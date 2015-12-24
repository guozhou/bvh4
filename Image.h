#pragma once
using namespace Eigen;

class Image
{
public:
	Image(int width, int height);
	~Image();

	const int width, height;

	void write(int row, int col, const Array3f& val);
	void write(const char* fileName);
private:
	std::vector<float> images[3]; // RGB
};

