#pragma once
using namespace Eigen;

struct Ray;

class Camera
{
public:
	Camera(int width, int height);
	~Camera();

	void lookAt(const Vector3f& eye, const Vector3f& center, const Vector3f& up);
	void perspective(float yFov = 30.f, float zNear = 1e-1f, float zFar = 1e4f);
	void test() const;
	Array3f sampleRay(Ray &ray, const Array2f& sample) const;
private:
	const float rwidth, rheight, aspect;
	Vector3f eye;
	float zNear, zFar;

	Matrix4f view;
	Matrix4f proj;
	Matrix4f rast;
	Matrix4f invView;
	Matrix4f invProj;
	Matrix4f invRast;
	Matrix4f invProjRast;
};

