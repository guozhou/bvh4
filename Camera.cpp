#include "stdafx.h"
#include "Camera.h"
#include "Ray.h"

Camera::Camera(int width, int height)
	:rwidth(1.f / width), rheight(1.f / height), aspect(float(width) / height)
{
	// assumes a lower - left origin for window coordinates and assumes pixel centers 
	// are located at half - pixel centers
	rast <<
		0.5f * width, 0, 0, 0.5f * width,
		0, 0.5f * height, 0, 0.5f * height,
		0, 0, 0.5f, 0.5f,
		0, 0, 0, 1;

	invRast << 
		2 * rwidth, 0, 0, -1, 
		0, 2 * rheight, 0, -1, 
		0, 0, 2, -1, 
		0, 0, 0, 1;

	std::cout << invRast * rast << std::endl;
}

Camera::~Camera()
{
}

void Camera::lookAt(const Vector3f& _eye, const Vector3f& center, const Vector3f& up)
{
	eye = _eye;

	const Vector3f f = (center - eye).normalized(), s = f.cross(up).normalized(), u = s.cross(f);
	view.setIdentity();
	view.block<1, 3>(0, 0) = s;
	view.block<1, 3>(1, 0) = u;
	view.block<1, 3>(2, 0) = -f;
	view.block<3, 1>(0, 3) = Vector3f(-s.dot(eye), -u.dot(eye), f.dot(eye));

	invView.setIdentity();
	invView.block<3, 3>(0, 0) = view.block<3, 3>(0, 0).transpose();
	invView.block<3, 1>(0, 3) = -Vector3f(view.col(0).dot(view.col(3)),
		view.col(1).dot(view.col(3)), view.col(2).dot(view.col(3)));

	std::cout << invView * view << std::endl;
}

void Camera::perspective(float yFov, float _zNear, float _zFar)
{
	zNear = _zNear;
	zFar = _zFar;

	const float yScale = std::tanf(float(M_PI) / 2 - yFov / 2 * float(M_PI) / 180.f),
		xScale = yScale / aspect;
	proj.setZero();
	proj(0, 0) = xScale;
	proj(1, 1) = yScale;
	proj(2, 2) = -(zFar + zNear) / (zFar - zNear);
	proj(3, 2) = -1;
	proj(2, 3) = -2 * zFar * zNear / (zFar - zNear);

	invProj.setZero();
	invProj(0, 0) = 1.f / proj(0, 0);
	invProj(1, 1) = 1.f / proj(1, 1);
	invProj(3, 2) = 1.f / proj(2, 3);
	invProj(0, 3) = invProj(0, 0) * proj(0, 2);
	invProj(1, 3) = invProj(1, 1) * proj(1, 2);
	invProj(2, 3) = -1.f;
	invProj(3, 3) = invProj(3, 2) * proj(2, 2);

	std::cout << invProj * proj << std::endl;

	invProjRast = invProj * invRast;
}

void Camera::test() const
{
	// when eye.z is -1 then clip.w (i.e. -eye.z) is 1 thereby 
	// the perspective division can be avoided
	Vector4f e = view * Vector4f(eye.x(), eye.y(), eye.z() - 1.f, 1.f); // eye coordinates
	Vector4f c = proj * e; // clip coordinates
	float rc = 1.f / c.w(); // perspective divide
	Vector4f n = c * rc; // normalized device coordinates
	Vector4f s = rast * n; // window coordinates
	s.w() = rc;
}

Array3f Camera::sampleRay(Ray &ray, const Array2f& sample) const
{
	// transform from window coordinates to world coordinates
	Vector3f v;
	v.block<2, 1>(0, 0) = invProjRast.block<2, 4>(0, 0) * 
		Vector4f(sample.x(), sample.y(), 0.f, 1.f);
	v.z() = -1.f;
	v.normalize(); // camera is at the origin in eye coordinates
	ray = Ray(eye, invView.block<3, 3>(0, 0) * v, zNear, zFar);

	return Array3f(1.f, 1.f, 1.f);
}
