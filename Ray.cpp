#include "stdafx.h"
#include "Ray.h"

const float Ray::eps = 1e-18f;

Matrix3f Ray::frame(const Vector3f& n)
{
	Matrix3f m;
	if (std::fabs(n.x()) > std::fabs(n.y())) {
		float invLen = 1.0f / std::sqrt(n.x() * n.x() + n.z() * n.z());
		m.col(1) = Vector3f(n.z() * invLen, 0.0f, -n.x() * invLen);
	}
	else {
		float invLen = 1.0f / std::sqrt(n.y() * n.y() + n.z() * n.z());
		m.col(1) = Vector3f(0.0f, n.z() * invLen, -n.y() * invLen);
	}
	m.col(0) = m.col(1).cross(n);
	m.col(2) = n;

	return m;
}
