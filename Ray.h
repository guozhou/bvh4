#pragma once
using namespace Eigen;

struct Ray
{
	Ray() {}
	Ray(const Vector3f& _org, const Vector3f& _dir,
		float _tnear = 1e-4f, float _tfar = finf)
		: org(_org), dir(_dir), tnear(_tnear), tfar(_tfar) {
		assert(tnear < tfar);
		rdir = reciprocal<Array3f>(dir.array()).matrix(); // explicit template argument!!!
	}

	Vector3f org;
	Vector3f dir;
	Vector3f rdir;
	float tnear;
	float tfar;

	Vector3f Ng;	// unnormalized geometry normal
	float u;		// barycentric u coordinate of hit
	float v;		// barycentric v coordinate of hit
	int geomID;		// geometry ID
	int primID;		// primitive ID

	template<typename T>
	static T reciprocal(const T& d)
	{
		T r = 1.f / ((d.abs() < eps).select(eps, d));
		// newton-raphson refinement https://en.wikipedia.org/wiki/Division_algorithm
		return r + r - d * r * r;
	}

	static float reciprocal(float d)
	{
		float r = 1.f / (std::fabs(d) < eps ? eps : d);
		return r + r - d * r * r;
	}

	static Matrix3f frame(const Vector3f& n);

	std::ostream& operator<<(std::ostream& cout) {
		return cout << "{ "
			<< "org = " << org
			<< ", dir = " << dir
			<< ", near = " << tnear
			<< ", far = " << tfar
			<< ", geomID = " << geomID
			<< ", primID = " << primID
			<< ", u = " << u
			<< ", v = " << v
			<< ", Ng = " << Ng
			<< " }";
	}

private:
	static const float eps;
};