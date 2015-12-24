#include "stdafx.h"
#include "Sampler.h"


Sampler::Sampler()
	:rng(pcg_extras::seed_seq_from<std::random_device>()), unif(0.f, 1.f)
{
}


Sampler::~Sampler()
{
}

Array2f Sampler::square2DiskConcentric(const Array2f& sample)
{
	float r1 = 2.0f*sample.x() - 1.0f;
	float r2 = 2.0f*sample.y() - 1.0f;

	/* Modified concentric map code with less branching (by Dave Cline), see
	http://psgraphics.blogspot.ch/2011/01/improved-code-for-concentric-map.html */
	float phi, r;
	if (r1 == 0 && r2 == 0) {
		r = phi = 0;
	}
	else if (r1*r1 > r2*r2) {
		r = r1;
		phi = (M_PI / 4.0f) * (r2 / r1);
	}
	else {
		r = r2;
		phi = (M_PI / 2.0f) - (r1 / r2) * (M_PI / 4.0f);
	}

	float cosPhi = std::cosf(phi), sinPhi = std::sinf(phi);

	return Array2f(r * cosPhi, r * sinPhi);
}

Vector3f Sampler::square2CosHemisphere(const Array2f& sample)
{
	Array2f p = square2DiskConcentric(sample);
	float z = std::sqrt(std::max(1e-10f, 1.0f - p.x()*p.x() - p.y()*p.y()));

	return Vector3f(p.x(), p.y(), z);
}
