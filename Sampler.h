#pragma once
using namespace Eigen;

class Sampler
{
public:
	Sampler();
	~Sampler();
	inline float roll(pcg32::state_type stream = 0) { 
		rng.set_stream(stream);
		return unif(rng); 
	}
	// sample transformation from unit square
	static Array2f square2DiskConcentric(const Array2f& sample); 
	static Vector3f square2CosHemisphere(const Array2f& sample);
private:
	pcg32 rng;
	const std::uniform_real_distribution<float> unif;
};

