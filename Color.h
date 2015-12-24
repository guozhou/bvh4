#pragma once
using namespace Eigen;

class Color
{
public:
	Color();
	~Color();
	// range of wavelength which is visible to human
	static constexpr float WAVELENGTH_MIN = 380.0f;
	static constexpr float WAVELENGTH_MAX = 700.0f;
	static constexpr float WAVELENGTH_LEN = WAVELENGTH_MAX - WAVELENGTH_MIN;
	static constexpr float RECP_INTY = 0.009358239977091027f;
	/**
	* Close analytic approximations of the CIE 1931 XYZ color curves.
	* From the paper "Simple Analytic Approximations to the CIE XYZ Color Matching
	* Functions" by Wyman et al.
	*
	* @param w The wavelength of light in nm.
	* @returns The sensitivity of the curve at that wavelength.
	*/
	static float X1931(float wavelength);
	static float Y1931(float wavelength);
	static float Z1931(float wavelength);

	static Array3f SPD2XYZ(float wavelength, float power);
	static Array3f XYZ2RGB(const Array3f& xyz); // CIE XYZ -> sRGB
	static Array3f SPD2RGB(float wavelength, float power) {
		return Color::XYZ2RGB(Color::SPD2XYZ(wavelength, power));
	}

	static inline float sampleWavelength(float sample) {
		return WAVELENGTH_MIN + WAVELENGTH_LEN * sample;
	}
	// analytical SPD of CIE illuminant A
	static float illuminantA(float wavelength);

//private:
	static const Matrix3f transMat; // matrix of color transformation
};

