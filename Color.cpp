#include "stdafx.h"
#include "Color.h"


Color::Color()
{
}


Color::~Color()
{
}

float Color::X1931(float wavelength)
{
	float t1 = (wavelength - 442.0f) * ((wavelength < 442.0f) ? 0.0624f : 0.0374f);
	float t2 = (wavelength - 599.8f) * ((wavelength < 599.8f) ? 0.0264f : 0.0323f);
	float t3 = (wavelength - 501.1f) * ((wavelength < 501.1f) ? 0.0490f : 0.0382f);
	return (0.362f * std::exp(-0.5f * t1 * t1)) + (1.056f * std::exp(-0.5f * t2 * t2)) - (0.065f * std::exp(-0.5f * t3 * t3));
}

float Color::Y1931(float wavelength)
{
	float t1 = (wavelength - 568.8f) * ((wavelength < 568.8f) ? 0.0213f : 0.0247f);
	float t2 = (wavelength - 530.9f) * ((wavelength < 530.9f) ? 0.0613f : 0.0322f);
	return (0.821f * std::exp(-0.5f * t1 * t1)) + (0.286f * std::exp(-0.5f * t2 * t2));
}

float Color::Z1931(float wavelength)
{
	float t1 = (wavelength - 437.0f) * ((wavelength < 437.0f) ? 0.0845f : 0.0278f);
	float t2 = (wavelength - 459.0f) * ((wavelength < 459.0f) ? 0.0385f : 0.0725f);
	return (1.217f * std::exp(-0.5f * t1 * t1)) + (0.681f * std::exp(-0.5f * t2 * t2));
}

Array3f Color::SPD2XYZ(float wavelength, float power)
{
	Array3f xyz;
	xyz[0] = X1931(wavelength);
	xyz[1] = Y1931(wavelength);
	xyz[2] = Z1931(wavelength);
	xyz *= (RECP_INTY * power);

	return xyz;
}

const Matrix3f Color::transMat((Eigen::Matrix3f() << 
	3.2406f, -0.9689f, 0.0557f,
	-1.5372f, 1.8758f, -0.2040f,
	-0.4986f, 0.0415f, 1.0570f).finished());

Array3f Color::XYZ2RGB(const Array3f& xyz)
{
	Array3f rgb = transMat * xyz.matrix();
	return rgb.max(0.f).min(1.f);
}

// https://en.wikipedia.org/wiki/Standard_illuminant
float Color::illuminantA(float wavelength)
{
	constexpr float c = 299792458.f; // speed of light
	constexpr float k = 1.3806488e-23f; // Boltzmann constant
	constexpr float h = 6.62606957e-34f; /* Planck constant */
	constexpr float c1 = 2 * M_PI * h*c*c;
	constexpr float c2 = (h / k)*c;
	constexpr float T = 2855.54f; // color temperature in Kelvin
	constexpr float M560 = 8.40753e+11f; // absolute SPD of illuminant A at 560nm
	wavelength *= 1e-9f;  // wavelength in meters
	// Watts per unit surface area (m^-2) per unit wavelength (nm^-1) per steradian (sr^-1)
	return c1 * std::pow(wavelength, -5.f) /
		(M560 * (std::exp(c2 / (wavelength*T)) - 1.f)); // relative SPD
}
