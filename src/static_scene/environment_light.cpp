#include "environment_light.h"

#include <algorithm>
#include <iostream>
#include <fstream>

namespace CGL { namespace StaticScene {

EnvironmentLight::EnvironmentLight(const HDRImageBuffer* envMap)
    : envMap(envMap) {
    	init();
}

EnvironmentLight::~EnvironmentLight() {
    delete[] pdf_envmap;
    delete[] conds_y;
    delete[] marginal_y;
}


void EnvironmentLight::init() {
	uint32_t w = envMap->w, h = envMap->h;
    pdf_envmap = new double[w * h];
	conds_y = new double[w * h];
	marginal_y = new double[h];

	std::cout << "[PathTracer] Initializing environment light...";

	// Store the environment map pdf to pdf_envmap
	// Store the marginal distribution for y to marginal_y
	// Store the conditional distribution for x given y to conds_y

	double sum = 0;
	for (int j = 0; j < h; ++j) {
		for (int i = 0; i < w; ++i) {
            pdf_envmap[w * j + i] = envMap->data[w * j + i].illum() * sin(M_PI * (j+.5) / h);
            sum += pdf_envmap[w * j + i];
		}
	}
	double I = sum;
	double msum = 0;
	for (int j = 0; j < h; ++j) {
		for (int i = 0; i < w; ++i) {
            pdf_envmap[w * j + i] = pdf_envmap[w * j + i] / I;
            msum += pdf_envmap[w * j + i];
		}
		marginal_y[j] = msum;
	}

	for (int j = 0; j < h; ++j) {
		double csum = 0;
		for (int i = 0; i < w; ++i) {
			csum += pdf_envmap[w * j + i] / (marginal_y[j] - (j>0?marginal_y[j-1]:0));
            conds_y[w * j + i] = csum;
		}
	}

	std::cout << "done." << std::endl;
}

// Helper functions

Vector2D EnvironmentLight::theta_phi_to_xy(const Vector2D &theta_phi) const {
    uint32_t w = envMap->w, h = envMap->h;
    double x = theta_phi.y / 2. / M_PI * w;
    double y = theta_phi.x / M_PI * h;
    return Vector2D(x, y);
}

Vector2D EnvironmentLight::xy_to_theta_phi(const Vector2D &xy) const {
    uint32_t w = envMap->w, h = envMap->h;
    double x = xy.x;
    double y = xy.y;
    double phi = x / w * 2.0 * M_PI;
    double theta = y / h * M_PI;
    return Vector2D(theta, phi);
}

Vector2D EnvironmentLight::dir_to_theta_phi(const Vector3D &dir) const {
    dir.unit();
    double theta = acos(dir.y);
    double phi = atan2(-dir.z, dir.x) + M_PI;
    return Vector2D(theta, phi);
}

Vector3D EnvironmentLight::theta_phi_to_dir(const Vector2D& theta_phi) const {
    double theta = theta_phi.x;
    double phi = theta_phi.y;

    double y = cos(theta);
    double x = cos(phi - M_PI) * sin(theta);
    double z = -sin(phi - M_PI) * sin(theta);

    return Vector3D(x, y, z);
}

Spectrum EnvironmentLight::bilerp(const Vector2D& xy) const {
	uint32_t w = envMap->w;
	const std::vector<Spectrum>& data = envMap->data;
	double x = xy.x, y = xy.y;
	Spectrum ret;
	for (int i = 0; i < 4; ++i)
		ret += (i%2 ? x-floor(x) : ceil(x)-x) * 
			   (i/2 ? y-floor(y) : ceil(y)-y) * 
			   data[w * (floor(y) + i/2) + floor(x) + i%2];
	return ret;
}


Spectrum EnvironmentLight::sample_L(const Vector3D& p, Vector3D* wi,
                                    float* distToLight,
                                    float* pdf) const {
  
	// First implement uniform sphere sampling for the environment light
	// Later implement full importance sampling

	// Uniform
//    *wi = sampler_uniform_sphere.get_sample();
//    *distToLight = INF_D;
//    *pdf = 1.0 / (4.0 * M_PI);

	// Importance
    Vector2D sample = sampler_uniform2d.get_sample();
    uint32_t w = envMap->w, h = envMap->h;
    uint32_t y = std::upper_bound(marginal_y, marginal_y + h, sample.y) - marginal_y;
    uint32_t x = std::upper_bound(conds_y + y * w, conds_y + (y+1) * w, sample.x) - (conds_y + y * w);
    x = clamp(x, 0, envMap->w - 1.5001);
    y = clamp(y, 0, envMap->h - 1.5001);
    Vector2D xy = Vector2D(x+.5, y+.5);

    Vector2D theta_phi = xy_to_theta_phi(xy);
    *wi = theta_phi_to_dir(theta_phi);
    *distToLight = INF_D;
    *pdf = w * h * pdf_envmap[y * w + x] / (2.0 * M_PI * M_PI * sin(theta_phi.x));

    return bilerp(xy);
}

Spectrum EnvironmentLight::sample_dir(const Ray& r) const {
  
	// Use the helper functions to convert r.d into (x,y)
	// then bilerp the return value 

	return bilerp(theta_phi_to_xy(dir_to_theta_phi(r.d)));

}

} // namespace StaticScene
} // namespace CGL
