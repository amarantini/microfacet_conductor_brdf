#include "bsdf.h"

#include <iostream>
#include <algorithm>
#include <utility>

using std::min;
using std::max;
using std::swap;
 
namespace CGL {

void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {

    Vector3D z = Vector3D(n.x, n.y, n.z);
    Vector3D h = z;
    if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z)) h.x = 1.0;
    else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z)) h.y = 1.0;
    else h.z = 1.0;

    z.normalize();
    Vector3D y = cross(h, z);
    y.normalize();
    Vector3D x = cross(z, y);
    x.normalize();

    o2w[0] = x;
    o2w[1] = y;
    o2w[2] = z;
}


// Diffuse BSDF //
Spectrum DiffuseBSDF::eval(const Vector3D& wo, const Vector3D& wi) {
  return reflectance * (1.0 / PI);
}

Spectrum DiffuseBSDF::sample(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *wi = sampler.get_sample(pdf);
  return reflectance * (1.0 / PI);
}


// Mirror BSDF //

Spectrum MirrorBSDF::eval(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum MirrorBSDF::sample(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1;
  reflect(wo,wi);
  return  reflectance * (1./wo[2]);
}





// Microfacet BSDF ///////////////////////

double MicrofacetBSDF::G(const Vector3D& wo, const Vector3D& wi) {
    // Shadowing-masking term
    return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D& h) {
    /*
	  TODO:
	  Compute Beckmann normal distribution function (NDF) here.
	  You will need the roughness "alpha".
	*/
    Vector3D h_normalized = h.unit();
    double cos_theta = h_normalized.z;
    if(cos_theta==0)
      return 0;
    double theta_h = acos(cos_theta); 
    return exp((-1)*pow(tan(theta_h),2.0)/pow(alpha,2.0))/(PI*pow(alpha,2.0)*pow(cos_theta,4.0));
}

Spectrum MicrofacetBSDF::F(const Vector3D& wi) {
    /*
      TODO:
      Compute Fresnel term for reflection on air-conductor interface.
      You will need both "eta" and "k", both of which are Spectrum.
    */
    Spectrum Rs = ((eta*eta+k*k)-2*eta*cos_theta(wi)+pow(cos_theta(wi),2.0))
                    /((eta*eta+k*k)+2*eta*cos_theta(wi)+pow(cos_theta(wi),2.0));
    Spectrum Rp = ((eta*eta+k*k)*pow(cos_theta(wi),2.0)-2*eta*cos_theta(wi)+Spectrum(1,1,1))
                      /((eta*eta+k*k)*pow(cos_theta(wi),2.0)+2*eta*cos_theta(wi)+Spectrum(1,1,1));
    return (Rs+Rp)/2;
}



Spectrum MicrofacetBSDF::eval(const Vector3D& wo, const Vector3D& wi) {
    /*
      TODO:
      Implement microfacet model here
      Note that you will return the BRDF only, without the cosine term
    */
    if(wo.z<=0 || wi.z<=0)
      return Spectrum();
    Spectrum f = F(wi);
    double g = G(wo,wi);
    double d = D(((wo+wi)/2).unit());
    return f*g*d / (4 * wo.z * wi.z);
}

Spectrum MicrofacetBSDF::sample(const Vector3D& wo, Vector3D* wi, float* pdf) {
    /*
	  TODO: 
	  *Importance* sample Beckmann normal distribution function (NDF) here.
	  Note: You should fill in the sampled direction *wi and the corresponding *pdf,
	   	    and return the sampled BRDF value.
	*/

    // TODO 1
    // - Sample theta_h and phi_h
    Vector2D sample = sampler.get_sample();
    double r1 = sample.x;
    double r2 = sample.y;

    double theta_h = atan(sqrt((-1)*alpha*alpha*log(1-r1)));
    double phi_h = 2 * PI * r2;
    Vector3D h = Vector3D(sin(theta_h)*cos(phi_h),sin(theta_h)*sin(phi_h),cos(theta_h));

    // TODO 2
    // - Calculate outgoing direction *wi from sampled h
    *wi = 2.0 * h * dot(h,wo)-wo;

    // This is left as a gift for you. Make sure to comment this out 
    // when you are using the cosineHemisphereSampler
    if (dot(wo, h) <= 0.0 || wo.z <= 0.0 || wi->z <= 0.0) {
        *pdf = 10.0;
        return Spectrum();
    }

    // TODO 3
    // - Calculate *pdf of sampling *wi 
    double sin_theta = sin(theta_h);
    double cos_theta = cos(theta_h);
    double pdf_theta;
    if(cos_theta > 0)
      pdf_theta = 2.0*sin_theta/(pow(alpha,2.0)*pow(cos_theta,3.0)) * exp(-pow(tan(theta_h),2.0)/pow(alpha,2.0));
    else
      pdf_theta = 0.0;
    double pdf_phi = 1.0 / (2.0 * PI);
    double pdf_w_h;
    if(sin_theta > 0) {
      pdf_w_h = pdf_theta * pdf_phi / sin_theta;
      *pdf = pdf_w_h / 4.0 / dot(*wi,h);
    } else {
      *pdf = 10.0;
    }
    

    // Comment this line once you have finished with the importance sampling 
    *wi = cosineHemisphereSampler.get_sample(pdf);
    
    return MicrofacetBSDF::eval(wo, *wi);
}

//////////////////////////////////////////




// Refraction BSDF //

Spectrum RefractionBSDF::eval(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum RefractionBSDF::sample(const Vector3D& wo, Vector3D* wi, float* pdf) {

  *pdf = 1.0f;

  float cosTheta = cos_theta(wo);
  bool entering = cosTheta > 0;
  float ei = 1.f, et = ior;
  if (!entering) {
      swap(ei, et);
      cosTheta = -cosTheta;
  }
  float inveta = et / ei;
  float inveta2 = inveta * inveta;

  if (refract(wo, wi, ior))
      return inveta2 / cosTheta * transmittance;
  else
      return Spectrum();
}

// Glass BSDF //

Spectrum GlassBSDF::eval(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlassBSDF::sample(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // compute Fresnel coefficient and use it as the probability of reflection
  float R0 = (ior-1.0)*(ior-1.0)/((ior + 1.0)*(ior + 1.0));
  float cosTheta = cos_theta(wo);
  float f = 1 - fabs(cosTheta);
  float g = ((f * f) * (f * f)) * f;
  float fresnel_coe = R0 + (1.0 - R0)*g;

  bool entering = cos_theta(wo) > 0;
  float ei = 1.f, et = ior;
  if (!entering) {
      swap(ei, et);
      cosTheta = -cosTheta;  // be careful here, want cosTheta to be
                             // positive for everything below
  }
  float inveta = et / ei;
  float inveta2 = inveta * inveta;

  if (!refract(wo, wi, ior)) {
    // total internal reflection; always reflect
    *pdf = 1.0;
    reflect(wo, wi);
    return (1 / cosTheta) * reflectance;
  }

  if ((double)(std::rand()) / RAND_MAX < fresnel_coe) {
    *pdf = fresnel_coe;
    reflect(wo, wi);
    return (fresnel_coe / cosTheta) * reflectance;
  } else {
    // refraction ray has already been computed
    float one_minus_fresnel = 1.0f - fresnel_coe;
    *pdf = one_minus_fresnel;
    return (one_minus_fresnel * inveta2 / cosTheta) * transmittance;
  }
}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {

  // Implement reflection of wo about normal (0,0,1) and store result in wi.
  *wi = Vector3D(-wo[0],-wo[1],wo[2]);

}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {

  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.

  bool entering = cos_theta(wo) > 0;

  float ei = 1.f, et = ior;
  if (!entering) swap(ei, et);

  float sini2 = sin_theta2(wo);
  float eta = ei / et;
  float sint2 = eta * eta * sini2;
  if (sint2 > 1.f) return false;
  float cost = sqrt(1.0f - sint2);

  if (entering) cost = -cost;
  float sint_over_sini = eta;

  *wi = Vector3D(-sint_over_sini * wo.x, -sint_over_sini * wo.y, cost);

  return true;

}

// Emission BSDF //

Spectrum EmissionBSDF::eval(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum EmissionBSDF::sample(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1.0 / PI;
  *wi  = sampler.get_sample(pdf);
  return Spectrum();
}

} // namespace CGL
