#include "sphere.h"

#include <cmath>

#include  "../bsdf.h"
#include "../misc/sphere_drawing.h"

namespace CGL { namespace StaticScene {

bool Sphere::test(const Ray& r, double& t1, double& t2) const {

  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.

  if (PART < 1)
    return false;
  else {
    double a = dot(r.d,r.d);
    double b = 2*dot(r.o-o,r.d);
    double c = dot(r.o-o,r.o-o)-r2;
    double disc = b*b-4.*a*c;
    if (disc < 0) return false;
    disc = sqrt(disc);
    t1 = (-b+disc)/(2.*a);
    t2 = (-b-disc)/(2.*a);;//2.*c/(-b+disc);
    if (t1 > t2) 
      std::swap(t1,t2);

    return true;
  }

}

bool Sphere::intersect(const Ray& r) const {

  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.
  if (PART < 1)
    return false;
  else {
    double t1, t2;
    if (!test(r,t1,t2))
      return false;

    if (t1 >= r.min_t && t1 <= r.max_t) {
      r.max_t = t1;
      return true;
    }
    if (t2 >= r.min_t && t2 <= r.max_t) {
      r.max_t = t2;
      return true;
    }
    return false;
  }

}

bool Sphere::intersect(const Ray& r, Intersection *i) const {

  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.
  if (PART < 1)
    return false;
  else {
    double t1, t2, t;
    if (!test(r,t1,t2))
      return false;

    if (t1 >= r.min_t && t1 <= r.max_t) 
      t = t1;
    else if (t2 >= r.min_t && t2 <= r.max_t) 
      t = t2;
    else 
      return false;

    r.max_t = t;
    i->t = t;
    i->bsdf = get_bsdf();
    i->primitive = this;
    i->n = (r.o+t*r.d - o).unit();

    return true;
  }

}

void Sphere::draw(const Color& c) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color& c) const {
    //Misc::draw_sphere_opengl(o, r, c);
}


} // namespace StaticScene
} // namespace CGL
