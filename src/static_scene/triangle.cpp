#include "triangle.h"

#include "CGL/CGL.h"
#include "GL/glew.h"

namespace CGL { namespace StaticScene {

Triangle::Triangle(const Mesh* mesh, size_t v1, size_t v2, size_t v3) :
    mesh(mesh), v1(v1), v2(v2), v3(v3) { }

BBox Triangle::get_bbox() const {

  Vector3D p1(mesh->positions[v1]), p2(mesh->positions[v2]), p3(mesh->positions[v3]);
  BBox bb(p1);
  bb.expand(p2); 
  bb.expand(p3);
  return bb;

}

bool Triangle::intersect(const Ray& r) const {
  
  Vector3D p1(mesh->positions[v1]), p2(mesh->positions[v2]), p3(mesh->positions[v3]);

  if (PART < 1)
    return false;
  else {
    Vector3D e1 = p2-p1, e2 = p3-p1, s = r.o - p1;
    Vector3D e1_d = cross(e1, r.d), s_e2 = cross(s, e2);
    double denom = dot(e1_d, e2);
    double u = -dot(s_e2, r.d) / denom;
    double v = dot(e1_d, s) / denom;

    if (!(u>=0 && v>=0 && u<=1 && v<=1 && 1-u-v>=0))
      return false;

    double t = -dot(s_e2, e1) / denom;
    return  !(t > r.max_t || t < r.min_t) ;
  }
}

bool Triangle::intersect(const Ray& r, Intersection *isect) const {
  
  // implement ray-triangle intersection. When an intersection takes
  // place, the Intersection data should be updated accordingly
  Vector3D p1(mesh->positions[v1]), p2(mesh->positions[v2]), p3(mesh->positions[v3]);
  Vector3D n1(mesh->normals[v1]), n2(mesh->normals[v2]), n3(mesh->normals[v3]);
  if (PART < 1)
    return false;
  else {

    Vector3D e1 = p2-p1, e2 = p3-p1, s = r.o - p1;
    Vector3D e1_d = cross(e1, r.d), s_e2 = cross(s, e2);
    double denom = dot(e1_d, e2);
    double u = -dot(s_e2, r.d) / denom;
    double v = dot(e1_d, s) / denom;

    if (!(u>=0 && v>=0 && u<=1 && v<=1 && 1-u-v>=0))
      return false;

    double t = -dot(s_e2, e1) / denom;
    if (t > r.max_t || t < r.min_t) 
      return false;

    isect->t = t;
    isect->primitive = this;
    isect->n = n2 * u + n3 * v + n1 * (1-u-v);
    isect->bsdf = get_bsdf();
    r.max_t = t;

    return true;
  }
}

void Triangle::draw(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_TRIANGLES);
  glVertex3d(mesh->positions[v1].x,
             mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x,
             mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x,
             mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}

void Triangle::drawOutline(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_LINE_LOOP);
  glVertex3d(mesh->positions[v1].x,
             mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x,
             mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x,
             mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}



} // namespace StaticScene
} // namespace CGL
