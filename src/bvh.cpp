#include "bvh.h"

#include "CGL/CGL.h"
#include "static_scene/triangle.h"

#include <iostream>
#include <stack>

using namespace std;

namespace CGL { namespace StaticScene {

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {

  root = construct_bvh(_primitives, max_leaf_size);

}

BVHAccel::~BVHAccel() {
  if (root) delete root;
}

BBox BVHAccel::get_bbox() const {
  return root->bb;
}

void BVHAccel::draw(BVHNode *node, const Color& c) const {
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims))
      p->draw(c);
  } else {
    draw(node->l, c);
    draw(node->r, c);
  }
}

void BVHAccel::drawOutline(BVHNode *node, const Color& c) const {
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims))
      p->drawOutline(c);
  } else {
    drawOutline(node->l, c);
    drawOutline(node->r, c);
  }
}

BVHNode *BVHAccel::construct_bvh(const std::vector<Primitive*>& prims, size_t max_leaf_size) {

  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.


  BBox centroid_box, bbox;

  for (Primitive *p : prims) {
    BBox bb = p->get_bbox();
    bbox.expand(bb);
    Vector3D c = bb.centroid();
    centroid_box.expand(c);
  }

  BVHNode *node = new BVHNode(bbox);

  if (PART < 2) {
    node->prims = new vector<Primitive *>(prims);
    return node;
  }

  // If it's a leaf, fill up its vector.
  if (prims.size() <= max_leaf_size) {
    //primitives.insert(primitives.end(),prims.begin(),prims.end());
    node->prims = new vector<Primitive *>(prims);
    return node;
  }

  // Take axis with maximum extent of centroids!!
  int axis; double extent = 0;
  for (int i = 0; i < 3; ++i)
    if (centroid_box.extent[i] >= extent)
      axis = i, extent = centroid_box.extent[i];

  double split = centroid_box.min[axis] + .5*centroid_box.extent[axis];

  // Do that split
  vector<Primitive *> left_prims, right_prims;
  for (Primitive *p : prims) 
    if (p->get_bbox().centroid()[axis] < split) 
      left_prims.push_back(p);
    else
      right_prims.push_back(p);

  if (left_prims.empty()) {

    delete node;
    return construct_bvh(right_prims,  right_prims.size());

  } else if (right_prims.empty()) {

    delete node;
    return construct_bvh(left_prims,  left_prims.size());

  } else {

    node->l = construct_bvh(left_prims,  max_leaf_size);
    node->r = construct_bvh(right_prims, max_leaf_size);

    return node;

  }

}


bool BVHAccel::intersect(const Ray& ray, BVHNode *node) const {
  if (PART < 2) {
    for (Primitive *p : *(root->prims)) {
      total_isects++;
      if (p->intersect(ray)) 
        return true;
    }
    return false;
  }

  double t_min, t_max;
  if (!node->bb.intersect(ray, t_min, t_max) ||
      t_min > ray.max_t || t_max < ray.min_t) 
    return false;

  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims)) {
        total_isects++;
        if (p->intersect(ray)) 
          return true;
      }
    return false;
  }
  return intersect(ray, node->l) || intersect(ray, node->r);

}

bool BVHAccel::intersect(const Ray& ray, Intersection* i, BVHNode *node) const {

  if (PART < 2) {
    bool hit = false;
    for (Primitive *p : *(root->prims)) {
      total_isects++;
      if (p->intersect(ray, i)) 
        hit = true;
    }
    return hit;
  }

  double t_min, t_max;
  if (!node->bb.intersect(ray, t_min, t_max) ||
      t_min > ray.max_t || t_max < ray.min_t) 
    return false;

  bool hit = false;
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims)) {
        total_isects++;
        if (p->intersect(ray, i)) 
          hit = true;
      }
    return hit;
  }

  bool left = intersect(ray, i, node->l);
  bool right= intersect(ray, i, node->r);
  return left || right;

}

}  // namespace StaticScene
}  // namespace CGL
