/*
 * Box.cpp
 *
 *  Created on: Feb 29, 2016
 *      Author: zxi
 */

#include "Box.h"

#include <float.h>

namespace mathtool {

Box3d::Box3d() {
}

Box3d::~Box3d() {
}

Box3d& Box3d::setFromPoints(const std::vector<Point3d>& points) {
  Vector3d min(FLT_MAX, FLT_MAX, FLT_MAX);
  Vector3d max(FLT_MIN, FLT_MIN, FLT_MIN);

  for (const auto& p : points) {
    for (int i = 0; i < 3; ++i) {
      if (p[i] < min[i])
        min[i] = p[i];
      if (p[i] > max[i])
        max[i] = p[i];
    }
  }

  return this->set(min, max);
}

Box3d& Box3d::setFromPoints(const std::vector<Vector3d>& points) {
  Vector3d min(FLT_MAX, FLT_MAX, FLT_MAX);
  Vector3d max(FLT_MIN, FLT_MIN, FLT_MIN);

  for (const auto& p : points) {
    for (int i = 0; i < 3; ++i) {
      if (p[i] < min[i])
        min[i] = p[i];
      if (p[i] > max[i])
        max[i] = p[i];
    }
  }

  return this->set(min, max);
}

} /* namespace mathtool */
