/*
 * Box.h
 *
 *  Created on: Feb 29, 2016
 *      Author: zxi
 */

#ifndef SRC_MATHTOOL_BOX_H_
#define SRC_MATHTOOL_BOX_H_

#include <vector>

#include "mathtool/Point.h"
#include "mathtool/Vector.h"


namespace mathtool {

class Box3d {
public:
  Box3d();
  virtual ~Box3d();

  Box3d& set(const Vector3d& min, const Vector3d& max) {
    this->m_min = min;
    this->m_max = max;
    return *this;
  }

  Box3d& setFromPoints(const std::vector<Point3d>& points);
  Box3d& setFromPoints(const std::vector<Vector3d>& points);

  const Vector3d& getMin() const {
    return m_min;
  }

  const Vector3d& getMax() const {
    return m_max;
  }

  const Vector3d getDim() const {
    return m_max - m_min;
  }

private:
  Vector3d m_min;
  Vector3d m_max;
};

} /* namespace mathtool */

#endif /* SRC_MATHTOOL_BOX_H_ */
