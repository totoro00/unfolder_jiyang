/*
 * OverlappingChecker.cpp
 *
 *  Created on: Feb 29, 2016
 *      Author: zxi
 */

#include "OverlappingChecker.h"

#include <float.h>

#include "unfolder.h"

#include "mathtool/Box.h"

#include "itree/RectKD.h"
#include "itree/MiddleStructure.h"

namespace masc {
namespace unfolding {

///////////////////////////////////////////////////////////////////////////////
// BruteForceChecker
///////////////////////////////////////////////////////////////////////////////

double BruteForceChecker::checkOverlapping(const MESH& mesh,
    const Config& config) {

  int count = 0;
  const int F = this->m_unfolder->getModel()->t_size;

  for (int i = 0; i < F; i++) {
    for (int j = i + 1; j < F; j++) {
      if (this->m_unfolder->checkOverlap(i, j)) {
        count++;
      }
    }
  }

  return count;
}

///////////////////////////////////////////////////////////////////////////////
// ITreeChecker
///////////////////////////////////////////////////////////////////////////////
double ITreeChecker::checkOverlapping(const MESH& mesh, const Config& config) {

  const int F = this->m_unfolder->getModel()->t_size;

  typedef Interval<EndPoint> Int;
  typedef RectKD<EndPoint, 2> Rect2D;
  typedef MiddleStructure<Rect2D, 2> MTree;
  typedef vector<Rect2D*>::iterator IT2;

  vector<Rect2D*> rects;

  for (int i = 0; i < F; i++) {

    float max_x = -FLT_MAX, min_x = FLT_MAX, max_y = -FLT_MAX, min_y =
    FLT_MAX;

    for (int p = 0; p < 3; p++) {
      const auto& vi1 = mesh[i][p].second;
      const double a[2] = { vi1[0], vi1[2] };
      max_x = (a[0] > max_x) ? a[0] : max_x;
      min_x = (a[0] < min_x) ? a[0] : min_x;
      max_y = (a[1] > max_y) ? a[1] : max_y;
      min_y = (a[1] < min_y) ? a[1] : min_y;
    }

    //
    Rect2D * p = new Rect2D(i);
    p->setInterval(Int(p, EndPoint(min_x), EndPoint(max_x)), 0);
    p->setInterval(Int(p, EndPoint(min_y), EndPoint(max_y)), 1);

    rects.push_back(p);
  }

  RectangleTree<MTree, 2> Tree;
  if (Tree.Build(rects) == false) {
    //failed, use brute force
    cerr
        << "! WARNING: Build interval tree failed. Resort to brute-force method"
        << endl;
    for (auto rect : rects)
      delete (rect);
    return UINT_MAX;
  }

  uint count = 0;
  for (const auto rect : rects) {
    uint i = rect->getVID();
    Tree.query(rect);
    const Rect2D::Intersect& intersections = rect->getIntersections();

    for (uint j : intersections) {
      if (i > j)
        continue; //this will be checked by (j,i) pair
      if (this->m_unfolder->checkOverlap(i, j)) {
        count++;
      }
    }
  }

  //free mem
  for (auto rect : rects)
    delete (rect);

  return count;
}

} /* namespace unfolding */
} /* namespace masc */
