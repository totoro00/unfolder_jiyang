/*
 * unfolder.h
 *
 *  Created on: Nov 14, 2014
 *      Author: zhonghua
 */

#ifndef UNFOLDER_H_
#define UNFOLDER_H_

#include <map>
#include <unordered_map>
#include <unordered_set>
#include <cfloat>
#include <cstdint>
using namespace std;

#include "config.h"
#include "model.h"
//#include "CD.h"

namespace masc {
// forward declaration
class Splitter;
}


using namespace masc;

struct DualGraphEdge {
  int fid1;
  int fid2;
  float weight;

  DualGraphEdge() {
    fid1 = -1;
    fid2 = -1;
    weight = FLT_MAX;
  }
  DualGraphEdge(const DualGraphEdge& rhs) {
    fid1 = rhs.fid1;
    fid2 = rhs.fid2;
    weight = rhs.weight;
  }
  DualGraphEdge(int fid1, int fid2, float weight) :
      fid1(fid1), fid2(fid2), weight(weight) {
  }

  bool operator <(const DualGraphEdge& rhs) const {
    return this->weight < rhs.weight;
  }

  bool operator >(const DualGraphEdge& rhs) const {
    return this->weight > rhs.weight;
  }
};

/// Crease
/// fold a face (fid) around an edge(vid1, vid2) by folding_angle rad according to parent face (pdi)
struct Crease {
  /// vertex id 1
  int vid1;

  /// vertex id 2
  int vid2;

  /// folding angle of the crease
  float folding_angle;

  /// parent face id
  int pid;

  /// face id
  int fid;

  Crease() :
      vid1(0), vid2(0), folding_angle(0), fid(0), pid(0) {
  }
  ;
  Crease(int vid1, int vid2, int fid, int pid, float folding_angle) :
      vid1(vid1), vid2(vid2), fid(fid), pid(pid), folding_angle(folding_angle) {
  }
};

class Unfolder {
public:
  /// constructor
  Unfolder(model* m, const Config& conifg);

  /// virtual deconstrutor
  virtual ~Unfolder();

  /// unfold the model
  void buildUnfolding();

  /// build the unfolding from weights file and return # overlaps
  int buildFromWeights(const string& path);

  /// build the unfolding from weights and return # overlaps
  int buildFromWeights(const vector<float>& weights, bool checkOverlaps = true);

  /// unfold the model to percentage
  void unfoldTo(double percentage);

  /// -------------------------------------------------------
  /// dumping
  /// -------------------------------------------------------

  //  type = 1: polygon boundary & single path creases
  //  type = 2: all
  //  type = 5: type 1 + spanning tree 4
  //  type = 9: type 1 + extra cuts 8
  void dumpSVG(const string& path, const int type = 2);


  /// -------------------------------------------------------
  /// access
  /// -------------------------------------------------------

  /// get model
  model* getModel() const;

  const Config& getConfig() const;

  /// get unfolded coordinates
  const MESH& getUnfolded() const;

  /// get input model's filename
  const string& getFilename() const;

  /// get rotation axis, for transforming the model back
  const Vector3d& getRotationAxis() const;

  /// get rotation angle in degree, for transforming the model back
  const double getRotationAngle() const;

  /// get translation, for transforming the model back
  const Vector3d& getTranslation() const;

  /// get model color
  const Vector3d& getColor() const;

  // set model color
  void setColor(double r, double g, double b);

  const bool isFlattened() const;

  /// check whether current unfolding has overlap or not
  /// return the number of pairs faces that overlapped
  int checkOverlap();

  /// check if a given pair of faces (indexed by i and j) overlaps.
  /// return true if they overlap
  bool checkOverlap(uint i, uint j);

  // count local overlaps due to insufficient count for hyperbolic vertex
  int checkLocalOverlaps();

  /// check whether current folded model has collision or not
  bool checkCollision();

  /// get parent face id
  const int getParentFaceId(const int fid) const {
    return this->m_parents[fid];
  }

  const set<pair<int, int>>& getSelectedEdges() const {
    return this->m_selected_edges;
  }

  const double getMaxPathLength() const {
    return this->m_max_path_len;
  }

  /// called once after instantiated
  void measureModel();

  /// build the topology of the unfolding
  void rebuildModel();

  /// get current base id
  uint getBaseFaceID() const {
    return m_base_face_id;
  }

  const vector<float>& getEdgeWeights() const {
    return this->m_weights;
  }

  /// get overlapping face pairs
  const vector<vector<uint> >& getOverlppingFacePairs() const {
    return m_overlap_pairs;
  }

  const uint64_t getCheckOverlappingCalls() const {
    return m_check_overlapping_calls;
  }

  const uint64_t getTotalCheckOverlappingTime() const {
    return m_total_check_overlapping_time;
  }

  const uint64_t getCheckCollisionCalls() const {
    return m_check_collision_calls;
  }

  const uint64_t getTotalCheckCollisionTime() const {
    return m_total_check_collison_time;
  }

private:

  /// initialize the unfolding
  void initUnfold();

  /// build the dual graph of the original mesh
  void buildDualGraph(GRAPH& g);

  /// compute the MST on the dual graph
  void buildMST(GRAPH& g);

  /// find best base face
  int findBestBaseFace();

  /// rebuild the tree rooted at base_face from selected edges
  /// return average path length
  void rebuildTree(int base_face, set<pair<int, int>> selected_edges);

  /// unfold to config: folding angles of each edge
  void unfoldTo(const vector<double>& folding_angles);

  /// unfold the model by linear interpolation
  /// p: percentage of folding
  /// count: number of faces to fold in order
  void linearUnfoldTo(double p, const int count = INT_MAX);

  /// unfold the model by order
  void orderedUnfoldTo(double p);

  /// compute the unfolding based the MST
  /// 	1. parent of each face
  ///		2. ordered face list for a fast implementation
  void computeUnfolding();

  /// align the base_face to XZ-plane (y=0) and with one vertex of face_base at origin (0,0,0)
  void alignModel();

  /// check whether two faces intersected or not
  bool hasIntersection(const int fid1, const int fid2);

  /// get the intersection of segment v1 -> v2 on face 1
  bool getIntersectionOfPlane(const int fid1, const Vector3d& v1,
      const Vector3d& v2, Vector3d& outPoint);

  // check whether a point is in face 1
  bool pointInTriangle(const int fid1, const Vector3d& p);

  // find the boundary of the unfolding
  void findBoundary();

  // find single path for all the creases
  void findSinglePath();

  // travel from vid
  void findSinglePath(int vid, unordered_set<int>& visited_vids,
      set<pair<int, int>>& visited_edges, vector<int>& path);

  /// get the coordinates of a vertex which has original vertex id org_vid in face[fid]
  /// note: one vertex in original mesh can be mapped to several different coordinates in the unfolding
  const Vector3d& getUnfoldedVertex(uint fid, uint org_vid) const;

  /// check whether an edge(vid1,vid2) is CCW in face[fid]
  bool isEdgeCCW(uint fid, uint vid1, uint vid2) const;

  /// shrink unfolded (flat) model's triangles
  void shrink();

  /////////////////////////////////////////////////////////////////////////
  /// add tabs on the net
  void addTabs();

  /// add tab for an original edge, which has two edges in the net...
  bool addTabForCutEdge(const int org_eid);

  /// get fid of the net by given vids
  int getNetFidByVids(const int vid1, const int vid2);

  /// add tab for a new edge of face fid in the net...
  bool addTabForNewEdge(const int vid1, const int vid2, const int fid);

  /// check whether an extended triangle from face fid, overlap with existing net
  bool isOverlapWithNet(const Vector3d& p1, const Vector3d& p2,
      const Vector3d& p3, const int fid);

  /// check whether current unfolding has overlap or not using interval tree
  /// return number of overlapping pairs
  uint checkOverlap_itree();

  //////////////////////////////////////////////////////////////////////////

  /// ordered face list
  vector<uint> m_ordered_face_list;

  /// shared edge <fid1, fid2> = <eid>
  unordered_map<uint, unordered_map<uint, uint> > m_shared_edge;

  /// model
  model* m_m;

  /// base face_id
  uint m_base_face_id;

  /// reference vector for unfolding
  Vector3d m_c;

  /// original mesh
  MESH m_org;

  /// unfolded mesh
  MESH m_unfolded;

  /// config model
  Config m_config;

  /// parent face id of edge face, -1 means no parent
  vector<int> m_parents;

  /// ---------------------------------------------------------------
  /// for visualization
  /// ---------------------------------------------------------------
  Vector3d m_color;

  /// ---------------------------------------------------------------
  /// for align back
  /// ---------------------------------------------------------------

  /// rotation axis for aligning
  Vector3d m_rotation_axis;

  /// rotation angle in rad
  double m_rotation_angle;

  /// translation
  Vector3d m_transilation;

  /// translation in ori file...
  Vector3d m_translation2;

  /// ---------------------------------------------------------------
  /// for output obj/ori
  /// ---------------------------------------------------------------

  /// output vertex array
  vector<Vector3d> m_vs;

  /// vid, adjacent face ids
  unordered_map<int, set<int>> m_vfmap;

  // org_vid -> new_vids
  unordered_map<int, set<int>> m_vvmap;

  // org_edge_id -> { <new_vid1,new_vid2> }
  unordered_map<int, set<pair<int, int>>> m_eemap;

  // <new_vid1, new_vid2> -> org_edge_id
  map<pair<int,int>, int> m_vemap;

  /// creases
  vector<Crease> m_cs;

  // vid1, vid2 is a crease, map to crease id
  map<pair<int,int>, int> m_crease_lines;

  /// output face array
  vector<vector<int> > m_fs;

  // polygonal boundary of the unfolding
  vector<int> m_boundary;

  // single path that go through all the crease lines
  vector<int> m_single_path;

  /// edge weight
  GRAPH m_edge_weight;

  // edge weights
  vector<float> m_weights;

  /// connected edges in the unfolding
  set<pair<int,int>> m_selected_edges;

  // eid -> true/false
  set<int> m_fold_edges;

  /// width, ?, height
  Vector3d m_net_size;

  // used in steepest edge rule, v+ should not be considered
  int m_top_vertex_id;

  float m_max_edge_length;
  float m_min_edge_length;

  // different heuristics
  Splitter* m_spliiter;

  // whether found a non-overlapping unfolding
  bool m_is_flattened;

  // count for flat edges, use for less cuts
  int m_flat_edges;

  // for collision detection
  //CD* m_cd;

  // unfolding properties
  double m_avg_path_len;
  double m_max_path_len;

  // overlapping triangles
  vector< vector<uint> > m_overlap_pairs;

  uint64_t m_check_collision_calls;
  uint64_t m_total_check_collison_time;

  uint64_t m_total_check_overlapping_time;
  uint64_t m_check_overlapping_calls;

};

#endif /* UNFOLDER_H_ */
