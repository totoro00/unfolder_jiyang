/*
 * unfolder.cpp
 *
 *  Created on: Nov 14, 2014
 *      Author: zhonghua
 */

#include "unfolder.h"

#include <cfloat>
#include <climits>
#include <ctime>
#include <queue>
#include <fstream>
#include <map>
#include <numeric>
#include <functional>
#include <algorithm>

#include "intersection.h"
//#include "CD.h"
#include "OverlappingChecker.h"
#include "Splitter.h"
using namespace masc;

Unfolder::Unfolder(model* m, const Config& config) {
  this->m_m = m;
  this->m_config = config;
  this->m_base_face_id = UINT_MAX;
  this->m_rotation_angle = 0.0;
  this->m_flat_edges = 0;
  this->m_color = Vector3d(0.8, 0.8, 0.8);
  this->m_is_flattened = false;
  this->m_max_edge_length = FLT_MIN;
  this->m_min_edge_length = FLT_MAX;
  this->m_top_vertex_id = INT_MIN;
  this->m_spliiter = nullptr;
  //this->m_cd = new RAPID_CD(this);

  this->m_max_path_len = -1;
  this->m_avg_path_len = -1;

  this->m_check_collision_calls = 0;
  this->m_check_overlapping_calls = 0;

  this->m_total_check_collison_time = 0;
  this->m_total_check_overlapping_time = 0;
}

Unfolder::~Unfolder() {
  //delete m_m;//don't delete this, this is from the client
  //delete m_cd;
}

void Unfolder::measureModel() {
  int count = 0;

  for (int i = 0; i < this->m_m->v_size; ++i)
    if (this->m_m->vertices[i].hyperbolic)
      ++count;

  cout << " - hyperbolic vertices = " << count << " ( "
      << (count * 100.0 / this->m_m->v_size) << " % )" << endl;

  count = 0;
  for (int i = 0; i < this->m_m->e_size; ++i)
    if (this->m_m->edges[i].type == 'r')
      ++count;

  cout << " - concave edges = " << count << " ( "
      << (count * 100.0 / this->m_m->e_size) << " % )" << endl;
}

void Unfolder::buildUnfolding() {
  cout << "- unfolding..." << endl;

  auto start = clock();

  // create a splitter
  this->m_spliiter = Splitter::createSplitter(m_config.heuristic);
  this->m_spliiter->measure(this->m_m);

  int max_tries = m_config.max_retries;

  int min_overlaps = INT_MAX;
  int max_overlaps = INT_MIN;
  int total_overlaps = 0;
  float avg_overlaps = 0.0f;

  for (int i = 0; i < max_tries; i++) {
    auto weights = this->m_spliiter->assignWeights(this->m_m, this->m_config);
    auto count = this->buildFromWeights(weights);
    cout << "- total overlaps = " << count << endl;

    min_overlaps = min(min_overlaps, count);
    max_overlaps = max(max_overlaps, count);

    total_overlaps += count;
    avg_overlaps = (float) total_overlaps / (i + 1);

    if (!count) {
      cout << "- no overlapping found in " << (i + 1) << " tries" << endl;
      cout << "- base face = " << this->m_base_face_id << endl;
      this->m_is_flattened = true;
      break;
    }

    cout << "--------------------------------------------" << endl;
  }

  if (this->m_is_flattened && this->m_config.find_best_base_face) {
    cout << "- finding best base face" << endl;

    this->m_base_face_id = this->findBestBaseFace();

    // update base face
    this->m_config.baseface = this->m_base_face_id;

    this->rebuildTree(m_base_face_id, this->m_selected_edges);
  }

  // redo everything....

  this->initUnfold();

  this->alignModel();

  this->computeUnfolding();

  this->rebuildModel();

  cout << "- total time " << (float) (clock() - start) / CLOCKS_PER_SEC << " s"
      << endl;

  cout << " - min/avg/max overlaps = " << min_overlaps << "/" << avg_overlaps
      << "/" << max_overlaps << endl;

  cout << " - total fat edges = " << m_flat_edges << endl;
}

int Unfolder::buildFromWeights(const string& path) {

  cout << " - Building Unfolding From Weights File = " << path << endl;

  this->m_weights.clear();

  ifstream fin(path, ifstream::in);

  if (!fin.good()) {
    cerr << " ! Error ! Can't open weights file " << path << endl;
    return INT_MAX;
  }

  float weight;

  while (!fin.eof()) {
    fin >> weight;
    m_weights.push_back(weight);
  }

  // last line is blank
  m_weights.pop_back();

  if (m_weights.size() != this->m_m->e_size) {
    cerr << " ! Error ! Weights size is incorrect! Excepted = " << m_m->e_size
        << " Actual = " << m_weights.size() << endl;
    return INT_MAX;
  }

  // build the model from weights
  return this->buildFromWeights(this->m_weights);
}

int Unfolder::buildFromWeights(const vector<float>& weights,
    bool checkOverlaps) {

  this->m_weights = weights;

  GRAPH g;
  vector<int> parents;

  auto start = clock();

  this->initUnfold();
  this->alignModel();

  // assign weights on dual graph
  for (auto i = 0; i < this->m_m->e_size; i++) {
    const auto& edge = m_m->edges[i];
    const auto fid1 = edge.fid[0];
    const auto fid2 = edge.fid[1];

    const auto concave_vetex = m_m->vertices[edge.vid[0]].concave
        || m_m->vertices[edge.vid[1]].concave;

    // assign weight, edge connect a concave vertex will be penalized
//    const float weight = weights[i] * (concave_vetex ? 0.2 : 1.0);
    // set the weight
    g[fid1][fid2] = g[fid2][fid1] = weights[i];
  }

  this->buildMST(g);
  this->computeUnfolding();

  int count = INT_MAX;
  if (checkOverlaps) {
    count = this->checkOverlap();
  }
  this->m_is_flattened = (count == 0);

  if (!m_config.quite)
    cout << "- total time " << (float) (clock() - start) / CLOCKS_PER_SEC
        << " s" << endl;

  return count;
}

void Unfolder::unfoldTo(const vector<double>& folding_angles) {

  this->m_unfolded = this->m_org;

  vector<Matrix4x4> unfolding_matrax(m_m->t_size);

  int c = 0;

  for (auto fid : m_ordered_face_list) {

    // no need to unfold base_face;
    if (fid == this->m_base_face_id)
      continue;

    // faces of tabs
    if (fid >= this->m_m->t_size)
      break;

    const auto pfid = m_parents[fid];				// get parent face id
    assert(pfid >= 0 && pfid < this->m_m->t_size);

    const auto& f = this->m_m->tris[fid];				// current_face
    const auto& pf = this->m_m->tris[pfid];				// parent_face
    const auto eid = this->m_shared_edge[fid][pfid];	// shared edge id
    const auto& e = this->m_m->edges[eid];				// shared edge

    this->m_m->edges[eid].cutted = false;

    // rotation axis order matters, find the correct order
    const auto tid = this->isEdgeCCW(pfid, e.vid[0], e.vid[1]) ? 0 : 1;

    const auto& ev1 = this->getUnfoldedVertex(pfid, e.vid[tid]);// get edge vertex1
    const auto& ev2 = this->getUnfoldedVertex(pfid, e.vid[1 - tid]);// get edge vertex2

    const auto folding_angle = folding_angles[eid];

    assert(!std::isnan(folding_angle));

    // compute unfolding matrix
    unfolding_matrax[fid] = Matrix4x4::getRotationMatrix(ev1, ev2,
        folding_angle) * unfolding_matrax[pfid];

    // unfold current face
    for (auto k = 0; k < 3; k++) {
      assert(!std::isnan(this->m_unfolded[fid][k].second[0]));
      this->m_unfolded[fid][k].second = (unfolding_matrax[fid]
          * this->m_unfolded[fid][k].second.from3dto4d()).from4dto3d();
      assert(!std::isnan(this->m_unfolded[fid][k].second[0]));
    }
  }
}

void Unfolder::linearUnfoldTo(double percentage, const int count) {

  vector<double> folding_angles(m_m->e_size);

  for (auto fid : m_ordered_face_list) {

    // no need to unfold base_face;
    if (fid == this->m_base_face_id)
      continue;

    // faces of tabs
    if (fid >= this->m_m->t_size)
      break;

    const auto pfid = m_parents[fid];				// get parent face id
    assert(pfid >= 0 && pfid < this->m_m->t_size);

    const auto& f = this->m_m->tris[fid];				// current_face
    const auto& pf = this->m_m->tris[pfid];				// parent_face
    const auto eid = this->m_shared_edge[fid][pfid];	// shared edge id
    const auto& e = this->m_m->edges[eid];				// shared edge

    const auto folding_angle = -e.folding_angle * percentage;

    assert(!std::isnan(folding_angle));

    folding_angles[eid] = folding_angle;
  }

  this->unfoldTo(folding_angles);
}

void Unfolder::orderedUnfoldTo(double percentage) {

  int count = this->m_m->t_size * percentage;

  double ppf = (1.0 / this->m_m->t_size);

  double left_percentage = (percentage - ppf * count) / ppf;

  vector<double> folding_angles(m_m->e_size);

  auto face_list = m_ordered_face_list;

  std::reverse(face_list.begin(), face_list.end());

  int c = 0;

  for (auto fid : face_list) {
    ++c;
    // no need to unfold base_face;
    if (fid == this->m_base_face_id)
      continue;

    // faces of tabs
    if (fid >= this->m_m->t_size)
      break;

    const auto pfid = m_parents[fid];				// get parent face id
    assert(pfid >= 0 && pfid < this->m_m->t_size);

    const auto& f = this->m_m->tris[fid];				// current_face
    const auto& pf = this->m_m->tris[pfid];				// parent_face
    const auto eid = this->m_shared_edge[fid][pfid];	// shared edge id
    const auto& e = this->m_m->edges[eid];				// shared edge

    auto folding_angle = -e.folding_angle;

    // fully unfolded expect the last one
    if (c > count)
      folding_angle *= left_percentage;

    assert(!std::isnan(folding_angle));

    folding_angles[eid] = folding_angle;

    if (c > count)
      break;
  }

  this->unfoldTo(folding_angles);

}

void Unfolder::unfoldTo(double percentage) {
  if (percentage < 0)
    percentage = 0.0f;
  if (percentage > 1)
    percentage = 1.0f;

  if (m_config.ordered_unfolding) {
    this->orderedUnfoldTo(percentage);
  } else {
    this->linearUnfoldTo(percentage);
  }
}

void Unfolder::shrink() {
  for (auto i = 0; i < this->m_m->t_size; ++i) {
    auto& p0 = this->m_unfolded[i][0].second;
    auto& p1 = this->m_unfolded[i][1].second;
    auto& p2 = this->m_unfolded[i][2].second;

    // TODO use circumcenter instead
    auto center = (p0 + p1 + p2) / 3;

    p0 = (p0 - center) * m_config.shrink_factor + center;
    p1 = (p1 - center) * m_config.shrink_factor + center;
    p2 = (p2 - center) * m_config.shrink_factor + center;
  }
}

void Unfolder::rebuildModel() {
  if (!m_config.quite) {
    cout << "- rebuilding the model...";
    cout.flush();
  }

  // <new_fid, org_vid>, new_vid>
  map<pair<uint, uint>, uint> mapping;

  auto start = clock();

  this->m_eemap.clear();
  this->m_vemap.clear();
  this->m_vvmap.clear();
  this->m_vfmap.clear();
  this->m_vs.clear();
  this->m_cs.clear();
  this->m_fs = vector<vector<int> >(m_m->t_size);

  // loop each original vertex
  for (auto org_vid = 0; org_vid < m_m->v_size; org_vid++) {
    const auto& org_v = m_m->vertices[org_vid];

    vector<Vector3d> pos_list;
    vector<uint> idx_list;

    for (const auto fid : org_v.m_f) {
      const auto& new_f = this->m_unfolded[fid];

      Vector3d new_v_pos;
      for (auto k = 0; k < 3; k++) {
        if (new_f[k].first == org_vid) {
          new_v_pos = new_f[k].second;
        }
      }

      auto v_idx = -1;

      for (auto idx = 0; idx < pos_list.size(); idx++) {
        // same coordinates
        if ((new_v_pos - pos_list[idx]).norm() < 1e-6) {
          v_idx = idx_list[idx];
          break;
        }
      }

      // not in the list, new vertex
      if (v_idx < 0) {
        pos_list.push_back(new_v_pos);
        v_idx = m_vs.size();
        idx_list.push_back(v_idx);
        m_vs.push_back(new_v_pos);
      }

      // mapping from original fid, vid to new_vid
      mapping[std::make_pair(fid, org_vid)] = v_idx;

      this->m_vvmap[org_vid].insert(v_idx);

      // push new vertex index
      m_fs[fid].push_back(v_idx);
    } // end for, loop face around one vertex
  }  // end for, loop original vertex

  // loop each face to check normal
  for (auto fid = 0; fid < m_fs.size(); fid++) {
    const auto& new_f = m_fs[fid];

    assert(new_f.size() == 3);

    auto v0 = m_vs[new_f[0]];
    auto v1 = m_vs[new_f[1]];
    auto v2 = m_vs[new_f[2]];

    auto e1 = v1 - v0;
    auto e2 = v2 - v1;
    auto n = e1 % e2;

    if (n[1] < 0) {
      m_fs[fid] = {new_f[2], new_f[1], new_f[0]};
    }

    this->m_vfmap[m_fs[fid][0]].insert(fid);
    this->m_vfmap[m_fs[fid][1]].insert(fid);
    this->m_vfmap[m_fs[fid][2]].insert(fid);
  }

  // build creases
  for (auto fid : m_ordered_face_list) {
    // no need to unfold base_face
    if (fid == m_base_face_id)
      continue;

    const auto pfid = m_parents[fid];				// get parent face id
    const auto& f = this->m_m->tris[fid];				// current_face
    const auto& pf = this->m_m->tris[pfid];				// parent_face
    const auto eid = this->m_shared_edge[fid][pfid];	// shared edge id
    const auto& e = this->m_m->edges[eid];				// shared edge

    // rotation axis order matters, find the correct order
    const auto tid = this->isEdgeCCW(pfid, e.vid[0], e.vid[1]) ? 0 : 1;

    auto cvid1 = mapping[std::make_pair(fid, e.vid[tid])];
    auto cvid2 = mapping[std::make_pair(fid, e.vid[1 - tid])];

    this->m_cs.push_back(Crease(cvid1, cvid2, fid, pfid, e.folding_angle));

    auto cid = m_cs.size() - 1;

    // cvid1 <-> cvid2 is a crease line
    m_crease_lines[make_pair(cvid1, cvid2)] = cid;
    m_crease_lines[make_pair(cvid2, cvid1)] = cid;
  }

  // build ve/ev map
  for (auto fid = 0; fid < this->m_m->t_size; ++fid) {
    const auto& face = this->m_m->tris[fid];
    for (auto i = 0; i < 3; ++i) {
      //const auto& e = this->m_m->edges[face.e[i]];
      auto vid1 = face.v[i];
      auto vid2 = face.v[(i + 1) % 3];

      //const auto pfid = m_parents[fid];                   // get parent face id

//	        // rotation axis order matters, find the correct order
//	        const auto tid = this->isEdgeCCW(pfid, e.vid[0], e.vid[1]) ? 0 : 1;

      int new_vid1 = mapping[std::make_pair(fid, vid1)];
      int new_vid2 = mapping[std::make_pair(fid, vid2)];

      this->m_vemap[std::make_pair(new_vid1, new_vid2)] = face.e[i];
      this->m_vemap[std::make_pair(new_vid2, new_vid1)] = face.e[i];

//      cout<<new_vid1<<"/"<<new_vid2<<"->"<<face.e[i]<<endl;

      this->m_eemap[face.e[i]].insert(std::make_pair(new_vid1, new_vid2));
    }
  }

  if (this->m_is_flattened && this->m_config.find_boundary) {
    // find polygonal boundary for output.
    //this->findBoundary();

    // find a single path that go though all the crease lines.
    //this->findSinglePath();
  }

  // measure the net size
  auto min_v = Vector3d(FLT_MAX, FLT_MAX, FLT_MAX);
  auto max_v = Vector3d(-FLT_MAX, -FLT_MAX, -FLT_MAX);

  // find the left top corner
  for (const auto& v : m_vs)
    for (int i = 0; i < 3; ++i) {
      min_v[i] = std::min(min_v[i], v[i]);
      max_v[i] = std::max(max_v[i], v[i]);
    }

  this->m_net_size = max_v - min_v;

  // add tabs only when the unfolding is a net
  if (this->m_is_flattened && this->m_config.add_tabs) {
    this->addTabs();
  }

  ///////////////////////////////////////////////////
  // translate and scale the net
  ///////////////////////////////////////////////////

  // shift and scale the vertices

  // translation2 ...
  this->m_translation2 = -min_v;

  // translate the scale
  for (auto& v : m_vs) {
    v = (v - min_v) * this->m_config.scale;
  }

  this->m_net_size = this->m_net_size * this->m_config.scale;

  if (!m_config.quite)
    cout << " done in " << (float) (clock() - start) / CLOCKS_PER_SEC << " s"
        << endl;
}

/// dump unfolded results to svg file to the given path
void Unfolder::dumpSVG(const string& path, const int type) {
  ofstream out(path);

  auto min_x = DBL_MAX;
  auto max_x = DBL_MIN;
  auto min_y = DBL_MAX;
  auto max_y = DBL_MIN;

  for (const auto& v : this->m_vs) {
    min_x = min(min_x, v[0]);
    max_x = max(max_x, v[0]);
    min_y = min(min_y, v[2]);
    max_y = max(max_y, v[2]);
  }

  const auto scale_factor = m_config.scale;

  const auto width = (max_x - min_x);
  const auto height = (max_y - min_y);
  const auto aspect = width / height;

  const int width_pixel = width;
  const int height_pixel = width_pixel / aspect;

  const auto dashed_length = width_pixel * 0.01;

  const auto stroke_width = 1.0 * max(width_pixel, height_pixel) / 800;

  const auto font_size = min(width_pixel, height_pixel) * 0.02
      * m_config.label_font_scale;

  const string dashed_line = "stroke-dasharray:" + std::to_string(dashed_length)
      + ", " + std::to_string(dashed_length) + ";";

  const string boundary_style = "stroke:rgb(0,0,0);stroke-width:"
      + std::to_string(stroke_width * 2) + ";fill:white";
  const string normal_style = "stroke:rgb(128,128,128);stroke-width:"
      + std::to_string(stroke_width);
  const string mountain_style = "stroke:rgb(255,0,0);stroke-width:"
      + std::to_string(stroke_width) + ";"; // + dashed_line;

  const string valley_style = "stroke:rgb(0,0,255);stroke-width:"
      + std::to_string(stroke_width) + ";";

  const string tree_style = "stroke:rgb(0,128,0);stroke-width:"
      + std::to_string(stroke_width) + ";";

  const string crease_style = "stroke:rgb(128,128,128);stroke-width:"
      + std::to_string(stroke_width) + ";" + dashed_line;

  const string cut_style = "stroke:rgb(0,0,0);stroke-width:"
      + std::to_string(stroke_width) + ";" + dashed_line;
  //
//      stroke-dasharray:"
//      + std::to_string(dashed_length) + ", " + std::to_string(dashed_length)
//      + ", " + std::to_string(dashed_length / 5) + ", "
//      + std::to_string(dashed_length) + ";";

  out << "<?xml version=\"1.0\"?>" << endl;
  out << "<svg width=\"" << width_pixel << "\" height=\"" << height_pixel
      << "\"" << " version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" >"
      << endl;

  out << "<defs>" << endl;
  out << "<style type=\"text/css\"><![CDATA[" << endl;

  out << ".m {" << mountain_style << "}" << endl;
  out << ".v {" << valley_style << "}" << endl;
  out << ".b {" << boundary_style << "}" << endl;
  out << ".n {" << normal_style << "}" << endl;
  out << ".t {" << tree_style << "}" << endl;
  out << ".c {" << crease_style << "}" << endl;

  out << "]]></style>" << endl;
  out << "</defs>" << endl;

  // dump cuts, in two colors
  if (type & 1) {
    // draw boundary in black color
    out << "<path fill=\"none\" style=\"stroke:black;stroke-width:"
        << std::to_string(stroke_width * 2) << ";\"";
    out << " d=\"";

    set<pair<int, int>> draw_edges;
    set<set<int>> processed_edges;

    for (auto& face : this->m_fs) {
      for (auto j = 1; j <= 3; ++j) {
        const auto vid0 = face[j - 1];
        const auto vid1 = face[j % 3];
        const auto& p0 = this->m_vs[vid0];
        const auto& p1 = this->m_vs[vid1];

        bool boundary = true;

        auto edge = make_pair(vid0, vid1);
        auto edge2 = make_pair(vid1, vid0);

        if (!this->m_crease_lines.count(edge)
            && !this->m_crease_lines.count(edge2)) {

          out << " M " << p0[0] << " " << p0[2];
          out << " L " << p1[0] << " " << p1[2];
        }
      }

      // draw extra cuts for tabs
      if (type & 8) {

        const int vids[3] = { face[0], face[1], face[2] };

        const Vector3d p[3] = { m_vs[vids[0]], m_vs[vids[1]], m_vs[vids[2]] };
        const Vector3d e[3] = { p[1] - p[0], p[2] - p[1], p[0] - p[2] };

        const auto shortest_e = std::min(std::min(e[0].norm(), e[1].norm()),
            e[2].norm());

        const auto x = shortest_e * m_config.extra_cuts_x;
        const auto y = shortest_e * m_config.extra_cuts_y;

        for (int i = 0; i < 3; ++i) {
          int vidc = vids[i];
          int vidn = vids[(i + 1) % 3];
          set<int> edge = { vidc, vidn };

          if (processed_edges.count(edge))
            continue;
          processed_edges.insert(edge);

          const auto ep = e[(i + 2) % 3];
          const auto ec = e[i];
          const auto en = e[(i + 1) % 3];

          const auto tp = acos(en.normalize() * ep.normalize());
          const auto tc = acos(ec.normalize() * ep.normalize());
          const auto tn = acos(ec.normalize() * en.normalize());

          const auto h = y * sin(tn);
          const auto b = y * sin(tp) / sin(tn);
          const auto a = ec.norm() - b - x;
          const auto z = y * sin(tn) / sin(tc);

          const auto pbx = p[(i + 1) % 3] - (ec.normalize()) * (b + x);
          const auto pbp = pbx + en.normalize() * y;
          const auto pz = p[i] - ep.normalize() * z;

          out << " M " << pbx[0] << " " << pbx[2];
          out << " L " << pbp[0] << " " << pbp[2];
          out << " L " << pz[0] << " " << pz[2];
          out << " L " << p[i][0] << " " << p[i][2];
          out << " L " << pbx[0] << " " << pbx[2];

          auto e01 = make_pair(face[i], face[(i + 1) % 3]);
          auto e10 = make_pair(face[(i + 1) % 3], face[i]);

          // check whether is boundary edge
          if (!this->m_crease_lines.count(e01)
              && !this->m_crease_lines.count(e10))
            continue;

          // draw mirrow part
          // height dir: perpendicular to current edge, ignore y...
          const auto h_dir = Vector3d(-ec[2], 0, ec[0]);
          const auto pbpp = pbp + h_dir.normalize() * 2 * h;
          const auto pbzp = pz + h_dir.normalize() * 2 * h;

          out << " L " << pbpp[0] << " " << pbpp[2];
          out << " L " << pbzp[0] << " " << pbzp[2];
          out << " L " << p[i][0] << " " << p[i][2];

        }
      } // end of if (type & 8) extra cuts
    }

    out << "\" />" << endl;

    {
      // draw creases in dotted grey as a single path
      out << "<path fill=\"none\" class=\"c\" d=\"";

      for (auto& c : this->m_cs) {
        const auto& p1 = this->m_vs[c.vid1];
        const auto& p2 = this->m_vs[c.vid2];

        // skip flat edges
        if (c.folding_angle == 0)
          continue;

        out << " M " << p1[0] << " " << p1[2];
        out << " L " << p2[0] << " " << p2[2];
      }

      out << "\" />" << endl;
    } // end of creases

    // draw spanning tree
    if (type & 4) {
      for (auto& c : this->m_cs) {
        int fid = c.fid;
        int pid = c.pid;
        if (pid == -1)
          continue;
        const auto& f1 = this->m_fs[fid];
        const auto& f2 = this->m_fs[pid];
        Vector3d c1;
        Vector3d c2;
        for (int i = 0; i < 3; ++i) {
          c1 += this->m_vs[f1[i]];
          c2 += this->m_vs[f2[i]];
        }
        c1 = c1 * (1.0 / 3);
        c2 = c2 * (1.0 / 3);

        char buf[1024];

        sprintf(buf,
            "<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" class=\"t\" />",
            c1[0], c1[2], c2[0], c2[2]);
        out << buf << endl;
      }
    } // end of spanning tree

  } // end of type 1

  // display mode
  // 1. draw crease and boundary in their own colors
  // 2. draw cuts with different thickness
  if (type == 2) {
    set<pair<int, int>> draw_edges;

    const double min_line_width = 0.2;
    const double max_line_width = 3.0;

    char buf[1024];

    // draw crease lines in color
    for (auto& face : this->m_fs) {
      for (auto j = 1; j <= 3; ++j) {
        int vid1 = face[j - 1];
        int vid2 = face[j % 3];
        const auto p0 = this->m_vs[vid1];
        const auto p1 = this->m_vs[vid2];

//        const auto edge = std::make_pair(min(vid1, vid2), max(vid1,vid2));
        const auto edge1 = std::make_pair(vid1, vid2);
        const auto edge2 = std::make_pair(vid2, vid1);

        int org_eid = this->m_vemap[edge1];
        const auto& org_e = this->m_m->edges[org_eid];

        string class_name = "b";
        bool boundary = true;

        // already drawn
        if (draw_edges.count(edge1) || draw_edges.count(edge2))
          continue;

        // const Crease& crease = this->m_cs[m_crease_lines.at(edge1)];
        if (org_e.folding_angle > 1e-3)
          class_name = "m";
        else if (org_e.folding_angle < -1e-3)
          class_name = "v";
        else
          class_name = "n";

        const double line_width = stroke_width
            * (min_line_width
                + fabs(org_e.folding_angle / PI)
                    * (max_line_width - min_line_width));

        const string style = "stroke-width:" + std::to_string(line_width) + ";";

        boundary = false;

        draw_edges.insert(edge1);
        draw_edges.insert(edge2);

        sprintf(buf,
            "<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" class=\"%s\" style=\"%s\" />",
            p0[0], p0[2], p1[0], p1[2], class_name.c_str(), style.c_str());
        out << buf << endl;
      }
    }

    // draw labels
    if (this->m_config.dump_labels) {
      char text_buf[1024];

      for (auto i = 0; i < this->m_m->e_size; ++i) {
        const auto& new_edges = this->m_eemap.at(i);

        // boundary edges
        if (new_edges.size() == 1u)
          continue;

        for (const auto& new_edge : new_edges) {
          if (this->m_crease_lines.count(new_edge))
            continue;

          const auto p0 = m_vs[new_edge.first];
          const auto p1 = m_vs[new_edge.second];

          const auto vec = (p1 - p0).normalize();
          const auto angle = RadToDeg(atan2(vec[2], vec[0])); // angle from x axis in degree

          const auto edge_length = (p1 - p0).norm();
          int label_size = std::to_string(i).size();
          double label_length = label_size * font_size * 0.4;
          double pp = min(label_length / 2 / edge_length, 0.5);
          const auto center = p0 + (p1 - p0) * (0.5 - pp);

          sprintf(text_buf,
              "<text x=\"%f\" y=\"%f\" transform=\"rotate(%f %f %f)\" fill=\"darkgreen\" font-weight=\"bold\" font-size=\"%f\">%d</text>",
              center[0], center[2],   // x, y
              angle, center[0], center[2], font_size, i);

          out << text_buf << endl;

        } // end for

      } // end for

    } // end if
  }

  out << "</svg>" << endl;

  out.close();

  cout << "- dumped svg to " << path << endl;
}

// =============================================================
// private
// =============================================================

void Unfolder::buildDualGraph(GRAPH& g) {
  auto start = clock();

  if (!m_config.quite) {
    cout << "- building dual graph...";
    cout.flush();
  }

  g.clear();

  auto weights = this->m_spliiter->assignWeights(this->m_m, this->m_config);

  if (!m_config.quite)
    cout << " done in " << (float) (clock() - start) / CLOCKS_PER_SEC << " s"
        << endl;
}

void Unfolder::buildMST(GRAPH& g) {
  auto start = clock();

  priority_queue<DualGraphEdge, vector<DualGraphEdge>, greater<DualGraphEdge>> q;

  if (!m_config.quite) {
    cout << "- building MST...";
    cout.flush();
  }

  m_parents.clear();
  m_selected_edges.clear();

  this->m_ordered_face_list.clear();

  m_parents.resize(m_m->t_size);

  auto in_the_tree = set<uint>();

  // base_face is in the tree
  in_the_tree.insert(m_base_face_id);
  m_ordered_face_list.push_back(m_base_face_id);
  m_parents[m_base_face_id] = -1;

  // insert all edges start from base face to the queue
  for (const auto& de : g[m_base_face_id]) {
    q.push(DualGraphEdge(m_base_face_id, de.first, de.second));
  }

  float total_weight = 0;

  const int edges_not_to_cut = this->m_flat_edges;

  // total t_size - 1 edges in the tree
  for (auto k = 1u; k < m_m->t_size; k++) {
    //auto e = q.top();
    DualGraphEdge e;

    // find the maximum weight edge with one node in the tree, but another not
    //while (true)
    while (q.empty() == false) {
      e = q.top();
      q.pop();

      // both node are already in the tree
      // not valid, remove that edge
      if (in_the_tree.count(e.fid1) && in_the_tree.count(e.fid2)) {
        continue;
      }

      break;
    }

    if (e.fid1 < 0 || e.fid2 < 0) {
      break;
    }
//    cout << "e = " << e.fid1 << "," << e.fid2 << " w = " << e.weight
//        << " selected" << endl;

    // add directed edge to selected edge set
    this->m_selected_edges.insert(make_pair(e.fid1, e.fid2));

    total_weight += e.weight;

    if (in_the_tree.count(e.fid1) + in_the_tree.count(e.fid2) == 2) {
      cerr << "!!!Warning!!!! both faces are already in the tree" << endl;
    }

    auto s = in_the_tree.count(e.fid1) ? e.fid1 : e.fid2;
    auto t = s == e.fid1 ? e.fid2 : e.fid1;

    m_parents[t] = s;

    in_the_tree.insert(t);
    this->m_ordered_face_list.push_back(t);

    // add new edges
    for (const auto& de : g[t]) {

      // if the other face is already in the tree
      if (in_the_tree.count(de.first))
        continue;

      q.push(DualGraphEdge(t, de.first, de.second));
    }
  }

  //if (this->m_ordered_face_list.size() != this->m_m->t_size)
  //{
  //	cout << "this->m_ordered_face_list.size()=" << this->m_ordered_face_list.size() << " this->m_m->t_size=" << this->m_m->t_size << endl;
  //}
  //assert(this->m_ordered_face_list.size() == this->m_m->t_size);

  if (!m_config.quite) {
    cout << " total weight = " << total_weight << endl;
    cout << " done in " << (float) (clock() - start) / CLOCKS_PER_SEC << " s"
        << endl;
  }
}

void Unfolder::rebuildTree(int base_face, set<pair<int, int>> selected_edges) {

  this->m_parents.clear();
  this->m_ordered_face_list.clear();

  this->m_parents.resize(m_m->t_size);

  auto in_the_tree = set<uint>();
  in_the_tree.insert(base_face);
  this->m_ordered_face_list.push_back(base_face);
  this->m_parents[base_face] = -1;
  this->m_m->tris[base_face].path_len = 0;
  this->m_max_path_len = -1;
  this->m_avg_path_len = -1;

  auto sum_path_len = 0.0;

  while (selected_edges.size()) {
    set<pair<int, int>> to_remove;
    for (const auto& e : selected_edges) {
      if (!in_the_tree.count(e.first) && !in_the_tree.count(e.second))
        continue;

      int pf = in_the_tree.count(e.first) ? e.first : e.second;
      int f = pf == e.first ? e.second : e.first;

      in_the_tree.insert(f);
      this->m_parents[f] = pf;
      this->m_ordered_face_list.push_back(f);

      this->m_m->tris[f].path_len = this->m_m->tris[pf].path_len + 1;

      this->m_max_path_len = max(this->m_max_path_len,
          this->m_m->tris[f].path_len);

      sum_path_len += this->m_m->tris[f].path_len;

      to_remove.insert(e);
    }

    for (const auto& e : to_remove) {
      selected_edges.erase(e);
    }
  }

  assert(this->m_ordered_face_list.size() == this->m_m->t_size);

  this->m_avg_path_len = sum_path_len / this->m_m->t_size;

  cout << "  - base face = " << base_face << " avg_path_length = "
      << m_avg_path_len << endl;
}

int Unfolder::findBestBaseFace() {
  int best_base_face = -1;
  double min_avg_path_len = FLT_MAX;

  for (int i = 0; i < m_m->t_size; ++i) {
    this->rebuildTree(i, this->m_selected_edges);
    if (m_avg_path_len < min_avg_path_len) {
      min_avg_path_len = m_avg_path_len;
      best_base_face = i;
    }
  }

  this->m_base_face_id = best_base_face;

  this->rebuildTree(m_base_face_id, this->m_selected_edges);

  this->alignModel();

  this->computeUnfolding();

  cout << "  - best base face = " << best_base_face << " avg_path_length = "
      << min_avg_path_len << endl;

  return best_base_face;
}

void Unfolder::initUnfold() {
  auto start = clock();

  if (!m_config.quite) {
    cout << "- Initializing unfolding" << endl;
    cout << "  - vertices = " << m_m->v_size		// number of vertices
        << " edges = " << (m_m->e_size - m_m->e_boundary_size)// number of non-bounday edges / total edges
        << "/" << m_m->e_size << " faces = " << m_m->t_size << endl;// number of total facet
  }

// use random base face
  if (this->m_config.random_baseface)
    m_base_face_id = (int) (mathtool::drand48() * m_m->t_size);
  else if (this->m_config.baseface >= 0) {
    if (this->m_config.baseface >= this->m_m->t_size) {
      std::cerr << "Error! base face out of range!" << endl;
      exit(-1);
    }
    m_base_face_id = this->m_config.baseface;
  } else {
    m_base_face_id = 0;
  }

  if (!m_config.quite)
    cout << "  - base face = " << m_base_face_id << endl;

  this->m_org.clear();

// copy vertex coordinates
  for (auto i = 0; i < this->m_m->t_size; i++) {
    const auto& f = this->m_m->tris[i];
    const auto vs = this->m_m->vertices;
    auto v0 = std::make_pair(f.v[0], Vector3d(vs[f.v[0]].p.get()));
    auto v1 = std::make_pair(f.v[1], Vector3d(vs[f.v[1]].p.get()));
    auto v2 = std::make_pair(f.v[2], Vector3d(vs[f.v[2]].p.get()));

    m_org.push_back( { v0, v1, v2 });
  }

// copy to unfolded
  m_unfolded = m_org;

// reset flat edges
  m_flat_edges = 0;

// shared edges
  m_shared_edge.clear();
  for (auto i = 0; i < m_m->e_size; i++) {
    const auto& edge = m_m->edges[i];
    const auto fid1 = edge.fid[0];
    const auto fid2 = edge.fid[1];
    m_shared_edge[fid1][fid2] = m_shared_edge[fid2][fid1] = i;

    if (edge.folding_angle == 0)
      m_flat_edges++;
  }

  if (!m_config.quite)
    cout << "  - Done in " << (float) (clock() - start) / CLOCKS_PER_SEC << " s"
        << endl;
}

void Unfolder::alignModel() {
// align based_face's normal with Y axis
  auto n0 = this->m_m->tris[m_base_face_id].n.normalize();
  // work around
  auto y_axis = Vector3d(0, 1, 0).normalize();
  this->m_rotation_angle = acos(n0 * y_axis);
  this->m_rotation_axis = n0 % y_axis;

// no rotation required
  if (this->m_rotation_angle == 0) {
    this->m_rotation_axis = Vector3d(0, 1, 0);
  } else if (fabs((fabs(this->m_rotation_angle) - PI)) < 1e-6) {
    this->m_rotation_axis = Vector3d(1, 0, 0);
    this->m_rotation_angle = PI;
  }

  auto r_matrix = Matrix4x4::getRotationMatrix(Vector3d(0, 0, 0),
      this->m_rotation_axis, this->m_rotation_angle);

  for (auto& f : this->m_org) {
    for (auto k = 0; k < 3; k++) {
      f[k].second = (r_matrix * f[k].second.from3dto4d()).from4dto3d();
    }
  }

  const auto v0 = this->m_org[m_base_face_id][0].second;

  this->m_transilation = -v0;

  for (auto& f : this->m_org) {
    for (auto k = 0; k < 3; k++) {
      f[k].second = f[k].second + m_transilation;
    }
  }

  this->m_unfolded = this->m_org;
}

void Unfolder::computeUnfolding() {
  if (!m_config.quite) {
    cout << "- Computing unfolding...";
    cout.flush();
  }

  auto start = clock();

  this->linearUnfoldTo(1.0);

  if (!m_config.quite)
    cout << " Done in " << (float) (clock() - start) / CLOCKS_PER_SEC << " s"
        << endl;
}

const Vector3d& Unfolder::getUnfoldedVertex(uint fid, uint org_vid) const {
  const auto& f = this->m_unfolded[fid];
  for (const auto& pair : f) {
    if (pair.first == org_vid)
      return pair.second;
  }

  assert(false);

// remove warning...
  return this->m_unfolded[0][0].second;
}

bool Unfolder::isEdgeCCW(uint fid, uint vid1, uint vid2) const {
  for (auto k = 0; k < 3; k++) {
    const auto fvid1 = m_m->tris[fid].v[k];
    const auto fvid2 = m_m->tris[fid].v[(k + 1) % 3];

    if (fvid1 == vid1 && fvid2 == vid2)
      return true;
  }

  return false;
}

// find the boundary of the net
void Unfolder::findBoundary() {
// clear the result
  this->m_boundary.clear();

  unordered_set<int> attached_vertices;

// an edge from vid1 - vid2
  int vid1 = INT_MAX;
  int vid2 = INT_MAX;

  for (int fid = 0; fid < this->m_m->t_size; ++fid) {
    for (int i = 0; i < 3; ++i) {
      const auto s = this->m_fs[fid][i];
      const auto t = this->m_fs[fid][(i + 1) % 3];

      const auto edge = make_pair(s, t);

      if (!this->m_crease_lines.count(edge)) {
        // we find the first edge
        vid1 = s;
        vid2 = t;
        break;
      }
    }

    // found an boundary edge
    if (vid1 != INT_MAX && vid2 != INT_MAX) {
      cout << " - initial boundary edge = " << vid1 << " -> " << vid2
          << " at face = " << fid << endl;
      break;
    }
  }

  m_boundary.push_back(vid1);
  m_boundary.push_back(vid2);

  attached_vertices.insert(vid1);
  attached_vertices.insert(vid2);

// we haven't find a loop
  while (vid2 != vid1) {
    bool found = false;
    // loop all adjacent faces of vid2
    for (auto fid : m_vfmap[vid2]) {
      for (int i = 0; i < 3; ++i) {
        auto s = this->m_fs[fid][i];
        auto t = this->m_fs[fid][(i + 1) % 3];

        // next edge must started with vid2
        if (s != vid2)
          continue;

        auto edge = make_pair(s, t);

        // it's the boundary
        if (!this->m_crease_lines.count(edge)) {
          found = true;

          // we find the edge
          m_boundary.push_back(t);
          attached_vertices.insert(t);

          vid2 = t;

          break;
        }
      }

      if (found)
        break;
    }
  }

  if (!m_config.quite) {
    cout << "boundary:" << endl;
    for (auto vid : m_boundary)
      cout << vid << " ";
    cout << endl;
  }
}

// find single path for all the creases
void Unfolder::findSinglePath() {
  this->m_single_path.clear();

  int vid1, vid2;

  bool found = false;

// find the start vertex and start edge
  for (auto i = 0; !found && i < this->m_vs.size(); ++i) {
    // loop all adjacent faces
    for (auto fid : this->m_vfmap[i]) {
      for (auto k = 0; k < 3; ++k) {
        auto s = this->m_fs[fid][k];
        auto t = this->m_fs[fid][(k + 1) % 3];

        // an edge must start with vid1
        if (s != i)
          continue;

        auto edge = make_pair(s, t);

        if (this->m_crease_lines.count(edge)) {
          vid1 = s;
          vid2 = t;
          found = true;
          break;
        }
      } // end for k
    } // end for fid
  } // end for i

  unordered_set<int> vistied_vs;
  set<pair<int, int>> visited_edges;

  this->findSinglePath(vid1, vistied_vs, visited_edges, this->m_single_path);

  cout << "single path = " << endl;
  for (auto vid : m_single_path)
    cout << vid << " ";
  cout << endl;
}

// travel from vid
void Unfolder::findSinglePath(int vid, unordered_set<int>& visited_vids,
    set<pair<int, int>>& visited_edges, vector<int>& path) {
// append current node to the path
  path.push_back(vid);

// mark vid as visited
  visited_vids.insert(vid);

// loop adjacent faces
  for (const auto fid : this->m_vfmap[vid]) {
    // loop all edges
    for (auto k = 0; k < 3; ++k) {
      const auto s = this->m_fs[fid][k];
      const auto t = this->m_fs[fid][(k + 1) % 3];

      // an edge must start with vid
      if (s != vid)
        continue;

      auto edge = make_pair(s, t);

      // is a crease line
      if (this->m_crease_lines.count(edge)) {
        // direct goto if edge not visited but vertex visited
        if (visited_vids.count(t) && !visited_edges.count(edge)) {
          visited_edges.insert(edge);
          path.push_back(t);
        }
        // recursive goto the vid
        else {
          visited_edges.insert(edge);
          this->findSinglePath(t, visited_vids, visited_edges, path);
        }

        // move back
        path.push_back(-vid - 1);
      }
    } // end for k
  } // end for fid
}

bool Unfolder::getIntersectionOfPlane(const int fid1, const Vector3d& v1,
    const Vector3d& v2, Vector3d& outPoint) {
  const double EP = 1e-3;

// the normal of the edge
  Vector3d vl = (v1 - v2).normalize();

// get 3 vertices of the face
  Vector3d p0 = this->m_unfolded[fid1][0].second;
  Vector3d p1 = this->m_unfolded[fid1][1].second;
  Vector3d p2 = this->m_unfolded[fid1][2].second;

// two vector on the face
  Vector3d vf0 = (p1 - p0);
  Vector3d vf1 = (p2 - p0);

// connect the end points of the edge to one point on the plane
  Vector3d vl1 = (v1 - p0);
  Vector3d vl2 = (v2 - p0);

// normal vector of the face
  Vector3d fno = (vf0 % vf1).normalize();

  double dp1 = vl1 * fno;
  double dp2 = vl2 * fno;

// line segment are on the same side from the plane
  if (sign(dp1) == sign(dp2) || (fabs(dp1) < EP && fabs(dp2) < EP))
    return false;

  double dot_prodcut = (vl * fno);

  if (fabs(dot_prodcut) < EP) {
    //the line is parallel to the plane
    //???
    return false;
  } else {
    double t = -vl1 * fno / dot_prodcut;
    if (fabs(t) < EP)
      return false;

    //calculate the intersection point
    outPoint[0] = v1[0] + vl[0] * t;
    outPoint[1] = v1[1] + vl[1] * t;
    outPoint[2] = v1[2] + vl[2] * t;

    Vector3d pOut = Vector3d(outPoint[0], outPoint[1], outPoint[2]);

    for (uint i = 0; i < 3; i++) {
      Vector3d p = this->m_unfolded[fid1][i].second;
      if ((p - pOut).norm() < EP) {
        return false;
      }
    }
    return true;
  }
}

bool Unfolder::pointInTriangle(const int fid1, const Vector3d& p) {
  Vector3d p0 = this->m_unfolded[fid1][0].second;
  Vector3d p1 = this->m_unfolded[fid1][1].second;
  Vector3d p2 = this->m_unfolded[fid1][2].second;

  Vector3d v0 = (p1 - p0);
  Vector3d v1 = (p2 - p0);
  Vector3d v2 = (p - p0);

  double dot00 = v0 * v0;
  double dot01 = v0 * v1;
  double dot02 = v0 * v2;
  double dot11 = v1 * v1;
  double dot12 = v1 * v2;

  double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
  double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
  double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

// Check if point is in triangle
  if ((u >= 0) && (v >= 0) && (u + v < 1)) {
    return true;
  }

  return false;
}

bool Unfolder::hasIntersection(const int fid1, const int fid2) {
  bool intersection = false;

  const auto face1 = this->m_m->tris[fid1];
  const auto face2 = this->m_m->tris[fid2];

  Vector3d p;

  for (uint i = 1; i <= 3; i++) {
    const Vector3d& v1 = this->m_unfolded[fid2][i - 1].second;
    const Vector3d& v2 = this->m_unfolded[fid2][i % 3].second;

    // get the intersection of the line contains the edge and the plane contains the face
    if (this->getIntersectionOfPlane(fid1, v1, v2, p)) {
      // check whether the intersection is in the face
      if (this->pointInTriangle(fid1, p)) {
        intersection = true;
      }
    }
  }

  return intersection;
}

#include "itree/RectKD.h"
#include "itree/MiddleStructure.h"

uint Unfolder::checkOverlap_itree() {
  const int F = this->m_unfolded.size();
  this->m_unfolded;

  typedef Interval<EndPoint> Int;
  typedef RectKD<EndPoint, 2> Rect2D;
  typedef MiddleStructure<Rect2D, 2> MTree;
  typedef vector<Rect2D*>::iterator IT2;

  vector<Rect2D*> rects;

  for (int i = 0; i < F; i++) {

    float max_x = -FLT_MAX, min_x = FLT_MAX, max_y = -FLT_MAX, min_y =
    FLT_MAX;

    for (int p = 0; p < 3; p++) {
      const auto& vi1 = m_unfolded[i][p].second;
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
  for (auto rect : rects) {
    uint i = rect->getVID();
    Tree.query(rect);
    Rect2D::Intersect& intersections = rect->getIntersections();

    for (uint j : intersections) {
      if (i > j)
        continue; //this will be checked by (j,i) pair
      if (checkOverlap(i, j)) {
        count++;
      }
    }
  }

  //free mem
  for (auto rect : rects)
    delete (rect);

  return count;
}

bool Unfolder::checkOverlap(uint i, uint j) {
  bool intersected = false;
  double pp[2] = { 0.0, 0.0 };

  // two faces shared an edge in the unfolding (not in the original mesh).
  if (this->m_selected_edges.count(make_pair(i, j))
      || this->m_selected_edges.count(make_pair(j, i)))
    return intersected;

  for (int p = 0; p < 3 && (!intersected); p++) {

    const auto& vi1 = m_unfolded[i][p].second;
    const auto& vi2 = m_unfolded[i][(p + 1) % 3].second;

    const double a[2] = { vi1[0], vi1[2] };
    const double b[2] = { vi2[0], vi2[2] };

    for (int q = 0; q < 3 && (!intersected); q++) {

      const auto& vj1 = m_unfolded[j][q].second;
      const auto& vj2 = m_unfolded[j][(q + 1) % 3].second;

      if ((vi1 == vj1 && vi2 == vj2) || (vi1 == vj2 && vi2 == vj1))
        continue;

      const double c[2] = { vj1[0], vj1[2] };
      const double d[2] = { vj2[0], vj2[2] };

      char r = SegSegInt<double>(a, b, c, d, pp);

      if (r == '1') {
        intersected = true;
      }
    }
  }

  if (intersected) {
    this->m_m->tris[i].overlapped = true;
    this->m_m->tris[j].overlapped = true;
    if (this->m_config.record_overlap) {
      this->m_overlap_pairs[i].push_back(j);
      this->m_overlap_pairs[j].push_back(i);
    }
  }

  return intersected;
}

int Unfolder::checkOverlap() {

  ++this->m_check_overlapping_calls;

  const int F = this->m_unfolded.size();

  if (!m_config.quite) {
    cout << "- checking overlapping... ";
    cout.flush();
  }

// backup the current status
  auto unfoled_back = this->m_unfolded;

// shrink the triangles only for checking overlapping...
  if (this->m_config.shrink) {
    this->shrink();
  }

  auto start = clock();

  uint count = 0;

#define USE_PIXEL_CHECKER 0

#if USE_PIXEL_CHECKER
  masc::unfolding::PixelChecker checker(this);
  double ratio = checker.checkOverlapping(this->m_unfolded, this->m_config);
  auto pixel_time = (clock() - start) * 1.0 / CLOCKS_PER_SEC;
#endif

  if (this->m_config.record_overlap) {
    this->m_overlap_pairs.clear();
    m_overlap_pairs.resize(F);
  }

  for (int i = 0; i < F; i++)
    this->m_m->tris[i].overlapped = false;

  //if (!m_config.use_rapid) {

    bool checked = false;
    if (F > 100) {
      count = checkOverlap_itree();
      if (count != UINT_MAX)
        checked = true; //good
    }

    if (checked == false) {
      for (int i = 0; i < F; i++) {
        for (int j = i + 1; j < F; j++) {
          if (checkOverlap(i, j)) {
            count++;
          }
        } //end for j
      } //end for i
    } //end if (checked == false)

  //} else {
  //  //use rapid
  //  //count = this->m_cd->hasCollision();
  //}

  const auto total_time = (clock() - start) * 1.0 / CLOCKS_PER_SEC;

#if USE_PIXEL_CHECKER
  cout << "ratio = " << ratio << " overlaps = " << count << " time = "
  << pixel_time << " s total time = " << total_time << " s" << endl;
#endif

  this->m_total_check_overlapping_time += (clock() - start);

  if (!m_config.quite)
    cout << " count = " << count << " Done in " << total_time << " s" << endl;

// recovery
  this->m_unfolded.swap(unfoled_back);

  return count;
}

// count local overlaps due to insufficient count for hyperbolic vertex
int Unfolder::checkLocalOverlaps() {
  int overlaps = 0;
  int concave_vertex = 0;

  for (auto i = 0; i < this->m_m->v_size; ++i) {
    const auto& v = this->m_m->vertices[i];

    // only consider hyperbolic vertices
    if (!v.hyperbolic)
      continue;

    concave_vertex++;

    int cuts = 0;
    for (const auto eid : v.m_e) {
      const auto& e = this->m_m->edges[eid];

      const auto pair1 = std::make_pair(e.fid[0], e.fid[1]);
      const auto pair2 = std::make_pair(e.fid[1], e.fid[0]);

      if (!this->m_selected_edges.count(pair1)
          && !this->m_selected_edges.count(pair2))
        cuts++;
    }

    if (cuts < 2)
      overlaps++;
  }

  return overlaps;
}

// collision is detected based on penetration
bool Unfolder::checkCollision() {
  const int F = this->m_unfolded.size();

  if (!m_config.quite) {
    cout << "- checking collision... ";
    cout.flush();
  }

  auto start = clock();

  int count = 0;

  for (int i = 0; i < F; i++)
    this->m_m->tris[i].overlapped = false;

  for (int i = 0; i < F; i++) {
    for (int j = 0; j < F; j++) {
      if (i == j)
        continue;

      // two faces shared an edge in the unfolding (not in the original mesh).
      if (this->m_selected_edges.count(make_pair(i, j))
          || this->m_selected_edges.count(make_pair(j, i)))
        continue;

      bool intersected = this->hasIntersection(i, j);

      if (intersected) {
        ++count;
        this->m_m->tris[i].overlapped = true;
        this->m_m->tris[j].overlapped = true;
      }
    }
  }

  if (!m_config.quite)
    cout << " count = " << count << " done in "
        << (float) (clock() - start) / CLOCKS_PER_SEC << " s" << endl;

  return count > 0;
}

// add tabs on the net
void Unfolder::addTabs() {
// add a tab for each cut edge in the original model
// tab width = edge width
// tab height = edge width * 15%
// tab angle = try [45, 30, 15]

// steps:
// 1. find cut/boundary edges to add tab
// 2. add two triangle for each edge

  int count_added = 0;
  int count_total = 0;

  for (const auto& ee : m_eemap) {
    set<int> vids;
    for (const auto pair : ee.second) {
      vids.insert(pair.first);
      vids.insert(pair.second);
    }

    bool boundary_edge = this->m_m->edges[ee.first].fid.size() == 1;

    if (boundary_edge) {
      cout << "boundary edge = " << ee.first << endl;
    }

    // a cut edge has at least 3 vertices in the net
    if (vids.size() >= 3 || boundary_edge) {
      // find an cut edge
      ++count_total;

      if (this->addTabForCutEdge(ee.first))
        ++count_added;
    }
  }

  cout << "\n\n";
  cout << "Tabs added = " << count_added << "/" << count_total << endl;
}

int Unfolder::getNetFidByVids(const int vid1, const int vid2) {
  set<int> vids = { vid1, vid2 };

  for (int i = 0; i < this->m_fs.size(); ++i) {
    const auto face = this->m_fs[i];
    int count = 0;

    for (int j = 0; j < 3; ++j) {
      if (vids.count(face[j]))
        ++count;
    }

    if (count == 2)
      return i;
  }

  return -1;
}

/// add one tab for each cut edge with original edge id
bool Unfolder::addTabForCutEdge(const int org_eid) {
  for (const auto& e : this->m_eemap[org_eid]) {
    const int vid1 = e.first;
    const int vid2 = e.second;
    const int fid = getNetFidByVids(vid1, vid2);

    if (this->addTabForNewEdge(vid1, vid2, fid)) {
      return true;
    }
  }

  cout << "can't add tab on to edge " << org_eid << endl;

  return false;
}

/// add one tab on an edge of the net
bool Unfolder::addTabForNewEdge(const int vid1, const int vid2, const int fid) {
// see: https://www.dropbox.com/s/16yrtfxnkomubwi/mesh_unfolder_tab.png?dl=0

// get the coordinates of the edge for adding tab
  const auto& p1 = this->m_vs[vid1];
  const auto& p2 = this->m_vs[vid2];
// vector of the edge
  auto edge = p2 - p1;
// length of the edge
  const auto edge_len = edge.norm();
// normalize the edge
  edge = edge.normalize();

//TODO, should also be based on the size of the entire net
  const auto min_net_size = min(m_net_size[0], m_net_size[2]);
  const auto tab_heights = { min_net_size * 0.05, min_net_size * 0.04,
      min_net_size * 0.03 };

  const vector<double> angles = { 45.0, 30.0, 20.0, 10.0 };

// try different combinations...
  for (const auto tab_height : tab_heights) {
    for (const auto angle3 : angles) {
      for (const auto angle4 : angles) {
        const auto m3 = Matrix4x4::getRotationMatrixY(DegToRad(-angle3));
        const auto m4 = Matrix4x4::getRotationMatrixY(DegToRad(angle4));
        const auto len3 = tab_height / sin(DegToRad(angle3));
        const auto len4 = tab_height / sin(DegToRad(angle4));

        // |p3p4|
        const auto len34 = edge_len
            - (tab_height / tan(DegToRad(angle3))
                + tab_height / tan(DegToRad(angle4)));

        if (len34 < 0 || len34 * 2 < edge_len)
          continue;

        // compute the vertices of the tab
        const auto p3 = (m3 * edge.from3dto4d()).from4dto3d() * len3 + p1;
        const auto p4 = (m4 * (-edge).from3dto4d()).from4dto3d() * len4 + p2;

        //check collision, latter
        if (this->isOverlapWithNet(p1, p3, p4, fid)
            || this->isOverlapWithNet(p1, p4, p2, fid)) {
          continue;
        }

        // add tab to the net...
        // 1. add vertices
        // 2. add creases
        // 3. add faces
        // 3. update ordered face list

        // create new vertices
        const int vid3 = this->m_vs.size();
        const int vid4 = this->m_vs.size() + 1;

        this->m_vs.push_back(p3);
        this->m_vs.push_back(p4);

        // create new faces
        int fid134 = this->m_fs.size();
        int fid142 = this->m_fs.size() + 1;

        // face 1 3 4
        this->m_fs.insert(this->m_fs.end(), { vid1, vid3, vid4 });

        // face 1 4 2
        this->m_fs.insert(this->m_fs.end(), { vid1, vid4, vid2 });

        // create new crease lines
        int cid12 = this->m_cs.size();
        int cid14 = this->m_cs.size() + 1;

        auto c12 = Crease(vid1, vid2, fid142, fid, PI);
        auto c14 = Crease(vid1, vid4, fid134, fid142, 0);
        this->m_cs.push_back(c12);
        this->m_cs.push_back(c14);

        this->m_crease_lines[make_pair(vid1, vid2)] = cid12;
        this->m_crease_lines[make_pair(vid1, vid4)] = cid14;

        this->m_ordered_face_list.push_back(fid142);
        this->m_ordered_face_list.push_back(fid134);

        return true;

      }   // for angle3
    }   // for angle4
  }   // for heights

  return false;
}

/// check whether a triangle overlap with existing net (included attached tabs)
bool Unfolder::isOverlapWithNet(const Vector3d& p1, const Vector3d& p2,
    const Vector3d& p3, const int fid) {
// output, intersection location
  double pp[2] = { 0.0, 0.0 };

  const auto tab = vector<Vector3d> { p1, p2, p3 };

  int cur_fid = 0;
  for (const auto& f : this->m_fs) {
    // do not test against parent face...
    if (cur_fid++ == fid)
      continue;

    for (int i = 1; i <= 3; ++i) {
      const auto& vi1 = this->m_vs[f[i - 1]];
      const auto& vi2 = this->m_vs[f[i % 3]];

      for (int j = 1; j <= 3; ++j) {
        const auto& vj1 = tab[j - 1];
        const auto& vj2 = tab[j % 3];

        if ((vi1 == vj1 && vi2 == vj2) || (vi1 == vj2 && vi2 == vj1))
          continue;

        const double a[2] = { vi1[0], vi1[2] };
        const double b[2] = { vi2[0], vi2[2] };

        const double c[2] = { vj1[0], vj1[2] };
        const double d[2] = { vj2[0], vj2[2] };

        char r = SegSegInt<double>(a, b, c, d, pp);

        if (r == '1') {
          cout << "has overlap!!!!" << endl;
          return true;
        }
      }
    }
  }

  return false;
}
// =================================================================
// public access
// =================================================================

model * Unfolder::getModel() const {
  return this->m_m;
}

const Config& Unfolder::getConfig() const {
  return this->m_config;
}

const MESH& Unfolder::getUnfolded() const {
  return this->m_unfolded;
}

const string& Unfolder::getFilename() const {
  return this->m_config.filename;
}

const Vector3d& Unfolder::getRotationAxis() const {
  return this->m_rotation_axis;
}

/// angle in degree
const double Unfolder::getRotationAngle() const {
  return this->m_rotation_angle * 180.0 / PI;
}

const Vector3d& Unfolder::getTranslation() const {
  return this->m_transilation;
}

/// get model color
const Vector3d& Unfolder::getColor() const {
  return this->m_color;
}

// set model color
void Unfolder::setColor(double r, double g, double b) {
  this->m_color.set(r, g, b);
}

const bool Unfolder::isFlattened() const {
  return this->m_is_flattened;
}
