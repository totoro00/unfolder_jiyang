//------------------------------------------------------------------------------
//  Copyright 2007-2009 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#include "model.h"
#include "ModelGraph.h"
#include "mathtool/Geometry.h"

#ifdef _WIN32
#pragma warning(disable:4244)
#endif

//-----------------------------------------------------------------------------
// two global scope models, P, Q
model P, Q;
model& getP() {
  return P;
}
model& getQ() {
  return Q;
}
//-----------------------------------------------------------------------------

//
// build the model from a file
// the file will contain an obj model
// once obj file is read, vertices and facets will be built
// other information associated with facets and vertices,
// as well as edges will be build also using model graph
//

bool model::build(const string & name) {
  MeshReader* reader = nullptr;

  if (name.find(".obj") == name.length() - 4) {
    reader = new objReader();
  } else if (name.find(".off") == name.length() - 4) {
    reader = new offReader();
  }

  assert(reader != nullptr);

  if (!reader->Read(name)) {
    delete reader;
    return false;
  }

  objModel& data = reader->getModel();

  bool result = build(data);

  delete reader;

	{//model name
		int start = name.find_last_of("/");
		if (start == string::npos) start = name.find_last_of("\\") + 1;
		if (start == string::npos) start = 0;
		int end = name.find_last_of(".");
		if (end == string::npos) end = name.length();
		this->name = name.substr(start, end - start);
	}

  return result;
}

bool model::build(const vector<Vector3d>& vertices,
    const vector<vector<int>>& faces) {
  objModel data(vertices, faces);

  return this->build(data);
}

bool model::build(const objModel& data) {

  //allocate memory
  v_size = data.pts.size();
  t_size = data.polys.size();

  vertices = new vertex[v_size];   //
  tris = new triangle[t_size];     //
  assert(vertices && tris);        //make sure enough memory

//copy vertices
  for (uint i = 0; i < v_size; i++) {
    vertices[i].p.set(&data.pts[i].x);
    vertices[i].bk_p = vertices[i].p;  //backup for modification
  }

  cout << "vertices copied!" << endl;

  //copy triangles
  int tid = 0;
  for (auto i = data.polys.begin(); i != data.polys.end(); i++) {
    const auto& ids = i->pts;
    //check if triangle
    if (ids.size() != 3) {
      cerr << "! Error: polygon " << tid << " is not a triangle, edge size is "
          << ids.size() << "." << endl;
      return false;
    }
    int vid = 0;
    for (auto j = ids.begin(); j != ids.end(); j++) {
      assert(*j >= 0 && *j < v_size);
      tris[tid].v[vid++] = *j;
      vertices[*j].m_f.push_back(tid);
    }
    tid++;
  }

  cout << "triangle copied!" << endl;

  {  //build model graph and collect informations
    this->e_boundary_size = 0;

    CModelGraph G;
    G.doInit(this);
    //create edges
    e_size = G.getEdgeSize();
    CModelEdge * ptrE = G.getEdges();
    edges = new edge[e_size];
    assert(edges);
    for (uint i = 0; i < e_size; i++, ptrE = ptrE->getNext()) {
      int v1 = edges[i].vid[0] = ptrE->getStartPt();
      int v2 = edges[i].vid[1] = ptrE->getEndPt();

      const vector<int>&tmp_fid = ptrE->getFacets();

      edges[i].fid.insert(edges[i].fid.end(), tmp_fid.begin(), tmp_fid.end());

      if (tmp_fid.size() < 2) { //check if boundary edge
        //edges[i].fid[1]=edges[i].fid[0];
        edges[i].fid.push_back(edges[i].fid[0]);
        edges[i].type = 'b';		//bd
        this->e_boundary_size++;
      }

      //compute parallel vector
      edges[i].v = edges[i].bk_v = ptrE->getV();
      //edges[i].in_n[0]=ptrE->getInNormal(0); edges[i].bk_in_n[0]=edges[i].in_n[0];
      //edges[i].in_n[1]=ptrE->getInNormal(1); edges[i].bk_in_n[1]=edges[i].in_n[1];
      vertices[v1].m_e.push_back(i);
      vertices[v2].m_e.push_back(i);
    }            //end i

//facets
    vector<CModelFacet>& facets = G.getFacets();
    for (uint i = 0; i < t_size; i++) {
      tris[i].n = tris[i].bk_n = facets[i].n;
      for (int j = 0; j < 3; j++) {
        tris[i].e[j] = facets[i].m_Edge[j]->getID();
      }            //end j

      // compute face area
      double ss = 0.0;
      double s[3]; // edge length
      for(int j = 1; j <= 3; ++j) {
        s[j-1] = (vertices[tris[i].v[j%3]].p - vertices[tris[i].v[j-1]].p).norm();
        ss += s[j-1];
        tris[i].center += Vector3d(vertices[tris[i].v[j%3]].p.get());
      } // end j

      // store center of the face
      tris[i].center = tris[i].center / 3.0;

      ss /= 2.0;

      tris[i].area = ss;
      for(int j = 0; j < 3; ++j) {
        tris[i].area *= (ss - s[j]);
      } // end j

      tris[i].area = std::sqrt(tris[i].area);


      this->surface_area += tris[i].area;
    }            //end i

//edge type
    ptrE = G.getEdges();
    for (uint i = 0; i < e_size; i++, ptrE = ptrE->getNext()) {
      if (edges[i].type == 'b')
        continue; //already know
      Vector3d& n1 = tris[edges[i].fid[0]].n;
      Vector3d& n2 = tris[edges[i].fid[1]].n;
      double d = n1 * n2;
      if (fabs(1 - d) < SMALLNUMBER) {
        edges[i].type = 'p'; //plane
        edges[i].folding_angle = 0;
      } else {
        double angle = acos(d);
        Vector3d vec = (n1 % n2).normalize();
        if (vec * edges[i].v > 0) {
          edges[i].type = 'c'; //convex
          edges[i].folding_angle = angle;
        } else {
          edges[i].type = 'r'; //reflex
          edges[i].folding_angle = -angle;
        }
      }

      // check whether the edge is diagonal edge or not
      const auto& t1 = tris[edges[i].fid[0]];

      const auto& p11 = vertices[t1.v[0]].p;
      const auto& p12 = vertices[t1.v[1]].p;
      const auto& p13 = vertices[t1.v[2]].p;

      if (!mathtool::isRightTriangle(p11, p12, p13))
        continue;

      const auto& t2 = tris[edges[i].fid[1]];

      const auto& p21 = vertices[t2.v[0]].p;
      const auto& p22 = vertices[t2.v[1]].p;
      const auto& p23 = vertices[t2.v[2]].p;

      if (!mathtool::isRightTriangle(p21, p22, p23))
        continue;

      double area1 = mathtool::triangleArea(p11, p12, p13);
      double area2 = mathtool::triangleArea(p21, p22, p23);

      if (fabs(area1 - area2) < SMALLNUMBER && edges[i].type == 'p') {
        edges[i].type = 'd';  // diagonal
        edges[i].diagonal = true;
      }

    }

//vertex type
    typedef list<uint>::iterator IT;
    for (uint i = 0; i < v_size; i++) {
      //int convex_e=-1;
      vertex& v = vertices[i];

      Vector3d evec;
      for (IT ie = v.m_e.begin(); ie != v.m_e.end(); ie++) {
        edge& e = edges[*ie];
        Vector3d vec = e.v;
        if (e.vid[1] == i)
          vec = -vec;
        evec = evec + vec;
      } //end ir

      Vector3d fvec;
      for (IT f = v.m_f.begin(); f != v.m_f.end(); f++) {
        triangle& t = tris[*f];
        fvec = fvec + t.n;
      } //end ir

      if (evec * fvec > 0)
        v.concave = true;

      double sum_angles = 0;

      for (const auto fid : v.m_f) {
        const auto& f = this->tris[fid];

        int vid1 = -1;
        int vid2 = -1;

        for (auto vid : f.v) {
          if (vid == i)
            continue;
          if (vid1 == -1)
            vid1 = vid;
          else
            vid2 = vid;
        }

        Vector3d e1 =
            (this->vertices[vid1].p - this->vertices[i].p).normalize();
        Vector3d e2 =
            (this->vertices[vid2].p - this->vertices[i].p).normalize();

        double angle = acos(e1 * e2);

        sum_angles += angle;
      }

      if (sum_angles > 2 * PI) {
        v.hyperbolic = true;
      }

    } //end i
  }

  cout << "- v_size = " << v_size << " t_size = " << t_size << " e_size = "<< e_size << endl;

  this->compute_COM_R();

  cout << "- COM = " << this->COM << " R = " << this->R << endl;

  return true;
}

//
//
//
// Transformation operations
// 
// Rotation
// Negation
// Reverse (turn inside out)
//

void model::perturb(double noise) {

  //build the matrix
  Vector3d r(noise * mathtool::drand48(), noise * mathtool::drand48(),
      noise * mathtool::drand48());
  for (int i = 0; i < 3; i++) {
    if (mathtool::drand48() > 0.5)
      r[i] = -r[i];
  }

  Quaternion q(r.get());
  Matrix3x3 M = q.getMatrix();

  //rotate edges
  for (uint i = 0; i < e_size; i++) {
//edges[i].in_n[0]=(M*edges[i].in_n[0]).normalize();
//edges[i].in_n[1]=(M*edges[i].in_n[1]).normalize();
    edges[i].v = (M * edges[i].v).normalize();
  }
  //rotate facets
  for (uint i = 0; i < t_size; i++)
    tris[i].n = (M * tris[i].n).normalize();
}

void model::unperturb() {
  //build the matrx
  Matrix3x3 M(current_rot[0][0], current_rot[0][1], current_rot[0][2],
      current_rot[1][0], current_rot[1][1], current_rot[1][2],
      current_rot[2][0], current_rot[2][1], current_rot[2][2]);

  //rotate edges
  for (uint i = 0; i < e_size; i++) {
//edges[i].in_n[0]=(M*edges[i].bk_in_n[0]).normalize();
//edges[i].in_n[1]=(M*edges[i].bk_in_n[1]).normalize();
    edges[i].v = (M * edges[i].bk_v).normalize();
  }
  //rotate facets
  for (uint i = 0; i < t_size; i++)
    tris[i].n = (M * tris[i].bk_n).normalize();
}

void model::rotate(const Matrix2x2& M) {
  Vector2d tmp;

  //rotate vertices
  for (uint i = 0; i < v_size; i++) {
    tmp.set(vertices[i].bk_p[0], vertices[i].bk_p[1]);
    tmp = M * tmp;
    vertices[i].p.set(tmp[0], tmp[1]);
  }

  //rotate edges
  for (uint i = 0; i < e_size; i++) {
//for(int j=0;j<2;j++){
//    tmp.set(edges[i].bk_in_n[j][0],edges[i].bk_in_n[j][1]);
//    tmp=M*tmp;
//    edges[i].in_n[j].set(tmp[0],tmp[1]);
//}

    tmp.set(edges[i].v[0], edges[i].v[1]);
    tmp = M * tmp;
    edges[i].v.set(tmp[0], tmp[1]);
  }

  //rotate facets
  for (uint i = 0; i < t_size; i++) {
    tmp.set(tris[i].n[0], tris[i].n[1]);
    tmp = M * tmp;
    tris[i].n.set(tmp[0], tmp[1]);
  }
}

void model::rotate(const Matrix3x3& M) {
  Vector3d tmp;

  //rotate vertices
  for (uint i = 0; i < v_size; i++) {
    tmp.set(vertices[i].bk_p.get());
    vertices[i].p = M * tmp;
  }

  //rotate edges
  for (uint i = 0; i < e_size; i++) {
//edges[i].in_n[0]=(M*edges[i].bk_in_n[0]).normalize();
//edges[i].in_n[1]=(M*edges[i].bk_in_n[1]).normalize();
    edges[i].v = (M * edges[i].bk_v).normalize();
  }
  //rotate facets
  for (uint i = 0; i < t_size; i++)
    tris[i].n = (M * tris[i].bk_n).normalize();
}

void model::scale(double s) {
  //rotate vertices
  Point3d tmp;
  for (uint i = 0; i < v_size; i++) {
    tmp.set(vertices[i].bk_p.get());
    vertices[i].p.set(tmp[0] * s, tmp[1] * s, tmp[2] * s);
    vertices[i].bk_p.set(tmp[0] * s, tmp[1] * s, tmp[2] * s);
  }
}

void model::negate() {

  for (uint i = 0; i < v_size; i++) {
    Point3d& pt = vertices[i].p;
    pt.set(-pt[0], -pt[1], -pt[2]);

    Point3d& pt2 = vertices[i].bk_p;
    pt2.set(-pt2[0], -pt2[1], -pt2[2]);
  }

  for (uint i = 0; i < t_size; i++) {
    tris[i].n = -tris[i].n;
    tris[i].bk_n = -tris[i].bk_n;
    swap(tris[i].v[1], tris[i].v[2]);
//swap(tris[i].e[1],tris[i].e[2]);
  }

  for (uint i = 0; i < e_size; i++) {
    edge& e = edges[i];
    e.v = -e.v;
//e.in_n[0]=-e.in_n[0];
//e.in_n[1]=-e.in_n[1];
    e.bk_v = -e.bk_v;
//e.bk_in_n[0]=-e.bk_in_n[0];
//e.bk_in_n[1]=-e.bk_in_n[1];
  }
}

void model::reverse() {
  for (uint i = 0; i < t_size; i++) {
    tris[i].n = -tris[i].n;
    tris[i].bk_n = -tris[i].bk_n;
    swap(tris[i].v[1], tris[i].v[2]);
  }

  for (uint i = 0; i < e_size; i++) {
    edge& e = edges[i];
    if (e.type == 'r')
      e.type = 'c';
    else if (e.type == 'c')
      e.type = 'r';
  }
}

//get the neighbors of triangle t
/*
 void model::get_neighbors(triangle * t, list<triangle *>& nei)
 {
 for(int i=0;i<3;i++){
 edge& e=edges[t->e[i]];
 uint o_t=(&(tris[e.fid[0]])==f)?e.fid[1]:e.fid[0];
 if(&tris[o_t]==t) continue;
 nei.push_back(&tris[o_t]);
 }
 }
 */

void model::compute_COM_R() {
  this->COM = Vector3d(0, 0, 0);

  for (int i = 0; i < this->v_size; ++i) {
    this->COM = this->COM + Vector3d(this->vertices[i].p.get());
  }

  this->COM = this->COM / this->v_size;

  double r = 0;

  for (int i = 0; i < this->v_size; ++i) {
    double d = (Vector3d(this->vertices[i].p.get()) - this->COM).normsqr();

    if (d > r)
      r = d;
  }

  this->R = std::sqrt(r);

}
