//------------------------------------------------------------------------------
//  Copyright 2007-2009 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _BF_MODEL_H_
#define _BF_MODEL_H_

#include <climits>

#include <Point.h>
#include <Vector.h>
#include <Matrix.h>
#include <Quaternion.h>
using namespace mathtool;

#include <string>
#include <cassert>
#include <set>
#include <unordered_map>
#include <utility>
using namespace std;

//#include "gmap.h"
#include "objReader.h"

typedef unsigned int uint;

//a triangle of the model
struct triangle {
	
	triangle() {
		cluster_id = 0;
		overlapped = false;
		parent_id = -1;
		source_fid = -1;
		path_len = 0.0;
		area = 0.0;
		score = 0.0;
	}

	uint v[3]; //vertex id
	uint e[3]; //edge id

	Vector3d n; //normal

	// source face id
	int source_fid;

	uint cluster_id; //id of the corresponding facet in kway-union.h

	//backups
	Vector3d bk_n;

	// for unfolding
	bool overlapped;

	// parent id, -1 means root face
	int parent_id;

	// path length from root
	double path_len;

	// area of the face
	double area;

	// center of mass
	Vector3d center;

	float score;
};

//a vertex of the model
struct vertex {
	vertex() {
		concave = false;
		hyperbolic = false;
		score = 0.0;
	}
	Point3d p;  //position
	list<uint> m_f;
	list<uint> m_e; //a list of edges

	//backups
	Point3d bk_p;

	//if concave, set to true
	bool concave;

	//whether the vertex is hyperbolic
	bool hyperbolic;

	// weighted score
	float score;
};

//an edge of the model
struct edge {
	edge() {
		type = 'x';
		vid[0] = vid[1] = UINT_MAX;
		folding_angle = 0.0;
		cutted = true;
		diagonal = false;
	}
	uint vid[2];
	vector<uint> fid;

	Vector3d v;       //parallel vector
	//Vector3d in_n[2]; //inface normals

	//backups
	Vector3d bk_v;       //parallel vector
	//Vector3d bk_in_n[2]; //inface normals

	//type, c-convex, r-reflex, p-plane, b-border
	char type;

	// added for origami
	double folding_angle;
	bool cutted;      // whether the edge is cut or not
	bool diagonal;    // whether the edge is a diagonal edge
};

struct model {
	//initialization
	model() {
		v_size = e_size = t_size = e_boundary_size = 0;
		vertices = NULL;
		edges = NULL;
		tris = NULL;
		//spheres=NULL;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				current_rot[i][j] = 0;
		current_rot[0][0] = current_rot[1][1] = current_rot[2][2] = 1;
		R = 1.0;

		surface_area = 0.0;
	}

	~model() {
	}

	void destroy() {
		delete[] tris;
		tris = NULL;
		delete[] edges;
		edges = NULL;
		delete[] vertices;
		vertices = NULL;
		v_size = e_size = t_size = 0;
	}

	//build from model file: OBJ/OFF
	bool build(const string& name);

	//build from obj_model
	bool build(const objModel& obj_model);

	//build from vertices and faces
	bool build(const vector<Vector3d>& vertices, const vector<vector<int>>& faces);

	//rotate points
	void rotate(const Matrix2x2& m);
	void rotate(const Matrix3x3& M);

	//scale the model
	void scale(double s);

	//negate point/facets ...
	void negate();

	//reverse facets ...
	void reverse();

	//perturb
	void perturb(double amount);
	void unperturb();

	// compute COM and R
	void compute_COM_R();

	//void get_neighbors(triangle * t, list<triangle *>& nei);

	//data
	vertex * vertices;  //vertices
	triangle * tris;      //triangles
	edge * edges;     //edges
	uint v_size;
	uint e_size;
	uint t_size;
	uint e_boundary_size; // size of boundary edges

	//current orientation
	double current_rot[3][3];

	Vector3d COM;
	double R;
	string name;

	// total surface area
	double surface_area;
};

typedef unordered_map<uint, unordered_map<uint, float> > GRAPH; // adjacency list of weight <v_i, <v_j, w_ij>>
typedef vector<vector<pair<uint, Vector3d> > > MESH; // coordinates faces[f_i][k] = <vid, coordinates> (k = 0..2)

#endif //_BF_MODEL_H_
