//------------------------------------------------------------------------------
//  Copyright 2007-2008 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _OBJ_READER_H_
#define _OBJ_READER_H_

#include <cctype>
#include <cmath>
#include <cassert>

#include <string>
#include <fstream>
#include <iostream>
#include <list>
#include <vector>

using namespace std;

#include <Point.h>
#include <Quaternion.h>
#include <Vector.h>
using namespace mathtool;

class Vpt {
public:
	Vpt() {
		x = y = z = nx = ny = nz = 0;
	}

	void normalize() {
		double norm = (double) sqrt(nx * nx + ny * ny + nz * nz);
		nx /= norm;
		ny /= norm;
		nz /= norm;
	}

	double x, y, z, nx, ny, nz;
};

class V {
public:
	V() {
		x = y = z = 0;
	}
	double x, y, z;
};

class polygon {
public:
	list<int> pts;
	list<int> normals;
};

class objModel {
public:
	objModel() {
	}

	// build from vertices and faces
	objModel(const vector<Vector3d>& vertices,
			const vector<vector<int>>& faces) {
		Vpt pt;
		//  vertices
		for (const auto& v : vertices) {
			pt.x = v[0];
			pt.y = v[1];
			pt.z = v[2];
			this->pts.push_back(pt);
		}

		// faces
		for (const auto& f : faces) {
			polygon ply;
			ply.pts.push_back(f[0]);
			ply.pts.push_back(f[1]);
			ply.pts.push_back(f[2]);
			this->polys.push_back(ply);
		}

		this->compute_v_normal();
	}

	void compute_v_normal() {
		//check if normal information is valid from the obj file
		if (normals.empty()) { //compute normal
			for (list<polygon>::iterator i = polys.begin(); i != polys.end();
					i++) {
				//get 3 points, compute normal and assign to all vertices
				list<int>::iterator pi = i->pts.begin();
				vector<Point3d> v3;
				for (; pi != i->pts.end(); pi++) {
					Vpt& pt = pts[*pi];
					Point3d pos(pt.x, pt.y, pt.z);
					v3.push_back(pos);
					if (v3.size() == 3)
						break; //we've collected 3 points
				}
				//compute normal
				Vector3d n = ((v3[1] - v3[0]) % (v3[2] - v3[0])).normalize();
				//copy normals
				pi = i->pts.begin();
				for (; pi != i->pts.end(); pi++) {
					Vpt& pt = pts[*pi];
					pt.nx += n[0];
					pt.ny += n[1];
					pt.nz += n[2];
				} //end copying normals
			} //end looping polygons
		} else { // use the information provided
			for (list<polygon>::iterator i = polys.begin(); i != polys.end();
					i++) {
				list<int>::iterator ni = i->normals.begin();
				list<int>::iterator pi = i->pts.begin();

				for (; pi != i->pts.end(); pi++, ni++) {
					V& n = normals[*ni];
					Vpt& pt = pts[*pi];
					pt.nx += n.x;
					pt.ny += n.y;
					pt.nz += n.z;
				} //end copying normals

			} //end looping polygons
		}

		//normalize
		for (vector<Vpt>::iterator i = pts.begin(); i != pts.end(); i++) {
			i->normalize();
		}
	}

	vector<Vpt> pts;
	vector<V> normals;
	list<polygon> polys;
};

class MeshReader {
public:
	virtual ~MeshReader() {
	}
	bool Read(const string& filename) {
		this->m_filename = filename;

		ifstream in(m_filename.c_str());
		if (!in.good()) {
			cerr << "Can't open file " << m_filename << endl;
			return false;
		}
		bool r = Read(in);
		in.close();
		return r;
	}
	const objModel& getModel() const {
		return m_data;
	}
	objModel& getModel() {
		return m_data;
	}

protected:
	MeshReader() {
	}
	virtual bool Read(istream& in)=0;

	string m_filename;
	objModel m_data;

};

class objReader: public MeshReader {
public:
	objReader() {
	}
	virtual ~objReader() {
	}
protected:
	virtual bool Read(istream& in) override {

		string tmp;

		//read pts
		while (true) {
			in >> tmp;
			if (tmp == "f")
				break;
			if (tmp == "v") {
				Vpt pt;
				in >> pt.x >> pt.y >> pt.z;
				m_data.pts.push_back(pt);
			} else if (tmp == "vn") {
				V pt;
				in >> pt.x >> pt.y >> pt.z;
				m_data.normals.push_back(pt);
			}
			getline(in, tmp);
		}

		//read faces
		polygon poly;
		do {

			in >> tmp;

			if (in.eof())
				break;

			if (isdigit(tmp[0])) { //this defines a vertex

				int pos1 = tmp.find('/');
				int pos2 = tmp.rfind('/');

				int id_v = atoi(tmp.substr(0, pos1).c_str()) - 1;
				int id_n = atoi(tmp.substr(pos2 + 1).c_str()) - 1;

				poly.pts.push_back(id_v);
				poly.normals.push_back(id_n);
			} else if (tmp == "f") {
				m_data.polys.push_back(poly);
				poly.pts.clear();
				poly.normals.clear();
			} else {
				getline(in, tmp);
			}

		} while (!in.eof());

		m_data.polys.push_back(poly);
		m_data.compute_v_normal();

		return true;
	}
};

class offReader: public MeshReader {
public:
	offReader() {
	}
	virtual ~offReader() {
	}
protected:
	virtual bool Read(istream& in) override {
		string tmp;
		in >> tmp;
		assert(tmp == "OFF");

		int vsize, fsize, tsize;

		// vertices, faces, unknown
		in >> vsize >> fsize >> tsize;

		// read vertices
		for (auto i = 0; i < vsize; ++i) {
			Vpt pt;
			in >> pt.x >> pt.y >> pt.z;
			m_data.pts.push_back(pt);
		}

		// read faces
		for (auto i = 0; i < fsize; ++i) {
			polygon ply;

			int vs, vid;
			in >> vs;
			assert(vs == 3);
			for (auto j = 0; j < 3; ++j) {
				in >> vid;
				ply.pts.push_back(vid);
			}

			m_data.polys.push_back(ply);
		}

		m_data.compute_v_normal();

		return true;
	}
};

#endif //_OBJ_READER_H_

