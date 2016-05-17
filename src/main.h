//------------------------------------------------------------------------------
//  Copyright 2007-2008 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _BF_MINKOWSKI_SUM_H_
#define _BF_MINKOWSKI_SUM_H_

//
//
//
// This is the main header file is shared by m+3d.cpp (single thread) and pm+3d.cpp
// (multi-thread). Functions defined in this file are mainly for text parsing and
// OpenGL rendering 
//
//
//
#include <list>
#include <cfloat>
using namespace std;

#include "objReader.h"
#include "model.h"
#include "unfolder.h"
#include "config.h"
#include "util/DataHelper.h"

//using namespace masc::unfolding;


//-----------------------------------------------------------------------------
// INPUTS
vector<string> filenames;
Config config;

//-----------------------------------------------------------------------------
// Intermediate data
vector<Unfolder*> unfolders;

double R = 0;       //radius
Point3d COM;     //center of mass

//-----------------------------------------------------------------------------
//for random rotation
double rot_theta;
Vector3d rot_vec;
Quaternion current_rot;

//-----------------------------------------------------------------------------
//read M+ from file
bool readfromfile();
void computeCOM_R();

//-----------------------------------------------------------------------------
//
//
//
//  Open GL stuff below
//
//
//-----------------------------------------------------------------------------

#include "draw.h"

//-----------------------------------------------------------------------------
bool parseArg(int argc, char ** argv) {
  if (argc == 1)
    return false;

  for (int i = 1; i < argc; i++) {
    auto const arg = string(argv[i]);
    if (arg == "-r") {
      config.max_retries = std::stoi(argv[++i]);
    } else if (arg == "-h") {
      auto const heuristic = string(argv[++i]);
      if (heuristic[0] == 's')
        config.heuristic = CutHeuristic::STEEPEST_EDGE;
      else if (heuristic[0] == 'f')
        config.heuristic = CutHeuristic::FLAT_TREE;
      else if (heuristic[0] == 'p')
        config.heuristic = CutHeuristic::MINIMUM_PERIMETER;
      else if (heuristic[0] == 'r')
        config.heuristic = CutHeuristic::RANDOM;
      else if (heuristic[0] == 'b')
        config.heuristic = CutHeuristic::BRUTE_FORCE;
	  else if (heuristic[0] == '1')
		  config.heuristic = CutHeuristic::MY_UNFOLDER_01;
	  else if (heuristic[0] == '2')
		  config.heuristic = CutHeuristic::MY_UNFOLDER_02;
      else {
        cerr << "!Error! Unknown heuristic type = " << heuristic << endl;
        return false;
      }
    } else if (arg == "-g") {
      config.disable_gui = true;
    } else if (arg == "-rb") {
      config.random_baseface = true;
    } else if (arg == "-bf" || arg == "--base-face") {
      config.baseface = std::stoi(argv[++i]);
    } else if (arg == "-lc") {
      config.less_cuts = true;
    } else if (arg == "-s") {
      config.seed = std::stoi(argv[++i]);
    } else if (arg == "-q") {
      config.quite = true;
    } else if (arg == "-ns") {
      config.shrink = false;
    } else if (arg == "-nfb") {
      config.find_boundary = false;
    } else if (arg == "-k") {
      config.k = std::stoi(argv[++i]);
    } else if (arg == "-i") {
      config.max_iterations = std::stoi(argv[++i]);
    } else if (arg == "-run") {
      config.run = std::stoi(argv[++i]);
    } else if (arg == "-weights") {
      config.weights_filename = argv[++i];
    } else if (arg == "-sf") {
      config.shrink_factor = stof(argv[++i]);
    } else if (arg == "-scale") {
      config.scale = stof(argv[++i]);
    } else if (arg == "-nl") {
      config.dump_labels = false;
	}
	else if (arg == "-lfs") {
		config.label_font_scale = stof(argv[++i]);
	}
	else if (arg == "-tab") {
		config.add_tabs = true;
	}
	else if (arg[0] == '-') {
      cerr << "!Error! Unknown arg = " << arg << endl;
      return false;
    } else {
      filenames.push_back(string(argv[i]));
    }
  }

  return true;
}

void printUsage(char * name) {

  Config default_config;

  //int offset=20;
  cerr << "Usage: " << name << " [options] *.obj/*.off \n\n";

  cerr << "Heuristic Methods\n";
  cerr << "  -h heuristic | use heuristic method\n";
  cerr << "      s        | STEEPEST_EDGE\n";
  cerr << "      f        | FLAT_TREE (default)\n";
  cerr << "      p        | MINIMUM_PERIMETER\n";
  cerr << "      r        | RANDOM\n";
  cerr << "      1        | Your Unfolder 01\n";
  cerr << "      2        | Your Unfolder 02\n";
  cerr << "\n";

  cerr << "Unfolding\n";
  cerr << "  -s seed      | specify random seed\n";

  cerr << "  -r times     | retry times, default is "
      << default_config.max_retries << "\n";

  cerr << "  -weights fn  | using the specify weights to unfold the mesh.\n";

  cerr << "  -q           | quite mode.\n";

  cerr << "  -bf          | specify the base face.\n";

  cerr << "  -rb          | random base face.\n";

  cerr << "  -lc          | less cuts / don't cut flat edges.\n";

  cerr << "  -g           | disable GUI. dump outputs.\n";

  cerr << "  -tab         | add tabs in the net.\n";

  cerr << "Dumping SVG\n";
  cerr
      << "  -scale       | scale factor applied. both *.svg and *.ori will be affected.\n";
  cerr << "  -nl          | do not dump labels in SVG file\n";
  cerr << "  -lfs         | label font size scale [default=1.0]\n";

  cerr << endl;
  cerr << "-- Complied on " << __DATE__ << endl;
  cerr << "-- Report bugs to: Jyh-Ming Lien jmlien@cs.gmu.edu" << endl;
}

//-----------------------------------------------------------------------------

bool readfromfiles() {
  long vsize = 0;
  long fsize = 0;

  int flattened = 0;

  uint id = 0;

  auto start = clock();

  for (const auto& filename : filenames) {

    // use the same seed for all the component
    mathtool::srand48(config.seed);
    srand(config.seed);

    cout << "- [" << ++id << "/" << filenames.size() << "] Start reading "
        << filename << endl;

    config.filename = filename;

    model* m = new model();
    if (!m->build(config.filename)) {
      delete m;
      return false;
    }
    cout << "- Done reading " << m->v_size << " vertices and " << m->t_size << " facets" << endl;

    Unfolder* unfolder = new Unfolder(m, config);
    unfolder->measureModel();

    // using basic heuristic
    unfolder->buildUnfolding();

    unfolder->unfoldTo(0.0);

    if (unfolder->isFlattened())
      flattened++;

    cout << string(40, '-') << endl;

    if (unfolder->getCheckOverlappingCalls() > 0) {
      cout << "- Total CO calls  = " << unfolder->getCheckOverlappingCalls()
          << endl;
      cout << "- Total CO time   = "
          << unfolder->getTotalCheckOverlappingTime() * 1.0 / CLOCKS_PER_SEC
          << " s" << endl;
      cout << "- Average CO time = "
          << unfolder->getTotalCheckOverlappingTime() * 1.0
              / unfolder->getCheckOverlappingCalls() / CLOCKS_PER_SEC << " s"
          << endl;
    }

    unfolders.push_back(unfolder);

    cout << string(40, '-') << endl;
  }

  auto time_cost = (clock() - start) * 1.0 / CLOCKS_PER_SEC;

  cout << "Total flattened = " << flattened << "/" << filenames.size() << endl;
  cout << "Total time = " << time_cost << " s" << endl;

  computeCOM_R();

  return true;
}

void computeCOM_R() {
//compute a bbox
  double box[6] = { FLT_MAX, -FLT_MAX, FLT_MAX, -FLT_MAX, FLT_MAX, -FLT_MAX };
//-------------------------------------------------------------------------
  for (auto& u : unfolders) {
    auto& m = *u->getModel();
    for (int j = 0; j < m.v_size; j++) {
      Point3d& p = m.vertices[j].p;
      if (p[0] < box[0])
        box[0] = p[0];
      if (p[0] > box[1])
        box[1] = p[0];
      if (p[1] < box[2])
        box[2] = p[1];
      if (p[1] > box[3])
        box[3] = p[1];
      if (p[2] < box[4])
        box[4] = p[2];
      if (p[2] > box[5])
        box[5] = p[2];
    }				//j
  }				//i

//-------------------------------------------------------------------------
// compute center of mass and R...
  COM.set((box[1] + box[0]) / 2, (box[3] + box[2]) / 2, (box[5] + box[4]) / 2);

//-------------------------------------------------------------------------
  R = 0;
  for (auto& u : unfolders) {
    auto& m = *u->getModel();
    for (int j = 0; j < m.v_size; j++) {
      Point3d& p = m.vertices[j].p;
      double d = (p - COM).normsqr();
      if (d > R)
        R = d;
    }    //j
  }    //i

  R = sqrt(R);
}

#endif //_BF_MINKOWSKI_SUM_H_

