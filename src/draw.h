#ifndef _MKSUM_DRAW_H_
#define _MKSUM_DRAW_H_

#include <ctime>
#include <map>
#include <algorithm>
#include <sstream>
using namespace std;

#include "GL/gli.h"

#include "Basic.h"
#include "model.h"
using namespace mathtool;

#include "util/ColorHelper.h"
using namespace masc::util;

//-----------------------------------------------------------------------------
//variables used in rendering

bool showUnfold = true;
bool showNumber = false;
bool showWire = false; //on/off wireframe
bool showWeights = false; //show edge weights
bool showSpanningTree = false;
bool showPathLength = false;
bool showSharpEdgeOnly = false;
bool showEdgeTypes = true;
bool showOverlapping = true;
bool saveImg = false;
bool printPDF = false;
bool animation = false;
bool randomColor = false;
bool background = false; //on/off background
bool showScore = false;
bool light = false; //on/off light
bool reverse_order = false;
uint current_imgID = 0; //id for dumping images

bool breakearly = false; //tmp
bool collision_dection = true;

map<model*, int> model_solid_gids;
map<model*, int> model_wire_gids;
map<model*, Vector3d> model_colors;

float unfolding_percentage = 0.0f;
double unfolding_speed = 0.02; //2%
float rot_z = 0.0f;

void updateStatus() {

  stringstream title;
  title << "Unfolder " << config.filename << " " << unfolding_percentage << " "
      << unfolding_speed;

  glutSetWindowTitle(title.str().c_str());
}

void DisplayPolygon(polygon& p);

void drawText(const char* text, const Vector3d& p) {
  glRasterPos3f(p[0], p[1], p[2]);
  for (int i = 0, size = strlen(text); i < size; i++) {
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, text[i]);
  }
}

inline void DisaplayUnfold(Unfolder& unfolder) {

  glPushMatrix();

  const auto& t = unfolder.getTranslation();
  const auto& rot_axis = unfolder.getRotationAxis().normalize();

  glRotatef(-unfolder.getRotationAngle(), rot_axis[0], rot_axis[1],
      rot_axis[2]);
  glTranslatef(-t[0], -t[1], -t[2]);

  const auto& unfolded = unfolder.getUnfolded();

  glEnable( GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(0.5f, 0.5f);
  glEnable(GL_CULL_FACE);

  // ----------------------------------------
  // draw facets
  auto fid = 0;

  for (int side = 0; side <= 1; ++side) {
    glCullFace(side == 0 ? GL_FRONT : GL_BACK);
    glBegin(GL_TRIANGLES);
    fid = 0;
    auto norm = Vector3d(0, 1, 0);
    for (const auto& f : unfolded) {
      if (showOverlapping && unfolder.getModel()->tris[fid].overlapped) {
        glColor3d(1.0, 0.0, 0.0);
      } else if (showPathLength) {
        double b = 1.0
            - unfolder.getModel()->tris[fid].path_len
                / unfolder.getMaxPathLength();
        glColor3d(b, b, b);
      } else {
        if (side == 0) {
          glColor3dv(unfolder.getColor().get());
        } else {
          glColor3d(1.0, 1.0, 0.8);
        }
      }

      glNormal3dv(norm.get());

      for (auto k = 0; k < 3; k++) {
        glVertex3dv(f[k].second.get());
      }

      fid++;
    }

    glEnd();
  }

  glDisable(GL_CULL_FACE);
  glDisable( GL_POLYGON_OFFSET_FILL);

  // -----------------------------------------

//  glColor3f(0, 0, 0);
//  glBegin(GL_LINES);
//  for (const auto f : unfolded) {
//    for (auto k = 0; k < 3; k++) {
//      glVertex3dv(f[k].second.get());
//      glVertex3dv(f[(k + 1) % 3].second.get());
//    }
//  }
//  glEnd();

  fid = 0;
  glBegin(GL_LINES);
  const auto m = unfolder.getModel();
  const auto edge_weithgs = unfolder.getEdgeWeights();
  const auto& se = unfolder.getSelectedEdges();
  for (const auto f : unfolded) {
    const auto& t = m->tris[fid];

    const auto& e1 = m->edges[t.e[1]];
    const auto& e2 = m->edges[t.e[2]];

    for (auto k = 0; k < 3; k++) {
      const auto& e = m->edges[t.e[k]];

      glColor3f(0.0f, 0.0f, 0.0f);

      if (showEdgeTypes) {
        if (e.diagonal) {
          glColor3f(0.9f, 0.9f, 0.9f);
        } else if (e.folding_angle == 0) {
          glColor3f(0.5f, 0.5f, 0.5f);
        }
      }

      if (showWeights) {
        if (se.count(std::make_pair((int) e.fid[0], (int) e.fid[1]))
            || se.count(std::make_pair((int) e.fid[1], (int) e.fid[0])))
        {
          // do nothing
        } else{
          glColor3f(1.0f, 0, 0);
        }
      }

      glVertex3dv(f[k].second.get());
      glVertex3dv(f[(k + 1) % 3].second.get());
    }

    fid++;
  }
  glEnd();

  if (showSpanningTree) {
    glColor3f(0.0f, 0.7f, 0.7f);
    glBegin(GL_LINES);
    const auto& se = unfolder.getSelectedEdges();
    for (auto& e : se) {
      const auto fid1 = e.first;
      const auto fid2 = e.second;

      const auto c1 = (unfolded[fid1][0].second + unfolded[fid1][1].second
          + unfolded[fid1][2].second) / 3;
      const auto c2 = (unfolded[fid2][0].second + unfolded[fid2][1].second
          + unfolded[fid2][2].second) / 3;

      glVertex3dv(c1.get());
      glVertex3dv(c2.get());
    }
    glEnd();
  }

  if (!showNumber) {

    glPopMatrix();

    return;
  }

  fid = 0;

  auto n = unfolded[0][0].second % unfolded[0][1].second;

  for (const auto f : unfolded) {
    Vector3d cof;
    for (auto k = 0; k < 3; k++) {
      cof += f[k].second;
    }
    cof = cof * 1.0 / 3 + n * 0.01;

    char buf[100];

    sprintf(buf, "%d", fid);
    drawText(buf, cof);

    fid++;
  }

  glPopMatrix();
}

inline void DisplayModel(Unfolder& unfolder, bool randcolor = false) {
  auto M = unfolder.getModel();

  if (model_solid_gids.find(M) == model_solid_gids.end()) {
    const Point3d& O = COM;
    model_solid_gids[M] = glGenLists(1);
    glNewList(model_solid_gids[M], GL_COMPILE);

    //Draw facets
    glEnable( GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(0.5f, 0.5f);
    //for(list<polygon>::iterator i=M.polys.begin();i!=M.polys.end();i++)
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < M->t_size; i++) {
      const triangle & t = M->tris[i];

      Point3d cluster_color;
      switch (t.cluster_id % 14) {
      case 0:
        cluster_color = Point3d(255, 255, 255);
        break;
      case 1:
        cluster_color = Point3d(255, 0, 255);
        break;
      case 2:
        cluster_color = Point3d(0, 255, 0);
        break;
      case 3:
        cluster_color = Point3d(0, 0, 255);
        break;
      case 4:
        cluster_color = Point3d(255, 255, 0);
        break;
      case 5:
        cluster_color = Point3d(0, 255, 255);
        break;
      case 6:
        cluster_color = Point3d(255, 0, 0);
        break;

      case 7:
        cluster_color = Point3d(127, 127, 127);
        break;
      case 8:
        cluster_color = Point3d(127, 0, 127);
        break;
      case 9:
        cluster_color = Point3d(0, 127, 0);
        break;
      case 10:
        cluster_color = Point3d(0, 0, 127);
        break;
      case 11:
        cluster_color = Point3d(127, 127, 0);
        break;
      case 12:
        cluster_color = Point3d(0, 127, 127);
        break;
      case 13:
        cluster_color = Point3d(127, 0, 0);
        break;
      }

      if (unfolder.getConfig().scalar_score || showScore) {
        auto color = masc::util::getColour(t.score, 0.0, 1.0);
        glColor3dv(color.get());
      } else {
        glColor3ub(cluster_color[0], cluster_color[1], cluster_color[2]);
      }
      glNormal3dv(M->tris[i].n.get());
      for (int k = 0; k < 3; k++) {
        const Point3d& pt = M->vertices[t.v[k]].p;
        glVertex3dv(pt.get());
      }
    }
    glEnd();
    glDisable( GL_POLYGON_OFFSET_FILL);
    glEndList();
  }

  //draw
  if (randcolor) {
    if (model_colors.find(M) == model_colors.end())
      model_colors[M] =
          Vector3d(mathtool::drand48() / 2 + 0.5, mathtool::drand48() / 2 + 0.5,
              mathtool::drand48() / 2 + 0.5).normalize(); // +Vector3d(0.25, 0.25, 0.25);
    glColor3dv(model_colors[M].get());
  }

  glCallList(model_solid_gids[M]);

  if (!showNumber)
    return;

  // draw fids
  for (int i = 0; i < M->t_size; i++) {
    const triangle & t = M->tris[i];
    Vector3d COF; // center of face

    for (int k = 0; k < 3; k++) {
      const Point3d& pt = M->vertices[t.v[k]].p;

      COF = COF + (pt - COM);
    }

    COF = COF / 3 + t.n * 0.05; //TEMP to fix

    char buf[100];

    sprintf(buf, "%d", i);
    drawText(buf, COF);
  }

}

inline void DisplayModelWireFrame(Unfolder& u, bool randcolor = false) {
  auto M = u.getModel();
  //Draw Edges
  if (showWire) {

    if (model_wire_gids.find(M) == model_wire_gids.end()) {
      const Point3d& O = COM;
      model_wire_gids[M] = glGenLists(1);
      glNewList(model_wire_gids[M], GL_COMPILE);
      glBegin(GL_LINES);
      for (uint i = 0; i < M->e_size; i++) {
        glColor3f(0, 0, 0);
        const edge & e = M->edges[i];
        if (e.fid.size() == 2) {	//normal case, check if e is sharp
          triangle& f1 = M->tris[e.fid.front()];
          triangle& f2 = M->tris[e.fid.back()];
          if (fabs(1 - f1.n * f2.n) < 1e-2) {
            if (showSharpEdgeOnly)
              continue; //not sharp
            else
              glColor3f(0.7f, 0.7f, 0.7f);
          }
          // show cut boundays
          glLineWidth(1.0);
          if (e.cutted) {
            glLineWidth(5.0);
            glColor3f(1.0f, 0, 0);
          }
        }

        Point3d& p1 = M->vertices[e.vid[0]].p;
        Point3d& p2 = M->vertices[e.vid[1]].p;
        glVertex3d(p1[0], p1[1], p1[2]);
        glVertex3d(p2[0], p2[1], p2[2]);
      }
      glEnd();
      glEndList();
    }

    glCallList(model_wire_gids[M]);
  }
}

//-----------------------------------------------------------------------------

//copied from meshlab
void DisplayBackground(void) {
  float topcolor[] = { 1, 1, 1 };
  float bottomcolor[] = { 1, 0, 0 };

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(-1, 1, -1, 1, -1, 1);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glPushAttrib(GL_ENABLE_BIT);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);
  glDisable(GL_TEXTURE_2D);
  glBegin(GL_TRIANGLE_STRIP);
  glColor3fv(topcolor);
  glVertex2f(-1, 1);
  glColor3fv(bottomcolor);
  glVertex2f(-1, -1);
  glColor3fv(topcolor);
  glVertex2f(1, 1);
  glColor3fv(bottomcolor);
  glVertex2f(1, -1);
  glEnd();

  glPopAttrib();
  glPopMatrix(); // restore modelview
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
}

void drawAll() {
  glEnable(GL_LIGHTING);

  //show the inputs
  glColor3f(1, 1, 1);
  if (light)
    glEnable(GL_LIGHTING);
  else
    glDisable(GL_LIGHTING);

  if (showUnfold) {

    for (auto unfolderPtr : unfolders) {
      DisaplayUnfold(*unfolderPtr);
    }

    return; // do not show model
  }

  glEnable(GL_LIGHTING);
  glColor3f(1, 1, 1);
  for (auto& u : unfolders) {
    DisplayModel(*u, randomColor);
    if (breakearly)
      break;
  }
  for (auto& u : unfolders) {
    DisplayModelWireFrame(*u);
    if (breakearly)
      break;
  }
}

//-----------------------------------------------------------------------------
void dumpUnfolding(bool flattened_only = false) {
  for (auto u : unfolders) {
    if (flattened_only && !u->isFlattened())
      continue;

    const auto base_name = u->getFilename().substr(0,
        u->getFilename().find_last_of("."));
    const auto prefix = base_name
        + (u->getConfig().no_tick ? "" : "_s" + std::to_string(config.seed));

    u->dumpSVG(prefix + ".svg");

    // dump cut svg
    u->dumpSVG(prefix + "_cut.svg", 1);
    u->dumpSVG(prefix + "_tree.svg", 5);
  }
}

//-----------------------------------------------------------------------------
void Display(void) {
  //Init Draw
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  if (background)
    DisplayBackground();

  glPushMatrix();
  glLoadIdentity();
  static GLfloat light_position1[] = { 100, 100, 100.0f, 1.0f };
  glLightfv(GL_LIGHT0, GL_POSITION, light_position1);
  static GLfloat light_position2[] = { -100, -100, 50.0f, 1.0f };
  glLightfv(GL_LIGHT1, GL_POSITION, light_position2);
  glPopMatrix();

  glRotatef(rot_z, 0, 0, 1); //rotate along Z-axis

  if (R > 0)
    glTranslatef(-COM[0], -COM[1], -COM[2]);

  drawAll();

  glDisable(GL_LIGHTING);
}

//-----------------------------------------------------------------------------
// regular openGL callback functions
bool InitGL() {
  // transparent
  glShadeModel(GL_SMOOTH);
  glEnable( GL_BLEND);
  glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glEnable( GL_LINE_SMOOTH);
  glHint( GL_LINE_SMOOTH_HINT, GL_NICEST);

  // others
  glEnable( GL_DEPTH_TEST);
  //glDepthMask(GL_TRUE);
  glEnable(GL_NORMALIZE);
  glClearColor(1.0, 1.0, 1.0, 0.0);

  //Let's have light!
  GLfloat Diffuse[] = { 0.9f, 0.9f, 0.9f, 1.0f };
  glMaterialfv(GL_FRONT, GL_DIFFUSE, Diffuse);
  glColorMaterial(GL_FRONT, GL_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);

  GLfloat WhiteLight[] = { 0.75f, 0.75f, 0.75f, 1.0f };
  glLightfv(GL_LIGHT0, GL_DIFFUSE, WhiteLight);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, WhiteLight);

  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);

  return true;
}

void Reshape(int w, int h) {
  if (w > h)
    glViewport(0, 0, (GLsizei) w, (GLsizei) w);
  else
    glViewport(0, 0, (GLsizei) h, (GLsizei) h);
  glMatrixMode( GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60, 1, R / 100, R * 100);
  glMatrixMode( GL_MODELVIEW);
  glLoadIdentity();
}

// Generate a random unfold using current heuristic
void randomUnfold() {
  for(auto u : unfolders) {
    u->buildUnfolding();
    u->unfoldTo(1.0);
  }
}

/*
 * Update unfolding to given percentage
 */
void updateUnfolding(float p) {
  showUnfold = true;
  p = std::min(std::max(p, 0.0f), 1.0f);

  unfolding_percentage = p;

  for (auto u : unfolders) {

    u->unfoldTo(unfolding_percentage);

    if (collision_dection) {
      if (unfolding_percentage == 1.0)
        u->checkOverlap();
      else
        u->checkCollision();
    }
  }
}

void TimerCallback(int value) {
  // not in animation
  if (!animation) {
    return;
  }

  unfolding_percentage += unfolding_speed * (reverse_order ? -1 : 1);

  updateUnfolding(unfolding_percentage);

  updateStatus();

  //in simulation state
  glutPostRedisplay();

  // fully folded/unfolded
  if (unfolding_percentage < -1 || unfolding_percentage > 1) {
    animation = false;
    return;
  }

  glutTimerFunc(30, TimerCallback, value);
}

void randomColors() {
  model_colors.clear();

  for (auto u : unfolders) {
    u->setColor(mathtool::drand48(), mathtool::drand48(), mathtool::drand48());
  }
}

void printGUIKeys() {

  const string mvmo = " [Model View Mode Only]";

  cout << "---------- GUI keys ----------" << endl;
  cout << "Displaying:" << endl;
  cout << "  b: show background" << endl;
  cout << "  c: toggle random colors" << endl;
  cout << "  e: show edge types" << endl;
  cout << "  n: show numbers" << endl;
  cout << "  o: show overlapping" << endl;
  cout << "  w: show wire frame" << mvmo << endl;
  cout << "  u: toggle showing unfolding" << endl;
  cout << "  C: toggle collision detection" << endl;
  cout << "  R: show animation in reverse order" << endl;
  cout << "  S: only show sharp edges" << mvmo << endl;
  cout << "  T: show spanning tree" << endl;
  cout << "z/Z: rotation model around z-axis" << endl;
  cout << endl;
  cout << "Control:" << endl;
  cout << "  r: fully folded" << endl;
  cout << "  t: fully unfolded" << endl;
  cout << " sp: toggle animation" << endl;
  cout << "Esc: exit" << endl;
  cout << endl;
  cout << "Dumping:" << endl;
  cout << "  D: dump unfolding (ori/svg)" << endl;
  cout << "  M: dump current model (obj)" << endl;
  cout << "  W: dump components (wrl)" << endl;
}

void Keyboard(unsigned char key, int x, int y) {
  // find closest colorPt3D if ctrl is pressed...
  switch (key) {
  case 27:
    exit(0);
  case 'a':
    randomUnfold();
    break;
  case 'e':
    showEdgeTypes = !showEdgeTypes;
    break;
  case 'w':
    showWire = !showWire;
    showWeights = !showWeights;
    break;
  case 'r':
    updateUnfolding(0.0f);
    break;
  case 'R':
    reverse_order = !reverse_order;
    break;
  case 't':
    updateUnfolding(1.0f);
    break;
  case 'c':
    randomColors();
    randomColor = !randomColor;
    break;
  case 's':
    saveImg = !saveImg;
    break;
  case 'L':
    light = !light;
    break;
  case 'n':
    showNumber = !showNumber;
    break;
  case 'p':
    printPDF = !printPDF;
    break;
  case 'b':
    background = !background;
    break;
  case 'S':
    showSharpEdgeOnly = !showSharpEdgeOnly;
    for (map<model*, int>::iterator i = model_wire_gids.begin();
        i != model_wire_gids.end(); i++)
      glDeleteLists(i->second, 1);
    model_wire_gids.clear();
    break;
  case 'u':
    showUnfold = !showUnfold;
    break;
  case 'T':
    showSpanningTree = !showSpanningTree;
    break;
  case 'P':
    showPathLength = !showPathLength;
    break;
  case 'D':
    dumpUnfolding();
    break;
  case 'C':
    collision_dection = !collision_dection;
    break;
  case 'o':
    showOverlapping = !showOverlapping;
    break;
  case ' ':
    animation = !animation; //rotate scene
    glutTimerFunc(10, TimerCallback, animation ? 1 : 0);
    break;
  case ',':
    updateUnfolding(unfolding_percentage + unfolding_speed);
    break;
  case '.':
    updateUnfolding(unfolding_percentage - unfolding_speed);
    break;
  case '-':
    unfolding_speed *= 0.9;
    unfolding_speed = max(unfolding_speed, 1e-4);
    break;
  case '=':
    unfolding_speed /= 0.9;
    unfolding_speed = min(unfolding_speed, 0.08);
    break;
  case 'z':
    rot_z += 0.618f;
    break;
  case 'Z':
    rot_z -= 0.618f;
    break;
  case '/':
    breakearly = !breakearly;
    break;
  case '?':
    printGUIKeys();
    break;
  }

  updateStatus();

  glutPostRedisplay();
}

#endif //_MKSUM_DRAW_H_

