//------------------------------------------------------------------------------
//  Copyright 2007-2008 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#include "main.h"

//
//
//
//
//   The MAIN function
//
//
//
//

int main(int argc, char ** argv) {
  if (!parseArg(argc, argv)) {
    printUsage(argv[0]);
    return 1;
  }

  cout << " - seed = " << config.seed << endl;
  mathtool::srand48(config.seed);

  //
  readfromfiles();

  if (config.disable_gui) 
  {
    // only dump if find an unfolding
    dumpUnfolding(true);
    return 0;
  }

  /////////////////////////////////////////////////////////////////
  //setup glut/gli
  glutInit(&argc, argv);
  glutInitDisplayMode(
  GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);
  glutInitWindowSize(880, 880);
  glutInitWindowPosition(50, 50);
  glutCreateWindow("unfolder");

  InitGL();
  gli::gliInit();
  gli::gliDisplayFunc(Display);
  glutReshapeFunc(Reshape);
  glutKeyboardFunc(Keyboard);

  //set camera position
  gli::setCameraPosZ(R * 2.1);
  /////////////////////////////////////////////////////////////
  gli::gliMainLoop();

  return 0;
}

