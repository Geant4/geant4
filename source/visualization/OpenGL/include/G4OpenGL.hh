//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id$
//
// G.Barrand.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

 #ifndef G4OpenGL_h
 #define G4OpenGL_h 

 #ifdef WIN32
 #include <windows.h>
 #undef min
 #undef max
 #endif


 #ifdef G4VIS_BUILD_OPENGLX_DRIVER
 #  include <GL/gl.h>
 #  include <GL/glu.h>
 #endif

 #ifdef G4VIS_BUILD_OPENGLXM_DRIVER
 #    include <GL/gl.h>
 #    include <GL/glu.h>
 #endif

 #ifdef G4VIS_BUILD_OPENGLWIN32_DRIVER
 #    include <GL/gl.h>
 #    include <GL/glu.h>
 #endif
//# Do NOT include glx Here ! It has to be done, after all <Qxx...> includes
//#  include <GL/glx.h>

#ifdef G4VIS_BUILD_OPENGLWT_DRIVER
 #include <Wt/WGLWidget>
 #include <qgl.h>
#endif
#ifdef  G4VIS_BUILD_OPENGLQT_DRIVER
  #ifndef G4VIS_BUILD_OPENGLX_DRIVER
    #ifdef __MACH__
      #include <OpenGL/gl.h>
      #include <OpenGL/glu.h>
    #endif
    #include <qgl.h>
  #endif
#endif

#define G4OPENGL_FLT_BIG 1.e20

 #endif

#endif
