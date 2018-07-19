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
// $Id: G4OpenGL.hh 101714 2016-11-22 08:53:13Z gcosmo $
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
   #ifndef  G4VIS_BUILD_OPENGLQT_DRIVER
     #  include <GL/gl.h>
   #else
     #ifdef __MACH__
       #include <OpenGL/gl.h>
     #else
       #include <GL/gl.h>
     #endif
   #endif
#endif

 #ifdef G4VIS_BUILD_OPENGLXM_DRIVER
  #ifndef  G4VIS_BUILD_OPENGLQT_DRIVER
    #  include <GL/gl.h>
  #else
    #ifdef __MACH__
      #include <OpenGL/gl.h>
    #else
      #include <GL/gl.h>
    #endif
  #endif
#endif

 #ifdef G4VIS_BUILD_OPENGLWIN32_DRIVER
 #    include <GL/gl.h>
 #endif
//# Do NOT include glx Here ! It has to be done, after all <Qxx...> includes
//#  include <GL/glx.h>

 #ifdef G4VIS_BUILD_OPENGLWT_DRIVER
 #  include <Wt/WGLWidget>
 #  define G4OPENGL_VERSION_2 1
 #endif
 #ifdef  G4VIS_BUILD_OPENGLQT_DRIVER
  #ifndef G4VIS_BUILD_OPENGLX_DRIVER
    #ifdef __MACH__
//#  define G4OPENGL_VERSION_2 1
      #include <OpenGL/gl.h>
    #else
      #include <GL/gl.h>
    #endif
    #include <qgl.h>
  #endif
#endif

#ifdef G4OPENGL_VERSION_2
#  undef G4VIS_BUILD_OPENGL_GL2PS
// include all redefinitions of openGl functions for Vertex Buffer Objects
#  include "G4OpenGLVboDrawer.hh"
#endif

#define G4OPENGL_FLT_BIG 1.e20

 #endif

#endif
