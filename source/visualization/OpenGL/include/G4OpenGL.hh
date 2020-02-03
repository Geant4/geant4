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
//
// G.Barrand.

#if defined (G4VIS_BUILD_OPENGL_DRIVER) || defined (G4VIS_USE_OPENGL)

 #ifndef G4OpenGL_h
 #define G4OpenGL_h 

 #ifdef WIN32
 #include <windows.h>
 #undef min
 #undef max
 #endif


 #if defined (G4VIS_BUILD_OPENGLX_DRIVER) || defined (G4VIS_USE_OPENGLX)
   #if defined (G4VIS_BUILD_OPENGLQT_DRIVER) || defined (G4VIS_USE_OPENGLQT)
     #ifdef __MACH__
       #include <OpenGL/gl.h>
     #else
       #include <GL/gl.h>
     #endif
   #else
     #include <GL/gl.h>
   #endif
 #endif

 #if defined (G4VIS_BUILD_OPENGLXM_DRIVER) || defined (G4VIS_USE_OPENGLXM)
   #if defined (G4VIS_BUILD_OPENGLQT_DRIVER) || defined (G4VIS_USE_OPENGLQT)
     #ifdef __MACH__
       #include <OpenGL/gl.h>
     #else
       #include <GL/gl.h>
     #endif
   #else
     #include <GL/gl.h>
   #endif
#endif

 #if defined (G4VIS_BUILD_OPENGLWIN32_DRIVER) || defined (G4VIS_USE_OPENGLWIN32)
 #    include <GL/gl.h>
 #endif
//# Do NOT include glx Here ! It has to be done, after all <Qxx...> includes
//#  include <GL/glx.h>

 #if defined (G4VIS_BUILD_OPENGLWT_DRIVER) || defined (G4VIS_USE_OPENGLWT)
 #  include <Wt/WGLWidget>
 #  define G4OPENGL_VERSION_2 1
 #endif
 #if defined (G4VIS_BUILD_OPENGLQT_DRIVER) || defined (G4VIS_USE_OPENGLQT)
   #if defined (G4VIS_BUILD_OPENGLX_DRIVER) || defined (G4VIS_USE_OPENGLX)
   #else
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
