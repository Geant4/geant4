//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4OpenGLTransform3D.hh,v 1.6 2004/04/07 15:17:12 gbarrand Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
// Andrew Walkden  24th October 1996
// G4OpenGLTransform3D provides OpenGL style transformation matrix
// from G4Transform3D.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#ifndef G4OPENGLTRANSFORM3D_HH
#define G4OPENGLTRANSFORM3D_HH

#include "G4Transform3D.hh"

#include "G4OpenGL.hh"

class G4OpenGLTransform3D : public G4Transform3D {
public:
  G4OpenGLTransform3D (const G4Transform3D &t);
  const GLdouble* GetGLMatrix ();
private:
  GLdouble m[16];
};

#endif

#endif
