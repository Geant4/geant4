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
// $Id: G4OpenGLTransform3D.cc,v 1.7 2004/09/13 12:12:08 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// 
// Andrew Walkden  24th October 1996
// G4OpenGLTransform3D provides OpenGL style transformation matrix
// from G4Transform3D.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#include "G4OpenGLTransform3D.hh"

G4OpenGLTransform3D::G4OpenGLTransform3D (const G4Transform3D &t):
  G4Transform3D (t) {}

const GLdouble* G4OpenGLTransform3D::GetGLMatrix () 
{
  GLdouble *p = m;
  for (size_t i=0; i<4; i++)
  { 
    for (size_t k=0; k<3; k++)
    {
      *p++ = operator()(k,i);
    }
    *p++ = 0.; 
  } 
  m[15] = 1.; 
  return m;
}

#endif
