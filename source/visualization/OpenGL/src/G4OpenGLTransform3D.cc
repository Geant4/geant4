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
// $Id: G4OpenGLTransform3D.cc,v 1.9 2009-01-19 16:53:42 lgarnier Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
