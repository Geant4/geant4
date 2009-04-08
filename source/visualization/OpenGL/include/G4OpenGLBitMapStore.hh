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
// $Id: G4OpenGLBitMapStore.hh,v 1.4 2009-04-08 15:15:07 lgarnier Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  6th January 2007
//
// Class description
//
// Keeps bit maps on byte boundaries suitable for drawing.  For
// example, in OpenGL:
//
//   const char* circle = G4OpenGLBitMapStore::GetCircle(size, true);
//   glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
//   glBitmap(size, size, size/2., size/2., 0., 0., circle);

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#ifndef G4OPENGLBITMAPSTORE_HH
#define G4OPENGLBITMAPSTORE_HH

#include "globals.hh"
#include <map>

#include "G4OpenGL.hh"
#include <GL/gl.h>
#include <GL/glu.h>

namespace G4OpenGLBitMapStore {

  enum Shape {circle, square};

  const GLubyte* GetBitMap(Shape, G4double& size, G4bool filled);
  // Size in pixels (gets changed to a rationalised value).

  struct Key{
    Key(Shape shape, G4int size, G4bool filled):
      fShape(shape), fSize(size), fFilled(filled) {}
    bool operator<(const Key& rhs) const {
      if (fShape < rhs.fShape) return true;
      else if (fSize < rhs.fSize) return true;
      else if (fFilled != rhs.fFilled) return true;
      else return false;
    }
    Shape fShape;
    G4int fSize;
    G4bool fFilled;
  };

  extern std::map<Key, GLubyte*> fStore;

}

#endif

#endif
