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
// $Id: G4OpenInventorTransform3D.cc 66373 2012-12-18 09:41:34Z gcosmo $
//
// 
// jck 17 Dec 1996
// G4OpenInventorTransform3D provides OpenGL style transformation matrix
// from G4Transform3D.
// gb 22 Nov 2004 : use SbMatrix instead of SoSFMatrix.

#ifdef G4VIS_BUILD_OI_DRIVER

// this :
#include "G4OpenInventorTransform3D.hh"

#include <Inventor/SbLinear.h>

G4OpenInventorTransform3D::G4OpenInventorTransform3D (const G4Transform3D &t) 
: G4Transform3D (t) {
#define elem(i,j) ((float)t(i,j))
  m[0]  = elem(0,0); //xx
  m[1]  = elem(1,0); //yx
  m[2]  = elem(2,0); //zx
  m[3]  = 0;
  m[4]  = elem(0,1); //xy
  m[5]  = elem(1,1); //yy
  m[6]  = elem(2,1); //zy
  m[7]  = 0;
  m[8]  = elem(0,2); //xz
  m[9]  = elem(1,2); //yz
  m[10] = elem(2,2); //zz
  m[11] = 0;
  m[12] = elem(0,3); //dx
  m[13] = elem(1,3); //dy
  m[14] = elem(2,3); //dz
  m[15] = 1;
#undef elem
}

SbMatrix* G4OpenInventorTransform3D::GetSbMatrix () const {
  SbMatrix* tm = new SbMatrix(m[0],m[1],m[2],m[3],
                              m[4],m[5],m[6],m[7],
                              m[8],m[9],m[10],m[11],
                              m[12],m[13],m[14],m[15]);
  return tm;
}

#endif
