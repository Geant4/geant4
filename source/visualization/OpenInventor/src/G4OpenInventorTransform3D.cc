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
// $Id: G4OpenInventorTransform3D.cc,v 1.6 2004/11/22 14:20:19 gbarrand Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
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
  m[0]  = xx; 
  m[1]  = yx; 
  m[2]  = zx; 
  m[3]  = 0;
  m[4]  = xy; 
  m[5]  = yy; 
  m[6]  = zy; 
  m[7]  = 0;
  m[8]  = xz; 
  m[9]  = yz; 
  m[10] = zz; 
  m[11] = 0;
  m[12] = dx; 
  m[13] = dy; 
  m[14] = dz; 
  m[15] = 1;
}

SbMatrix* G4OpenInventorTransform3D::GetSbMatrix () const {
  SbMatrix* tm = new SbMatrix(m[0],m[1],m[2],m[3],
                              m[4],m[5],m[6],m[7],
                              m[8],m[9],m[10],m[11],
                              m[12],m[13],m[14],m[15]);
  return tm;
}

#endif
