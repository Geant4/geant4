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
// $Id: G4OpenInventorTransform3D.hh,v 1.4 2001-07-11 10:09:00 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// jck  17 Dec 1996
// G4OpenGLInventorTransform3D provides SoSFMatrix transformation matrix
// from G4Transform3D.

#ifndef G4OpenInventorTRANSFORM3D_HH
#define G4OpenInventorTRANSFORM3D_HH

#ifdef G4VIS_BUILD_OI_DRIVER

#include "G4Transform3D.hh"

class SoSFMatrix;

class G4OpenInventorTransform3D : public G4Transform3D {
public:
  G4OpenInventorTransform3D (const G4Transform3D &t);
  SoSFMatrix* GetOIMatrix () const;

private:
  G4float m[16];
};

#endif

#endif
