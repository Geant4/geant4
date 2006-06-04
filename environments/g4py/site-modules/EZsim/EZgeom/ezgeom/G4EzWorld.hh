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
// $Id: G4EzWorld.hh,v 1.2 2006-06-04 21:36:34 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   G4EzWorld.hh
//
//   Utility class for supporting creating user geometry
//
//   management of a world volume
//
//                                         2005 Q
// ====================================================================
#ifndef G4_EZ_WORLD_H
#define G4_EZ_WORLD_H

#include "globals.hh"

// ====================================================================
//
// class definition
//
// ====================================================================
class G4Material;
class G4VPhysicalVolume;

class G4EzWorld {
protected:
  // world volume is automatically presented.
  static G4VPhysicalVolume* world;
  static G4VPhysicalVolume* CreateWorld(G4double dx=1.*m,
					G4double dy=1.*m,
					G4double dz=1.*m);
public:
  G4EzWorld();
  ~G4EzWorld();
  
  static G4VPhysicalVolume* GetWorldVolume();

  static void Reset(G4double dx, G4double dy, G4double dz);
  static void Resize(G4double dx, G4double dy, G4double dz);
  static void SetMaterial(G4Material* amaterial);
  static void SetVisibility(G4bool qvis);

};

// ====================================================================
//   inline functions
// ====================================================================

inline G4VPhysicalVolume* G4EzWorld::GetWorldVolume() { return world; }

#endif
