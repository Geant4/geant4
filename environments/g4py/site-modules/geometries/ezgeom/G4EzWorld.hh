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
// $Id: G4EzWorld.hh,v 1.1 2008-12-01 07:07:34 kmura Exp $
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
#include "G4SystemOfUnits.hh"

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
