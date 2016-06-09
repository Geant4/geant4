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
// $Id: G4VPrimaryGenerator.cc,v 1.4 2003/08/02 00:18:30 asaim Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//

// G4VPrimaryGenerator
#include "G4VPrimaryGenerator.hh"

G4VPrimaryGenerator::G4VPrimaryGenerator()
{;}

G4VPrimaryGenerator::~G4VPrimaryGenerator()
{;}

#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"

G4bool G4VPrimaryGenerator::CheckVertexInsideWorld
                         (const G4ThreeVector& pos)
{
  G4Navigator* navigator= G4TransportationManager::GetTransportationManager()
                                                 -> GetNavigatorForTracking();
  
  G4VPhysicalVolume* world= navigator-> GetWorldVolume();
  G4VSolid* solid= world-> GetLogicalVolume()-> GetSolid();
  EInside qinside= solid-> Inside(pos);
  
  if( qinside != kInside) return false;
  else return true;
}

