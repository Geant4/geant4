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
//

#include "globals.hh"

#include "G4GeometryManager.hh"

#include "F04GlobalField.hh"

#include "F04FocusSolenoid.hh"
#include "F04SimpleSolenoid.hh"

F04FocusSolenoid::F04FocusSolenoid(G4double Ba, G4double Bb, G4double fz,
                                G4LogicalVolume* lv, 
                                G4ThreeVector c) : F04SimpleSolenoid(B1, fz, lv, c) 

{  
   B1 = Ba;
   B2 = Bb;
   Half = false;
}

void F04FocusSolenoid::addFieldValue(const G4double point[4],
                                         G4double field[6]) const
{
   G4ThreeVector global(point[0],point[1],point[2]);

   G4ThreeVector local = global2local.TransformPoint(global);

   if (isOutside(local)) return;

   G4double length = ((F04SimpleSolenoid*)this)->getLength();

   G4double Bz = (B2-B1) * std::abs(local.z())/(length/2.) + B1;

   if (Half) { if (local.z() >= 0.) Bz = B1; }

   G4ThreeVector B(0.0,0.0,Bz);

   B = global2local.Inverse().TransformAxis(B);

   field[0] += B[0];
   field[1] += B[1];
   field[2] += B[2];
}
