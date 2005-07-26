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
#ifndef G4HumanPhantomConstruction_H
#define G4HumanPhantomConstruction_H 1

#include "G4VUserDetectorConstruction.hh"

#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VisAttributes;

class G4HumanPhantomSD;
class G4VPhysicalVolume;

class G4HumanPhantomConstruction : public G4VUserDetectorConstruction
{
  public:
     G4HumanPhantomConstruction();
    ~G4HumanPhantomConstruction();
     G4VPhysicalVolume* Construct();

  void SetBodyPartSensitivity(G4String, G4bool);
  G4VPhysicalVolume* GetMotherVolume(){return mother;};
 
 private:
    G4HumanPhantomSD* userPhantomSD; 
  G4VPhysicalVolume* mother;

};

#endif

