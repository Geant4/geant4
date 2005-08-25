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
#include "G4HumanPhantomMessenger.hh"

#include "globals.hh"
#include <map>

class G4VPhysicalVolume;
class G4LogicalVolume;

class G4HumanPhantomSD;

class G4HumanPhantomConstruction : public G4VUserDetectorConstruction
{
  public:
     G4HumanPhantomConstruction();
    ~G4HumanPhantomConstruction();
     G4VPhysicalVolume* Construct();

  std::map<std::string,G4bool> sensitivities;

  void SetBodyPartSensitivity(G4String, G4bool);
  void CleanPhantom();
  void UpdatePhantom();
  void SetPhantomSex(G4String);
  void SetPhantomModel(G4String);

  G4VPhysicalVolume* GetMotherVolume(){return mother;};
 
 private:
  
  G4HumanPhantomMessenger* messenger;

  G4HumanPhantomSD*        userPhantomSD; 
  G4VPhysicalVolume*       mother;

  G4String                 model;
  G4String                 sex;
};

#endif

