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
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
// 
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
//
#ifndef G4HumanPhantomConstruction_H
#define G4HumanPhantomConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "G4HumanPhantomMessenger.hh"

#include "globals.hh"
#include <map>

class G4VPhysicalVolume;
//class G4LogicalVolume;
//class G4HumanPhantomSD;
class G4HumanPhantomMaterial;
//class G4PhantomBuilder;
class G4HumanPhantomEnergyDeposit;
class G4HumanPhantomConstruction : public G4VUserDetectorConstruction
{
  public:
     G4HumanPhantomConstruction(G4HumanPhantomEnergyDeposit*);
    ~G4HumanPhantomConstruction();
     G4VPhysicalVolume* Construct();

  void SetBodyPartSensitivity(G4String, G4bool);
  void CleanPhantom();
  void UpdatePhantom();
  void SetPhantomSex(G4String);
  void SetPhantomModel(G4String);

  G4VPhysicalVolume* GetMotherVolume(){return mother;};
 
 private:
  G4HumanPhantomEnergyDeposit* edepTot;
  G4HumanPhantomMessenger* messenger;

  G4VPhysicalVolume*       mother;

  G4String                 model;
  G4String                 sex;
  G4HumanPhantomMaterial* material;
  std::map<std::string,G4bool> sensitivities;
};

#endif

