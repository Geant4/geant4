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
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli,University of Wollongong, Australia
// Contributions by F. Ambroglini INFN Perugia, Italy
//
#ifndef G4HumanPhantomConstruction_H
#define G4HumanPhantomConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "G4HumanPhantomMessenger.hh"

#include "globals.hh"
#include <map>

class G4VPhysicalVolume;
class G4HumanPhantomMaterial;
class G4HumanPhantomConstruction : public G4VUserDetectorConstruction
{
  public:
     G4HumanPhantomConstruction();
    ~G4HumanPhantomConstruction();
     G4VPhysicalVolume* Construct();

  void SetBodyPartSensitivity(G4String, G4bool);

  void SetPhantomSex(G4String);
  void SetPhantomModel(G4String);
  void ConstructSDandField();
  //G4VPhysicalVolume* GetMotherVolume(){return mother;};
 
 private:
  G4VPhysicalVolume* ConstructWorld();
  G4HumanPhantomMaterial* fMaterial;
  G4HumanPhantomMessenger* fMessenger;
  std::map<std::string,G4bool> fSensitivities;
  G4String                 fModel;
  G4String                 fSex;
};
#endif

