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
/// \file exoticphysics/dmparticle/include/PhysicsList.hh
/// \brief Definition of the PhysicsList class
//
// $Id: PhysicsList.hh 92047 2015-08-14 07:23:37Z gcosmo $
//
/////////////////////////////////////////////////
//
// ClassName:   PhysicsList
//
// Authors:  01.06.17 V.Ivanchenko 
//
//
///////////////////////////////////////////

#ifndef PhysicsList_h
#define PhysicsList_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4VPhysicsConstructor;
class PhysicsListMessenger;
class StepMax;

class PhysicsList: public G4VModularPhysicsList
{
public:

  PhysicsList();
  virtual ~PhysicsList();

  virtual void ConstructParticle();
    
  virtual void SetCuts();
        
  virtual void ConstructProcess();

  void AddPhysicsList(const G4String& name);

  void SetLDMPhotonMass(G4double val);

  void SetLDMHiMass(G4double val);

private:

  void AddDarkMatter(); 

  G4VPhysicsConstructor*  fEmPhysicsList;
  G4VPhysicsConstructor*  fDecayPhysicsList;
  G4String                fEmName;
    
  PhysicsListMessenger* fMessenger;

  G4double fLDMPhotonMass; 
  G4double fLDMHiMass;

  G4bool fLDMPhoton;
  G4bool fLDMHi;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

