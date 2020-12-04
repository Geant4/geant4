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

#ifndef GammaKnifePhysicsList_h
#define GammaKnifePhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class GammaKnifePhysicsListMessenger;
class G4VPhysicsConstructor;

class GammaKnifePhysicsList: public G4VModularPhysicsList
{
public:
  GammaKnifePhysicsList();
  virtual ~GammaKnifePhysicsList();

  void ConstructParticle();

  void SetCuts();
  void SetCutForGamma(G4double);
  void SetCutForElectron(G4double);
  void SetCutForPositron(G4double);
  void AddPhysicsList(const G4String& name);  
  void AddPackage(const G4String& pack);
  void ConstructProcess();

private:
  G4bool radioactiveDecayIsRegistered;

  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;

  G4String                             emName;
  G4VPhysicsConstructor*               emPhysicsList;
  G4VPhysicsConstructor*               decPhysicsList;
  std::vector<G4VPhysicsConstructor*>  hadronPhys;

  GammaKnifePhysicsListMessenger* messenger;
};

#endif



