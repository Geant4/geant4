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
/// \file electromagnetic/TestEm6/include/PhysicsList.hh
/// \brief Definition of the PhysicsList class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "G4EmConfigurator.hh"


class G4VPhysicsConstructor;
class StepMax;
class PhysicsListMessenger;
class G4GammaConversionToMuons;
class G4AnnihiToMuPair;
class G4eeToHadrons;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysicsList: public G4VModularPhysicsList
{
  public:
    PhysicsList();
    virtual ~PhysicsList();

    // Construct particles
    virtual void ConstructParticle();
    virtual void ConstructProcess();

    void AddPhysicsList(const G4String& name);
    void ConstructHighEnergy();

    void AddStepMax();

    // Construct processes and register them

    void SetGammaToMuPairFac(G4double);
    void SetAnnihiToMuPairFac(G4double);
    void SetAnnihiToHadronFac(G4double);

  private:
 
    G4VPhysicsConstructor* fEmPhysicsList;
    G4VPhysicsConstructor* fDecayPhysicsList;
    G4String fEmName;

    StepMax* fStepMaxProcess;

    PhysicsListMessenger*  fMes;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

