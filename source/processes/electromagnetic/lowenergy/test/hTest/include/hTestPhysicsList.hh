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
#ifndef hTestPhysicsList_h
#define hTestPhysicsList_h 1

//---------------------------------------------------------------------------
//
// ClassName:   hTestPhysicsList
//  
// Description: hTest PhysicsList 
//
// Authors:    07.04.01 V.Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VUserPhysicsList.hh"
#include "hTestDetectorConstruction.hh"
#include "hTestVEMPhysicsList.hh"
#include "hTestVHadronPhysicsList.hh"
#include "globals.hh"

class hTestPhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestPhysicsList: public G4VUserPhysicsList
{
public: // Without description
    hTestPhysicsList(const hTestDetectorConstruction*);
   ~hTestPhysicsList();

public: // Without description
    void SetGammaCut(G4double);
    void SetElectronCut(G4double);
    void SetProtonCut(G4double);
    void SetElectronCutByEnergy(G4double);
    void SetLowEnergyLimit(G4double);
    void SetHighEnergyLimit(G4double);
    void SetMaxStep(G4double);
    void SetEMPhysicsList(const G4String&);  
    void SetHadronPhysicsList(const G4String&);  
    inline void SetVerbose(G4int val) {verbose = val;};    
    inline G4int GetVerbose() const {return verbose;};    

protected:
    // Construct particle and physics
    void ConstructParticle();
    void ConstructProcess();

private:
    void InitializeMe();
    void SetCuts();

    // these methods Construct particles 
    void ConstructMyBosons();
    void ConstructMyLeptons();
    void ConstructMyMesons();
    void ConstructMyBarions();
    void ConstructMyIons();

  // these methods Construct physics processes and register them
    void ConstructDecay();
    
  private:

    const hTestDetectorConstruction* pDet;
    hTestPhysicsListMessenger* theMessenger;
    hTestVEMPhysicsList* theEMList;
    hTestVHadronPhysicsList* theHadList;

    G4double cutForGamma;
    G4double cutForElectron;
    G4double cutForProton;
    G4double maxChargedStep;    
    G4double lowEnergyLimit;
    G4double highEnergyLimit;

    G4String emPhysics;
    G4String hadronPhysics;

    G4int    verbose;
    G4bool   physicsIsDefined;
};

#endif



