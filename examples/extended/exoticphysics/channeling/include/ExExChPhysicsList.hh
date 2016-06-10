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
// --------------------------------------------------------------
//

#ifndef ExExChPhysicsList_h
#define ExExChPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "ExExChPhysicsListMessenger.hh"

class ExExChPhysicsList: public G4VModularPhysicsList
{
private:
    G4VPhysicsConstructor*  fParticleList;
    G4VPhysicsConstructor*  fEmPhysicsList;
    G4VPhysicsConstructor*  fDecayList;
    G4VPhysicsConstructor*  fRaddecayList;
    G4VPhysicsConstructor*  fHadronElasticPhysicsList;
    G4VPhysicsConstructor*  fHadronInelasticPhysicsList;
    G4VPhysicsConstructor*  fStoppingPhysics;
    G4VPhysicsConstructor*  fIonPhysics;
    G4VPhysicsConstructor*  fNeutronTrackingCut;
    G4VPhysicsConstructor*  fEmExtraPhysics;
    
    G4String fFilePotentialName;
    ExExChPhysicsListMessenger *fMessenger;
    
public:

    G4double GetTransverseVariationMax() {return fTransverseVariationMax;};
    void SetTransverseVariationMax(G4double aDouble) {fTransverseVariationMax = aDouble;};

    G4double GetTimeStepMin() {return fTimeStepMin;};
    void SetTimeStepMin(G4double aDouble) {fTimeStepMin = aDouble;};
public:
    ExExChPhysicsList();
    ~ExExChPhysicsList();
    
    //Add processes
    void AddChanneling();
    void AddInelaticProcesses();

    // Construct particle and physics
    void ConstructParticle();
    void ConstructProcess();
    void SetCuts();
    
    void SetFilePotentialName(const G4String&);
    G4String GetFilePotentialName();
    G4double fTimeStepMin;
    G4double fTransverseVariationMax;
};

#endif
