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
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef F04PhysicsList_h
#define F04PhysicsList_h 1

#include "globals.hh"
#include "G4VModularPhysicsList.hh"

class G4VPhysicsConstructor;
class F04PhysicsListMessenger;

class F04StepMax;

class F04PhysicsList: public G4VModularPhysicsList
{
public:

    F04PhysicsList(G4String);
    virtual ~F04PhysicsList();

    void ConstructParticle();
    
    void SetCuts();
    void SetCutForGamma(G4double);
    void SetCutForElectron(G4double);
    void SetCutForPositron(G4double);

    void SetStepMax(G4double);
    F04StepMax* GetStepMaxProcess();
    void AddStepMax();
    
    /// Add physics to the Physics List    
    void AddPhysicsList(const G4String& name);

    /// Remove specific EM physics from EM physics list.
    void RemoveFromEMPhysicsList(const G4String&);

    /// Remove specific Hadron physics from Hadron physics list.
    void RemoveFromHadronPhysicsList(const G4String&);

    /// Make sure that the EM physics list is empty.
    void ClearEMPhysics();

    /// Make sure that the hadron physics list is empty.
    void ClearHadronPhysics();

    void ConstructProcess();
    void List();
  
private:

    typedef std::vector<G4VPhysicsConstructor*>  PhysicsListVector;

    void SetStandardList(G4bool flagHP = false, G4bool glauber = false);

    G4double fCutForGamma;
    G4double fCutForElectron;
    G4double fCutForPositron;

    G4VPhysicsConstructor*  fParticleList;

    PhysicsListVector* fEMPhysics;
    PhysicsListVector* fHadronPhysics;

    G4double MaxChargedStep;
    F04StepMax* stepMaxProcess;
    
    F04PhysicsListMessenger* fMessenger;

    G4bool fDump;
};
#endif
