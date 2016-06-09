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
// $Id: Em5PhysicsList.hh,v 1.9 2003/04/30 14:12:34 maire Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// 14.10.02 (V.Ivanchenko) provide modular list on base of old Em5PhysicsList
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Em5PhysicsList_h
#define Em5PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class G4VPhysicsConstructor;
class Em5DetectorConstruction;
class Em5StepMax;
class Em5PhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Em5PhysicsList: public G4VModularPhysicsList
{
  public:
    Em5PhysicsList(Em5DetectorConstruction*);
   ~Em5PhysicsList();

    void ConstructParticle();
    
    void SetCuts();
    void SetCutForGamma   (G4double);
    void SetCutForElectron(G4double);
    void SetCutForPositron(G4double);
    
    void AddPhysicsList(const G4String& name);    
    void ConstructProcess();
    
    void AddStepMax();
    Em5StepMax* GetStepMaxProcess() {return stepMaxProcess;}
    
    G4double GetRange(G4double);
    

  private:
    G4double cutForGamma;
    G4double cutForElectron;
    G4double cutForPositron;
    
    Em5DetectorConstruction* pDet;        
    Em5PhysicsListMessenger* pMessenger;
        
    G4VPhysicsConstructor* particleList;
    G4VPhysicsConstructor* generalPhysicsList;    
    G4VPhysicsConstructor* emPhysicsList;
    G4String               emName;
    Em5StepMax*            stepMaxProcess;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

