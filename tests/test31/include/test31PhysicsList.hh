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
// $Id: test31PhysicsList.hh,v 1.2 2002-10-28 09:57:36 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//
// ------------------------------------------------------------
//	GEANT 4 class header file 
// Class Description:
//      This class is an derived class of G4VPhysicsConstructor
//
// --------------------------------------------------------------------------- 
//	History
//        Created:       14.10.02  V.Ivanchenko provide modular list on base 
//                                 of old test31PhysicsList
//
//        Modified:      
// 
// ---------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef test31PhysicsList_h
#define test31PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "test31StepCut.hh"
#include "globals.hh"

class test31PhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class test31PhysicsList: public G4VModularPhysicsList
{
  public:

    test31PhysicsList();
   ~test31PhysicsList();

    void AddPhysicsList(const G4String& name);
    void SetCuts();

    void SetCutForGamma(G4double);
    void SetCutForElectron(G4double);
    void SetCutForPositron(G4double);
    void SetCutForAll(G4double val) {SetCutForGamma(val);
                                     SetCutForElectron(val);
                                     SetCutForPositron(val);}

    void SetMaxStep(G4double val)  {theStepCut.SetMaxStep(val);};

    void SetNuclearStopping(G4String st) {if(st == "on") nuclStop = true;
                                          else           nuclStop = false;}
    void SetBarkas(G4String st)          {if(st == "on") barkas = true;
                                          else           barkas = false;}
    void SetStoppingTable(G4String tab)  {table = tab;}

    void SetVerbose(G4int val) {verbose = val;};    
    G4int GetVerbose() const {return verbose;};    
       
  private:
    G4double cutForGamma;
    G4double cutForElectron; 
    G4double cutForPositron;
    G4double currentDefaultCut;
    G4double theMaxStep;
    G4bool   nuclStop;
    G4bool   barkas;
    G4String table;
    G4int    verbose;

    G4bool   emPhysicsListIsRegistered;
    G4bool   hadPhysicsListIsRegistered;
    
    test31PhysicsListMessenger* pMessenger;         

    test31StepCut  theStepCut;    

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

