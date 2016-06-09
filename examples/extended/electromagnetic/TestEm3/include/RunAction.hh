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
// $Id: RunAction.hh,v 1.13 2004/06/15 11:39:57 maire Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "DetectorConstruction.hh"

#include "G4UserRunAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction;
class RunActionMessenger;
class HistoManager;

class G4Run;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
  public:

    RunAction(DetectorConstruction*, PrimaryGeneratorAction*, HistoManager*);
   ~RunAction();

    void BeginOfRunAction(const G4Run*);
    void   EndOfRunAction(const G4Run*);

    void fillPerEvent(G4int,G4double,G4double,G4double);
    
    void PrintDedxTables();
    
     // Acceptance parameters
     void     SetEdepAndRMS(G4int, G4ThreeVector);
     G4double GetAverageEdep(G4int i) const    {return edeptrue[i];};
     G4double GetRMSEdep(G4int i) const        {return rmstrue[i];};
     G4double GetLimitEdep(G4int i) const      {return limittrue[i];};
         
  private:

    G4double sumEAbs [MaxAbsor], sum2EAbs [MaxAbsor]; 
    G4double sumLAbs [MaxAbsor], sum2LAbs [MaxAbsor];
    G4double sumEleav[MaxAbsor], sum2Eleav[MaxAbsor];           

    DetectorConstruction*   Detector;
    PrimaryGeneratorAction* Primary;    
    RunActionMessenger*     runMessenger;
    HistoManager*           histoManager;

    G4double edeptrue [MaxAbsor];
    G4double rmstrue  [MaxAbsor];
    G4double limittrue[MaxAbsor];                
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

