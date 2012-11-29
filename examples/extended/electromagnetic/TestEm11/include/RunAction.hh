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
/// \file electromagnetic/TestEm11/include/RunAction.hh
/// \brief Definition of the RunAction class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "DetectorConstruction.hh"

class G4Run;
class PhysicsList;
class PrimaryGeneratorAction;
class HistoManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
  public:
    RunAction(DetectorConstruction*, PhysicsList*,PrimaryGeneratorAction*);
   ~RunAction();

  public:
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);
    
    void AddEdep (G4double e)        { fEdeposit  += e; fEdeposit2  += e*e;};
    void AddTrackLength (G4double t) { fTrackLen  += t; fTrackLen2  += t*t;};
    void AddProjRange   (G4double x) { fProjRange += x; fProjRange2 += x*x;};
    void AddStepSize    (G4int nb, G4double s)
                                     { fNbOfSteps += nb; fNbOfSteps2 += nb*nb;
                                       fStepSize  += s ; fStepSize2  += s*s;  };
    void AddTrackStatus (G4int i)    { fStatus[i]++ ;};
                                       
    G4double GetCsdaRange (G4int i) {return fCsdaRange[i];};
    G4double GetXfrontNorm(G4int i) {return fXfrontNorm[i];};    
           
  private:
    DetectorConstruction*   fDetector;
    PhysicsList*            fPhysics;
    PrimaryGeneratorAction* fKinematic;
    HistoManager*           fHistoManager;
    
    G4double                fEdeposit,  fEdeposit2;
    G4double                fTrackLen,  fTrackLen2;
    G4double                fProjRange, fProjRange2;
    G4int                   fNbOfSteps, fNbOfSteps2;
    G4double                fStepSize,  fStepSize2;
    G4int                   fStatus[3];
    
    G4double                fCsdaRange[MaxAbsor];
    G4double                fXfrontNorm[MaxAbsor];    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

