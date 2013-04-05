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
/// \file electromagnetic/TestEm1/include/RunAction.hh
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
#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;
class DetectorConstruction;
class PrimaryGeneratorAction;
class HistoManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
  public:
    RunAction(DetectorConstruction*, PrimaryGeneratorAction*);
   ~RunAction();

  public:
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);
    
    void CountTraks0(G4int nt) { fNbOfTraks0 += nt;}
    void CountTraks1(G4int nt) { fNbOfTraks1 += nt;}
    void CountSteps0(G4int ns) { fNbOfSteps0 += ns;}
    void CountSteps1(G4int ns) { fNbOfSteps1 += ns;}
    void CountProcesses(G4String procName) { fProcCounter[procName]++;};
    
    void AddEdep(G4double val)     { fEdep += val;}
    void AddTrueRange (G4double l) { fTrueRange += l; fTrueRange2 += l*l;};
    void AddProjRange (G4double x) { fProjRange += x; fProjRange2 += x*x;};
    void AddTransvDev (G4double y) { fTransvDev += y; fTransvDev2 += y*y;};    
                                  
  private:
    DetectorConstruction*   fDetector;
    PrimaryGeneratorAction* fPrimary;
    HistoManager*           fHistoManager;
      
    G4int           fNbOfTraks0, fNbOfTraks1;
    G4int           fNbOfSteps0, fNbOfSteps1;
    G4double        fEdep;
    G4double        fTrueRange, fTrueRange2;             
    G4double        fProjRange, fProjRange2;
    G4double        fTransvDev, fTransvDev2;
    std::map<G4String,G4int>    fProcCounter;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

