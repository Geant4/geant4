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
/// \file radioactivedecay/rdecay01/include/RunAction.hh
/// \brief Definition of the RunAction class
//
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

class G4Run;
class HistoManager;
class PrimaryGeneratorAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
  public:
    RunAction(PrimaryGeneratorAction*);
   ~RunAction();
   
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);
    
    void ParticleCount(G4String, G4double);
    void Balance(G4double,G4double);
    void EventTiming(G4double);
    void PrimaryTiming(G4double);
    
  private:
    HistoManager*           fHistoManager;
    PrimaryGeneratorAction* fPrimary;
    
    std::map<G4String,G4int> fParticleCount;
    std::map<G4String,G4double> fEmean;
    std::map<G4String,G4double> fEmin;
    std::map<G4String,G4double> fEmax;
    G4int    fDecayCount, fTimeCount;
    G4double fEkinTot[3];
    G4double fPbalance[3];
    G4double fEventTime[3];
    G4double fPrimaryTime;                        
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

