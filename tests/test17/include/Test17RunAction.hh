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
// Class Description:
// The user has a possibility to define and to fill his histograms in this class.
// Class Description - end
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Test17RunAction_h
#define Test17RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "G4Run.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Test17RunAction : public G4UserRunAction
{
public: // Without description

  Test17RunAction();
  ~Test17RunAction();

public: // With description
 
  void BeginOfRunAction(const G4Run*);

  void   EndOfRunAction(const G4Run*);

public: // Without description

  G4int RunID() const {return run->GetRunID();};
  G4int EventNo() const {return nEvents;};
  void CountEvent()  {nEvents++;};
  void CountParticles(G4int nc, G4int nn)
                       {nCharged += nc; nNeutral += nn;};
  void AddEdep(G4double val) {edepTot += val;}; 
                       
  void AddTrackLength(G4double val)
                       {length += val; length2 += val*val;}; 
  void EndOfTrackCharged(G4double val)
                       {xend += val; xend2 += val*val;}; 

  void FillEn(G4double e) {kinEnergy0 = e;};
  void FillDef(const G4ParticleDefinition* p) {part0 = p;};


private:

  const G4Run* run;

  G4double edepTot;
  G4double length;
  G4double length2;
  G4double xend;
  G4double xend2;
  G4int nEvents;
  G4int nCharged;
  G4int nNeutral;

  G4ParticleDefinition* theProton ;
  const G4ParticleDefinition* theElectron ;
  const G4ParticleDefinition* part0;

  G4double kinEnergy0;
};

#endif

