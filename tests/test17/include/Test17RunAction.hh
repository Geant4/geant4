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
#include "G4Proton.hh"
#include "G4Alpha.hh"
#include "G4IonC12.hh"
#include "G4IonAr40.hh"
#include "G4IonFe56.hh"
#include "G4Electron.hh"
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
    const G4Electron* theElectron ;

    G4double kinEnergy0;
    const G4ParticleDefinition* part0;
};

#endif

