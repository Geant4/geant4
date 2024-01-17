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
/// \file Run.hh
/// \brief Definition of the Run class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Run_h
#define Run_h 1

#include "DetectorConstruction.hh"

#include "G4Run.hh"
#include "G4VProcess.hh"
#include "globals.hh"
#include <map>

class DetectorConstruction;
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run
{
  public:
    Run(DetectorConstruction* detector);
   ~Run() override = default;

  public:
    void SetPrimary(G4ParticleDefinition* particle, G4double energy);
    void CountProcesses(const G4VProcess* process);
    void ParticleCount(G4int, G4String, G4double, G4double);
    void AddEdep (G4int, G4double);
    void AddTotEdep     (G4double);
    void AddEleak       (G4double);
    void AddEtotal      (G4double);    
    void AddTrackStatus (G4int);

    void Merge(const G4Run*) override;
    void EndOfRun();
    
private:
  struct ParticleData {
   ParticleData()
     : fCount(0), fEmean(0.), fEmin(0.), fEmax(0.), fTmean(-1.) {}
   ParticleData(G4int count, G4double ekin, G4double emin, G4double emax,
                 G4double meanLife)
     : fCount(count),fEmean(ekin),fEmin(emin),fEmax(emax),fTmean(meanLife) {}
   G4int     fCount;
   G4double  fEmean;
   G4double  fEmin;
   G4double  fEmax;
   G4double  fTmean;
  };
  
  private:
    DetectorConstruction*  fDetector = nullptr;
    G4ParticleDefinition*  fParticle = nullptr;
    G4double  fEkin = 0.; 

    G4int      fStatus[3] = {0, 0, 0};
    
    G4double   fEdeposit[kMaxAbsor], fEmin[kMaxAbsor], fEmax[kMaxAbsor];
    G4double   fTotEdep[3] = {0.,0.,0.};
    G4double   fEleak[3]   = {0.,0.,0.};
    G4double   fEtotal[3]  = {0.,0.,0.};
   
    std::map<G4String,G4int>        fProcCounter;
    std::map<G4String,ParticleData> fParticleDataMap[kMaxAbsor];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

