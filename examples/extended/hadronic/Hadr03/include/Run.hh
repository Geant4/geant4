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

#include "G4Run.hh"
#include "G4VProcess.hh"
#include "globals.hh"
#include <map>

class DetectorConstruction;
class G4ParticleDefinition;
class G4HadronicProcessStore;
class G4Material;
class G4Element;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run
{
  public:
    Run(DetectorConstruction*);
   ~Run() override = default;

  public:
    void SetPrimary(G4ParticleDefinition* particle, G4double energy);
    void SetTargetXXX(G4bool);        
    void CountProcesses(G4VProcess* process);
    void SumTrack (G4double);    
    void CountNuclearChannel(G4String, G4double);                          
    void ParticleCount(G4String, G4double);
    void Balance(G4double);
    void CountGamma(G4int);
        
    void Merge(const G4Run*) override;      
    void EndOfRun(G4bool); 
   
  private:

    void PrintXS(const G4VProcess*, const G4Material*, const G4Element*,
                 G4HadronicProcessStore*, G4double dens,
                 G4double& sum1, G4double& sum2);
           
    struct ParticleData {
     ParticleData()
       : fCount(0), fEmean(0.), fEmin(0.), fEmax(0.) {}
     ParticleData(G4int count, G4double ekin, G4double emin, G4double emax)
       : fCount(count), fEmean(ekin), fEmin(emin), fEmax(emax) {}
     G4int     fCount;
     G4double  fEmean;
     G4double  fEmin;
     G4double  fEmax;
    };
    
    struct NuclChannel {
     NuclChannel()
       : fCount(0), fQ(0.) {}
     NuclChannel(G4int count, G4double Q)
       : fCount(count), fQ(Q) {}
     G4int     fCount;
     G4double  fQ;
    };
         
  private:
    DetectorConstruction* fDetector = nullptr;
    G4ParticleDefinition* fParticle = nullptr;
    G4double              fEkin = 0.;
        
    std::map<G4String,G4int> fProcCounter;            
    
    G4int fTotalCount = 0;       //all processes counter
    G4int fGammaCount = 0;       //nb of events with gamma
    G4double fSumTrack  = 0.;    //sum of trackLength
    G4double fSumTrack2 = 0.;    //sum of trackLength*trackLength
         
    std::map<G4String,NuclChannel>  fNuclChannelMap;    
    std::map<G4String,ParticleData> fParticleDataMap;

    G4bool   fTargetXXX = false;
    G4double fPbalance[3];
    G4int    fNbGamma[3];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

