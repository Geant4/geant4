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
/// \file electromagnetic/TestEm3/include/Run.hh
/// \brief Definition of the Run class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Run_h
#define Run_h 1

#include "DetectorConstruction.hh"

#include "G4Run.hh"
#include "globals.hh"
#include <map>

class DetectorConstruction;
class G4ParticleDefinition;
class G4Track;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run
{
  public:
    Run(DetectorConstruction*);
   ~Run();

  public:
    void SetPrimary(G4ParticleDefinition* particle, G4double energy);
      
    void FillPerEvent(G4int,G4double,G4double);
    
    void SumEnergyFlow (G4int plane, G4double Eflow);
    void SumLateralEleak(G4int cell, G4double Eflow);
    void AddChargedStep();
    void AddNeutralStep();
    void AddSecondaryTrack(const G4Track*);
    
    void SetEdepAndRMS(G4int, G4double, G4double, G4double);
    void SetApplyLimit(G4bool);
                
    virtual void Merge(const G4Run*);
    void EndOfRun();
     
  private:
    DetectorConstruction*  fDetector;
    G4ParticleDefinition*  fParticle;
    G4double  fEkin;
                           
    G4double fSumEAbs [kMaxAbsor], fSum2EAbs [kMaxAbsor]; 
    G4double fSumLAbs [kMaxAbsor], fSum2LAbs [kMaxAbsor];
    
    std::vector<G4double> fEnergyFlow;
    std::vector<G4double> fLateralEleak;
    std::vector<G4double> fEnergyDeposit[kMaxAbsor];
    
    G4double fChargedStep;
    G4double fNeutralStep;

    G4int  fN_gamma;
    G4int  fN_elec;
    G4int  fN_pos;
    
    G4double fEdeptrue [kMaxAbsor];
    G4double fRmstrue  [kMaxAbsor];
    G4double fLimittrue[kMaxAbsor];        
    G4bool fApplyLimit;  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

