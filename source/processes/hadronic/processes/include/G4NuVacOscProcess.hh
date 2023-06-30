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
//
// Geant4 muon neutrino nucleus scattering Process -- header file
// 
// Created  from G4HadronElasticProcess 1.3.19 V. Grichine
//  
// Modified:
//
// 5.4.23 V.Grichine first implementation

// Class Description
// General process for neutrino nucleus 2->X scattering  
// Class Description - End

#ifndef G4NuVacOscProcess_h
#define G4NuVacOscProcess_h 1
 
#include "globals.hh"
#include "G4HadronicProcess.hh"

class G4ParticleDefinition;

class G4NuVacOscProcess : public G4HadronicProcess
{
public:

  G4NuVacOscProcess(const G4String& anEnvelopeName,
                    const G4String& procName = "nu-vacuum-oscillation");

  ~G4NuVacOscProcess() override = default;

  void InitParameters();
  
  G4double GetMeanFreePath(const G4Track &aTrack, G4double, G4ForceCondition *) override;
 
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep) override;

  G4int NuVacProbability( G4int aa, G4double Enu, G4double Lnu );

  void ProcessDescription(std::ostream& outFile) const override;

  // set internal limit
  void SetLowestEnergy(G4double e) { fLowestEnergy = e; }

  void SetBiasingFactor(G4double bf) {
    if(bf > 1.0) { fNuNuclTotXscBias = bf; fBiased = true; }
  }
  
  // hide assignment operator as private 
  G4NuVacOscProcess& operator=(const G4NuVacOscProcess &right)=delete;
  G4NuVacOscProcess(const G4NuVacOscProcess& )=delete;

private:

  G4bool fBiased{false};
  G4bool fAnti{false}; 
  G4bool fNormOrd{true};

  G4String fEnvelopeName;

  G4double fLowestEnergy;
  G4double fNuNuclTotXscBias {1.0};

  G4double fSin2t12, fSin2t23, fSin2t13; 
  G4double fDsm21, fDsm32, fdcp;

  G4complex fUdcp[3][3]; // vacuum PMNS matrix, Dirac nu, CP delta13
  G4double  fDms[3][3];  // dm2 matrix, 
  // involved oscillationg paricles
  G4ParticleDefinition* theNuE; 
  G4ParticleDefinition* theAntiNuE; 
  G4ParticleDefinition* theNuMu; 
  G4ParticleDefinition* theAntiNuMu; 
  G4ParticleDefinition* theNuTau; 
  G4ParticleDefinition* theAntiNuTau; 
};

#endif
