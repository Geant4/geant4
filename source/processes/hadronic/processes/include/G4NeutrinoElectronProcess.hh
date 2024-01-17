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
// Geant4 Neutrino Electron Scattering Process -- header file
// 
// Created  from G4HadronElasticProcess 15.12.17 V. Grichine
//  
// Modified:
//
// 02.02.18 V.Grichine PostStepDoIt implementation

// Class Description
// General process for neutrino electron 2->2 scattering  
// Class Description - End

#ifndef G4NeutrinoElectronProcess_h
#define G4NeutrinoElectronProcess_h 1
 
#include "globals.hh"
#include "G4HadronicProcess.hh"

class G4ParticleDefinition;
class G4CrossSectionDataStore;
class G4NeutrinoElectronTotXsc;
class G4SafetyHelper;

class G4NeutrinoElectronProcess : public G4HadronicProcess
{
public:

  G4NeutrinoElectronProcess(const G4String& anEnvelopeName,
                            const G4String& procName = "nuElectron");

  ~G4NeutrinoElectronProcess() override = default;
  
  G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double previousStepSize,
                             G4ForceCondition* condition
                            ) override;

  G4double GetMeanFreePath(const G4Track &aTrack,
                           G4double, G4ForceCondition*) override;
  
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, 
       				  const G4Step& aStep) override;

  void ProcessDescription(std::ostream& outFile) const override;

  // set internal parameters
  void SetLowestEnergy(G4double);
  void SetBiasingFactors(G4double bfCc, G4double bfNc);
  void SetBiasingFactor(G4double bf);
  
  // hide assignment operator as private 
  G4NeutrinoElectronProcess& operator=
  (const G4NeutrinoElectronProcess &right) = delete;
  G4NeutrinoElectronProcess(const G4NeutrinoElectronProcess&) = delete;

private:

  G4NeutrinoElectronTotXsc* fTotXsc;
  G4SafetyHelper* safetyHelper;
  G4double lowestEnergy;
  G4double fNuEleCcBias{1.0};
  G4double fNuEleNcBias{1.0};
  G4double fNuEleTotXscBias{1.0};
  G4String fEnvelopeName;
};

#endif
