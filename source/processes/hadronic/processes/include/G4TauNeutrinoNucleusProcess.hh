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
// Created  from G4MuNeutrinoNucleusProcess 1.11.22 V. Grichine
//  

// Class Description
// Hadronic inelastic process for tau neutrino nucleus 2->X scattering  
// Class Description - End

#ifndef G4TauNeutrinoNucleusProcess_h
#define G4TauNeutrinoNucleusProcess_h 1
 
#include "globals.hh"
#include "G4HadronicProcess.hh"

class G4ParticleDefinition;
class G4CrossSectionDataStore;
class G4LogicalVolume;
class G4TauNeutrinoNucleusTotXsc;
class G4SafetyHelper;

class G4TauNeutrinoNucleusProcess : public G4HadronicProcess
{
public:

  G4TauNeutrinoNucleusProcess(G4String anEnvelopeName , const G4String& procName = "tau-neutrino-nucleus");

  virtual ~G4TauNeutrinoNucleusProcess();
 
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, 
					  const G4Step& aStep) override;

  // initialise thresholds
  void PreparePhysicsTable(const G4ParticleDefinition&) override;

  // set internal limit
  virtual void SetLowestEnergy(G4double);

  void ProcessDescription(std::ostream& outFile) const override;

  void SetBiasingFactors(G4double bfCc, G4double bfNc);
  void SetBiasingFactor(G4double bf);
  G4double GetMeanFreePath(const G4Track &aTrack, G4double, G4ForceCondition *) override;
  
private:

  // hide assignment operator as private 
  G4TauNeutrinoNucleusProcess& operator=(const G4TauNeutrinoNucleusProcess &right);
  G4TauNeutrinoNucleusProcess(const G4TauNeutrinoNucleusProcess& );

  G4double lowestEnergy;
  G4bool   isInitialised, fBiased;
  G4LogicalVolume* fEnvelope;
  G4String fEnvelopeName;
  G4TauNeutrinoNucleusTotXsc* fTotXsc;
  G4double fNuNuclCcBias, fNuNuclNcBias, fNuNuclTotXscBias;
  G4SafetyHelper* safetyHelper;
};

#endif
