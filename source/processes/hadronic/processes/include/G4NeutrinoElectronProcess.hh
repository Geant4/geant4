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
// Geant4 Neutrino Electron Scattering Process -- header file
// 
// Created  from G4HadronElasticProcess 15.12.17 V. Grichine
//  
// Modified:
//
// 2.2.18 V.Grichine PostStepDoIt implementation

// Class Description
// General process for neutrino electron 2->2 scattering  
// Class Description - End

#ifndef G4NeutrinoElectronProcess_h
#define G4NeutrinoElectronProcess_h 1
 
#include "globals.hh"
#include "G4HadronicProcess.hh"

class G4ParticleDefinition;
class G4CrossSectionDataStore;
class G4LogicalVolume;
class G4NeutrinoElectronTotXsc;
class G4SafetyHelper;

class G4NeutrinoElectronProcess : public G4HadronicProcess
{
public:

  G4NeutrinoElectronProcess(G4String anEnvelopeName , const G4String& procName = "neutrino-electron");

  virtual ~G4NeutrinoElectronProcess();
 
  virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack, 
					  const G4Step& aStep);

  // initialise thresholds
  virtual void PreparePhysicsTable(const G4ParticleDefinition&);

  // set internal limit
  virtual void SetLowestEnergy(G4double);

  virtual void ProcessDescription(std::ostream& outFile) const;

  void SetBiasingFactors(G4double bfCc, G4double bfNc);
  void SetBiasingFactor(G4double bf);
  G4double GetMeanFreePath(const G4Track &aTrack, G4double, G4ForceCondition *);
private:

  // hide assignment operator as private 
  G4NeutrinoElectronProcess& operator=(const G4NeutrinoElectronProcess &right);
  G4NeutrinoElectronProcess(const G4NeutrinoElectronProcess& );

  G4double lowestEnergy;
  G4bool   isInitialised, fBiased;
  G4LogicalVolume* fEnvelope;
  G4String fEnvelopeName;
  G4NeutrinoElectronTotXsc* fTotXsc;
  G4double fNuEleCcBias, fNuEleNcBias, fNuEleTotXscBias;
  G4SafetyHelper* safetyHelper;
};

#endif
