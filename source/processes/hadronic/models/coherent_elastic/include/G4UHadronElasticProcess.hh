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
//
// Geant4 Hadron Elastic Scattering Process -- header file
// 
// Created 21 April 2006 V.Ivanchenko
//  
// Modified:
//

// Class Description
// Process for hadron nuclear elastic scattering using optimal 
// combination of Geant4 models
// Class Description - End

#ifndef G4UHadronElasticProcess_h
#define G4UHadronElasticProcess_h 1
 
#include "globals.hh"
#include "G4HadronicProcess.hh"
#include "G4Nucleus.hh"

class G4VQCrossSection;
class G4ParticleDefinition;
class G4CrossSectionDataStore;

class G4UHadronElasticProcess : public G4HadronicProcess
{
public:

  G4UHadronElasticProcess(const G4String& processName = "hElastic", G4bool fl = true);

  virtual ~G4UHadronElasticProcess();
 
  virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack, 
					  const G4Step& aStep);

  virtual G4bool IsApplicable(const G4ParticleDefinition& aParticleType);

  virtual void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);

  virtual void DumpPhysicsTable(const G4ParticleDefinition& aParticleType);

  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

  virtual G4double GetMicroscopicCrossSection(const G4DynamicParticle* aParticle,
					      const G4Element* anElement,
					      G4double aTemp);

private:

  G4VQCrossSection*           qCManager;
  const G4ParticleDefinition* theProton;
  const G4ParticleDefinition* theNeutron;
  const G4ParticleDefinition* theParticle;

  G4CrossSectionDataStore* store;
  G4Nucleus                targetNucleus;

  G4double        xsec[40];
  G4double        xsecH[2];
  G4double        cross;
  G4double        thEnergy;

  G4int    pPDG;
  G4bool   flagHP;
};

#endif
