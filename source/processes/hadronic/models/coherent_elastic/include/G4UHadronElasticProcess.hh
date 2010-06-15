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
// $Id: G4UHadronElasticProcess.hh,v 1.12 2010-06-15 15:24:34 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Geant4 Hadron Elastic Scattering Process -- header file
// 
// Created 21 April 2006 V.Ivanchenko
//  
// Modified:
// 26.09.06 V.Ivanchenko add lowestEnergy
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
#include "G4StableIsotopes.hh"

class G4VQCrossSection;
class G4ParticleDefinition;
class G4CrossSectionDataStore;

class G4UHadronElasticProcess : public G4HadronicProcess
{
public:

  G4UHadronElasticProcess(const G4String& procName = "hElastic", 
			  G4double elow = 19.*MeV);

  virtual ~G4UHadronElasticProcess();
 
  virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack, 
					  const G4Step& aStep);

  virtual G4bool IsApplicable(const G4ParticleDefinition& aParticleType);

  virtual void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);

  virtual void DumpPhysicsTable(const G4ParticleDefinition& aParticleType);

  virtual G4double GetMeanFreePath(const G4Track&, G4double, 
				   G4ForceCondition*);

  virtual G4double GetMicroscopicCrossSection(const G4DynamicParticle*,
					      const G4Element*,
					      G4double aTemp);

  void SetQElasticCrossSection(G4VQCrossSection*);

private:

  G4StableIsotopes            theDefaultIsotopes;
  G4VQCrossSection*           pCManager;
  G4VQCrossSection*           nCManager;
  const G4ParticleDefinition* theProton;
  const G4ParticleDefinition* theNeutron;
  const G4ParticleDefinition* theParticle;

  G4CrossSectionDataStore* store;
  G4Nucleus                targetNucleus;

  G4double        xsec[40];
  G4double        xsecH[4];
  G4double        cross;
  G4double        thEnergy;
  G4double        lowestEnergy;

  G4int    pPDG;
  G4bool   first;
};

#endif
