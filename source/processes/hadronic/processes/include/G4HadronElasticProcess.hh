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
// $Id: G4HadronElasticProcess.hh 90262 2015-05-22 09:03:22Z gcosmo $
//
// Geant4 Hadron Elastic Scattering Process -- header file
// 
// Created 26 July 2012 V.Ivanchenko from G4WHadronElasticProcess
//  
// Modified:
//

// Class Description
// General process for hadron nuclear elastic scattering  
// Class Description - End

#ifndef G4HadronElasticProcess_h
#define G4HadronElasticProcess_h 1
 
#include "globals.hh"
#include "G4HadronicProcess.hh"

class G4ParticleDefinition;
class G4CrossSectionDataStore;
class G4VCrossSectionRatio;

class G4HadronElasticProcess : public G4HadronicProcess
{
public:

  G4HadronElasticProcess(const G4String& procName = "hadElastic");

  virtual ~G4HadronElasticProcess();
 
  virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack, 
					  const G4Step& aStep);

  // initialise thresholds
  virtual void PreparePhysicsTable(const G4ParticleDefinition&);

  // set internal limit
  virtual void SetLowestEnergy(G4double);

  // obsolete method - will be removed
  virtual void SetLowestEnergyNeutron(G4double);

  virtual void ProcessDescription(std::ostream& outFile) const;

  // enable sampling of low-mass diffraction process
  void SetDiffraction(G4HadronicInteraction*, G4VCrossSectionRatio*);

private:

  // hide assignment operator as private 
  G4HadronElasticProcess& operator=(const G4HadronElasticProcess &right);
  G4HadronElasticProcess(const G4HadronElasticProcess& );

  G4double lowestEnergy;
  G4bool   isInitialised;
  G4HadronicInteraction* fDiffraction;
  G4VCrossSectionRatio*  fDiffractionRatio;
};

#endif
