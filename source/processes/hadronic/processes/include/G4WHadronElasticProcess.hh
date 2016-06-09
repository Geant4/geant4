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
// $Id$
//
// Geant4 Hadron Elastic Scattering Process -- header file
// 
// Created 21 April 2006 V.Ivanchenko
//  
// Modified:
// 26.09.06 V.Ivanchenko add lowestEnergy
// 17.10.11 V.Ivanchenko move to process sub-library and cleaned up
//

// Class Description
// Process for hadron nuclear elastic scattering using optimal 
// combination of Geant4 models
// Class Description - End

#ifndef G4WHadronElasticProcess_h
#define G4WHadronElasticProcess_h 1
 
#include "globals.hh"
#include "G4HadronicProcess.hh"

class G4ParticleDefinition;
class G4CrossSectionDataStore;

class G4WHadronElasticProcess : public G4HadronicProcess
{
public:

  G4WHadronElasticProcess(const G4String& procName = "hadElastic");

  virtual ~G4WHadronElasticProcess();
 
  virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack, 
					  const G4Step& aStep);

  inline void SetLowestEnergy(G4double);

  inline void SetLowestEnergyNeutron(G4double);

  virtual void Description() const;

private:

  // hide assignment operator as private 
  G4WHadronElasticProcess& operator=(const G4WHadronElasticProcess &right);
  G4WHadronElasticProcess(const G4WHadronElasticProcess& );

  const G4ParticleDefinition* theNeutron;

  G4double        lowestEnergy;
  G4double        lowestEnergyNeutron;
};

inline void 
G4WHadronElasticProcess::SetLowestEnergy(G4double val)
{
  lowestEnergy = val;
}

inline void 
G4WHadronElasticProcess::SetLowestEnergyNeutron(G4double val)
{
  lowestEnergyNeutron = val;
}


#endif
