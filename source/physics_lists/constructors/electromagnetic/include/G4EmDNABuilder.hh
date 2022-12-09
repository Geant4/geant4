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
// Geant4 header G4EmDNABuilder
//
// Author V.Ivanchenko 18.02.2022
//
// Utilities to build DNA physics
//

#ifndef G4EmDNABuilder_h
#define G4EmDNABuilder_h 1

#include "globals.hh"
// DNA processes
#include "G4DNAElectronSolvation.hh"
#include "G4DNAElastic.hh"
#include "G4DNAIonisation.hh"
#include "G4DNAExcitation.hh"
#include "G4DNAVibExcitation.hh"
#include "G4DNAAttachment.hh"
#include "G4DNAChargeDecrease.hh"
#include "G4DNAChargeIncrease.hh"
#include "G4LowECapture.hh"

enum G4EmDNAMscModelType
{
  dnaUrban = 0,
  dnaWVI,
  dnaGS
};

class G4ParticleDefinition;
class G4Region;

class G4EmDNABuilder
{
public:

  static void ConstructDNAParticles();

  static void ConstructStandardEmPhysics(const G4double emin_electron,
                                         const G4double emin_proton,
                                         const G4double emin_alpha,
                                         const G4double emin_ion,
                                         const G4EmDNAMscModelType mscType,
                                         const G4bool fast);

  static void ConstructDNAElectronPhysics(const G4double emaxDNA,
                                          const G4int opt,
                                          const G4bool fast,
                                          const G4bool stationary,
                                          const G4Region* reg = nullptr);

  static void ConstructDNAProtonPhysics(const G4double e1DNA,
                                        const G4double emaxDNA,
                                        const G4int opt,
                                        const G4bool fast,
                                        const G4bool stationary,
                                        const G4Region* reg = nullptr);

  static void ConstructDNAIonPhysics(const G4double emax,
                                     const G4bool stationary,
                                     const G4Region* reg = nullptr);

  static void ConstructDNALightIonPhysics(G4ParticleDefinition* part,
                                          const G4int charge,
                                          const G4int opt,
                                          const G4double emax,
                                          const G4bool fast,
                                          const G4bool stationary,
                                          const G4Region* reg = nullptr);

  static G4DNAElectronSolvation* FindOrBuildElectronSolvation();

  static G4DNAElastic* 
  FindOrBuildElastic(G4ParticleDefinition* part, const G4String& name);

  static G4DNAExcitation*
  FindOrBuildExcitation(G4ParticleDefinition* part, const G4String& name);

  static G4DNAVibExcitation*
  FindOrBuildVibExcitation(G4ParticleDefinition* part, const G4String& name);

  static G4DNAIonisation*
  FindOrBuildIonisation(G4ParticleDefinition* part, const G4String& name);

  static G4DNAAttachment*
  FindOrBuildAttachment(G4ParticleDefinition* part, const G4String& name);

  static G4DNAChargeDecrease*
  FindOrBuildChargeDecrease(G4ParticleDefinition* part, const G4String& name);

  static G4DNAChargeIncrease*
  FindOrBuildChargeIncrease(G4ParticleDefinition* part, const G4String& name);

  static G4LowECapture*
  FindOrBuildCapture(const G4double elim, G4ParticleDefinition* part);

private:

  static void StandardHadronPhysics(G4ParticleDefinition*,
			            const G4double lowELimitForMSC,
			            const G4double lowELimitForIoni,
			            const G4double maxEnergy,
				    const G4EmDNAMscModelType mscType,
                                    const G4bool isIon);
};

#endif
