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
/// \file HadronicGenerator.hh
/// \brief Definition of the HadronicGenerator class
// 
//------------------------------------------------------------------------
// Class: HadronicGenerator
// Author: Alberto Ribon (CERN EP/SFT)
// Date: May 2020
//
// This class shows how to use Geant4 as a generator for simulating
// inelastic hadron-nuclear interactions.
// Some of the most used hadronic models are currently supported in
// this class:
// - the hadronic string models Fritiof (FTF) and Quark-Gluon-String (QGS)
//   coupled with Precompound/de-excitation
// - the intranuclear cascade models: Bertini (BERT), Binary Cascade (BIC),
//                                    and Liege (INCL)
// Combinations of two models - in a transition energy interval, with a
// linear probability as a function of the energy - are also available to
// "mimic" the transition between hadronic models as in the most common
// Geant4 reference physics lists.
//
// The current version of this class does NOT support:
// -  hadron elastic interactions
// -  neutron capture and fission
// -  precise low-energy inelastic interactions of neutrons and
//    charged particles (i.e. ParticleHP)
// -  gamma/lepton-nuclear inelastic interactions
//
// This class does NOT use the Geant4 run-manager, and therefore should
// be usable in a multi-threaded application, with one instance of this
// class in each thread.
// 
// This class has been inspired by test30 (whose author is Vladimir
// Ivanchenko), with various simplifications and restricted to hadronic
// inelastic interactions.
//------------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HadronicGenerator_h
#define HadronicGenerator_h 1

#include <iomanip>
#include "globals.hh"
#include "G4ios.hh"
#include "G4ThreeVector.hh"
#include <map>
#include "G4HadronicProcess.hh"

class G4ParticleDefinition;
class G4VParticleChange;
class G4ParticleTable;
class G4Material;
class G4HadronicInteraction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HadronicGenerator {
  // This class provides the functionality of a "hadronic generator"
  // for Geant4 final-state inelastic hadronic collisions.
  // Only a few of the available Geant4 final-state hadronic inelastic
  // "physics cases" are currently available in this class - but it can
  // be extended to other cases if needed.
  // It is important to notice that this class does NOT use the Geant4
  // run-manager, so it should work fine in a multi-threaded environment,
  // with a separate instance of this class in each thread.
  public:

    explicit HadronicGenerator( const G4String physicsCase = "FTFP_BERT_ATL" );
    // Currently supported final-state hadronic inelastic "physics cases":
    // -  Hadronic models :        BERT, BIC, IonBIC, INCL, FTFP, QGSP
    // -  "Physics-list proxies" : FTFP_BERT_ATL (default), FTFP_BERT,
    //                             QGSP_BERT, QGSP_BIC, FTFP_INCLXX
    //    (i.e. they are not real, complete physics lists - for instance
    //     they do not have: transportation, electromagnetic physics,
    //     hadron elastic scattering, neutron fission and capture, etc. -
    //     however, they cover all hadron types and all energies by
    //     combining different hadronic models, i.e. there are transitions
    //     between two hadronic models in well-defined energy intervals,
    //     e.g. "FTFP_BERT" has the transition between BERT and FTFP
    //     hadronic models; moreover, the transition intervals used in
    //     our "physics cases"might not be the same as in the corresponding
    //     physics lists).

    ~HadronicGenerator();

    inline G4bool IsPhysicsCaseSupported() const;
    // Returns "true" if the physicsCase is supported; "false" otherwise. 
  
    G4bool IsApplicable( const G4String &nameProjectile, const G4double projectileEnergy ) const;
    G4bool IsApplicable( G4ParticleDefinition* projectileDefinition,
                         const G4double projectileEnergy ) const;
    // Returns "true" if the specified projectile (either by name or particle definition)
    // of given energy is applicable, "false" otherwise.

    G4VParticleChange* GenerateInteraction( const G4String &nameProjectile,
                                            const G4double projectileEnergy,
                                            const G4ThreeVector &projectileDirection ,
                                            G4Material* targetMaterial );
    G4VParticleChange* GenerateInteraction( G4ParticleDefinition* projectileDefinition,
                                            const G4double projectileEnergy,
                                            const G4ThreeVector &projectileDirection ,
                                            G4Material* targetMaterial );
    // This is the main method provided by the class:
    // in input it receives the projectile (either by name or particle definition),
    // its energy, its direction and the target material, and it returns one sampled
    // final-state of the inelastic hadron-nuclear collision as modelled by the
    // final-state hadronic inelastic "physics case" specified in the constructor.
    // If the required hadronic collision is not possible, then the method returns
    // immediately an empty "G4VParticleChange", i.e. without secondaries produced.

    inline G4HadronicProcess* GetHadronicProcess() const;
    inline G4HadronicInteraction* GetHadronicInteraction() const;
    // Returns the hadronic process and the hadronic interaction, respectively,
    // that handled the last call of "GenerateInteraction".

    G4double GetImpactParameter() const;
    G4int GetNumberOfTargetSpectatorNucleons() const;
    G4int GetNumberOfProjectileSpectatorNucleons() const;
    G4int GetNumberOfNNcollisions() const;
    // In the case of hadronic interactions handled by the FTF model, returns,
    // respectively, the impact parameter, the number of target/projectile
    // spectator nucleons, and the number of nucleon-nucleon collisions,
    // else, returns a negative value (-999).

  private:

    G4String fPhysicsCase;
    G4bool fPhysicsCaseIsSupported;
    G4HadronicProcess* fLastHadronicProcess;
    G4ParticleTable* fPartTable;
    std::map< G4ParticleDefinition*, G4HadronicProcess* > fProcessMap;  
};


inline G4bool HadronicGenerator::IsPhysicsCaseSupported() const {
  return fPhysicsCaseIsSupported;
}


inline G4HadronicProcess* HadronicGenerator::GetHadronicProcess() const {
  return fLastHadronicProcess;
}


inline G4HadronicInteraction* HadronicGenerator::GetHadronicInteraction() const {
  return fLastHadronicProcess == nullptr ? nullptr
                                         : fLastHadronicProcess->GetHadronicInteraction();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
