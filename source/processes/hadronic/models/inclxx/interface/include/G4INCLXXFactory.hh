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
// INCL++ intra-nuclear cascade model
// Pekka Kaitaniemi, CEA and Helsinki Institute of Physics
// Davide Mancusi, CEA
// Alain Boudard, CEA
// Sylvie Leray, CEA
// Joseph Cugnon, University of Liege
//
// INCL++ revision: v5.0_rc3
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#ifndef G4INCLXXFactory_hh
#define G4INCLXXFactory_hh 1

#include "G4Nucleus.hh"
#include "G4HadProjectile.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4ParticleTable.hh"

#include "G4INCLConfig.hh"
#include "G4INCLCascade.hh"

/**
 * Build a new INCL++ instance with correct configuration and input information.
 */
class G4INCLXXFactory {
public:
  /**
   * Convert G4ParticleDefinition to corresponding INCL particle type
   */
  static G4INCL::ParticleType toINCLParticleType(const G4ParticleDefinition*);

  /**
   * Convert INCL particle type to corresponding G4ParticleDefinition
   */
  static const G4ParticleDefinition* fromINCLParticleType(G4INCL::ParticleType);

  /**
   * Create INCL projectile particle
   */
  static G4INCL::Particle* createProjectile(const G4HadProjectile &);

  /**
   * Create the INCL model initialized with the target information
   */
  static G4INCL::INCL* createModel(const G4Nucleus &);

  static G4DynamicParticle* toG4Particle(G4int A, G4int Z , G4double kinE, G4double px, G4double py, G4double pz);

  static G4ParticleDefinition* toG4ParticleDefinition (G4int A,
						       G4int Z);

  static G4double remnant4MomentumScaling(G4double mass,
					  G4double kineticE,
					  G4double px, G4double py, G4double pz);
protected:
  G4INCLXXFactory() {};
  ~G4INCLXXFactory() {};

};

#endif
