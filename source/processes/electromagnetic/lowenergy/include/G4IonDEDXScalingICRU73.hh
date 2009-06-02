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
// ===========================================================================
// GEANT4 class header file
//
// Class:                G4IonDEDXScalingICRU73
//
// Base class:           G4VIonDEDXScalingAlgorithm
//
// Author:               Anton Lechner (Anton.Lechner@cern.ch)
//
// First implementation: 10. 05. 2009
//
// Modifications: 
//
//
// Class description:
//    dE/dx scaling algorithm applied on top of ICRU 73 data (for ions not
//    covered by the ICRU 73 report) 
//
// Comments:
//
// =========================================================================== 

#ifndef G4IONDEDXSCALINGICRU73_HH
#define G4IONDEDXSCALINGICRU73_HH

#include "globals.hh"
#include "G4VIonDEDXScalingAlgorithm.hh"
#include <vector>


class G4IonDEDXScalingICRU73 : public G4VIonDEDXScalingAlgorithm {

 public:
   G4IonDEDXScalingICRU73(G4int minAtomicNumberIon,
                          G4int maxAtomicNumberIon,
                          G4int atomicNumberReference, 
                          G4int massNumberReference);
   ~G4IonDEDXScalingICRU73();

   // Function for scaling the kinetic energy (no scaling by default).
   // Returns scaling factor for a given ion.
   G4double ScalingFactorEnergy(
           const G4ParticleDefinition* particle);      // Projectile (ion) 
                                                         

   // Function for scaling the dE/dx value (no scaling by default).
   // Returns scaling factor for a given ion-material couple and
   // a given kinetic energy.
   G4double ScalingFactorDEDX(
             const G4ParticleDefinition* particle,     // Projectile (ion) 
             const G4Material*,                        // Target material
             G4double kineticEnergy);                  // Kinetic energy


   // Function for defining a base particle for dE/dx calculation.
   // (no base particle by default). Returns atomic number of base
   // particle.
   G4int AtomicNumberBaseIon(
             G4int atomicNumberIon,           // Atomic number of ion 
             const G4Material*);              // Target material

   void AddException(G4int atomicNumberIon);

 private:
   void UpdateCache(
             const G4ParticleDefinition* particle);    // Projectile (ion) 

   void CreateReferenceParticle();
 
   G4double EquilibriumCharge(
             G4double mass,                            // Ion mass
             G4double charge,                          // Ion charge
             G4double atomicNumberPow,                 // Power of atomic nmb  
             G4double kineticEnergy);                  // Kinetic energy

   // Scaling is only applied for ions with atomic numbers in the range
   // defined by the following parameters:
   G4int minAtomicNumber;
   G4int maxAtomicNumber;

   // Vector with atomic numbers excepted from scaling procedure
   std::vector<G4int> excludedAtomicNumbers;
   G4bool excludedIon;

   // Some properties of reference particles are stored for faster access
   G4ParticleDefinition* reference; 
   G4int atomicNumberRef;
   G4int massNumberRef;
   G4double atomicNumberRefPow23;
   G4double chargeRef;
   G4double massRef;

   // Some properties of projectiles are stored for faster access
   const G4ParticleDefinition* cacheParticle;
   G4int cacheMassNumber;
   G4int cacheAtomicNumber;
   G4double cacheAtomicNumberPow23;
   G4double cacheCharge;
   G4double cacheMass;
};

// ###########################################################################

inline void G4IonDEDXScalingICRU73::UpdateCache (
            const G4ParticleDefinition* particle) {    // Projectile (ion) 

  if(particle != cacheParticle) {

     cacheParticle = particle;
     cacheAtomicNumber = particle -> GetAtomicNumber();
     cacheMassNumber = particle -> GetAtomicMass();
     cacheCharge = particle -> GetPDGCharge();
     cacheMass = particle -> GetPDGMass();

     cacheAtomicNumberPow23 = std::pow(G4double(cacheAtomicNumber), 2./3.);

     excludedIon = false;
     size_t nmb = excludedAtomicNumbers.size();
     for(size_t i = 0; i < nmb; i++) {
        
        if(cacheAtomicNumber == excludedAtomicNumbers[i]) 
           excludedIon = true;
     }
  }
}

// ###########################################################################

inline G4double G4IonDEDXScalingICRU73::EquilibriumCharge(
                                    G4double mass, 
                                    G4double charge,
                                    G4double atomicNumberPow, 
                                    G4double kineticEnergy) {

  G4double totalEnergy  = kineticEnergy + mass;
  G4double betaSquared  = kineticEnergy * 
                  (totalEnergy + mass) / (totalEnergy * totalEnergy);

  G4double beta = std::sqrt( betaSquared );

  G4double velOverBohrVel = beta / CLHEP::fine_structure_const;

  G4double q1 = 1.0 - std::exp(-velOverBohrVel / atomicNumberPow);
 
  return q1 * charge;
}

// ###########################################################################

inline void G4IonDEDXScalingICRU73::AddException(G4int atomicNumber) {

  if(atomicNumber >= minAtomicNumber &&
     atomicNumber <= maxAtomicNumber) {

     excludedAtomicNumbers.push_back( atomicNumber );
  }
}

// ###########################################################################

#endif
