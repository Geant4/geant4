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
// GEANT4 class source file
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

#include "G4IonDEDXScalingICRU73.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"


// ###########################################################################

G4IonDEDXScalingICRU73::G4IonDEDXScalingICRU73(
                          G4int minAtomicNumberIon,
                          G4int maxAtomicNumberIon,
                          G4int atomicNumberReference, 
                          G4int massNumberReference) :
     minAtomicNumber( minAtomicNumberIon ),
     maxAtomicNumber( maxAtomicNumberIon ),
     excludedIon( false ),
     reference( 0 ),
     atomicNumberRef( atomicNumberReference ),
     massNumberRef( massNumberReference ) { 

}

// ###########################################################################

G4IonDEDXScalingICRU73::~G4IonDEDXScalingICRU73() {

}

// ###########################################################################

void G4IonDEDXScalingICRU73::CreateReferenceParticle() {

   
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4double excitationEnergy = 0.0;

  reference = 
    particleTable -> GetIon(atomicNumberRef, massNumberRef, excitationEnergy); 
  
  massRef = reference -> GetPDGMass();
  chargeRef = reference -> GetPDGCharge();

  atomicNumberRefPow23 = std::pow(G4double(atomicNumberRef), 2./3.);
}


// ###########################################################################


G4double G4IonDEDXScalingICRU73::ScalingFactorEnergy (
            const G4ParticleDefinition* particle) {    // Projectile (ion) 
                                                         
  G4double factor = 1.0;
 
  UpdateCache(particle);

  if(! excludedIon &&
     cacheAtomicNumber >= minAtomicNumber &&
     cacheAtomicNumber <= maxAtomicNumber) {

     if(reference == 0) CreateReferenceParticle();

     factor = cacheMassNumber * (massRef / cacheMass) / massNumberRef;
  }

  return factor;
}

// ###########################################################################

G4double G4IonDEDXScalingICRU73::ScalingFactorDEDX(
             const G4ParticleDefinition* particle,     // Projectile (ion) 
             const G4Material*,                        // Target material
             G4double kineticEnergy) {                 // Kinetic energy

  G4double factor = 1.0;

  UpdateCache(particle);

  if(! excludedIon &&
     cacheAtomicNumber >= minAtomicNumber &&
     cacheAtomicNumber <= maxAtomicNumber) {
      
      if(reference == 0) CreateReferenceParticle();

      G4double equilibriumCharge = EquilibriumCharge(cacheMass,
                                                     cacheCharge,
                                                     cacheAtomicNumberPow23,
                                                     kineticEnergy);

      G4double scaledKineticEnergy = kineticEnergy * (massRef / cacheMass);
      
      G4double equilibriumChargeRef = EquilibriumCharge(massRef,
                                                        chargeRef,
                                                        atomicNumberRefPow23,
                                                        scaledKineticEnergy);

      factor = equilibriumCharge * equilibriumCharge/ 
                ( equilibriumChargeRef * equilibriumChargeRef );
 }  

  return factor;
}

// ###########################################################################

G4int G4IonDEDXScalingICRU73::AtomicNumberBaseIon(
             G4int atomicNumberIon,           // Atomic number of ion 
             const G4Material*) {             // Target material

  G4int atomicNumber = atomicNumberIon;

  if(atomicNumberIon >= minAtomicNumber &&
     atomicNumberIon <= maxAtomicNumber) {
 
     G4bool ionIsExcluded = false;
     size_t nmb = excludedAtomicNumbers.size();
 
     for(size_t i = 0; i < nmb; i++) {
        
        if(atomicNumberIon == excludedAtomicNumbers[i]) 
           ionIsExcluded = true;
     }

     if(! ionIsExcluded) {
        if(reference == 0) CreateReferenceParticle();

        atomicNumber = atomicNumberRef;
     }
  }

  return atomicNumber;
}

// ###########################################################################
