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
// Modifications: 06. 08. 2009 - Minor bug fix (initialization of cache) AL
//                12. 11. 2009 - Moved all decision logic concerning ICRU 73
//                               scaling for heavy ions into this class.
//                               Adapting ScalingFactorEnergy class according
//                               to changes in base class (AL).
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
#include "G4IonTable.hh"
#include "G4Material.hh"


// ###########################################################################

G4IonDEDXScalingICRU73::G4IonDEDXScalingICRU73(
                          G4int minAtomicNumberIon,
                          G4int maxAtomicNumberIon) :
     minAtomicNumber( minAtomicNumberIon ),
     maxAtomicNumber( maxAtomicNumberIon ),
     referencePrepared( false ),
     atomicNumberRefFe( 26 ),
     massNumberRefFe( 56 ),
     atomicNumberRefPow23Fe( 0 ),
     chargeRefFe( 0 ),
     massRefFe( 0 ),
     atomicNumberRefAr( 18 ),
     massNumberRefAr( 40 ),
     atomicNumberRefPow23Ar( 0 ),
     chargeRefAr( 0 ),
     massRefAr( 0 ),
     useFe( true ),
     cacheParticle( 0 ),
     cacheMassNumber( 0 ),
     cacheAtomicNumber( 0 ),
     cacheAtomicNumberPow23( 0 ),
     cacheCharge( 0 ),
     cacheMass( 0 ),
     cacheMaterial( 0 ) { 

}

// ###########################################################################

G4IonDEDXScalingICRU73::~G4IonDEDXScalingICRU73() {
}

// ###########################################################################

void G4IonDEDXScalingICRU73::CreateReferenceParticles() {
  
  G4IonTable* ionTable = G4IonTable::GetIonTable(); 

  massRefFe = ionTable->GetIonMass(atomicNumberRefFe,massNumberRefFe);
  massRefAr = ionTable->GetIonMass(atomicNumberRefAr,massNumberRefAr);

  chargeRefFe = G4double(atomicNumberRefFe)*CLHEP::eplus;
  chargeRefAr = G4double(atomicNumberRefAr)*CLHEP::eplus;

  atomicNumberRefPow23Fe = std::pow(G4double(atomicNumberRefFe), 2./3.);
  atomicNumberRefPow23Ar = std::pow(G4double(atomicNumberRefAr), 2./3.);

  referencePrepared = true;
}


// ###########################################################################


G4double G4IonDEDXScalingICRU73::ScalingFactorEnergy (
            const G4ParticleDefinition* particle,     // Projectile (ion) 
            const G4Material* material) {             // Target material
                                                         
  G4double factor = 1.0;
 
  UpdateCacheParticle(particle);
  UpdateCacheMaterial(material);

  if(cacheAtomicNumber >= minAtomicNumber &&
     cacheAtomicNumber <= maxAtomicNumber &&
     cacheAtomicNumber != atomicNumberRefFe &&
     cacheAtomicNumber != atomicNumberRefAr) {

     if(!referencePrepared) CreateReferenceParticles();

     if( useFe )
         factor = cacheMassNumber * (massRefFe / cacheMass) / massNumberRefFe;
     else
         factor = cacheMassNumber * (massRefAr / cacheMass) / massNumberRefAr;
  }

  return factor;
}

// ###########################################################################

G4double G4IonDEDXScalingICRU73::ScalingFactorDEDX(
             const G4ParticleDefinition* particle,     // Projectile (ion) 
             const G4Material* material,               // Target material
             G4double kineticEnergy) {                 // Kinetic energy

  G4double factor = 1.0;

  UpdateCacheParticle(particle);
  UpdateCacheMaterial(material);

  if(cacheAtomicNumber >= minAtomicNumber &&
     cacheAtomicNumber <= maxAtomicNumber &&
     cacheAtomicNumber != atomicNumberRefFe &&
     cacheAtomicNumber != atomicNumberRefAr) {
      
      if(!referencePrepared) CreateReferenceParticles();

      if( useFe ) {

         G4double equilibriumCharge = EquilibriumCharge(cacheMass,
                                                     cacheCharge,
                                                     cacheAtomicNumberPow23,
                                                     kineticEnergy);

         G4double scaledKineticEnergy = kineticEnergy * (massRefFe / cacheMass);
      
         G4double equilibriumChargeRefFe = EquilibriumCharge(massRefFe,
                                                        chargeRefFe,
                                                        atomicNumberRefPow23Fe,
                                                        scaledKineticEnergy);

         factor = equilibriumCharge * equilibriumCharge/ 
                   ( equilibriumChargeRefFe * equilibriumChargeRefFe );

      }
      else {

         G4double equilibriumCharge = EquilibriumCharge(cacheMass,
                                                     cacheCharge,
                                                     cacheAtomicNumberPow23,
                                                     kineticEnergy);

         G4double scaledKineticEnergy = kineticEnergy * (massRefAr / cacheMass);
      
         G4double equilibriumChargeRefAr = EquilibriumCharge(massRefAr,
                                                        chargeRefAr,
                                                        atomicNumberRefPow23Ar,
                                                        scaledKineticEnergy);

         factor = equilibriumCharge * equilibriumCharge/ 
                   ( equilibriumChargeRefAr * equilibriumChargeRefAr );

      }
  }  

  return factor;
}

// ###########################################################################

G4int G4IonDEDXScalingICRU73::AtomicNumberBaseIon(
             G4int atomicNumberIon,           // Atomic number of ion 
             const G4Material* material) {    // Target material

  UpdateCacheMaterial(material);

  G4int atomicNumber = atomicNumberIon;

  if(atomicNumberIon >= minAtomicNumber &&
     atomicNumberIon <= maxAtomicNumber &&
     atomicNumberIon != atomicNumberRefFe &&
     atomicNumberIon != atomicNumberRefAr) {

     if(!referencePrepared) CreateReferenceParticles();

     if( useFe ) atomicNumber = atomicNumberRefFe;
     else atomicNumber = atomicNumberRefAr;     
  }

  return atomicNumber;
}

// ###########################################################################
