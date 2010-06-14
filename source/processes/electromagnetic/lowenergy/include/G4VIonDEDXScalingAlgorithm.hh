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
// Class:                G4VIonDEDXScalingAlgorithm
//
// Author:               Anton Lechner (Anton.Lechner@cern.ch)
//
// First implementation: 11. 03. 2009
//
// Modifications: 12. 11. 2009 - Added second argument (material) to energy
//                               scaling function (AL)
//
// Class description:
//    Base class for dE/dx scaling algorithms, used by G4IonDEDXHandler
//
// Comments:
//
// =========================================================================== 


#ifndef G4VIONDEDXSCALINGALGORITHM_HH
#define G4VIONDEDXSCALINGALGORITHM_HH

#include "globals.hh"
#include "G4ParticleDefinition.hh"

class G4Material;


class G4VIonDEDXScalingAlgorithm {

 public:
   G4VIonDEDXScalingAlgorithm();
   virtual ~G4VIonDEDXScalingAlgorithm();

   // Function for scaling the kinetic energy (no scaling by default).
   // Returns scaling factor for a given ion.
   virtual G4double ScalingFactorEnergy(
             const G4ParticleDefinition*,     // Projectile (ion) 
             const G4Material*)               // Target material
                 { return 1.0; }

   // Function for scaling the dE/dx value (no scaling by default).
   // Returns scaling factor for a given ion-material couple and
   // a given kinetic energy.
   virtual G4double ScalingFactorDEDX(
             const G4ParticleDefinition*,     // Projectile (ion) 
             const G4Material*,               // Target material
             G4double)                        // Kinetic energy of projectile 
                 { return 1.0; }

   // Function for defining a base particle for dE/dx calculation.
   // (no base particle by default). Returns atomic number of base
   // particle.
   virtual G4int AtomicNumberBaseIon(
             G4int atomicNumberIon,           // Atomic number of ion 
             const G4Material*)               // Target material
                 { return atomicNumberIon; }
				
};

#endif
