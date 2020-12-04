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

#ifndef G4FASTHIT_HH
#define G4FASTHIT_HH

#include "G4ThreeVector.hh"

/**
 * @brief Minimal hit created in the fast simulation
 *
 * Minimal hit containing energy and position, for use in the fast simulation
 * classes.
 * Hits of G4FastHit type can be created in user implementation of fast
 * simulation model and then deposited in the detector using G4FastSimHitMaker 
 * helper class. The helper will locate the sensitive volume and check if it
 * inherits from both base classes:
 * - G4VSensitiveDetector: for processing of detailed/non-fast simulation hits;
 * - G4VFastSimSensitiveDetector: for processing of fast sim (G4FastSim) hits;
 * An extended example extended/parameterisations/Par03 demonstrates how to use
 * G4FastHit to create multiple deposits from the fast simulation model.
 */

class G4FastHit
{
 public:
  G4FastHit();
  G4FastHit(const G4ThreeVector& aPosition, G4double aEnergy);
  G4FastHit(const G4ThreeVector& aPosition, G4double aEnergy, G4bool aDebug);
  virtual ~G4FastHit(){};

  /// Set energy
  inline void SetEnergy(const G4double& aEnergy) { fEnergy = aEnergy; }
  /// Get energy
  inline G4double GetEnergy() const { return fEnergy; }
  /// Set position
  inline void SetPosition(const G4ThreeVector& aPosition)
  {
    fPosition = aPosition;
  }
  /// Get position
  inline G4ThreeVector GetPosition() const { return fPosition; }
 private:
  /// energy
  G4double fEnergy = 0;
  /// position
  G4ThreeVector fPosition = G4ThreeVector();
};

#endif /* G4FASTHIT_HH */
