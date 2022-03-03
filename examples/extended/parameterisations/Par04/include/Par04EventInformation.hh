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
#ifndef PAR04EVENTINFORMATION_HH
#define PAR04EVENTINFORMATION_HH

#include <G4Types.hh>                    // for G4bool
#include "G4ThreeVector.hh"              // for G4ThreeVector
#include "G4VUserEventInformation.hh"    // for G4VUserEventInformation

/**
 * @brief Event information
 *
 * Additional data associated to the primary particle.
 * Contains information on the direction of the primary particle used to determine
 * how energy is scored (mesh around the particle direction).
 *
 */

class Par04EventInformation : public G4VUserEventInformation
{
 public:
  Par04EventInformation();
  virtual ~Par04EventInformation();

  /// Set particle direction
  inline void SetDirection(const G4ThreeVector& aDirection) { fDirection = aDirection; };
  /// Get particle direction
  inline G4ThreeVector GetDirection() const { return fDirection; };
  /// Set particle position
  inline void SetPosition(const G4ThreeVector& aPosition) { fPosition = aPosition; };
  /// Get particle position
  inline G4ThreeVector GetPosition() const { return fPosition; };
  /// Set flag
  inline void SetFlag(G4bool aFlag) { fIfSet = aFlag; };
  /// Get flag
  inline G4bool GetFlag() const { return fIfSet; };
  /// Print
  void Print() const final;

 private:
  /// Particle direction. By default equal to the default particle gun direction.
  G4ThreeVector fDirection = { 0, 1, 0 };
  /// Particle position. By default equal to the default inner radius.
  G4ThreeVector fPosition = { 0, 800, 0 };
  /// Flag
  G4bool fIfSet = false;
};

#endif /* PAR04EVENTINFORMATION_HH */
