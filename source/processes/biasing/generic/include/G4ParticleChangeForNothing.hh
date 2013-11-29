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
// $Id: $
//
// ---------------------------------------------------------------
//
// G4ParticleChangeForNothing
//
// Class Description:
//     A G4VParticleChange used in case of no interaction.
//
// ---------------------------------------------------------------
//   Initial version                         Nov. 2013 M. Verderi

#ifndef G4ParticleChangeForNothing_hh
#define G4ParticleChangeForNothing_hh 1

#include "G4VParticleChange.hh"

class G4ParticleChangeForNothing : public G4VParticleChange {
public:
  G4ParticleChangeForNothing() : G4VParticleChange() {}
  ~G4ParticleChangeForNothing() {}

public:
  // -- from base class G4VParticleChange:
  virtual void Initialize(const G4Track &track)
  {
    theStatusChange = track.GetTrackStatus();
    theNumberOfSecondaries = 0;
  }
  virtual G4Step* UpdateStepForAtRest   (G4Step* step) {return step;}
  virtual G4Step* UpdateStepForAlongStep(G4Step* step) {return step;}
  virtual G4Step* UpdateStepForPostStep (G4Step* step) {return step;}
};

#endif
