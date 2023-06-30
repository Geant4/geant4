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
// G4VTrackingManager
//
// Class description:
//
// Interface class for implementing a custom tracking manager that is
// specialized for stepping one or a small number of particle types.
//
// Original author: Jonas Hahnfeld, 2021

#ifndef G4VTrackingManager_hh
#define G4VTrackingManager_hh 1

class G4ParticleDefinition;
class G4Track;

////////////////////////
class G4VTrackingManager
////////////////////////
{
 public:
  virtual ~G4VTrackingManager() = default;

  // Messaged by the Particle definition whenever cross-section tables have
  // to be rebuilt (i.e. if new materials have been defined).
  virtual void BuildPhysicsTable(const G4ParticleDefinition&) {}

  // Messaged by the Particle definition whenever cross-section tables have
  // to be prepared for rebuild (i.e. if new materials have been defined).
  virtual void PreparePhysicsTable(const G4ParticleDefinition&) {}

  // Invoking this function, a G4Track given by the argument will be
  // handed over to this tracking manager. It may be tracked immediately
  // or processing may be deferred to a later time, at the latest when
  // calling FlushEvent().
  virtual void HandOverOneTrack(G4Track* aTrack) = 0;

  // Signal that all tracks in the current event have been finished and
  // this manager should process all tracks that may have been deferred.
  // When called via this method, the tracking manager may stack new
  // secondaries which will be tracked afterwards.
  virtual void FlushEvent() {}
};

#endif
