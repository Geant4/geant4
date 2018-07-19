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
/// \file processes/phonon/include/G4VPhononProcess.hh
/// \brief Definition of the G4VPhononProcess base class
//
// $Id: G4VPhononProcess.hh 75725 2013-11-05 16:52:30Z mkelsey $
//
#ifndef G4VPhononProcess_h
#define G4VPhononProcess_h 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4ThreeVector.hh"

class G4PhononTrackMap;
class G4LatticePhysical;


class G4VPhononProcess : public G4VDiscreteProcess {
public:
  G4VPhononProcess(const G4String& processName);
  virtual ~G4VPhononProcess();

  virtual G4bool IsApplicable(const G4ParticleDefinition& aPD);

  // Initialize wave vectors for currently active track(s)
  // NOTE:  These functions must call back to base class implementations!
  virtual void StartTracking(G4Track* track);
  virtual void EndTracking();

protected:
  // For convenience, map phonon type to polarization code
  virtual G4int GetPolarization(const G4Track& track) const;
  virtual G4int GetPolarization(const G4Track* track) const {
    return GetPolarization(*track);
  }

  // For convenience, generate random polarization from density of states
  // Values passed may be zero to suppress particular states
  virtual G4int ChoosePolarization(G4double Ldos, G4double STdos,
				   G4double FTdos) const;

  // Construct new track with correct momentum, position, etc.
  virtual G4Track* CreateSecondary(G4int polarization, const G4ThreeVector& K,
				   G4double energy) const;

protected:
  G4PhononTrackMap* trackKmap;		// For convenient access by processes
  const G4LatticePhysical* theLattice;

private:
  const G4Track* currentTrack;		// For use by Start/EndTracking

  // hide assignment operators as private 
  G4VPhononProcess(G4VPhononProcess&);
  G4VPhononProcess& operator=(const G4VPhononProcess& right);
};

#endif	/* G4VPhononProcess_h */
