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
// $Id: G4TrajectoryOriginVolumeFilter.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
// Filter trajectories according to volume name. Only registered 
// volumes will pass the filter.
//
// Jane Tinslay, May 2006
//
#ifndef G4TRAJECTORYORIGINVOLUMEFILTER_HH
#define G4TRAJECTORYORIGINVOLUMEFILTER_HH

#include "G4SmartFilter.hh"
#include "G4VTrajectory.hh"
#include <vector>

class G4TrajectoryOriginVolumeFilter : public G4SmartFilter<G4VTrajectory> {

public: // With description
 
  // Construct with filter name
  G4TrajectoryOriginVolumeFilter(const G4String& name = "Unspecified");
  
  virtual ~G4TrajectoryOriginVolumeFilter();

  // Evaluate this trajectory
  virtual bool Evaluate(const G4VTrajectory&) const;

  // Print configuration
  virtual void Print(std::ostream& ostr) const;

  // Clear filter
  virtual void Clear();

  // Configuration function
  void Add(const G4String& volume);

private:

  // Data member
  std::vector<G4String> fVolumes;

};

#endif
