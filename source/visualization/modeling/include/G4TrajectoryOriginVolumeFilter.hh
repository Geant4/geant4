// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4TrajectoryOriginVolumeFilter.hh,v 1.1 2006-05-30 18:44:36 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
  virtual bool Evaluate(const G4VTrajectory&);

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
