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
// $Id: G4TrajectoryDrawByOriginVolume.hh,v 1.1 2006-03-17 03:24:02 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Class Description:
// Trajectory model which colours a trajectory according to
// the origin volume
// Class Description - End:
// Jane Tinslay March 2006

#ifndef G4TRAJECTORYDRAWBYORIGINVOLUME
#define G4TRAJECTORYDRAWBYORIGINVOLUME

#include "G4VTrajectoryModel.hh"
#include "G4Colour.hh"
#include "G4ModelColourMap.hh"
#include "G4String.hh"
#include <map>

class G4TrajectoryDrawByOriginVolume : public G4VTrajectoryModel {

public: // With description
 
  G4TrajectoryDrawByOriginVolume(const G4String& name = "Unspecified");
  
  virtual ~G4TrajectoryDrawByOriginVolume();

  virtual void Draw(const G4VTrajectory&, G4int) const;
  // Draw the trajectory with optional i_mode parameter

  virtual void Print(std::ostream& ostr) const;
  // Print configuration

  void SetDefault(const G4String&);
  void SetDefault(const G4Colour&);

  void Set(const G4String& particle, const G4String& colour);
  void Set(const G4String& particle, const G4Colour& colour);
  // Configuration functions

private:

  // Data members
  G4ModelColourMap<G4String> fMap;
  G4Colour fDefault;

};

#endif
