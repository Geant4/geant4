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
// $Id: G4TrajectoryDrawByParticleID.hh,v 1.3 2005/11/23 05:19:23 tinslay Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// Jane Tinslay, John Allison, Joseph Perl November 2005
//
// Class Description:
// Trajectory model which colours a trajectory according to
// particle type.
// Class Description - End:

#ifndef G4TRAJECTORYDRAWBYPARTICLEID
#define G4TRAJECTORYDRAWBYPARTICLEID

#include "G4VTrajectoryModel.hh"
#include "G4Colour.hh"
#include "G4String.hh"
#include <map>

class G4TrajectoryDrawByParticleID : public G4VTrajectoryModel {

public: // With description
 
  G4TrajectoryDrawByParticleID(const G4String& name = "Unspecified");
  
  virtual ~G4TrajectoryDrawByParticleID();

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
  std::map<G4String, G4Colour> fMap;
  G4Colour fDefault;

};

#endif
