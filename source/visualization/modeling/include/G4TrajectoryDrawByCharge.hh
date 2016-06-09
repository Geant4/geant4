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
// $Id: G4TrajectoryDrawByCharge.hh,v 1.3 2005/11/23 05:19:23 tinslay Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// Jane Tinslay, John Allison, Joseph Perl November 2005
//
// Class Description:
// Trajectory model which colours a trajectory according to  
// charge. Guts taken from G4VTrajectory::DrawTrajectory method.
// Class Description - End:

#ifndef G4TRAJECTORYDRAWBYCHARGE_HH
#define G4TRAJECTORYDRAWBYCHARGE_HH

#include "G4Colour.hh"
#include "G4VTrajectoryModel.hh"
#include <map>

class G4TrajectoryDrawByCharge : public G4VTrajectoryModel {

public: // With description
 
  enum Charge {Negative=-1, Neutral=0, Positive=1}; 

  G4TrajectoryDrawByCharge(const G4String& name = "Unspecified");

  G4TrajectoryDrawByCharge(const G4String& name,
			   const G4Colour& positive,
			   const G4Colour& negative,
			   const G4Colour& neutral);
  
  virtual ~G4TrajectoryDrawByCharge();

  virtual void Draw(const G4VTrajectory& trajectory, G4int i_mode = 0) const;
  // Draw the trajectory with optional i_mode parameter

  virtual void Print(std::ostream& ostr) const;
  // Print configuration

  void Set(Charge charge, const G4String& colour);
  void Set(Charge charge, const G4Colour& colour);
  // Configuration functions 

private:
  
  // Data member
  std::map<Charge, G4Colour> fMap;
  
};

#endif
