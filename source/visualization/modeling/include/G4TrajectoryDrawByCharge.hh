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
// $Id: G4TrajectoryDrawByCharge.hh,v 1.1 2005-11-02 00:41:13 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Jane Tinslay October 2005
//
// Class Description:
// Trajectory drawer which colours trajectories according to particle charge.
// Guts taken from G4VTrajectory::DrawTrajectory method.
// Class Description - End:

#ifndef G4TRAJECTORYDRAWBYCHARGE_HH
#define G4TRAJECTORYDRAWBYCHARGE_HH

#include "G4VTrajectoryDrawer.hh"
#include "G4Colour.hh"
#include "G4String.hh"

class G4TrajectoryDrawByCharge : public G4VTrajectoryDrawer {

public: // With description
 
  G4TrajectoryDrawByCharge(const G4String& name = "Unspecified");
  // Use default charge colour scheme

  G4TrajectoryDrawByCharge(const G4Colour& positive,
			   const G4Colour& negative,
			   const G4Colour& neutral);

  virtual ~G4TrajectoryDrawByCharge();

  virtual void Draw(const G4VTrajectory&, G4int);
  virtual void Print() const;

private:
  
  //Data members
  G4Colour fPositive;
  G4Colour fNegative;
  G4Colour fNeutral;
  
};

#endif
