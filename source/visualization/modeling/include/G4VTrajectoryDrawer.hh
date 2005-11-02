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
// $Id: G4VTrajectoryDrawer.hh,v 1.1 2005-11-02 16:56:05 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Jane Tinslay, John Allison, Joseph Perl October 2005
//
// Class Description:
// Abstract base class for trajectory drawers. Drawers are responsible
// for drawing trajectories according to a particular style.
// Class Description - End:

#ifndef G4VTRAJECTORYDRAWER_HH
#define G4VTRAJECTORYDRAWER_HH

#include "globals.hh"
#include "G4String.hh"

class G4VTrajectory;

class G4VTrajectoryDrawer {

public: // With description

  G4VTrajectoryDrawer(const G4String& name = "Unspecified");

  virtual ~G4VTrajectoryDrawer();

  const G4String& GetName() const;
   
  virtual void Draw(const G4VTrajectory&, G4int i_mode = 0) = 0;
  // Draw the trajectory with optional i_mode parameter

  virtual void Print() const = 0;
  // Print drawer configuration

private:
  
  G4String fName; // Drawer name

};

#endif
