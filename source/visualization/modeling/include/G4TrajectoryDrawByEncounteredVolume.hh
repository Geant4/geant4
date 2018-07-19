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
// $Id: G4TrajectoryDrawByEncounteredVolume.hh 66264 2012-12-14 10:17:44Z allison $
//
// Class Description:
// Trajectory model which colours a trajectory according to
// the encountered volume
// Class Description - End:
// John Allison February 2016, based on
// G4TrajectoryDrawByVolume.hh  Jane Tinslay March 2006

#ifndef G4TRAJECTORYDRAWBYENCOUNTEREDVOLUME
#define G4TRAJECTORYDRAWBYENCOUNTEREDVOLUME

#include "G4VTrajectoryModel.hh"
#include "G4Colour.hh"
#include "G4ModelColourMap.hh"
#include "G4String.hh"
#include <map>

class G4TrajectoryDrawByEncounteredVolume : public G4VTrajectoryModel {

public: // With description
 
  G4TrajectoryDrawByEncounteredVolume(const G4String& name = "Unspecified", G4VisTrajContext* context=0);
  
  virtual ~G4TrajectoryDrawByEncounteredVolume();

  // Draw method
  virtual void Draw(const G4VTrajectory& trajectory, 
		    const G4bool& visible = true) const;
  
  virtual void Print(std::ostream& ostr) const;
  // Print configuration

  void SetDefault(const G4String&);
  void SetDefault(const G4Colour&);

  void Set(const G4String& pvname, const G4String& colour);
  void Set(const G4String& pvname, const G4Colour& colour);
  // Configuration functions

private:

  // Data members
  G4ModelColourMap<G4String> fMap;
  G4Colour fDefault;

};

#endif
