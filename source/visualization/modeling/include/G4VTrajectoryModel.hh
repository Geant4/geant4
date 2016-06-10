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
// $Id: G4VTrajectoryModel.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
// Jane Tinslay, John Allison, Joseph Perl October 2005
//
// Class Description:
// Abstract base class for trajectory drawing models. Trajectory drawing
// models are responsible for drawing individual trajectories according 
// to a particular style.
// Class Description - End:

#ifndef G4VTRAJECTORYMODEL_HH
#define G4VTRAJECTORYMODEL_HH

#include "G4String.hh"
#include "G4VTrajectory.hh"

class G4VisTrajContext;

class G4VTrajectoryModel {

public:

  // Construct with context object
  G4VTrajectoryModel(const G4String& name, G4VisTrajContext* fpContext=0);

  // Destructor
  virtual ~G4VTrajectoryModel();
  
  // Draw method
  virtual void Draw(const G4VTrajectory& trajectory, 
		    const G4bool& visible = true) const = 0;
  
  // Print configuration
  virtual void Print(std::ostream& ostr) const = 0;
  
  // Accessors
  G4String Name() const ;
  const G4VisTrajContext& GetContext() const;
  
  // Set verbosity
  void SetVerbose(const G4bool&);
  G4bool GetVerbose() const;

private:

  // Private copy constructor and assigment operator - copying and
  // assignment not allowed.  Keeps Coverity happy.
  G4VTrajectoryModel (const G4VTrajectoryModel&);
  G4VTrajectoryModel& operator = (const G4VTrajectoryModel&);

  G4String fName;
  G4bool fVerbose;
  G4VisTrajContext* fpContext;
  
};

#endif

