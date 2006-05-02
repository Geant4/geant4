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
// $Id: G4VTrajectoryModel.hh,v 1.5 2006-05-02 20:47:40 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
  virtual void Draw(const G4VTrajectory& model, const G4int& i_mode = 0, 
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

  G4String fName;
  G4bool fVerbose;
  G4VisTrajContext* fpContext;
  
};

#endif

