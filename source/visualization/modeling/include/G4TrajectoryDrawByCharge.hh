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
// $Id: G4TrajectoryDrawByCharge.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
// Jane Tinslay, John Allison, Joseph Perl November 2005
//
// Class Description:
// Trajectory model which colours a trajectory according to  
// charge. 
// Class Description - End:

#ifndef G4TRAJECTORYDRAWBYCHARGE_HH
#define G4TRAJECTORYDRAWBYCHARGE_HH

#include "G4Colour.hh"
#include "G4ModelColourMap.hh"
#include "G4VTrajectoryModel.hh"
#include <map>

class G4TrajectoryDrawByCharge : public G4VTrajectoryModel {

public: // With description
 
  enum Charge {Negative=-1, Neutral=0, Positive=1}; 

  G4TrajectoryDrawByCharge(const G4String& name = "Unspecified", G4VisTrajContext* context=0);

  G4TrajectoryDrawByCharge(const G4String& name,
			   const G4Colour& positive,
			   const G4Colour& negative,
			   const G4Colour& neutral);
  
  virtual ~G4TrajectoryDrawByCharge();

  // Draw method
  virtual void Draw(const G4VTrajectory& trajectory, 
		    const G4bool& visible = true) const;

  // Print configuration
  virtual void Print(std::ostream& ostr) const;

  void Set(const Charge& charge, const G4Colour& colour);
  void Set(const Charge& charge, const G4String& colour);

  void Set(const G4String& charge, const G4Colour& colour);
  void Set(const G4String& charge, const G4String& colour);
  // Configuration functions 

private:
  
  G4bool ConvertToCharge(const G4String&, Charge&);

  // Data member
  G4ModelColourMap<Charge> fMap;
  
};

#endif
