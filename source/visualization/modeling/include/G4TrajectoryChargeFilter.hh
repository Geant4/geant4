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
// $Id: G4TrajectoryChargeFilter.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
// Filter trajectories according to charge. Only registered 
// charges will pass the filter.
//
// Jane Tinslay, May 2006
//
#ifndef G4TRAJECTORYCHARGEFILTER_HH
#define G4TRAJECTORYCHARGEFILTER_HH

#include "G4SmartFilter.hh"
#include "G4VTrajectory.hh"
#include <vector>

class G4TrajectoryChargeFilter : public G4SmartFilter<G4VTrajectory> {

public: // With description
 
  // Construct with filter name
  G4TrajectoryChargeFilter(const G4String& name = "Unspecified");
  
  virtual ~G4TrajectoryChargeFilter();

  // Evaluate this trajectory
  virtual bool Evaluate(const G4VTrajectory&) const;

  // Print configuration
  virtual void Print(std::ostream& ostr) const;

  // Clear filter
  virtual void Clear();

  // Configuration function
  void Add(const G4String& particle);

private:

  enum MyCharge {Negative=-1, Neutral=0, Positive=1}; 

  G4bool ConvertToCharge(const G4String&, MyCharge&);

  void Add(const MyCharge& chgear);  

  // Data member
  std::vector<MyCharge> fCharges;

};

#endif
