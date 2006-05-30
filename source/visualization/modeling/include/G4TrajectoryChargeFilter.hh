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
// $Id: G4TrajectoryChargeFilter.hh,v 1.1 2006-05-30 18:44:36 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
  virtual bool Evaluate(const G4VTrajectory&);

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
