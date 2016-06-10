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
// $Id: G4TrajectoryChargeFilter.cc 66373 2012-12-18 09:41:34Z gcosmo $
//
// Filter trajectories according to charge. Only registered 
// charges will pass the filter.
//
// Jane Tinslay May 2006
//
#include "G4TrajectoryChargeFilter.hh"
#include <sstream>

G4TrajectoryChargeFilter::G4TrajectoryChargeFilter(const G4String& name)
  :G4SmartFilter<G4VTrajectory>(name)
{}

G4TrajectoryChargeFilter::~G4TrajectoryChargeFilter() {}

bool
G4TrajectoryChargeFilter::Evaluate(const G4VTrajectory& traj) const
{
  G4double charge = traj.GetCharge();

  if (GetVerbose()) G4cout<<"G4TrajectoryChargeFilter processing trajectory with charge: "<<charge<<G4endl;

  MyCharge myCharge;

  if(charge>0.)      myCharge = Positive; 
  else if(charge<0.) myCharge = Negative; 
  else               myCharge = Neutral; 

  std::vector<MyCharge>::const_iterator iter = std::find(fCharges.begin(), fCharges.end(), myCharge);

  // Fail if charge not registered 
  if (iter == fCharges.end()) return false;

  return true;
}

void
G4TrajectoryChargeFilter::Add(const G4String& charge) 
{
  MyCharge myCharge;
  
  if (!ConvertToCharge(charge, myCharge)) {
    G4ExceptionDescription ed;
    ed << "Invalid charge "<<charge;
    G4Exception   
      ("G4TrajectoryChargeFilter::Add(const G4String& charge)",
       "modeling0115", JustWarning, ed);
    return;
  }
  
  return Add(myCharge);
}

void
G4TrajectoryChargeFilter::Add(const MyCharge& charge)
{
  fCharges.push_back(charge);
}

void
G4TrajectoryChargeFilter::Print(std::ostream& ostr) const
{
  ostr<<"Charges registered: "<<G4endl;
  std::vector<MyCharge>::const_iterator iter = fCharges.begin();
  
  while (iter != fCharges.end()) {
    ostr<<*iter<<G4endl;    
    iter++;
  }
}

void 
G4TrajectoryChargeFilter::Clear()
{
  // Clear registered charge vector
  fCharges.clear();
}

G4bool
G4TrajectoryChargeFilter::ConvertToCharge(const G4String& string, MyCharge& myCharge)
{
  bool result(true);
 
  G4int charge;
  std::istringstream is(string.c_str());
  is >> charge;

  switch (charge) {
  case 1:
    myCharge = G4TrajectoryChargeFilter::Positive;
    break;
  case 0:
    myCharge = G4TrajectoryChargeFilter::Neutral;  
    break;
  case -1:
    myCharge = G4TrajectoryChargeFilter::Negative;   
    break;
  default:
    result = false;
  }
  
  return result;
}
