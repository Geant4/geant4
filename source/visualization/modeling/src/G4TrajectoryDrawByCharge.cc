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
// $Id: G4TrajectoryDrawByCharge.cc 66373 2012-12-18 09:41:34Z gcosmo $
//
// Jane Tinslay, John Allison, Joseph Perl November 2005
#include "G4TrajectoryDrawByCharge.hh"
#include "G4TrajectoryDrawerUtils.hh"
#include "G4VisTrajContext.hh"
#include "G4VTrajectory.hh"
#include <sstream>

G4TrajectoryDrawByCharge::G4TrajectoryDrawByCharge(const G4String& name, G4VisTrajContext* context)
  :G4VTrajectoryModel(name, context)
{
  // Default configuration
  fMap[Positive] = G4Colour::Blue();
  fMap[Negative] = G4Colour::Red();
  fMap[Neutral] = G4Colour::Green();
}

G4TrajectoryDrawByCharge::G4TrajectoryDrawByCharge(const G4String& name,
						   const G4Colour& positive,
						   const G4Colour& negative,
						   const G4Colour& neutral)
  :G4VTrajectoryModel(name)
{
  fMap[Positive] = positive;
  fMap[Negative] = negative;
  fMap[Neutral] = neutral;
}

G4TrajectoryDrawByCharge::~G4TrajectoryDrawByCharge() {}

void
G4TrajectoryDrawByCharge::Draw(const G4VTrajectory& traj, const G4bool& visible) const
{
  G4Colour colour;

  const G4double charge = traj.GetCharge();

  if(charge>0.)      fMap.GetColour(Positive, colour); 
  else if(charge<0.) fMap.GetColour(Negative, colour); 
  else               fMap.GetColour(Neutral, colour); 

  G4VisTrajContext myContext(GetContext());
  
  myContext.SetLineColour(colour);
  myContext.SetVisible(visible);
  
  if (GetVerbose()) {
    G4cout<<"G4TrajectoryDrawByCharge drawer named "<<Name();
    G4cout<<", drawing trajectory with charge, "<<charge<<G4endl;
    G4cout<<", with configuration:"<<G4endl;
    myContext.Print(G4cout);
  }

  G4TrajectoryDrawerUtils::DrawLineAndPoints(traj, myContext);
}

void
G4TrajectoryDrawByCharge::Print(std::ostream& ostr) const
{
  ostr<<"G4TrajectoryDrawByCharge model "<< Name() <<" colour scheme: "<<std::endl;
  fMap.Print(ostr);
  ostr<<"Default configuration:"<<G4endl;
  GetContext().Print(G4cout);
}

void
G4TrajectoryDrawByCharge::Set(const Charge& charge, const G4String& colour)
{
  fMap.Set(charge, colour);
}

void
G4TrajectoryDrawByCharge::Set(const Charge& charge, const G4Colour& colour)
{
  fMap[charge] = colour;
}

void
G4TrajectoryDrawByCharge::Set(const G4String& charge, const G4String& colour)
{  
  Charge myCharge;
  
  if (!ConvertToCharge(charge, myCharge)) {
    G4ExceptionDescription ed;
    ed << "Invalid charge "<<charge;
    G4Exception   
      ("G4TrajectoryDrawByCharge::Set(const G4int& charge, const G4String& colour)", "modeling0121", JustWarning, ed);
    return;
  }

  return Set(myCharge, colour);
}

void
G4TrajectoryDrawByCharge::Set(const G4String& charge, const G4Colour& colour)
{  
  Charge myCharge;
  
  if (!ConvertToCharge(charge, myCharge)) {
    G4ExceptionDescription ed;
    ed << "Invalid charge "<<charge;
    G4Exception   
      ("G4TrajectoryDrawByCharge::Set(const G4int& charge, const G4Colour& colour)", "modeling0122", JustWarning, ed);
  }

  return Set(myCharge, colour);
}

G4bool
G4TrajectoryDrawByCharge::ConvertToCharge(const G4String& string, Charge& myCharge)
{
  bool result(true);
 
  G4int charge;
  std::istringstream is(string.c_str());
  is >> charge;

  switch (charge) {
  case 1:
    myCharge = G4TrajectoryDrawByCharge::Positive;
    break;
  case 0:
    myCharge = G4TrajectoryDrawByCharge::Neutral;  
    break;
  case -1:
    myCharge = G4TrajectoryDrawByCharge::Negative;   
    break;
  default:
    result = false;
  }
  
  return result;
}
