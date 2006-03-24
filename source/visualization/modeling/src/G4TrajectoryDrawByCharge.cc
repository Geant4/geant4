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
// $Id: G4TrajectoryDrawByCharge.cc,v 1.5 2006-03-24 20:22:43 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Jane Tinslay, John Allison, Joseph Perl November 2005
#include "G4TrajectoryDrawByCharge.hh"
#include "G4TrajectoryDrawerUtils.hh"
#include "G4VTrajectory.hh"
#include <sstream>

G4TrajectoryDrawByCharge::G4TrajectoryDrawByCharge(const G4String& name)
  :G4VTrajectoryModel(name)
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
G4TrajectoryDrawByCharge::Draw(const G4VTrajectory& traj, const G4int& i_mode, const G4bool& visible) const
{
  G4Colour colour;

  const G4double charge = traj.GetCharge();

  if(charge>0.)      fMap.GetColour(Positive, colour); 
  else if(charge<0.) fMap.GetColour(Negative, colour); 
  else               fMap.GetColour(Neutral, colour); 

  G4TrajectoryDrawerUtils::DrawLineAndPoints(traj, i_mode, colour, visible);
}

void
G4TrajectoryDrawByCharge::Print(std::ostream& ostr) const
{
  ostr<<"G4TrajectoryDrawByCharge model "<< Name() <<" colour scheme: "<<std::endl;
  fMap.Print(ostr);
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
    std::ostringstream o;
    o << "Invalid charge "<<charge;
    G4Exception   
      ("G4TrajectoryDrawByCharge::Set(const G4int& charge, const G4String& colour)", "InvalidCharge", JustWarning, o.str().c_str());
  }

  return Set(myCharge, colour);
}

void
G4TrajectoryDrawByCharge::Set(const G4String& charge, const G4Colour& colour)
{  
  Charge myCharge;
  
  if (!ConvertToCharge(charge, myCharge)) {
    std::ostringstream o;
    o << "Invalid charge "<<charge;
    G4Exception   
      ("G4TrajectoryDrawByCharge::Set(const G4int& charge, const G4Colour& colour)", "InvalidCharge", JustWarning, o.str().c_str());
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
