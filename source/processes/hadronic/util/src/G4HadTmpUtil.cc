//
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
#include "G4HadTmpUtil.hh"
#include <sstream>

G4int G4lrint(double ad)
  {
    return (ad>0) ? static_cast<int>(ad+.5) : static_cast<int>(ad-.5);
  }

G4int G4lint(double ad)
  {
    return (ad>0) ? static_cast<int>(ad) : static_cast<int>(ad-1.);
  }

G4int G4rint(double ad)
  {
    return (ad>0) ? static_cast<int>(ad+1) : static_cast<int>(ad);
  }

G4String G4inttostring(int ai)
{
  std::ostringstream ost;
  ost << ai;
  G4String result = ost.str();
  return result;
}
