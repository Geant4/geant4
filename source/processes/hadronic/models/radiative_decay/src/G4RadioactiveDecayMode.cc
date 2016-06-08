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

#include "G4RadioactiveDecayMode.hh"

G4std::istream &operator >> (G4std::istream &s, G4RadioactiveDecayMode &q)
{
  G4String a;
  s >> a;
  if (a == "IT")
    {q = IT;}
  else if (a == "BetaMinus")
    {q = BetaMinus;}
  else if (a == "BetaPlus")
    {q = BetaPlus;}
  else if (a == "KshellEC")
    {q = KshellEC;}
  else if (a == "LshellEC")
    {q = LshellEC;}
  else if (a == "MshellEC")
    {q = MshellEC;}
  else if (a == "Alpha")
    {q = Alpha;}
  else
    {q = ERROR;}
  return s;
}
////////////////////////////////////////////////////////////////////////////////

