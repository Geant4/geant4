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

#include "G4UppTrack.hh"


void G4UppTrack::dump() const
{
  cout << "Name: " << GetDefinition()->GetParticleName() << endl;
  cout << "  at " << Get4Position() << endl;
  cout << "  p= " << Get4Momentum() << endl;
  cout << "  nColl=" << numberOfCollisions << G4endl;
}


G4bool G4UppTrack::isLastInteractionPartner(const G4UppTrack* PartPtr) const
{
  for (G4int i=0; i<lastInteractionPartners.size(); i++) {
    if (PartPtr==lastInteractionPartners[i]) return true;
  }
  return false;
}
