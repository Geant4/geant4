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


#include "G4UppSimpleFieldtransport.hh"
#include "G4LorentzVector.hh"
#include "G4UppGlobal.hh"


void G4UppSimpleFieldtransport::propagate(G4UppTrackVector& allTracks, 
					  const G4double dTime) const
{
  if (uppdebug) G4cout << "(debug) fieldtransport dt=" << dTime*c_light/fermi << endl;
  if (uppdebug) {
    cout << "(debug)   from " << allTracks.getGlobalTime()*c_light/fermi;
    cout << " to " << (allTracks.getGlobalTime()+dTime)*c_light/fermi << endl;
  }
  for (G4int i=0; i<allTracks.size(); i++) {
    G4LorentzVector aMomentum = allTracks[i]->Get4Momentum();
    G4double E = aMomentum.e();
    G4LorentzVector aPosition = allTracks[i]->Get4Position();
    aPosition.setT(aPosition.t() + dTime);
    aPosition.setVect(aPosition.vect() + aMomentum.vect() * (dTime*c_light / E));
    allTracks[i]->Set4Momentum(aMomentum);
    allTracks[i]->Set4Position(aPosition);
  }
  allTracks.setGlobalTime(allTracks.getGlobalTime() + dTime);
} 
