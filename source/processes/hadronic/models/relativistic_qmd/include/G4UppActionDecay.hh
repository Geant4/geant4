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

#ifndef G4UPPACTIONDECAY_H
#define G4UPPACTIONDECAY_H


#include "G4VUppAction.hh"


class G4UppActionDecay : public G4VUppAction
{
public:

  G4UppActionDecay(const G4double decayTime, 
		   const G4UppTrack& decayingParticle);

  G4bool isValid() const;

  G4UppTrackChange* perform(const G4UppTrackVector& allTracks) const;
  
  void G4UppActionDecay::dump() const;

private:

  G4UppTrack* particlePtr;

};


#endif // G4UPPACTIONDECAY_H
