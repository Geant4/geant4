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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//

#ifndef G4DummyCascade_h
#define G4DummyCascade_h 1

#include "G4Fancy3DNucleus.hh"
#include "G4Nucleon.hh"
#include "G4Nucleus.hh"
#include "G4VIntraNuclearTransportModel.hh"
#include "G4KineticTrackVector.hh"
#include "G4FragmentVector.hh"
#include "G4ParticleChange.hh"
#include "G4DynamicParticleVector.hh"

class G4DummyCascade : public G4VIntraNuclearTransportModel 
{
public:
   G4DummyCascade(){}      
   G4DummyCascade(G4double anEnergy)
   {
     theEnergy = anEnergy;
   }      
   ~G4DummyCascade(){}

private:
   G4int operator==(G4DummyCascade& right) {return (this == &right);}
   G4int operator!=(G4DummyCascade& right) {return (this != &right);}
   
   G4double theEnergy;
      
public:
   G4VParticleChange* ApplyYourself(const G4Track& aTrack, G4Nucleus& theNucleus)
   { return new G4ParticleChange;}

   G4ReactionProductVector* Propagate(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus)
   {
     // decay the strong resonances
     for (G4int aRes=0; aRes < theSecondaries->size(); aRes++)
     {
       delete theSecondaries->operator[](aRes);
     }
     delete theSecondaries;
     return new G4ReactionProductVector;
   }


private:   
};

#endif // G4DummyCascade_h


