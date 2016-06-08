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

#ifndef G4GeneratorPrecompoundInterface_h
#define G4GeneratorPrecompoundInterface_h 1

#include "G4Fancy3DNucleus.hh"
#include "G4Nucleon.hh"
#include "G4Nucleus.hh"
#include "G4VIntraNuclearTransportModel.hh"
#include "G4KineticTrackVector.hh"
#include "G4FragmentVector.hh"
#include "G4ParticleChange.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"

// Class Description
// Trivial implementation of an intra-nuclear transport. It pworvides coupling
// of high energy generators with pre equilibrium decay models. 
// To be used in your physics list in case you need this physics.
// Class Description - End

class G4GeneratorPrecompoundInterface : public G4VIntraNuclearTransportModel 
{
public:
   G4GeneratorPrecompoundInterface(){}      
   ~G4GeneratorPrecompoundInterface(){}

private:
   G4int operator==(G4GeneratorPrecompoundInterface& right) {return (this == &right);}
   G4int operator!=(G4GeneratorPrecompoundInterface& right) {return (this != &right);}
      
public:
   G4VParticleChange* ApplyYourself(const G4Track& aTrack, G4Nucleus& theNucleus);
   G4ReactionProductVector* Propagate(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus);


private:   
};

#endif // G4GeneratorPrecompoundInterface_h


