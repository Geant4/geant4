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

#include "globals.hh"
#include "G4XAnnihilationChannel.hh"
#include "G4ConcreteMesonBaryonToResonance.hh"

G4BaryonWidth & G4ConcreteMesonBaryonToResonance::theBaryonWidth()
{
  static G4BaryonWidth theWidth;
  return theWidth;
}

G4BaryonPartialWidth & G4ConcreteMesonBaryonToResonance::theBaryonPartialWidth()
{
  static G4BaryonPartialWidth theWidth;
  return theWidth;
}

G4ParticleTypeConverter & G4ConcreteMesonBaryonToResonance::myConv() 
{
    static G4ParticleTypeConverter theConv;
    return theConv;
}

G4ConcreteMesonBaryonToResonance::G4ConcreteMesonBaryonToResonance(const G4ParticleDefinition* aPrimary,
								   const G4ParticleDefinition* bPrimary,
								   const G4ParticleDefinition* aSecondary,
								   const G4String& partWidthLabel)
  : thePrimary1(aPrimary), thePrimary2(bPrimary), theSecondary(aSecondary)
{
  crossSectionSource = new G4XAnnihilationChannel(aSecondary, 
						  theBaryonWidth(),
						  theBaryonPartialWidth(),
						  partWidthLabel);
}


G4ConcreteMesonBaryonToResonance::~G4ConcreteMesonBaryonToResonance()
{ 
  delete crossSectionSource;
}


G4bool G4ConcreteMesonBaryonToResonance::IsInCharge(const G4KineticTrack& trk1, 
						    const G4KineticTrack& trk2) const
{
  if (myConv().GetGenericType(trk1)==myConv().GetGenericType(thePrimary1) && 
      myConv().GetGenericType(trk2)==myConv().GetGenericType(thePrimary2)) return true;
  if (myConv().GetGenericType(trk1)==myConv().GetGenericType(thePrimary2) && 
      myConv().GetGenericType(trk2)==myConv().GetGenericType(thePrimary1)) return true;
  return false;
}

const G4ParticleDefinition* G4ConcreteMesonBaryonToResonance::GetOutgoingParticle(const G4KineticTrack& trk1, 
										  const G4KineticTrack& trk2) const
{
  G4int secondaryIso3 = trk1.GetDefinition()->GetPDGiIsospin3() + trk2.GetDefinition()->GetPDGiIsospin3();
  const G4ParticleDefinition* state;
  if ( (state = myConv().FindIso3State(myConv().GetGenericType(theSecondary), secondaryIso3)) == NULL) 
  {
    G4cerr << "for "<<static_cast<G4int>(myConv().GetGenericType(theSecondary))<<" "<<secondaryIso3<<G4endl;
    G4Exception("G4ConcreteMesonBaryonToResonance: Can't find iso3 state!");
  }
  return state;
}
