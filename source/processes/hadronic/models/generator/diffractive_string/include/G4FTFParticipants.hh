// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FTFParticipants.hh,v 1.2 1998/08/26 17:12:00 gcosmo Exp $
// GEANT4 tag $Name: geant4-00 $
//

#ifndef G4FTFParticipants_h
#define G4FTFParticipants_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4FTFParticipants----------------
//             by Gunter Folger, June 1998.
//       class finding colliding particles in FTFPartonStringModel
// ------------------------------------------------------------

#include "G4VParticipants.hh"
#include <rw/tpordvec.h>
#include "G4Nucleon.hh"
#include "G4V3DNucleus.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4ReactionProduct.hh"
#include "G4InteractionContent.hh"

class G4FTFParticipants : public G4VParticipants
{

  public:
      G4FTFParticipants();
      G4FTFParticipants(const G4FTFParticipants &right);
      const G4FTFParticipants & operator=(const G4FTFParticipants &right);
      ~G4FTFParticipants();

      int operator==(const G4FTFParticipants &right) const;
      int operator!=(const G4FTFParticipants &right) const;

      void BuildInteractions(const G4ReactionProduct  &thePrimary);
      G4bool Next();
      const G4InteractionContent & GetInteraction() const;

      void StartLoop();
      
  private:

      RWTPtrOrderedVector<G4InteractionContent> theInteractions;
  
      G4int currentInteraction;

};


inline
void G4FTFParticipants::StartLoop()
{
	currentInteraction=-1;
}

inline
G4bool G4FTFParticipants::Next()
{
	return ++currentInteraction < theInteractions.entries();
}


inline
const G4InteractionContent & G4FTFParticipants::GetInteraction() const
{
	return *theInteractions[currentInteraction];
}

#endif


