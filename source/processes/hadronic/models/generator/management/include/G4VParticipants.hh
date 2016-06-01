// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VParticipants.hh,v 1.1 1998/08/22 08:55:46 hpw Exp $
// GEANT4 tag $Name: geant4-00 $
//

#ifndef G4VParticipants_h
#define G4VParticipants_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4VParticipants ----------------
//             by Gunter Folger, May 1998.
//      abstract class finding participants in a hadron Nucleus collision
//       in Parton String Models.
// ------------------------------------------------------------
#include "globals.hh"

class G4V3DNucleus;

#include "G4Fancy3DNucleus.hh"


class G4VParticipants 
{

  public:
      G4VParticipants();
      G4VParticipants(const G4VParticipants &right);
      virtual ~G4VParticipants();

      const G4VParticipants & operator=(const G4VParticipants &right);
      int operator==(const G4VParticipants &right) const;
      int operator!=(const G4VParticipants &right) const;

      void Init(G4double theZ, G4double theA);
      
      void SetNucleus(G4V3DNucleus * aNucleus);
      G4V3DNucleus * GetWoundedNucleus() const;


  protected:

  
      G4V3DNucleus *theNucleus;
      
  private:
  
};

// Class G4VParticipants 

inline G4V3DNucleus * G4VParticipants::GetWoundedNucleus() const
{
  return theNucleus;
}

inline void G4VParticipants::SetNucleus(G4V3DNucleus * aNucleus)
{
  theNucleus = aNucleus;
}

inline void G4VParticipants::Init(G4double theA, G4double theZ)
{
	if ( theNucleus == NULL ) theNucleus = new G4Fancy3DNucleus();
	theNucleus->Init(theA, theZ);
}


#endif


