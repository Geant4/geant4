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
//
// $Id: G4FTFParticipants.hh,v 1.2 2003/10/08 13:48:47 hpw Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//

#ifndef G4FTFParticipants_h
#define G4FTFParticipants_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4FTFParticipants----------------
//             by Gunter Folger, June 1998.
//       class finding colliding particles in FTFPartonStringModel
// ------------------------------------------------------------

#include "G4VParticipants.hh"
#include <vector>
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

      std::vector<G4InteractionContent *> theInteractions;
  
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
	return ++currentInteraction < static_cast<G4int>(theInteractions.size());
}


inline
const G4InteractionContent & G4FTFParticipants::GetInteraction() const
{
	return *theInteractions[currentInteraction];
}

#endif


