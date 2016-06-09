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
// $Id: G4VParticipants.hh,v 1.1 2003/10/07 11:26:00 hpw Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//

#ifndef G4VParticipants_h
#define G4VParticipants_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
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


