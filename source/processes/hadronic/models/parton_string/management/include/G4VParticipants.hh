//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4VParticipants.hh,v 1.7 2010-09-08 16:58:04 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

      void Init(G4int theZ, G4int theA);
      
      void SetNucleus(G4V3DNucleus * aNucleus);
      G4V3DNucleus * GetWoundedNucleus() const;


//  protected:   // Uzhi 26 July 09

  
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

inline void G4VParticipants::Init(G4int theA, G4int theZ)
{
	if ( theNucleus == NULL ) theNucleus = new G4Fancy3DNucleus();
	theNucleus->Init(theA, theZ);
        theNucleus->SortNucleonsIncZ();
}


#endif


