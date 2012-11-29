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
// $Id$
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
#include "G4FTFParameters.hh"
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
      const G4FTFParticipants & operator=(const G4FTFParticipants &right);
      ~G4FTFParticipants();

      int operator==(const G4FTFParticipants &right) const;
      int operator!=(const G4FTFParticipants &right) const;
//---------------------------------------------------
      void InitProjectileNucleus(G4int theZ, G4int theA);
      void SetProjectileNucleus(G4V3DNucleus * aNucleus);
      G4V3DNucleus * GetProjectileNucleus();
//---------------------------------------------------
      void GetList(const G4ReactionProduct  &thePrimary, 
                         G4FTFParameters    *theParameters);

      void StartLoop();
      G4bool Next();
//Vova      const G4InteractionContent & GetInteraction() const;
      G4InteractionContent & GetInteraction();
      
      std::vector<G4InteractionContent *> theInteractions;
      G4V3DNucleus *theProjectileNucleus;
  private:

      //A.R. 25-Jul-2012 Coverity fix : copy constructor becomes private.
      G4FTFParticipants(const G4FTFParticipants &right);

//      std::vector<G4InteractionContent *> theInteractions;
  
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


//inline
//const G4InteractionContent & G4FTFParticipants::GetInteraction() const
inline
G4InteractionContent & G4FTFParticipants::GetInteraction()
{
	return *theInteractions[currentInteraction];
}
#endif
