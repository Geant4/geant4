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
// $Id: G4VParticipants.cc 100828 2016-11-02 15:25:59Z gcosmo $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4VParticipants ----------------
//             by Gunter Folger, May 1998.
//      abstract class finding participants in a hadron Nucleus collision
//       in Parton String Models.
// ------------------------------------------------------------
// 20110805  M. Kelsey -- Reduce external rebuilds: move #include, Init()
//		and SetNucleus() here.

#include "G4VParticipants.hh"
#include "G4Fancy3DNucleus.hh"


G4VParticipants::G4VParticipants() : theNucleus(NULL), 
                                     theProjectileNucleus(NULL)
{}


G4VParticipants::~G4VParticipants()
{
  // G4cout << "G4VParticipants::~G4VParticipants()" << G4endl;
  if ( theNucleus != NULL ) delete theNucleus;
  if ( theProjectileNucleus != NULL ) delete theProjectileNucleus;
}


void G4VParticipants::Init(G4int theA, G4int theZ)
{
  if ( theNucleus == NULL ) theNucleus = new G4Fancy3DNucleus();
  theNucleus->Init(theA, theZ);
  theNucleus->SortNucleonsIncZ();
}


void G4VParticipants::SetNucleus(G4V3DNucleus * aNucleus)
{
  if (theNucleus != NULL) delete theNucleus;
  theNucleus = aNucleus;
}

void G4VParticipants::InitProjectileNucleus(G4int theA, G4int theZ)
{
  if ( theProjectileNucleus == NULL ) theProjectileNucleus = new G4Fancy3DNucleus();
  theProjectileNucleus->Init(theA, theZ);
  theProjectileNucleus->SortNucleonsDecZ();
}


void G4VParticipants::SetProjectileNucleus(G4V3DNucleus * aNucleus)
{
  if (theProjectileNucleus != NULL) delete theProjectileNucleus;
  theProjectileNucleus = aNucleus;
}

