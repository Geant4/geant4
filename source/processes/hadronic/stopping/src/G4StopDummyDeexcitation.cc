// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4StopDummyDeexcitation.cc,v 1.4 1999-12-15 14:53:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      File name:     G4StopDummyDeexcitation
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 8 May 1998
//
//      Modifications: 
// -------------------------------------------------------------------

#include "G4ios.hh"

#include "G4StopDummyDeexcitation.hh"

#include "g4rw/tpordvec.h"
#include "g4rw/tvordvec.h"

#include "globals.hh"
#include "G4DynamicParticleVector.hh"
#include "G4ThreeVector.hh"

// Constructor

G4StopDummyDeexcitation::G4StopDummyDeexcitation()
  
{}


// Destructor

G4StopDummyDeexcitation::~G4StopDummyDeexcitation()
{
}

G4ReactionProductVector* G4StopDummyDeexcitation::BreakUp(G4double A, G4double Z, 
							 G4double excitation, 
							 const G4ThreeVector& p)
{
  return 0;
}



