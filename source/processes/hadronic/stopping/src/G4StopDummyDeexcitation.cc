// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4StopDummyDeexcitation.cc,v 1.1 1999-01-07 16:13:47 gunter Exp $
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

#include <rw/tpordvec.h>
#include <rw/tvordvec.h>

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

G4DynamicParticleVector* G4StopDummyDeexcitation::BreakUp(G4double A, G4double Z, 
							 G4double excitation, 
							 const G4ThreeVector& p)
{
  return 0;
}



