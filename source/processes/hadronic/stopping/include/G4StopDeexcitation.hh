// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4StopDeexcitation.hh,v 1.1 1999-01-07 16:13:41 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      File name:     G4StopDeexcitation.hh 
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 12 May 1998
//
//      Modifications: 
// -------------------------------------------------------------------

#ifndef G4STOPDEEXCITATION_HH
#define G4STOPDEEXCITATION_HH

#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4DynamicParticleVector.hh"
#include "G4ThreeVector.hh"
#include "G4StopDeexcitationAlgorithm.hh"

class G4StopDeexcitation
{  

public:

  // Constructor
  G4StopDeexcitation(G4StopDeexcitationAlgorithm* algorithm);

  // Destructor
  ~G4StopDeexcitation();

  // Return final absorption products
  G4DynamicParticleVector* DoBreakUp(G4double A, G4double Z, 
				     G4double excitation, const G4ThreeVector& p) const;


private:

  // Hide assignment operator as private 
  G4StopDeexcitation& operator=(const G4StopDeexcitation &right);

  // Copy constructor
  G4StopDeexcitation(const G4StopDeexcitation& );

  G4StopDeexcitationAlgorithm* _algorithm; // owned pointer

};
 
#endif




















