// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4StopDummyDeexcitation.hh,v 1.1 1999-01-07 16:13:42 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      File name:     G4StopDummyDeexcitation.hh 
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 18 May 1998
//
//      Modifications: 
// -------------------------------------------------------------------

#ifndef G4STOPDUMMYDEEXCITATION_HH
#define G4STOPDUMMYDEEXCITATION_HH

#include "G4StopDeexcitationAlgorithm.hh"

#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4DynamicParticleVector.hh"
#include "G4ThreeVector.hh"

class G4StopDummyDeexcitation: public G4StopDeexcitationAlgorithm
{  

private:

  // Hide assignment operator as private 
  G4StopDummyDeexcitation& operator=(const G4StopDummyDeexcitation &right);

  // Copy constructor
  G4StopDummyDeexcitation(const G4StopDummyDeexcitation& );

public:

  // Constructor
  G4StopDummyDeexcitation();

  // Destructor
  virtual ~G4StopDummyDeexcitation();

  // Products
  virtual G4DynamicParticleVector* BreakUp(G4double A, G4double Z, 
				    G4double excitation, const G4ThreeVector& p);

protected:
   
private:

  G4DynamicParticleVector* _products;

};
 

#endif
