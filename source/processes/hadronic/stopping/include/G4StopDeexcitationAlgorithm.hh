// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4StopDeexcitationAlgorithm.hh,v 1.2 1999-04-18 11:28:56 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      File name:     G4StopDeexcitationAlgorithm.hh 
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 18 May 1998
//
//      Modifications: 
// -------------------------------------------------------------------

#ifndef G4STOPDEEXCITATIONALGORITHM_HH
#define G4STOPDEEXCITATIONALGORITHM_HH

#include "globals.hh"
#include "G4ReactionProduct.hh"
#include "G4ReactionProductVector.hh"
#include "G4ThreeVector.hh"

class G4StopDeexcitationAlgorithm 
{  

private:

  // Hide assignment operator as private 
  G4StopDeexcitationAlgorithm& operator=(const G4StopDeexcitationAlgorithm &right);

  // Copy constructor
  G4StopDeexcitationAlgorithm(const G4StopDeexcitationAlgorithm& );

public:

  // Constructor
  G4StopDeexcitationAlgorithm() {};

  // Destructor
  virtual ~G4StopDeexcitationAlgorithm() {};

  // Products
  virtual G4ReactionProductVector* BreakUp(G4double A, G4double Z, 
					   G4double excitation, const G4ThreeVector& p) =0;

protected:
   
private:

  
};
 

#endif
