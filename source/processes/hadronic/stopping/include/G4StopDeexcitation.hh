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
// $Id: G4StopDeexcitation.hh,v 1.4 2001-07-11 10:08:07 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
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
#include "G4ReactionProduct.hh"
#include "G4ReactionProductVector.hh"
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
  G4ReactionProductVector* DoBreakUp(G4double A, G4double Z, 
				     G4double excitation, const G4ThreeVector& p) const;


private:

  // Hide assignment operator as private 
  G4StopDeexcitation& operator=(const G4StopDeexcitation &right);

  // Copy constructor
  G4StopDeexcitation(const G4StopDeexcitation& );

  G4StopDeexcitationAlgorithm* _algorithm; // owned pointer

};
 
#endif




















