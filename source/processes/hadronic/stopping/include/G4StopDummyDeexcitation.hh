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
//      File name:     G4StopDummyDeexcitation.hh 
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 18 May 1998
//
// -------------------------------------------------------------------

#ifndef G4STOPDUMMYDEEXCITATION_HH
#define G4STOPDUMMYDEEXCITATION_HH

#include "G4StopDeexcitationAlgorithm.hh"

#include "globals.hh"
#include "G4ReactionProduct.hh"
#include "G4ReactionProductVector.hh"
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
  virtual G4ReactionProductVector* BreakUp(G4double A, G4double Z, 
				    G4double excitation, const G4ThreeVector& p);

protected:
   
private:

  G4ReactionProductVector* _products;

};
 

#endif
