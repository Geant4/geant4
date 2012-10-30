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
//      File name:     G4StopDeexcitation
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 8 May 1998
//
// -------------------------------------------------------------------


#include "G4StopDeexcitation.hh"
#include <vector>

#include "globals.hh"
#include "Randomize.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "G4HadronicDeprecate.hh"


// Constructor

G4StopDeexcitation::G4StopDeexcitation(G4StopDeexcitationAlgorithm* algorithm)
{
  G4HadronicDeprecate("G4StopDeexcitation");
  _algorithm = algorithm;
}


// Destructor

G4StopDeexcitation::~G4StopDeexcitation()
{
  delete _algorithm;
}

G4ReactionProductVector* G4StopDeexcitation::DoBreakUp(G4double A, G4double Z, 
						       G4double excitation, 
						       const G4ThreeVector& p) const
{
  G4ReactionProductVector* v = 0;
  if (_algorithm != 0) 
    {
      v = _algorithm->BreakUp(A,Z,excitation,p);
    }
  return v;
}
