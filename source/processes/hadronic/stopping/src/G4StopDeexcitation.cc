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
// $Id: G4StopDeexcitation.cc,v 1.8 2002-12-12 19:18:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     G4StopDeexcitation
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 8 May 1998
//
//      Modifications: 
// -------------------------------------------------------------------


#include "G4StopDeexcitation.hh"
#include "g4std/vector"

#include "globals.hh"
#include "Randomize.hh"
#include "G4ParticleTypes.hh"
#include "G4NucleiPropertiesTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"


// Constructor

G4StopDeexcitation::G4StopDeexcitation(G4StopDeexcitationAlgorithm* algorithm)
  
{
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
  if (_algorithm != 0) 
    {
      return _algorithm->BreakUp(A,Z,excitation,p);
    }
  else 
    return 0;
}
