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
// Short description: The CHIPS model provides the G4QHadronVector
// output, which is converted to the G4 particle-singletons
// --------------------------------------------------------------------
//
// Created: J.P. Wellisch, 2000/08/18 
// 01.09.2008 V.Ivanchenko 
//

#include "G4GammaNuclearReaction.hh"
#include "G4Gamma.hh"
#include "G4Nucleus.hh"
#include "G4HadFinalState.hh"
#include "G4HadProjectile.hh"

G4GammaNuclearReaction::G4GammaNuclearReaction(): 
  G4HadronicInteraction("CHIPS")
{}

G4GammaNuclearReaction::~G4GammaNuclearReaction()
{}

G4HadFinalState * G4GammaNuclearReaction::ApplyYourself(
	   const G4HadProjectile& aTrack, 
	   G4Nucleus& aTargetNucleus)
{
  if(aTrack.GetDefinition() != G4Gamma::GammaDefinition())
  {
    throw G4HadronicException(__FILE__, __LINE__, 
			      "Called G4GammaNuclearReaction for particle other than gamma");
  }
  return theModel.ApplyYourself(aTrack, aTargetNucleus);
}

