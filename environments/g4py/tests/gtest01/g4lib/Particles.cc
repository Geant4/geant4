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
// $Id: Particles.cc,v 1.3 2006-06-29 15:37:02 gunter Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   Particles.cc
//
//   Physics list for defining particles
//
// ====================================================================
#include "Particles.hh"

#include "G4LeptonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4IonConstructor.hh"

// ====================================================================
//
// class description
//
// ====================================================================

//////////////////////////////////////
Particles::Particles()
  : G4VPhysicsConstructor("Particles")
//////////////////////////////////////
{
}


///////////////////////
Particles::~Particles()
///////////////////////
{
}


///////////////////////////////////
void Particles::ConstructParticle()
///////////////////////////////////
{
  G4LeptonConstructor::ConstructParticle();
  G4BosonConstructor::ConstructParticle();
  G4MesonConstructor::ConstructParticle();
  G4BaryonConstructor::ConstructParticle();
  G4ShortLivedConstructor::ConstructParticle();
  G4IonConstructor::ConstructParticle(); 
}


//////////////////////////////////
void Particles::ConstructProcess()
//////////////////////////////////
{
}

