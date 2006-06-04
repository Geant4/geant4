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
// $Id: Particles.cc,v 1.2 2006-06-04 21:37:25 kmura Exp $
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

