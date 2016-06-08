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
// $Id: G4ConcreteNNToNDeltaStar.cc,v 1.4 2002/12/12 19:17:50 gunter Exp $ //

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ConcreteNNToNDeltaStar.hh"
#include "G4NDeltastarBuilder.hh"

G4XNDeltastarTable G4ConcreteNNToNDeltaStar::theSigmaTable;

G4ConcreteNNToNDeltaStar::G4ConcreteNNToNDeltaStar(const G4ParticleDefinition* aPrimary,
					   const G4ParticleDefinition* bPrimary,
					   const G4ParticleDefinition* aSecondary,
					   const G4ParticleDefinition* bSecondary)
  : G4ConcreteNNTwoBodyResonance(aPrimary, bPrimary, aSecondary, bSecondary,
                                 G4NDeltastarBuilder(bSecondary->GetParticleName(), theSigmaTable))
{
}

G4ConcreteNNToNDeltaStar::~G4ConcreteNNToNDeltaStar()
{ 
}
