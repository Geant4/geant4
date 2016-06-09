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
// $Id: G4ConcreteNStarNToNN.cc,v 1.1 2003/10/07 12:37:38 hpw Exp $ //

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ConcreteNStarNToNN.hh"
#include "G4NNstarBuilder.hh"

G4XNNstarTable G4ConcreteNStarNToNN::theSigmaTable;

G4ConcreteNStarNToNN::G4ConcreteNStarNToNN(const G4ParticleDefinition* aPrimary,
					   const G4ParticleDefinition* bPrimary,
					   const G4ParticleDefinition* aSecondary,
					   const G4ParticleDefinition* bSecondary)
  : G4ConcreteNNTwoBodyResonance(aPrimary, bPrimary, aSecondary, bSecondary,
                                 G4NNstarBuilder(aPrimary->GetParticleName(), theSigmaTable))
{
}

G4ConcreteNStarNToNN::~G4ConcreteNStarNToNN()
{ 
}
