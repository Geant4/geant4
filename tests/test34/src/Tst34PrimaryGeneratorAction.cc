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
#include "Tst34PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4ios.hh"

Tst34PrimaryGeneratorAction::Tst34PrimaryGeneratorAction()
{
	particleGun=new G4GeneralParticleSource;
}

Tst34PrimaryGeneratorAction::~Tst34PrimaryGeneratorAction()
{
	delete particleGun;
}

void Tst34PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 
	particleGun->GeneratePrimaryVertex(anEvent);
}
