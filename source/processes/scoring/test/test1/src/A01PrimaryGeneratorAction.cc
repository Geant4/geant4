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
// $Id: A01PrimaryGeneratorAction.cc,v 1.1 2006-07-14 14:43:33 asaim Exp $
// --------------------------------------------------------------
//

#include "A01PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"

A01PrimaryGeneratorAction::A01PrimaryGeneratorAction()
{
    G4int n_particle = 1;
    particleGun  = new G4ParticleGun(n_particle);
    particleGun->SetParticlePosition(G4ThreeVector(0.,0.,-0.65*m));
    particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
}

A01PrimaryGeneratorAction::~A01PrimaryGeneratorAction() 
{
    delete particleGun;
}

void A01PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) 
{
    particleGun->GeneratePrimaryVertex(anEvent);
}

