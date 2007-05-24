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
// $Id: PrimaryGeneratorAction.cc,v 1.1 2007-05-24 21:57:03 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Jane Tinslay, May 2007. Creation.
//
#include "PrimaryGeneratorAction.hh"

#include "ConfigData.hh"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  particleGun  = new G4ParticleGun(1);

  particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("e-"));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,-10.*cm));
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  particleGun->SetParticleEnergy(ConfigData::GetPrimaryEnergy());
  particleGun->GeneratePrimaryVertex(anEvent);
}
