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
// $Id: RemSimPrimaryGeneratorAction.cc,v 1.1 2004-01-30 12:25:44 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "RemSimPrimaryGeneratorAction.hh"
#include "RemSimBasicGenerator.hh"
#include "RemSimVPrimaryGeneratorFactory.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

RemSimPrimaryGeneratorAction::RemSimPrimaryGeneratorAction()
{
  primaryFactory = new RemSimBasicGenerator();
}

RemSimPrimaryGeneratorAction::~RemSimPrimaryGeneratorAction()
{
  delete primaryFactory;
}

void RemSimPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  primaryFactory -> GeneratePrimaries(anEvent);
}


