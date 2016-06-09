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
// ====================================================================
//
//   ExN04PrimaryGeneratorAction.cc
//   $Id: ExN04PrimaryGeneratorAction.cc,v 1.3 2006/06/29 17:08:48 gunter Exp $
//
// ====================================================================
#include "ExN04PrimaryGeneratorAction.hh"
#include "ExN04PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4HepMCAsciiReader.hh"
#include "H01PythiaInterface.hh"


//////////////////////////////////////////////////////////
ExN04PrimaryGeneratorAction::ExN04PrimaryGeneratorAction()
//////////////////////////////////////////////////////////
{
  // default generator is particle gun.
  currentGenerator= particleGun= new G4ParticleGun();
  currentGeneratorName= "particleGun";
  hepmcAscii= new G4HepMCAsciiReader();
  pythiaGen= new H01PythiaInterface();

  gentypeMap["particleGun"]= particleGun;
  gentypeMap["hepmcAscii"]= hepmcAscii;
  gentypeMap["pythia"]= pythiaGen;

  messenger= new ExN04PrimaryGeneratorMessenger(this);  
}

///////////////////////////////////////////////////////////
ExN04PrimaryGeneratorAction::~ExN04PrimaryGeneratorAction()
///////////////////////////////////////////////////////////
{
  delete messenger;
}

/////////////////////////////////////////////////////////////////////
void ExN04PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
/////////////////////////////////////////////////////////////////////
{
  if(currentGenerator)
    currentGenerator-> GeneratePrimaryVertex(anEvent);
  else 
    G4Exception("generator is not instanciated.");
}

