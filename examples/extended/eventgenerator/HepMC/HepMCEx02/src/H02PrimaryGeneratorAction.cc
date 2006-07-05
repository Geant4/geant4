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
//   H02PrimaryGeneratorAction.cc
//   $Id: H02PrimaryGeneratorAction.cc,v 1.4 2006-07-05 12:04:13 gcosmo Exp $
//
// ====================================================================
#include "H02PrimaryGeneratorAction.hh"
#include "H02PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "HepMCG4AsciiReader.hh"
#include "HepMCG4PythiaInterface.hh"

////////////////////////////////////////////////////
H02PrimaryGeneratorAction::H02PrimaryGeneratorAction()
////////////////////////////////////////////////////
{
  // default generator is particle gun.
  currentGenerator= particleGun= new G4ParticleGun();
  currentGeneratorName= "particleGun";
  hepmcAscii= new HepMCG4AsciiReader();
#ifdef G4LIB_USE_PYTHIA
  pythiaGen= new HepMCG4PythiaInterface();
#else
  pythiaGen= 0;
#endif

  gentypeMap["particleGun"]= particleGun;
  gentypeMap["hepmcAscii"]= hepmcAscii;
  gentypeMap["pythia"]= pythiaGen;

  messenger= new H02PrimaryGeneratorMessenger(this);  
}

/////////////////////////////////////////////////////
H02PrimaryGeneratorAction::~H02PrimaryGeneratorAction()
/////////////////////////////////////////////////////
{
  delete messenger;
}

//////////////////////////////////////////////////////////////////
void H02PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
//////////////////////////////////////////////////////////////////
{
  if(currentGenerator)
    currentGenerator-> GeneratePrimaryVertex(anEvent);
  else 
    G4Exception("generator is not instanciated.");
}

