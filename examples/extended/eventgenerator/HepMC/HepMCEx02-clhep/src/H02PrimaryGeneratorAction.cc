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

// ====================================================================
//
//   H02PrimaryGeneratorAction.cc
//   $Id: H02PrimaryGeneratorAction.cc,v 1.1 2002-11-19 10:36:20 murakami Exp $
//
// ====================================================================
#include "H02PrimaryGeneratorAction.hh"
#include "H02PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4HepMCAsciiReader.hh"
#include "H02PythiaInterface.hh"

////////////////////////////////////////////////////
H02PrimaryGeneratorAction::H02PrimaryGeneratorAction()
////////////////////////////////////////////////////
{
  // default generator is particle gun.
  currentGenerator= particleGun= new G4ParticleGun();
  currentGeneratorName= "particleGun";
  hepmcAscii= new G4HepMCAsciiReader();
  pythiaGen= new H02PythiaInterface();

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

