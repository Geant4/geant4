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
/// \file persistency/gdml/G02/include/G02PrimaryGeneratorAction.hh
/// \brief Definition of the G02PrimaryGeneratorAction class
//
//
//
// Class G02PrimaryGeneratorAction
//
// Simple primary-generator action class.
//
// ----------------------------------------------------------------------------

#ifndef G02PrimaryGeneratorAction_h
#define G02PrimaryGeneratorAction_h 1

#include "G4ParticleTable.hh"
#include "G4ParticleGun.hh"
#include "G4Event.hh"

#include "G4VUserPrimaryGeneratorAction.hh"

// ----------------------------------------------------------------------------

/// Primary generator action used in GDML read/write example

class G02PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:

    // Constructor and destructor
    //
    G02PrimaryGeneratorAction();
   ~G02PrimaryGeneratorAction();

    // Used by Geant4 to generate the primary particles of the event
    //
    virtual void GeneratePrimaries(G4Event* anEvent);

  private:

    G4ParticleGun* fParticleGun;
    G4ParticleTable * fParticleTable;
};

// ----------------------------------------------------------------------------

#endif
