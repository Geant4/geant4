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
//
/// \file Par02PrimaryGeneratorAction.hh
/// \brief Definition of the Par02PrimaryGeneratorAction class

#ifndef PAR02PrimaryGeneratorAction_h
#define PAR02PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

/// Construction of a primary generation action.
///
/// For simplicity, we use here the particle gun, but in the original
/// application for FCC (developed by Anna Zaborowska), the Monte Carlo
/// event generator Pythia8 is used as generator and it is interfaced
/// to Geant4 via HepMC.

class Par02PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
  public:
    Par02PrimaryGeneratorAction();
    ~Par02PrimaryGeneratorAction();

    virtual void GeneratePrimaries( G4Event* anEvent );
    G4ParticleGun* GetParticleGun();

  private:
    G4ParticleGun* fParticleGun;
};

#endif

