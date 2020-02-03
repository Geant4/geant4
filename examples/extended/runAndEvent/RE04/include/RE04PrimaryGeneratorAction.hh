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
/// \file runAndEvent/RE04/include/RE04PrimaryGeneratorAction.hh
/// \brief Definition of the RE04PrimaryGeneratorAction class
//
//
#ifndef RE04PrimaryGeneratorAction_h
#define RE04PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

//
/// User primary particle generator class
///
/// - void GeneratePrimaries(G4Event*)
///     an incident particle is mu- with 10 GeV energy at the position 
///     (-75 cm,y,0) toward the (1,0,0) direction. The y position is uniformly
///     varied from 95.5 cm to 96.5 cm.
//
class RE04PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    RE04PrimaryGeneratorAction();    
    virtual ~RE04PrimaryGeneratorAction();

  public:
    virtual void GeneratePrimaries(G4Event*);

  private:
    G4ParticleGun* fParticleGun;

};

#endif


