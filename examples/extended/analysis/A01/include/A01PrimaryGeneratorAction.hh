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
/// \file analysis/A01/include/A01PrimaryGeneratorAction.hh
/// \brief Definition of the A01PrimaryGeneratorAction class
//
// $Id$
// --------------------------------------------------------------
//

#ifndef A01PrimaryGeneratorAction_h
#define A01PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class G4ParticleDefinition;
class A01PrimaryGeneratorMessenger;

class A01PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    A01PrimaryGeneratorAction();
    virtual ~A01PrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event*);

  private:
    G4ParticleGun* fParticleGun;
    A01PrimaryGeneratorMessenger* fGunMessenger;
    G4ParticleDefinition* fPositron;
    G4ParticleDefinition* fMuon;
    G4ParticleDefinition* fPion;
    G4ParticleDefinition* fKaon;
    G4ParticleDefinition* fProton;
    G4double fMomentum;
    G4double fSigmaMomentum;
    G4double fSigmaAngle;
    G4bool fRandomizePrimary;

  public:
    inline void SetMomentum(G4double val) { fMomentum = val; }
    inline G4double GetMomentum() const { return fMomentum; }
    inline void SetSigmaMomentum(G4double val) { fSigmaMomentum = val; }
    inline G4double GetSigmaMomentum() const { return fSigmaMomentum; }
    inline void SetSigmaAngle(G4double val) { fSigmaAngle = val; }
    inline G4double GetSigmaAngle() const { return fSigmaAngle; }
    inline void SetRandomize(G4bool val) { fRandomizePrimary = val; }
    inline G4bool GetRandomize() const { return fRandomizePrimary; }
};

#endif
