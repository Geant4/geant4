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
// $Id: A01PrimaryGeneratorAction.hh,v 1.4 2006-06-29 16:31:35 gunter Exp $
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
    G4ParticleGun* particleGun;
    A01PrimaryGeneratorMessenger* gunMessenger;
    G4ParticleDefinition* positron;
    G4ParticleDefinition* muon;
    G4ParticleDefinition* pion;
    G4ParticleDefinition* kaon;
    G4ParticleDefinition* proton;
    G4double momentum;
    G4double sigmaMomentum;
    G4double sigmaAngle;
    G4bool randomizePrimary;

  public:
    inline void SetMomentum(G4double val) { momentum = val; }
    inline G4double GetMomentum() const { return momentum; }
    inline void SetSigmaMomentum(G4double val) { sigmaMomentum = val; }
    inline G4double GetSigmaMomentum() const { return sigmaMomentum; }
    inline void SetSigmaAngle(G4double val) { sigmaAngle = val; }
    inline G4double GetSigmaAngle() const { return sigmaAngle; }
    inline void SetRandomize(G4bool val) { randomizePrimary = val; }
    inline G4bool GetRandomize() const { return randomizePrimary; }
};

#endif


