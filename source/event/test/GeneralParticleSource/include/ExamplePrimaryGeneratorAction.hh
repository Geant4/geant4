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
#ifndef ExamplePrimaryGeneratorAction_h
#define ExamplePrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4GeneralParticleSource;
class G4Event;

class ExamplePrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    ExamplePrimaryGeneratorAction();
    ~ExamplePrimaryGeneratorAction();

  // DEBUG arrays
  G4int debugx[100], debugy[100], debugz[100];
  G4int debugpx[100], debugpy[100], debugpz[100];
  G4int debugtheta[100], debugphi[100];
  G4int debugenergy[100];
  G4double debug_energy_step;
  G4double DebugXmin, DebugXmax, DebugYmin, DebugYmax;
  G4double DebugZmin, DebugZmax, DebugXStep, DebugYStep;
  G4double DebugZStep;

  G4String SourceType;
  G4String SourceShape;
  G4double radius;
  G4double halfx;
  G4double halfy;
  G4double halfz;
  G4ThreeVector centre;

    // for use with energy
  G4String EneDisType ;
  G4String InterpolationType;
  G4double emin;
  G4double emax;

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    G4GeneralParticleSource* particleGun;
  
};

#endif



