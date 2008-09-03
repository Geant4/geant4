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
// $Id: AXPETPrimaryGeneratorAction.hh,v 1.1 2008-09-03 13:34:03 gcosmo Exp $
// ------------------------------------------------------------
// Geant4 class header file
//
// 03/09/2008, by T.Nikitina
// ------------------------------------------------------------

#ifndef AXPETPrimaryGeneratorAction_h
#define AXPETPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

#include "globals.hh"

class G4ParticleGun;
class G4Event;
class AXPETPrimaryGeneratorMessenger;

class AXPETPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    AXPETPrimaryGeneratorAction();
   ~AXPETPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event*);

    void SetOptPhotonPolar();
    void SetOptPhotonPolar(G4double);

  private:
    G4ParticleGun* particleGun;
    AXPETPrimaryGeneratorMessenger* gunMessenger;
};

#endif 
