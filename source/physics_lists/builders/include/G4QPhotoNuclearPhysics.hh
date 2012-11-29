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
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   G4QPhotoNuclearPhysics
//
// Author: 2009 M. V. Kosov
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef G4QPhotoNuclearPhysics_h
#define G4QPhotoNuclearPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

#include "G4QMessenger.hh"

#include "G4QSynchRad.hh"
#include "G4QInelastic.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4TauPlus.hh"
#include "G4TauMinus.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ProcessManager.hh"


class G4QPhotoNuclearPhysics : public G4VPhysicsConstructor
{
public:
  G4QPhotoNuclearPhysics(G4int verbose =1);
  G4QPhotoNuclearPhysics(const G4String& name);
  virtual ~G4QPhotoNuclearPhysics();

  void ConstructParticle();
  void ConstructProcess();

  G4String GetSynchRadOnOff()     {return synchrOn ? "on" : "off";}
  G4String GetGammaNuclearOnOff() {return gamNucOn ? "on" : "off";}
  G4String GetElPosNuclearOnOff() {return eleNucOn ? "on" : "off";}
  G4String GetMuonNuclearOnOff()  {return muoNucOn ? "on" : "off";}
  G4String GetTauNuclearOnOff()   {return tauNucOn ? "on" : "off";}

  void SetSynchRadOnOff(G4String& aSwitch);
  void SetGammaNuclearOnOff(G4String& aSwitch);
  void SetElPosNuclearOnOff(G4String& aSwitch);
  void SetMuonNuclearOnOff(G4String& aSwitch);
  void SetTauNuclearOnOff(G4String& aSwitch);
  void SetMinGammaSR(G4double newValue);
  void SetPhotoNucBias(G4double newValue);

private:

  void BuildSynchRad();
  void BuildGammaNuclear();
  void BuildElectroNuclear();
  void BuildMuonNuclear();
  void BuildTauNuclear();

  G4bool   wasBuilt;                    // Flag of forbidden reactivation of processes
  G4bool   SynchRActivated;             // Flag of finished activation of SynchroRadiation
  G4bool   GamNucActivated;             // Flag of finished activation of gamma-nuclear
  G4bool   EleNucActivated;             // Flag of finished activation of electron-nuclear
  G4bool   MuoNucActivated;             // Flag of finished activation of muon-nuclear
  G4bool   TauNucActivated;             // Flag of finished activation of tau-nuclear
  G4bool   synchrOn;                    // Switch flag for Synchrotron Radiation process
  G4double synchrMinGam;                // MinimumGamma for SynchrotronRadiation activation
  G4bool   gamNucOn;                    // Switch flag for the gamma-nuclear process
  G4bool   eleNucOn;                    // Switch flag for the electron-nuclear process
  G4bool   muoNucOn;                    // Switch flag for the electron-nuclear process
  G4bool   tauNucOn;                    // Switch flag for the electron-nuclear process
  G4double photoNucBias;                // Biasing factor for photo-nuclear processes

  G4QInelastic*           inelastic;
  G4QSynchRad*            synchrad;
  G4QMessenger*           theMessenger;
};

#endif
