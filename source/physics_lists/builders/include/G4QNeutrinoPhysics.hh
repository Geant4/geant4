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
// $Id: G4QNeutrinoPhysics.hh,v 1.3 2010-06-03 14:37:24 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4QNeutrinoPhysics
//
// Author: 2009 M. V. Kosov 
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef G4QNeutrinoPhysics_h
#define G4QNeutrinoPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

#include "G4QMessenger.hh"

#include "G4QInelastic.hh"

#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4NeutrinoTau.hh"
#include "G4AntiNeutrinoTau.hh"

#include "G4ProcessManager.hh"

class G4QNeutrinoPhysics : public G4VPhysicsConstructor
{
public:
  G4QNeutrinoPhysics(G4int verbose =1);
  G4QNeutrinoPhysics(const G4String& name);
  virtual ~G4QNeutrinoPhysics();

  void ConstructParticle();
  void ConstructProcess();

  G4String GetNuElNuclearOnOff()  {return nuEleOn ? "on" : "off";}
  G4String GetNuMuNuclearOnOff()  {return nuMuoOn ? "on" : "off";}
  G4String GetNuTauNuclearOnOff() {return nuTauOn ? "on" : "off";}

  void SetNuElNuclearOnOff(G4String& aState);
  void SetNuMuNuclearOnOff(G4String& aState);
  void SetNuTauNuclearOnOff(G4String& aState);
  void SetNuNuclearBias(G4double newValue);

private:

  void BuildNuEleNuclear();
  void BuildNuMuoNuclear();
  void BuildNuTauNuclear();

  G4bool wasBuilt;
  G4bool nuEleActivated;
  G4bool nuMuoActivated;
  G4bool nuTauActivated;
  G4bool nuEleOn;
  G4bool nuMuoOn;
  G4bool nuTauOn;
  G4double nuNucBias;                // Biasing factor for neutrino-nuclear processes

  G4QInelastic*            inelastic;
  G4QMessenger*            theMessenger;
};

#endif
