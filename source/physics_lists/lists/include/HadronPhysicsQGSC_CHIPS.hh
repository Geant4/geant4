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
// $Id: HadronPhysicsQGSC_CHIPS.hh,v 1.7 2010-06-03 10:42:44 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   HadronPhysicsQGSC_CHIPS
//
// Author: 2009  M. Kosov
//
// Modified:
//
//----------------------------------------------------------------------------
// Short description: In fact this is the definition of the Hadronic Inelastic
// physics. The definition of the Hadronic Elastic physics one can find in the
// G4HadronQElasticPhysics, which is stable (the same for all physics lists).
// The only "unstable" part of the physics is the Hadronic Inelastic physics,
// which is usually composed of the wixing of the High Energy Inelastic Model
// (HEIM) and the Low Energy Inelastic Model (LEIM), which are applied only for
// some hadrons (mostly nucleons and pi-mesons), above the LHEP model, which
// usually covers all particles (but for Sigma_0 ?) and sometimes covers the
// "hole" between the LEIM and HIME at intermediate energies. The name of the
// Physics list is usually have a form HEIM_LEIM and the inelastic interactions
// are defined in the HadronicPhysicsHEIM_LEIM class. So in this particular
// physics list the low energy model is CHIPS (G4QCollision process) and the
// high energy model is QGSC (QGS with the Energy Flow interface to CHIPS),
// which are in terms of the energy boundary are mixed not on the model level,
// but on the process level (G4DiscProcessMixer class). The LHEP is completely
// excluded from this physics list, because the MiscLHEP is substituted by the
// MiscQGSC class (QGS with the Energy Flow interface to CHIPS), covering all
// particles, which are not N, pi, or K, defined by the separate builders. 
//---------------------------------------------------------------------------

#ifndef HadronPhysicsQGSC_CHIPS_h
#define HadronPhysicsQGSC_CHIPS_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4MiscQGSCBuilder.hh"

#include "G4PiKBuilder.hh"
#include "G4QGSC_CHIPSPiKBuilder.hh"

#include "G4QProtonBuilder.hh"
#include "G4QGSC_CHIPSProtonBuilder.hh"

#include "G4QNeutronBuilder.hh"
#include "G4QGSC_CHIPSNeutronBuilder.hh"

class HadronPhysicsQGSC_CHIPS : public G4VPhysicsConstructor
{
public: 
  HadronPhysicsQGSC_CHIPS(G4int verbose =1);
  HadronPhysicsQGSC_CHIPS(const G4String& name, G4bool quasiElastic=true);
  virtual ~HadronPhysicsQGSC_CHIPS();

public: 
  virtual void ConstructParticle();
  virtual void ConstructProcess();

private:
  void CreateModels();

  G4QNeutronBuilder* theNeut;
  G4QGSC_CHIPSNeutronBuilder* theQGSCNeut;
    
  G4PiKBuilder* thePiK;
  G4QGSC_CHIPSPiKBuilder* theQGSCPiK;
    
  G4QProtonBuilder* theProt;
  G4QGSC_CHIPSProtonBuilder * theQGSCProt;    
    
  G4MiscQGSCBuilder* theMiscQGSC;
    
  G4bool QuasiElastic;
};

#endif

