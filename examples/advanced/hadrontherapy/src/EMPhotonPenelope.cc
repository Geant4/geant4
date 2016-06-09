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
// $Id: HadrontherapyPhotonPenelope.cc; Version 4.0 May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#include "EMPhotonPenelope.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4PenelopeCompton.hh"
#include "G4PenelopeGammaConversion.hh"
#include "G4PenelopePhotoElectric.hh"
#include "G4PenelopeRayleigh.hh"
#include "G4StepLimiter.hh"


EMPhotonPenelope::EMPhotonPenelope(const G4String& name): 
   G4VPhysicsConstructor(name)
{ 
  G4cout<< "ELECTROMAGNETIC PROCESS(ES): G4PenelopePhotoElectric (photon)" 
        << G4endl
        << "                             G4PenelopeCompton (photon)" 
        << G4endl
        << "                             G4PenelopeGammaConversion (photon)" 
        << G4endl
        << "                             G4PenelopeRayleigh (photon)" 
        << G4endl
        << "APPLIED MODEL(S): -" 
        << G4endl;
}

EMPhotonPenelope::~EMPhotonPenelope()
{ }

void EMPhotonPenelope::ConstructProcess()
{

  // **************
  // *** Photon ***
  // **************
  
  G4PenelopePhotoElectric* photonPhotoElectricProcess = new G4PenelopePhotoElectric();
  G4PenelopeCompton* photonComptonProcess = new G4PenelopeCompton;
  G4PenelopeGammaConversion* photonGammaConvProcess = new G4PenelopeGammaConversion;
  G4PenelopeRayleigh* photonRayleighProcess = new G4PenelopeRayleigh;

  G4StepLimiter* photonStepLimiter = new G4StepLimiter();

  G4ParticleDefinition* particle = G4Gamma::Gamma(); 
  G4ProcessManager* processManager = particle -> GetProcessManager();
  processManager -> AddDiscreteProcess(photonPhotoElectricProcess);
  processManager -> AddDiscreteProcess(photonComptonProcess);
  processManager -> AddDiscreteProcess(photonGammaConvProcess);
  processManager -> AddDiscreteProcess(photonRayleighProcess);
  processManager -> AddProcess(photonStepLimiter, -1, -1,  3);

}


