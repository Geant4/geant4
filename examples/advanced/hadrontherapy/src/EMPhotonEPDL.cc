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
// $Id: HadrontherapyPhotonEPDL.cc; May 2005
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

#include "EMPhotonEPDL.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4LowEnergyCompton.hh"
#include "G4LowEnergyGammaConversion.hh"
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyRayleigh.hh"
#include "G4StepLimiter.hh"


EMPhotonEPDL::EMPhotonEPDL(const G4String& name): 
   G4VPhysicsConstructor(name)
{ 
  G4cout<< "ELECTROMAGNETIC PROCESS(ES): G4LowEnergyPhotoElectric (photon)" 
        << G4endl
        << "                             G4LowEnergyCompton (photon)" 
        << G4endl
        << "                             G4LowEnergyGammaConversion (photon)" 
        << G4endl
        << "                             G4LowEnergyRayleigh (photon)" 
        << G4endl
        << "APPLIED MODEL(S): -" 
        << G4endl;
}

EMPhotonEPDL::~EMPhotonEPDL()
{ }

void EMPhotonEPDL::ConstructProcess()
{

  // **************
  // *** Photon ***
  // **************
  
  G4LowEnergyPhotoElectric* photonPhotoElectricProcess = new G4LowEnergyPhotoElectric();
  G4LowEnergyCompton* photonComptonProcess = new G4LowEnergyCompton;
  G4LowEnergyGammaConversion* photonGammaConvProcess = new G4LowEnergyGammaConversion;
  G4LowEnergyRayleigh* photonRayleighProcess = new G4LowEnergyRayleigh;

  G4StepLimiter* photonStepLimiter = new G4StepLimiter();

  G4ParticleDefinition* particle = G4Gamma::Gamma(); 
  G4ProcessManager* processManager = particle -> GetProcessManager();
  processManager -> AddDiscreteProcess(photonPhotoElectricProcess);
  processManager -> AddDiscreteProcess(photonComptonProcess);
  processManager -> AddDiscreteProcess(photonGammaConvProcess);
  processManager -> AddDiscreteProcess(photonRayleighProcess);
  processManager -> AddProcess(photonStepLimiter, -1, -1,  3);

}
