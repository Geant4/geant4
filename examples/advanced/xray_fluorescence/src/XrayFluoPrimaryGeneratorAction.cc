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
//
// $Id: XrayFluoPrimaryGeneratorAction.cc
// GEANT4 tag $Name: xray_fluo-V03-02-00
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
//
// -------------------------------------------------------------------

#include "XrayFluoPrimaryGeneratorAction.hh"
//#include "G4DataVector.hh"
#include "XrayFluoDetectorConstruction.hh"
#include "XrayFluoPrimaryGeneratorMessenger.hh"
#include "XrayFluoRunAction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "XrayFluoAnalysisManager.hh"
#include "XrayFluoDataSet.hh"
#ifdef G4ANALYSIS_USE  
#include "AIDA/AIDA.h"
#endif
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoPrimaryGeneratorAction::XrayFluoPrimaryGeneratorAction(XrayFluoDetectorConstruction* XrayFluoDC)
  :rndmFlag("off"),beam("off"),spectrum("off"),isoVert("off"),phaseSpaceGunFlag(false), 
   rayleighFlag(true), particleEnergies(0),  particleTypes(0), detectorPosition(0) 
{

  XrayFluoDetector = XrayFluoDC;

  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  
  //create a messenger for this class
  gunMessenger = new XrayFluoPrimaryGeneratorMessenger(this);
  runManager = new XrayFluoRunAction();
  
  // default particle kinematic
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="gamma");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  particleGun->SetParticleEnergy(10. * keV);

  G4double position = -0.5*(XrayFluoDetector->GetWorldSizeZ());
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,position));


  
  G4cout << "XrayFluoPrimaryGeneratorAction created" << G4endl;
  
}


void XrayFluoPrimaryGeneratorAction::ActivatePhaseSpace(G4String fileName) {

  // load phase-space
#ifdef G4ANALYSIS_USE     
  phaseSpaceGunFlag = true;

  // reads the data stored on disk form previous runs 
  // and get these data to data members

  XrayFluoAnalysisManager* analysis =  XrayFluoAnalysisManager::getInstance();
  analysis->LoadGunData(fileName, rayleighFlag);
  particleEnergies = analysis->GetEmittedParticleEnergies();
  particleTypes = analysis->GetEmittedParticleTypes();
  detectorPosition = XrayFluoDetector->GetDetectorPosition();
  detectorPosition.setR(detectorPosition.r()-(5.*cm)); // 5 cm before the detector, so in front of it.
#endif

}


void XrayFluoPrimaryGeneratorAction::SetRayleighFlag (G4bool value)
{

  G4cout <<  "value: " << value << G4endl; 
  rayleighFlag = value;
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoPrimaryGeneratorAction::~XrayFluoPrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
  delete runManager;

  G4cout << "XrayFluoPrimaryGeneratorAction deleted" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  // 
  G4double z0 = -0.5*(XrayFluoDetector->GetWorldSizeZ());
  G4double y0 = 0.*cm, x0 = 0.*cm;
  if (rndmFlag == "on")
    {
      y0 = (XrayFluoDetector->GetDia3SizeXY())/std::sqrt(2.)*(G4UniformRand()-0.5); // it was GetSampleSizeXY(), 
      x0 = (XrayFluoDetector->GetDia3SizeXY())/std::sqrt(2.)*(G4UniformRand()-0.5); // not divided by std::sqrt(2.)
    } 
  particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  
  //randomize starting point
  if (beam == "on")
    {
      G4double radius = 0.5 * mm;
      G4double rho = radius*std::sqrt(G4UniformRand());
      G4double theta = 2*pi*G4UniformRand()*rad;
      G4double position = -0.5*(XrayFluoDetector->GetWorldSizeZ());
      
      G4double y = rho * std::sin(theta);
      G4double x = rho * std::cos(theta);
      
      particleGun->SetParticlePosition(G4ThreeVector(x,y,position));
    }
  //shoot particles according to a certain spectrum
  if (spectrum =="on")
    {
      G4String particle =  particleGun->GetParticleDefinition()
	->GetParticleName();
      if(particle == "proton"|| particle == "alpha")
	{
	  G4DataVector* energies =  runManager->GetEnergies();
	  G4DataVector* data =  runManager->GetData();
	 
	  G4double sum = runManager->GetDataSum();
	  G4double partSum = 0;
	  G4int j = 0;
	  G4double random= sum*G4UniformRand();
	  while (partSum<random)
	    {
	      partSum += (*data)[j];
	      j++;
	    }
	 
	  particleGun->SetParticleEnergy((*energies)[j]);
	
	}
      else if (particle == "gamma")
	{
	  const XrayFluoDataSet* dataSet = runManager->GetGammaSet();
	  
	  G4int i = 0;
	  G4int id = 0;
	  G4double minEnergy = 0. * keV;
	  G4double particleEnergy= 0.;
	  G4double maxEnergy = 10. * keV;
	  G4double energyRange = maxEnergy - minEnergy;

	   while ( i == 0)
	    {
	      G4double random = G4UniformRand();
	      
	      G4double randomNum = G4UniformRand(); //*5.0E6;
	      
	      particleEnergy = (random*energyRange) + minEnergy;
	      
	      if ((dataSet->FindValue(particleEnergy,id)) > randomNum)
		{
		  i = 1;
		  
		}
	    }
	   particleGun->SetParticleEnergy(particleEnergy);
	}
    }

  // Randomize starting point and direction
  
  if (isoVert == "on")
    {
      G4double rho = 1. *m;
      //theta in [0;pi/2]
      G4double theta = (pi/2)*G4UniformRand();
      //phi in [-pi;pi]
      G4double phi = (G4UniformRand()*2*pi)- pi;
      G4double x = rho*std::sin(theta)*std::sin(phi);
      G4double y = rho*std::sin(theta)*std::cos(phi);
      G4double z = -(rho*std::cos(theta));
      particleGun->SetParticlePosition(G4ThreeVector(x,y,z));
      
      G4double Xdim = XrayFluoDetector->GetSampleSizeXY();
      G4double Ydim = XrayFluoDetector->GetSampleSizeXY();
      
      G4double Dx = Xdim*(G4UniformRand()-0.5);
      
      G4double Dy = Ydim*(G4UniformRand()-0.5);
      
      particleGun->SetParticleMomentumDirection(G4ThreeVector(-x+Dx,-y+Dy,-z));
      
    }

  // using prevoiously genereated emissions from sample.....

  if (phaseSpaceGunFlag){

    size_t i = anEvent->GetEventID();

    particleGun->SetParticlePosition(detectorPosition); 
    particleGun->SetParticleMomentumDirection(detectorPosition);
    G4double energy;
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle;
    if (i < particleEnergies->size()) {      
      energy = (*particleEnergies)[i];
      particle = particleTable->FindParticle((*particleTypes)[i]);
    }

    else {
      energy = 0.;
      particle = particleTable->FindParticle("gamma");
    }

    particleGun->SetParticleEnergy(energy);
    particleGun->SetParticleDefinition(particle);


  }

#ifdef G4ANALYSIS_USE 

  G4double partEnergy = particleGun->GetParticleEnergy();
  XrayFluoAnalysisManager* analysis =  XrayFluoAnalysisManager::getInstance();
  analysis->analysePrimaryGenerator(partEnergy/keV);

#endif
 
  particleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....









