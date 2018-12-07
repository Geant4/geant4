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
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
//
// -------------------------------------------------------------------

#include "XrayFluoPrimaryGeneratorAction.hh"
#include "XrayFluoDetectorConstruction.hh"
#include "XrayFluoPrimaryGeneratorMessenger.hh"
#include "XrayFluoRunAction.hh"
#include "G4Event.hh"
#include "G4Gamma.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4MTRunManager.hh"
#include "Randomize.hh"
#include "XrayFluoAnalysisManager.hh"
#include "XrayFluoDataSet.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoPrimaryGeneratorAction::XrayFluoPrimaryGeneratorAction(const 
							       XrayFluoDetectorConstruction* XrayFluoDC)
  :rndmFlag("off"),beam("off"),spectrum("off"),isoVert("off"),phaseSpaceGunFlag(false), 
   rayleighFlag(true), detectorPosition(0) 
{
  runAction = 0;
  XrayFluoDetector = XrayFluoDC;

  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  
  //create a messenger for this class
  gunMessenger = new XrayFluoPrimaryGeneratorMessenger(this); 

  // default particle kinematic
  G4ParticleDefinition* particle
    = G4Gamma::Definition();
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  particleGun->SetParticleEnergy(10. * keV);

  G4double position = -0.5*(XrayFluoDetector->GetWorldSizeZ());
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,position));

  G4cout << "XrayFluoPrimaryGeneratorAction created" << G4endl;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoPrimaryGeneratorAction::ActivatePhaseSpace(G4String fileName) {

  // load phase-space
  phaseSpaceGunFlag = true;

  // reads the data stored on disk form previous runs 
  // and get these data to data members

  XrayFluoAnalysisManager* analysis =  XrayFluoAnalysisManager::getInstance();
  analysis->LoadGunData(fileName, rayleighFlag);
  detectorPosition = XrayFluoDetector->GetDetectorPosition();
  detectorPosition.setR(detectorPosition.r()-(5.*cm)); // 5 cm before the detector, so in front of it.

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoPrimaryGeneratorAction::SetRayleighFlag (G4bool value)
{
  rayleighFlag = value; 
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoPrimaryGeneratorAction::~XrayFluoPrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //retrieve runAction, if not done
  if (!runAction)
    {
      //Sequential runaction
      if (G4RunManager::GetRunManager()->GetRunManagerType() == 
	  G4RunManager::sequentialRM)
	runAction = static_cast<const XrayFluoRunAction*>
	  (G4RunManager::GetRunManager()->GetUserRunAction());  
      else //MT master runaction
	runAction = static_cast<const XrayFluoRunAction*>
	  (G4MTRunManager::GetMasterRunManager()->GetUserRunAction());  
      if (!runAction)
	G4cout << "Something wrong here!" << G4endl;
    }
 
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
	  G4DataVector* energies =  runAction->GetEnergies();
	  G4DataVector* data =  runAction->GetData();
	 
	  G4double sum = runAction->GetDataSum();
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
	  const XrayFluoDataSet* dataSet = runAction->GetGammaSet();
	  
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

    particleGun->SetParticlePosition(detectorPosition); 
    particleGun->SetParticleMomentumDirection(detectorPosition);
   
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
   
    const std::pair<G4double,G4String> kine = 
      XrayFluoAnalysisManager::getInstance()->GetEmittedParticleEnergyAndType();

    G4double energy = kine.first;
    G4ParticleDefinition* particle = particleTable->FindParticle(kine.second);

    particleGun->SetParticleEnergy(energy);
    particleGun->SetParticleDefinition(particle);


  }

  G4double partEnergy = particleGun->GetParticleEnergy();
  XrayFluoAnalysisManager* analysis =  XrayFluoAnalysisManager::getInstance();
  analysis->analysePrimaryGenerator(partEnergy/keV);

 
  particleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....









