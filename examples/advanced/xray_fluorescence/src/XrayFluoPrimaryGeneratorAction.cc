//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: XrayFluoPrimaryGeneratorAction.cc
// GEANT4 tag $Name: xray_fluo-V04-01-03
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
//
// -------------------------------------------------------------------

#include "XrayFluoPrimaryGeneratorAction.hh"
#include "G4DataVector.hh"
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoPrimaryGeneratorAction::XrayFluoPrimaryGeneratorAction(
							       XrayFluoDetectorConstruction* XrayFluoDC)
  :XrayFluoDetector(XrayFluoDC),rndmFlag("off"),
   beam("off"),spectrum("off"),isoVert("off")
{
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
  
  particleGun->SetParticleEnergy(6.*keV);
  G4double position = -0.5*(XrayFluoDetector->GetWorldSizeZ());
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,position));
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoPrimaryGeneratorAction::~XrayFluoPrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  // 
  G4double z0 = -0.5*(XrayFluoDetector->GetWorldSizeZ());
  G4double y0 = 0.*cm, x0 = 0.*cm;
  if (rndmFlag == "on")
    {y0 = (XrayFluoDetector->GetSampleSizeXY())*(G4UniformRand()-0.5);
    x0 = (XrayFluoDetector->GetSampleSizeXY())*(G4UniformRand()-0.5);
    } 
  particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  
  //randomize starting point
  if (beam == "on")
    {
      G4double radius = 0.5 * mm;
      G4double rho = radius*sqrt(G4UniformRand());
      G4double theta = 2*pi*G4UniformRand()*rad;
      G4double position = -0.5*(XrayFluoDetector->GetWorldSizeZ());
      
      G4double y = rho * sin(theta);
      G4double x = rho * cos(theta);
      
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
	      
	      G4double randomNum = G4UniformRand()*5.0E6;
	      
	      particleEnergy = random *  energyRange;
	      
	      if ((dataSet->FindValue(particleEnergy,id)) > randomNum)
		{
		  i = 1;
		  
		}
	    }
	   particleGun->SetParticleEnergy(particleEnergy);
	}
    }
  
  if (isoVert == "on")
    {
      G4double rho = 1. *m;
      //theta in [0;pi/2]
      G4double theta = (pi/2)*G4UniformRand();
      //phi in [-pi;pi]
      G4double phi = (G4UniformRand()*2*pi)- pi;
      G4double x = rho*sin(theta)*sin(phi);
      G4double y = rho*sin(theta)*cos(phi);
      G4double z = -(rho*cos(theta));
      particleGun->SetParticlePosition(G4ThreeVector(x,y,z));
      
      G4double Xdim = XrayFluoDetector->GetSampleSizeXY();
      G4double Ydim = XrayFluoDetector->GetSampleSizeXY();
      
      G4double Dx = Xdim*(G4UniformRand()-0.5);
      
      G4double Dy = Ydim*(G4UniformRand()-0.5);
      
      particleGun->SetParticleMomentumDirection(G4ThreeVector(-x+Dx,-y+Dy,-z));
      
    }
#ifdef G4ANALYSIS_USE 

  G4double partEnergy = particleGun->GetParticleEnergy();
  XrayFluoAnalysisManager* analysis =  XrayFluoAnalysisManager::getInstance();
  analysis->analysePrimaryGenerator(partEnergy/keV);

#endif
 
  particleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

