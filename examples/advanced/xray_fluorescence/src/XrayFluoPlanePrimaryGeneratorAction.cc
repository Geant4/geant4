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
// Author: Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// History:
// -----------
// 02 Sep 2003 Alfonso Mantero created
//
// -------------------------------------------------------------------

#include "XrayFluoPlanePrimaryGeneratorAction.hh"
#include "XrayFluoPlaneDetectorConstruction.hh"
#include "XrayFluoPlanePrimaryGeneratorMessenger.hh"
#include "XrayFluoRunAction.hh"
#include "XrayFluoAnalysisManager.hh"
#include "XrayFluoDataSet.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4DataVector.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoPlanePrimaryGeneratorAction::XrayFluoPlanePrimaryGeneratorAction(const XrayFluoPlaneDetectorConstruction* XrayFluoDC)
  :rndmFlag("on"),beam("off"),spectrum("off"),isoVert("off")
{

  XrayFluoDetector = XrayFluoDC;

  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  
  //create a messenger for this class
  gunMessenger = new XrayFluoPlanePrimaryGeneratorMessenger(this);
  runManager = new XrayFluoRunAction();
  
  // default particle kinematic
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="gamma");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  

  particleGun->SetParticleEnergy(10.*keV);
  G4double position = -0.5*(XrayFluoDetector->GetWorldSizeZ());
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,position));

  G4cout << "XrayFluoPlanePrimaryGeneratorAction created  UUUUUUUUUUAAAAAAAAAAAAAAAAAAAAAAAaa" << G4endl;
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoPlanePrimaryGeneratorAction::~XrayFluoPlanePrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
  delete runManager;

  G4cout << "XrayFluoPlanePrimaryGeneratorAction deleted" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoPlanePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  // 
  G4double z0 = -0.5*(XrayFluoDetector->GetWorldSizeZ());
  G4double y0 = 0.*m, x0 = 0.*m;
  G4double dX = 0.5*(XrayFluoDetector->GetWorldSizeXY())-(XrayFluoDetector->GetPlaneSizeXY());
  if (rndmFlag == "on")

    {y0 = (XrayFluoDetector->GetPlaneSizeXY())*(G4UniformRand()-0.5); 
    x0 = (XrayFluoDetector->GetPlaneSizeXY())*(G4UniformRand()-0.5) + dX; 
    } 

  z0 = -1 * dX;

  particleGun->SetParticleMomentumDirection(G4ThreeVector(-1.,0.,1.));

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
      
      G4double Xdim = XrayFluoDetector->GetPlaneSizeXY();
      G4double Ydim = XrayFluoDetector->GetPlaneSizeXY();
      
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









