// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: FluoTestPrimaryGeneratorAction.cc,v 1.10 2001-11-15 13:04:12 guardi Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "FluoTestPrimaryGeneratorAction.hh"

#include "FluoTestDetectorConstruction.hh"
#include "FluoTestPrimaryGeneratorMessenger.hh"
#include "FluoTestRunAction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "FluoTestAnalysisManager.hh"
#include "FluoTestDataSet.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4ANALYSIS_USE
FluoTestPrimaryGeneratorAction::FluoTestPrimaryGeneratorAction( FluoTestDetectorConstruction* FluoTestDC,FluoTestAnalysisManager* analysisMgr)
  :FluoTestDetector(FluoTestDC),analysisManager(analysisMgr),rndmFlag("off"),
   rndmPart("off"),beam("off"),spectrum("off"),isoVert("off")
{ G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  
  //create a messenger for this class
  gunMessenger = new FluoTestPrimaryGeneratorMessenger(this);
  runManager = new FluoTestRunAction();

  // default particle kinematic

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
                    = particleTable->FindParticle(particleName="e-");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  particleGun->SetParticleEnergy(50.*MeV);
  G4double position = -0.5*(FluoTestDetector->GetWorldSizeZ());
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,position));
}
#else
FluoTestPrimaryGeneratorAction::FluoTestPrimaryGeneratorAction(
                                               FluoTestDetectorConstruction* FluoTestDC)
  :FluoTestDetector(FluoTestDC),rndmFlag("off"),
 rndmPart("off"),beam("off"),spectrum("off"),isoVert("off")
{
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  
  //create a messenger for this class
  gunMessenger = new FluoTestPrimaryGeneratorMessenger(this);
  runManager = new FluoTestRunAction();

  // default particle kinematic

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
                    = particleTable->FindParticle(particleName="e-");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
 
  particleGun->SetParticleEnergy(50.*MeV);
  G4double position = -0.5*(FluoTestDetector->GetWorldSizeZ());
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,position));

}
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestPrimaryGeneratorAction::~FluoTestPrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FluoTestPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  // 
  G4double z0 = -0.5*(FluoTestDetector->GetWorldSizeZ());
  G4double y0 = 0.*cm, x0 = 0.*cm;
  if (rndmFlag == "on")
     {y0 = (FluoTestDetector->GetSampleSizeXY())*(G4UniformRand()-0.5);
      x0 = (FluoTestDetector->GetSampleSizeXY())*(G4UniformRand()-0.5);
     } 
  particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

 //randomize particles

  G4double random = G4UniformRand();
  if (rndmPart == "on")
    {  
      G4String particleName;
      G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
      G4ParticleDefinition* particle;
      G4double electronAbundance = 0.5;
      G4double photonAbundance = 0.5;
      if ( random < electronAbundance)
	{
	  particle = particleTable->FindParticle(particleName="e-");
	}
      else{
	particle = particleTable->FindParticle(particleName="gamma");
      }
      particleGun->SetParticleDefinition(particle);
      
    }
  //randomize starting point
  if (beam == "on")
    {
      G4double radius = 0.5 * mm;
      G4double rho = radius*sqrt(G4UniformRand());
      G4double theta = 2*pi*G4UniformRand()*rad;
      G4double position = -0.5*(FluoTestDetector->GetWorldSizeZ());

      G4double y = rho * sin(theta);
      G4double x = rho * cos(theta);
      
      particleGun->SetParticlePosition(G4ThreeVector(x,y,position));
    }
//shoot particles according to a certain spectrum
  if (spectrum =="on")
    {
      const FluoTestDataSet* dataSet = runManager->GetSet();
      
      G4int i = 0;
      G4int id = 0;
      G4double minEnergy = 0. * keV;
      G4double particleEnergy= 0.;
      G4double maxEnergy = 10. * keV;
      G4double energyRange = maxEnergy - minEnergy;
      
    
      while ( i == 0)
	{
	  G4double random = G4UniformRand();
	  
	  G4double randomNum = G4UniformRand();
	  
	  particleEnergy = random *  energyRange;
	  
	  if ((dataSet->FindValue(particleEnergy,id)) > randomNum)
	    {
	      i = 1;
	      
	    }
	}
      /*      
	      #ifdef G4ANALYSIS_USE 
	      analysisManager->InsSpectrum(particleEnergy/keV);
	      #endif
      */
 particleGun->SetParticleEnergy(particleEnergy);

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
      G4cout<<"x = "<<x/m <<" m"<<G4endl;
      G4cout<<"y = "<<y/m <<" m"<<G4endl;
      G4cout<<"z = "<<z/m <<" m"<<G4endl;
      G4double Xdim = FluoTestDetector->GetSampleSizeXY();
      G4double Ydim = FluoTestDetector->GetSampleSizeXY();
     
      G4double Dx = Xdim*(G4UniformRand()-0.5);
      G4cout<<"Dx = "<<Dx/m<<" m"<<G4endl;
      G4double Dy = Ydim*(G4UniformRand()-0.5);
      G4cout<<"Dy = "<<Dy/m<<" m"<<G4endl;
      particleGun->SetParticleMomentumDirection(G4ThreeVector(-x+Dx,-y+Dy,-z));
      
    }
#ifdef G4ANALYSIS_USE 
  G4double partEnergy = particleGun->GetParticleEnergy();
  analysisManager->InsSpectrum(partEnergy/keV);
 
#endif
 
  particleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

