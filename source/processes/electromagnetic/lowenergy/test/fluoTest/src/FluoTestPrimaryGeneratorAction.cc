// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: FluoTestPrimaryGeneratorAction.cc,v 1.9 2001-10-31 12:33:46 elena Exp $
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
  :FluoTestDetector(FluoTestDC),analysisManager(analysisMgr),rndmFlag("off")
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
:FluoTestDetector(FluoTestDC),rndmFlag("off")
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
      const FluoTestDataSet* dataSet = runManager->GetSet();
     
      G4int i = 0;
      G4int id = 0;
      G4double minEnergy = 0. * keV;

      G4double maxEnergy = 10. * keV;
      G4double energyRange = maxEnergy - minEnergy;
     
      G4double randomEnergy  = 0;
      while ( i == 0)
	{
	  G4double random = G4UniformRand();
	 
	  G4double randomNum = G4UniformRand();
	 
	  randomEnergy = random *  energyRange;
	 
	  if ((dataSet->FindValue(randomEnergy,id)) > randomNum)
	    {
	      i = 1;
	     
	    }
	}
     
#ifdef G4ANALYSIS_USE 
      analysisManager->InsSpectrum(randomEnergy/keV);
#endif
     
      particleGun->SetParticleEnergy(randomEnergy);

    } 

  particleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

