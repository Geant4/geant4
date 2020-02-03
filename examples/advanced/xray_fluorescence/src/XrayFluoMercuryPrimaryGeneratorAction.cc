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

#include "XrayFluoMercuryPrimaryGeneratorAction.hh"
#include "XrayFluoMercuryDetectorConstruction.hh"
#include "XrayFluoMercuryPrimaryGeneratorMessenger.hh"
#include "XrayFluoRunAction.hh"
#include "XrayFluoAnalysisManager.hh"
#include "XrayFluoDataSet.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DataVector.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoMercuryPrimaryGeneratorAction::XrayFluoMercuryPrimaryGeneratorAction(const XrayFluoMercuryDetectorConstruction* XrayFluoDC)
  :globalFlag(false),spectrum("off")
{

  XrayFluoDetector = XrayFluoDC;

  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  
  //create a messenger for this class
  gunMessenger = new XrayFluoMercuryPrimaryGeneratorMessenger(this);
  runManager = new XrayFluoRunAction();
  
  // default particle kinematic
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="gamma");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
  

  particleGun->SetParticleEnergy(10.*keV);
  G4double position = -0.5*(XrayFluoDetector->GetWorldSizeZ());
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,position));

  G4cout << "XrayFluoMercuryPrimaryGeneratorAction created" << G4endl;
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoMercuryPrimaryGeneratorAction::~XrayFluoMercuryPrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
  delete runManager;

  G4cout << "XrayFluoMercuryPrimaryGeneratorAction deleted" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoMercuryPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  // 

  // Conidering the sunas a Poin-like source.

  G4double z0 = -0.5*(XrayFluoDetector->GetWorldSizeZ());
  G4double y0 = 0.*m, x0 = 0.*m;


  // Let's try to illuminate only the prtion of Mercury surface that can be seen by the detector.

  G4double spacecraftLatitude = XrayFluoDetector->GetOrbitInclination();
  G4double mercuryDia = XrayFluoDetector->GetMercuryDia();
  G4double sunDia = XrayFluoDetector->GetSunDia();
  G4double opticField = XrayFluoDetector->GetOpticAperture();
  
 
  G4double a = 2*std::tan(opticField/2);
  
  //  if (!pointLikeFlag) {

    // let's decide from wich point of the sun surface the particle is coming:

  G4double theta = std::acos(2.*G4UniformRand() - 1.0);
    G4double phi = 2. * pi * G4UniformRand();
    G4double rho = sunDia/2;                         

    G4double sunPosX = x0 + rho * std::sin(theta) * std::cos(phi);
    G4double sunPosY = y0 + rho * std::sin(theta) * std::sin(phi);
    G4double sunPosZ = z0 + rho * std::cos(theta);           

    particleGun->SetParticlePosition(G4ThreeVector(sunPosX,sunPosY,sunPosZ));

    // the angle at the center of Mercury subtending the area seen by the optics:
    G4double alpha = 2 * a/mercuryDia; 
    
    if(!globalFlag){
      theta = alpha * G4UniformRand() + (180.*deg - spacecraftLatitude)-alpha/2.;
      phi = alpha * G4UniformRand() + 90. * deg - alpha/2.;     
    }    
    
    else if(globalFlag){
      theta = pi/2. * rad * G4UniformRand() + 90.*deg ; //was 900., probably an error
      phi = 2*pi*rad * G4UniformRand() ;
    }
    
    rho = mercuryDia/2.;                        
    
    G4double mercuryPosX = rho * std::sin(theta) * std::cos(phi);
    G4double mercuryPosY = rho * std::sin(theta) * std::sin(phi);
    G4double mercuryPosZ = rho * std::cos(theta);           
    
    particleGun->SetParticleMomentumDirection(
			    G4ThreeVector(mercuryPosX-sunPosX ,mercuryPosY-sunPosY,mercuryPosZ-sunPosZ));

    //  }
//   if (pointLikeFlag) {

//   // theta is the angle that the mean direction of the incident light (on the desired 
//   // point of the surface of Mercury) makes with  the Z-axis
//   G4double theta = std::asin( mercuryDia/2. * std::sin(spacecraftLatitude) / 
// 			 std::sqrt(std::pow(z0,2)+std::pow(mercuryDia/2.,2)-2*mercuryDia/2.*z0*std::cos(spacecraftLatitude)) );	 

//   // on the y axis, the light emitted from the Sun must be in [theta-phi;theta+phi]
//   G4double phi = std::asin( mercuryDia/2.*std::sin(spacecraftLatitude) + a*std::cos(spacecraftLatitude) /
// 		       std::sqrt( std::pow(mercuryDia/2.*std::sin(spacecraftLatitude) + a*std::cos(spacecraftLatitude) , 2) +
// 			     std::pow(z0 - mercuryDia/2.*std::cos(spacecraftLatitude) - a*std::sin(spacecraftLatitude) , 2)) ) 
//     - theta;  
  
//   // on the x axis, the light emitted from the Sun must be in [-zeta;zeta]
//   G4double zeta = std::atan( a/std::sqrt(std::pow(z0,2)+std::pow(mercuryDia,2)-2*mercuryDia*z0*std::cos(spacecraftLatitude)) );
  
  
  
//   //alpha in [-zeta;zeta]
//   G4double alpha = (2*zeta)*G4UniformRand() - zeta;
//   //beta in [theta-phi;theta+phi]
//   G4double beta = (G4UniformRand()*2*phi) - phi + theta;
  
//   G4double dirY = std::sin(beta);
//   G4double dirX = std::sin(alpha);
  
//   particleGun->SetParticleMomentumDirection(G4ThreeVector(dirX.,dirY,1.));
  
//   particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

//   }


  
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
  

#ifdef G4ANALYSIS_USE 

  G4double partEnergy = particleGun->GetParticleEnergy();
  XrayFluoAnalysisManager* analysis =  XrayFluoAnalysisManager::getInstance();
  analysis->analysePrimaryGenerator(partEnergy/keV);

#endif

  particleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....









