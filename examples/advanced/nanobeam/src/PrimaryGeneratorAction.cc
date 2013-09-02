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
// Please cite the following paper if you use this software
// Nucl.Instrum.Meth.B260:20-27, 2007

#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* DC)
:fDetector(DC)
{
  fAngleMax = 0.09;
  
  // Default
  fEmission = 1;

  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
  
  G4ParticleDefinition* particle =
    G4ParticleTable::GetParticleTable()->FindParticle("proton");
      fParticleGun->SetParticleDefinition(particle);
    
  fParticleGun->SetParticleEnergy(3*MeV);

  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  
  fGunMessenger = new PrimaryGeneratorMessenger(this);  
  
  // Matrix
  
  fBeamMatrix = CLHEP::HepMatrix(32,32);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fGunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4int numEvent;
  
  numEvent=anEvent->GetEventID()+1;
  
  G4double x0,y0,z0,theta,phi,xMom0,yMom0,zMom0,e0,de;
  
  fShoot=false;

  x0=0; y0=0; z0=0; theta=0; phi=0; e0=0;

// Coefficient computation

if (fEmission==1)
{
	fDetector->SetCoef(1);
	fShoot=true;
	de = 0;

	if (numEvent==	1)  {theta =	-0.00500; phi =	-0.00500; de = -0.0050; }
	if (numEvent==	2)  {theta =	-0.005/3; phi =	-0.00500; de = -0.0050; }
	if (numEvent==	3)  {theta =	 0.005/3; phi =	-0.00500; de = -0.0050; }
	if (numEvent==	4)  {theta =	 0.00500; phi =	-0.00500; de = -0.0050; }
	if (numEvent==	5)  {theta =	-0.00500; phi =	-0.005/3; de = -0.0050; }
	if (numEvent==	6)  {theta =	-0.005/3; phi =	-0.005/3; de = -0.0050; }
	if (numEvent==	7)  {theta =	 0.005/3; phi =	-0.005/3; de = -0.0050; }
	if (numEvent==	8)  {theta =	 0.00500; phi =	-0.005/3; de = -0.0050; }
	if (numEvent==	9)  {theta =	-0.00500; phi =	 0.005/3; de = -0.0050; }
	if (numEvent==	10) {theta =	-0.005/3; phi =	 0.005/3; de = -0.0050; }
	if (numEvent==	11) {theta =	 0.005/3; phi =	 0.005/3; de = -0.0050; }
	if (numEvent==	12) {theta =	 0.00500; phi =	 0.005/3; de = -0.0050; }
	if (numEvent==	13) {theta =	-0.00500; phi =	 0.00500; de = -0.0050; }
	if (numEvent==	14) {theta =	-0.005/3; phi =	 0.00500; de = -0.0050; }
	if (numEvent==	15) {theta =	 0.005/3; phi =	 0.00500; de = -0.0050; }
	if (numEvent==	16) {theta =	 0.00500; phi =	 0.00500; de = -0.0050; }
	if (numEvent==	17) {theta =	-0.00500; phi =	-0.00500; de =  0.0050; }
	if (numEvent==	18) {theta =	-0.005/3; phi =	-0.00500; de =  0.0050; }
	if (numEvent==	19) {theta =	 0.005/3; phi =	-0.00500; de =  0.0050; }
	if (numEvent==	20) {theta =	 0.00500; phi =	-0.00500; de =  0.0050; }
	if (numEvent==	21) {theta =	-0.00500; phi =	-0.005/3; de =  0.0050; }
	if (numEvent==	22) {theta =	-0.005/3; phi =	-0.005/3; de =  0.0050; }
	if (numEvent==	23) {theta =	 0.005/3; phi =	-0.005/3; de =  0.0050; }
	if (numEvent==	24) {theta =	 0.00500; phi =	-0.005/3; de =  0.0050; }
	if (numEvent==	25) {theta =	-0.00500; phi =	 0.005/3; de =  0.0050; }
	if (numEvent==	26) {theta =	-0.005/3; phi =	 0.005/3; de =  0.0050; }
	if (numEvent==	27) {theta =	 0.005/3; phi =	 0.005/3; de =  0.0050; }
	if (numEvent==	28) {theta =	 0.00500; phi =	 0.005/3; de =  0.0050; }
	if (numEvent==	29) {theta =	-0.00500; phi =	 0.00500; de =  0.0050; }
	if (numEvent==	30) {theta =	-0.005/3; phi =	 0.00500; de =  0.0050; }
	if (numEvent==	31) {theta =	 0.005/3; phi =	 0.00500; de =  0.0050; }
	if (numEvent==	32) {theta =	 0.00500; phi =	 0.00500; de =  0.0050; }

	theta=theta*200*fAngleMax*1e-3;
	phi=phi*200*fAngleMax*1e-3;

	de=de/100;
	e0 = 3*MeV + 2*de*3*MeV;

	x0 = 0;
	y0 = 0;
	z0 = -8870*mm;
	
	// triplet only
	//z0 = -3230*mm;
	
        G4float cte = 1 ;
        G4float DE = 100*de;

        phi = phi * 1000;
        theta = theta * 1000;

        // Matrix filling up
       
        fBeamMatrix(numEvent,1) = cte;
      
        fBeamMatrix(numEvent,2) = theta;
        fBeamMatrix(numEvent,3) = phi;
        fBeamMatrix(numEvent,4) = DE;

        fBeamMatrix(numEvent,5) = theta * theta;	
        fBeamMatrix(numEvent,6) = theta * phi;	
        fBeamMatrix(numEvent,7) = phi * phi;	
        fBeamMatrix(numEvent,8) = theta * DE;
        fBeamMatrix(numEvent,9) = phi * DE;

        fBeamMatrix(numEvent,10) = theta * theta * theta;	
        fBeamMatrix(numEvent,11) = theta * theta * phi;	
        fBeamMatrix(numEvent,12) = theta * phi * phi;	
        fBeamMatrix(numEvent,13) = phi *  phi *  phi;
        fBeamMatrix(numEvent,14) = theta * theta * DE;
        fBeamMatrix(numEvent,15) = theta * phi * DE;		
        fBeamMatrix(numEvent,16) = phi *  phi * DE;

        fBeamMatrix(numEvent,17) = theta * theta * theta * phi;	
        fBeamMatrix(numEvent,18) = theta * theta * phi * phi;	
        fBeamMatrix(numEvent,19) = theta * phi * phi * phi;	
        fBeamMatrix(numEvent,20) = theta * theta * theta * DE;
        fBeamMatrix(numEvent,21) = theta * theta * phi *  DE;
        fBeamMatrix(numEvent,22) = theta * phi * phi * DE;
        fBeamMatrix(numEvent,23) = phi * phi * phi * DE;

        fBeamMatrix(numEvent,24) = theta * theta * theta * phi * phi;	
        fBeamMatrix(numEvent,25) = theta * theta * phi * phi * phi;	
        fBeamMatrix(numEvent,26) = theta * theta * theta * phi * DE;	
        fBeamMatrix(numEvent,27) = theta * theta * phi * phi * DE;	
        fBeamMatrix(numEvent,28) = theta * phi * phi * phi * DE;

        fBeamMatrix(numEvent,29) = theta * theta * theta * phi * phi * phi;	
        fBeamMatrix(numEvent,30) = theta * theta * theta * phi * phi * DE;
        fBeamMatrix(numEvent,31) = theta * theta * phi * phi * phi * DE;	

        fBeamMatrix(numEvent,32) = theta * theta * theta * phi * phi * phi * DE;	

       //	            
		            
        phi = phi / 1000;
        theta = theta / 1000;

        /*
        G4cout << fDetector->GetG1() << G4endl;
        G4cout << fDetector->GetG2() << G4endl;
        G4cout << fDetector->GetG3() << G4endl;
        G4cout << fDetector->GetG4() << G4endl;
        G4cout << fDetector->GetModel() << G4endl;
        G4cout << fDetector->GetCoef() << G4endl;
        G4cout << fDetector->GetGrid() << G4endl;
        */

} // end coefficient

// Full beam

if (fEmission==2)
{
	fDetector->SetCoef(0);
	fShoot=false;
	
	G4double aR, angle, rR;
	aR = -1;
	
	e0= G4RandGauss::shoot(3*MeV,5.0955e-5*MeV);  // AIFIRA ENERGY RESOLUTION
	
	while (aR < 0) aR = G4RandGauss::shoot(0.10e-3 , 0.06e-3/2.35) * rad; // old =0.08e-3 displacement
	
	angle = G4UniformRand() * 2 * CLHEP::pi *rad;
	
	theta = aR * std::cos(angle);
	
	phi = aR * std::sin(angle);
	
	rR = XYofAngle(aR);
	
	x0 = rR*std::cos(angle);
	
	y0 = rR*std::sin(angle);
	
	x0 = G4RandGauss::shoot(x0,220/2.35) *micrometer;
	
	theta = G4RandGauss::shoot(theta,0.03e-3/2.35);
	
	y0 = G4RandGauss::shoot(y0,220/2.35) *micrometer;
	
	phi = G4RandGauss::shoot(phi,0.03e-3/2.35);
	
	z0 = (-9120+250)*mm;  //position C0
	
	if (std::sqrt(x0*x0+y0*y0)/micrometer<2000) fShoot=true;

	/*
	xMom0 = std::sin(theta);
  	yMom0 = std::sin(phi);
   	zMom0 = std::cos(theta);  
	*/
	
} // end full beam

// shoot

if (fShoot)
{
  
  G4cout << "-> Event= " << numEvent<< ": X0(mm)= " << x0/mm << " X0(mm) = " << y0/mm << " Z0(m) = " << z0/m  
         << " THETA(mrad)= " << theta/mrad << " PHI(mrad)= " <<	phi/mrad << " E0(MeV)= " << e0/MeV << G4endl;
  
  xMom0 = std::sin(theta);
  
  yMom0 = std::sin(phi);
  
  zMom0 = std::sqrt(1.-xMom0*xMom0-yMom0*yMom0);  

  fParticleGun->SetParticleEnergy(e0);

  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(xMom0,yMom0,zMom0));

  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  
  G4ParticleDefinition* particle=
    G4ParticleTable::GetParticleTable()->FindParticle("proton");
  
  fParticleGun->SetParticleDefinition(particle);
  
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

// end shoot

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PrimaryGeneratorAction::SetEmission(G4int value)
{
  fEmission = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double PrimaryGeneratorAction::XYofAngle(G4double angle)//  returns position in micrometers
{
  return std::pow(20000*angle, 13);
}

