///////////////////////////////////////////////////////////////////////////////
// File: CMSPrimaryGeneratorAction.cc
// Author: I. Gonzalez
// Last modification: 11/98 I.G.
//		      06/08/99 V.L.
//                    08/09/99 I.G. -> Add gunMessenger. Pythia file clean up.
//                    18/04/00 P.A., S.B. Added functionality
//                    10/01 P.Arce use COBRA GeneratorInterface
///////////////////////////////////////////////////////////////////////////////

#include "CMSPrimaryGeneratorAction.hh"
#include "CMSPrimaryGeneratorMessenger.hh"
#include "G4HEPEvtInterface.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "CLHEP/Random/RandFlat.h"
#include "G4HEPEvtInterface.hh"
#include "G4RunManager.hh"
#include "SystemOfUnits.h"

#define debug

CMSPrimaryGeneratorAction::CMSPrimaryGeneratorAction(): particleGun(0),
  generatorInput(singleFixed), isInitialized(0),
  verboseLevel(0), scanSteps(0), n_particle(1), particleName("pi-"),
  particleEnergy(100*GeV), particlePosition(0.,0.,0.), particleDir(1.,1.,0.1) {
  
  //Initialise the messenger
  gunMessenger = new CMSPrimaryGeneratorMessenger(this);
    
  //Default settings:
  SetMinimumEnergy(1.*GeV);
  SetMaximumEnergy(1.*TeV);
  SetMinimumEta(0.);
  SetMaximumEta(3.5);
  SetMinimumPhi(0.*degree);
  SetMaximumPhi(360.*degree);
  SetStepsPhi(1);
  SetStepsEta(1);
    
  // Number of particles per gun
  particleGun = new G4ParticleGun(n_particle);
    
  // Getting particle definition
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName);
  particleGun->SetParticleDefinition(particle);
    
  // Setting particle properties
  particleGun->SetParticleEnergy(particleEnergy);
  particleGun->SetParticleMomentumDirection(particleDir);
  particleGun->SetParticlePosition(particlePosition);
  print(0);
}

CMSPrimaryGeneratorAction::~CMSPrimaryGeneratorAction() {
  if (gunMessenger)
    delete gunMessenger;
  if (particleGun)
    delete particleGun;
}

void CMSPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

  if (isInitialized == 0) initialize();

  if (generatorInput == singleRandom) {
    particleEnergy = RandFlat::shoot(energyMin,energyMax);
    particleGun->SetParticleEnergy(particleEnergy);

    G4double eta = RandFlat::shoot(etaMin,etaMax);
    G4double phi = RandFlat::shoot(phiMin,phiMax);
    G4double theta = atan(exp(-eta))*2.;
    G4double randomX = sin(theta)*cos(phi);
    G4double randomY = sin(theta)*sin(phi);
    G4double randomZ = cos(theta);
  
    particleDir = G4ThreeVector(randomX,randomY,randomZ);
    particleGun->SetParticleMomentumDirection(particleDir);
#ifdef debug
    if (verboseLevel >= 2 ) {
      G4cout << "Energy " << particleEnergy/GeV << " GeV; Theta " 
	     << theta/deg << " degree; Phi " << phi/deg << " degree" << endl;
      G4cout << "Shooting in " << particleDir << " direction "<< endl;
    }
#endif
  } else if (generatorInput == singleScan) {
    G4double scanEtaStep, scanPhiStep;
    if (scanSteps == 0) { 
      scanPhiStep = (phiMax - phiMin) / phiSteps;
      phiValue = phiMin - scanPhiStep; //first step is going to change scanCurrentPhiValue
      etaValue = etaMin;
    }

    scanEtaStep = (etaMax - etaMin) / etaSteps;
    scanPhiStep = (phiMax - phiMin) / phiSteps;
#ifdef debug
    if (verboseLevel > 2 ) {
      G4cout << " scanEtaStep " << scanEtaStep << " # of Steps " << etaSteps << endl;
      G4cout << " scanPhiStep " << scanPhiStep << " # of Steps " << phiSteps << endl;
    }
#endif

    //----- First scan in phi, then in eta
    if (phiMax - phiValue < 1.E-6 * scanPhiStep) { // !only <= 1.E6 steps allowed
      if (etaMax - etaValue < 1.E-6 * scanEtaStep) { // !only <= 1.E6 steps allowed
	G4cout << " Scan completed!" << endl;
	return;
      } else {
	etaValue += scanEtaStep; 
	phiValue  = phiMin;
      }
    } else {
      phiValue += scanPhiStep;
    }    
    G4double theta = atan(exp(-etaValue))*2.;

    G4double scanX = sin(theta)*cos(phiValue);
    G4double scanY = sin(theta)*sin(phiValue);
    G4double scanZ = cos(theta);
    if (verboseLevel >= 2 ) {
      G4cout << "Scan eta " << etaValue << " Phi " << phiValue/deg
	     << " theta " << theta/deg << endl;
    }
    particleDir = G4ThreeVector(scanX,scanY,scanZ);
    particleGun->SetParticleMomentumDirection(particleDir);
#ifdef debug
    if (verboseLevel > 2 ) {
      G4cout  << "Shooting in " << particleDir << " direction "<< endl;
    }
#endif
    scanSteps++;
  }
  
  // Generate GEANT4 primary vertex
  particleGun->GeneratePrimaryVertex(anEvent);
} 


void CMSPrimaryGeneratorAction::SetVerboseLevel(G4int val){
  verboseLevel = val;
}


void CMSPrimaryGeneratorAction::SetRandom(G4String val) { 

  if (val=="on") {
    generatorInput = singleRandom;
    print (1);
  } else {
    generatorInput = singleFixed;
    print (1);
  }
}


void CMSPrimaryGeneratorAction::SetScan(G4String val) { 

  if (val=="on") {
    generatorInput = singleScan;
    scanSteps = 0;
    print (1);
  } else {
    generatorInput = singleFixed;
    print (1);
  }  
}


void CMSPrimaryGeneratorAction::SetMinimumEnergy(G4double p){

  if (p <= 0.) {
    G4cerr<<"CMSPrimaryGeneratorAction::SetMinimumEnergy: value " << p/GeV << "GeV is out of bounds, it will not be used"<<endl;
    G4cerr<<" Should be  >0.  Please check"<<endl; 
  } else {
    energyMin = p;
#ifdef debug
    if (verboseLevel >= 1 ) {
      G4cout << " CMSPrimaryGeneratorAction: setting min. value of energy to "
	     << energyMin/GeV << " GeV " << endl;
    }
#endif
  }
}


void CMSPrimaryGeneratorAction::SetMaximumEnergy(G4double p){

  if (p <= 0.) {
    G4cerr<<"CMSPrimaryGeneratorAction::SetMaximumEnergy: value " << p/GeV << "GeV is out of bounds, it will not be used"<<endl;
    G4cerr<<" Should be  >0.  Please check"<<endl; 
  } else {
    energyMax = p;
#ifdef debug
    if (verboseLevel >= 1 ) {
      G4cout << " CMSPrimaryGeneratorAction: setting max. value of energy to "
	     << energyMax/GeV << " GeV " << endl;
    }
#endif
  }
}


void CMSPrimaryGeneratorAction::SetMinimumPhi(G4double p){

  if (fabs(p)>2.*pi) {
    G4cerr<<"CMSPrimaryGeneratorAction::SetMinimumPhi: setting value quite large"<<endl;
    G4cerr<<" Should be given in radians - Please check"<<endl;
  } else {
    phiMin = fabs(p);
#ifdef debug
    if (verboseLevel >= 1 ) {
      G4cout << " CMSPrimaryGeneratorAction: setting min. value of phi to "
	     << phiMin << endl;
    }
#endif
  }
}


void CMSPrimaryGeneratorAction::SetMaximumPhi(G4double p){

  if (fabs(p)>2.*pi) {
    G4cerr<<"CMSPrimaryGeneratorAction::SetMaximumPhi: setting value quite large"<<endl;
    G4cerr<<" Should be given in radians - Please check"<<endl;
  } else {
    phiMax = fabs(p);
#ifdef debug
    if (verboseLevel >= 1 ) {
      G4cout << " CMSPrimaryGeneratorAction: setting max. value of phi to "
	     << phiMax << endl;
    }
#endif
  }
}


void CMSPrimaryGeneratorAction::SetStepsPhi(G4int val){

  if (val <= 0) {
    G4cerr<<"CMSPrimaryGeneratorAction::SetStepsPhi: value " << val << " is out of bounds, it will not be used"<<endl;
    G4cerr<<" Should be  > 0  Please check"<<endl; 
  } else {
    phiSteps = val;
#ifdef debug
    if (verboseLevel >= 1 ) {
      G4cout << " CMSPrimaryGeneratorAction: setting no. of steps in phi to "
	     << phiSteps << endl;
    }
#endif
  }
}


void CMSPrimaryGeneratorAction::SetMinimumEta(G4double p){

  etaMin = p;
#ifdef debug
  if (verboseLevel >= 1 ) {
    G4cout << " CMSPrimaryGeneratorAction: setting min. value of eta to "
	   << etaMin << endl;
  }
#endif
}


void CMSPrimaryGeneratorAction::SetMaximumEta(G4double p){

  etaMax = p;
#ifdef debug
  if (verboseLevel >= 1 ) {
    G4cout << " CMSPrimaryGeneratorAction: setting max. value of eta to "
	   << etaMax << endl;
  }
#endif
}


void CMSPrimaryGeneratorAction::SetStepsEta(G4int val){

  if (val <= 0) {
    G4cerr<<"CMSPrimaryGeneratorAction::SetStepsEta: value " << val << " is out of bounds, it will not be used"<<endl;
    G4cerr<<" Should be  > 0  Please check"<<endl; 
  } else {
    etaSteps = val;
#ifdef debug
    if (verboseLevel >= 1 ) {
      G4cout << " CMSPrimaryGeneratorAction: setting no. of steps in eta to "
	     << etaSteps << endl;
    }
#endif
  }
}

void CMSPrimaryGeneratorAction::SetGunPosition(const G4ThreeVector & pos) const {

  particleGun->SetParticlePosition(pos);
}

void CMSPrimaryGeneratorAction::SetRunNo(G4int val){
  G4RunManager::GetRunManager()->SetRunIDCounter( val );
}

void CMSPrimaryGeneratorAction::initialize(){

  isInitialized = 1;

  print (1);
}


void CMSPrimaryGeneratorAction::print(G4int val){

#ifdef debug
  if (verboseLevel >= val) {

    if (generatorInput == singleRandom) {
      G4cout << endl
      	     << "**********************************************************************" << endl
	     << "*                                                                    *" << endl  
      	     << "* CMSPrimaryGeneratorAction DEFAULT Random Energy/Direction settings:*" << endl
	     << "*                                                                    *" << endl  
	     << "*                                                                    *" << endl  
	     << "*   Energy in    [ "<< energyMin/GeV   << " - " << energyMax/GeV   << "] (GeV) "<< endl
	     << "*   Phi angle in [ "<< phiMin          << " - " << phiMax          << "] (rad) "<< endl
	     << "*                [ "<< phiMin/degree   << " - " << phiMax/degree   << "] (deg) "<< endl 
	     << "*   Eta in       [ "<< etaMin          << " - " << etaMax          << "]       "<< endl
	     << "*                                                                    *" << endl  
	     << "*                                                                    *" << endl  
      	     << "**********************************************************************" << endl;
    } else if (generatorInput == singleScan) {
      G4cout << endl
	     << "**********************************************************************" << endl
	     << "*                                                                    *" << endl  
      	     << "* CMSPrimaryGeneratorAction DEFAULT Scan Direction settings :        *" << endl
	     << "*                                                                    *" << endl  
	     << "*                                                                    *" << endl  
	     << "*   Phi angle in [ " << phiMin/degree   << " - " << phiMax/degree << "] (deg) " << endl
	     << "*   Eta in       [ " << etaMin          << " - " << etaMax        << "]       " << endl
	     << "*   Steps along eta " << etaSteps << " and along phi " << phiSteps << endl
	     << "*                                                                    *" << endl  
	     << "*                                                                    *" << endl  
      	     << "**********************************************************************" << endl;
    } else  if (generatorInput == singleFixed) {
      G4cout << endl
	     << "*******************************************************************" << endl
	     << "*                                                                 *" << endl  
	     << "* CMSPrimaryGeneratorAction: Current settings :                   *" << endl
	     << "*                                                                 *" << endl  
	     << "* " << particleGun->GetNumberOfParticles() 
	     << "  " << particleGun->GetParticleDefinition()->GetParticleName() 
	     << " of " << particleGun->GetParticleEnergy()/GeV << " GeV" << endl
	     << "* will be shot from " << particleGun->GetParticlePosition() << endl;
      G4cout << "* in direction " << particleGun->GetParticleMomentumDirection() << endl;
      G4cout  << "*                                                                 *" << endl  
	      << "*                                                                 *" << endl  
	      << "*******************************************************************" << endl;
    }
  }
#endif
 
}

