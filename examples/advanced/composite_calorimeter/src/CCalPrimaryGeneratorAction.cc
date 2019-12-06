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
///////////////////////////////////////////////////////////////////////////////
// File: CCalPrimaryGeneratorAction.cc
// Description: CCalPrimaryGeneratorAction Sets up particle beam
///////////////////////////////////////////////////////////////////////////////

#include <CLHEP/Random/RandFlat.h>

#include "CCalPrimaryGeneratorAction.hh"
#include "CCalPrimaryGeneratorMessenger.hh"
#include "G4HEPEvtInterface.hh"

#include "G4Exp.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4HEPEvtInterface.hh"
#include "G4RunManager.hh"

//#define debug

CCalPrimaryGeneratorAction::CCalPrimaryGeneratorAction(): particleGun(0),
  generatorInput(singleFixed),  verboseLevel(0), n_particle(1), 
  particleName("pi-"), particleEnergy(100*GeV), particlePosition(0.,0.,0.),
  particleDir(1.,1.,0.1), isInitialized(0), scanSteps(0) {
  
  //Initialise the messenger
  gunMessenger = new CCalPrimaryGeneratorMessenger(this);
    
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

CCalPrimaryGeneratorAction::~CCalPrimaryGeneratorAction() {
  if (gunMessenger)
    delete gunMessenger;
  if (particleGun)
    delete particleGun;
}

void CCalPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

  if (isInitialized == 0) initialize();

  if (generatorInput == singleRandom) {
    particleEnergy = CLHEP::RandFlat::shoot(energyMin,energyMax);
    particleGun->SetParticleEnergy(particleEnergy);

    G4double eta = CLHEP::RandFlat::shoot(etaMin,etaMax);
    G4double phi = CLHEP::RandFlat::shoot(phiMin,phiMax);
    G4double theta = std::atan(G4Exp(-eta))*2.;
    G4double randomX = std::sin(theta)*std::cos(phi);
    G4double randomY = std::sin(theta)*std::sin(phi);
    G4double randomZ = std::cos(theta);
  
    particleDir = G4ThreeVector(randomX,randomY,randomZ);
    particleGun->SetParticleMomentumDirection(particleDir);
    if (verboseLevel >= 2 ) {
      G4cout << "Energy " << particleEnergy/GeV << " GeV; Theta " 
             << theta/deg << " degree; Phi " << phi/deg << " degree" << G4endl;
      G4cout << "Shooting in " << particleDir << " direction "<< G4endl;
    }
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
      G4cout << " scanEtaStep " << scanEtaStep << " # of Steps " << etaSteps 
             << G4endl;
      G4cout << " scanPhiStep " << scanPhiStep << " # of Steps " << phiSteps 
             << G4endl;
    }
#endif

    //----- First scan in phi, then in eta
    if (phiMax - phiValue < 1.E-6 * scanPhiStep) { // !only <= 1.E6 steps allowed
      if (etaMax - etaValue < 1.E-6 * scanEtaStep) { // !only <= 1.E6 steps allowed
        G4cout << " Scan completed!" << G4endl;
        return;
      } else {
        etaValue += scanEtaStep; 
        phiValue  = phiMin;
      }
    } else {
      phiValue += scanPhiStep;
    }    
    G4double theta = std::atan(G4Exp(-etaValue))*2.;

    G4double scanX = std::sin(theta)*std::cos(phiValue);
    G4double scanY = std::sin(theta)*std::sin(phiValue);
    G4double scanZ = std::cos(theta);
    if (verboseLevel >= 2 ) {
      G4cout << "Scan eta " << etaValue << " Phi " << phiValue/deg
             << " theta " << theta/deg << G4endl;
    }
    particleDir = G4ThreeVector(scanX,scanY,scanZ);
    particleGun->SetParticleMomentumDirection(particleDir);
#ifdef debug
    if (verboseLevel > 2 ) {
      G4cout  << "Shooting in " << particleDir << " direction "<< G4endl;
    }
#endif
    scanSteps++;
  }
  
  // Generate GEANT4 primary vertex
  particleGun->GeneratePrimaryVertex(anEvent);
} 


void CCalPrimaryGeneratorAction::SetVerboseLevel(G4int val){
  verboseLevel = val;
}


void CCalPrimaryGeneratorAction::SetRandom(G4String val) { 

  if (val=="on") {
    generatorInput = singleRandom;
    print (1);
  } else {
    generatorInput = singleFixed;
    print (1);
  }
}


void CCalPrimaryGeneratorAction::SetScan(G4String val) { 

  if (val=="on") {
    generatorInput = singleScan;
    scanSteps = 0;
    print (1);
  } else {
    generatorInput = singleFixed;
    print (1);
  }  
}


void CCalPrimaryGeneratorAction::SetMinimumEnergy(G4double p){

  if (p <= 0.) {
    G4cerr << "CCalPrimaryGeneratorAction::SetMinimumEnergy: value " << p/GeV 
           << "GeV is out of bounds, it will not be used" << G4endl;
    G4cerr << " Should be  >0.  Please check" << G4endl; 
  } else {
    energyMin = p;
#ifdef debug
    if (verboseLevel >= 1 ) {
      G4cout << " CCalPrimaryGeneratorAction: setting min. value of energy to "
             << energyMin/GeV << " GeV " << G4endl;
    }
#endif
  }
}


void CCalPrimaryGeneratorAction::SetMaximumEnergy(G4double p){

  if (p <= 0.) {
    G4cerr << "CCalPrimaryGeneratorAction::SetMaximumEnergy: value " << p/GeV 
           << "GeV is out of bounds, it will not be used" << G4endl;
    G4cerr << " Should be  >0.  Please check" << G4endl; 
  } else {
    energyMax = p;
#ifdef debug
    if (verboseLevel >= 1 ) {
      G4cout << " CCalPrimaryGeneratorAction: setting max. value of energy to "
             << energyMax/GeV << " GeV " << G4endl;
    }
#endif
  }
}


void CCalPrimaryGeneratorAction::SetMinimumPhi(G4double p){

  if (std::fabs(p)>2.*pi) {
    G4cerr << "CCalPrimaryGeneratorAction::SetMinimumPhi: setting value quite "
           << "large" << G4endl;
    G4cerr << " Should be given in radians - Please check" << G4endl;
  } else {
    phiMin = std::fabs(p);
#ifdef debug
    if (verboseLevel >= 1 ) {
      G4cout << " CCalPrimaryGeneratorAction: setting min. value of phi to "
             << phiMin << G4endl;
    }
#endif
  }
}


void CCalPrimaryGeneratorAction::SetMaximumPhi(G4double p){

  if (std::fabs(p)>2.*pi) {
    G4cerr << "CCalPrimaryGeneratorAction::SetMaximumPhi: setting value quite "
           << "large" << G4endl;
    G4cerr << " Should be given in radians - Please check" << G4endl;
  } else {
    phiMax = std::fabs(p);
#ifdef debug
    if (verboseLevel >= 1 ) {
      G4cout << " CCalPrimaryGeneratorAction: setting max. value of phi to "
             << phiMax << G4endl;
    }
#endif
  }
}


void CCalPrimaryGeneratorAction::SetStepsPhi(G4int val){

  if (val <= 0) {
    G4cerr << "CCalPrimaryGeneratorAction::SetStepsPhi: value " << val 
           << " is out of bounds, it will not be used" << G4endl;
    G4cerr << " Should be  > 0  Please check" << G4endl; 
  } else {
    phiSteps = val;
#ifdef debug
    if (verboseLevel >= 1 ) {
      G4cout << " CCalPrimaryGeneratorAction: setting no. of steps in phi to "
             << phiSteps << G4endl;
    }
#endif
  }
}


void CCalPrimaryGeneratorAction::SetMinimumEta(G4double p){

  etaMin = p;
#ifdef debug
  if (verboseLevel >= 1 ) {
    G4cout << " CCalPrimaryGeneratorAction: setting min. value of eta to "
           << etaMin << G4endl;
  }
#endif
}


void CCalPrimaryGeneratorAction::SetMaximumEta(G4double p){

  etaMax = p;
#ifdef debug
  if (verboseLevel >= 1 ) {
    G4cout << " CCalPrimaryGeneratorAction: setting max. value of eta to "
           << etaMax << G4endl;
  }
#endif
}


void CCalPrimaryGeneratorAction::SetStepsEta(G4int val){

  if (val <= 0) {
    G4cerr<<"CCalPrimaryGeneratorAction::SetStepsEta: value " << val << " is out of bounds, it will not be used"<<G4endl;
    G4cerr<<" Should be  > 0  Please check"<<G4endl; 
  } else {
    etaSteps = val;
#ifdef debug
    if (verboseLevel >= 1 ) {
      G4cout << " CCalPrimaryGeneratorAction: setting no. of steps in eta to "
             << etaSteps << G4endl;
    }
#endif
  }
}

void CCalPrimaryGeneratorAction::SetGunPosition(const G4ThreeVector & pos) const {

  particleGun->SetParticlePosition(pos);
}

void CCalPrimaryGeneratorAction::SetRunNo(G4int val){
  G4RunManager::GetRunManager()->SetRunIDCounter( val );
}

void CCalPrimaryGeneratorAction::initialize(){

  isInitialized = 1;

  print (1);
}


void CCalPrimaryGeneratorAction::print(G4int val){

  if (verboseLevel >= val) {

    if (generatorInput == singleRandom) {
      G4cout << G4endl
             << "**********************************************************************" << G4endl
             << "*                                                                    *" << G4endl  
             << "* CCalPrimaryGeneratorAction DEFAULT Random Energy/Direction setting:*" << G4endl
             << "*                                                                    *" << G4endl  
             << "*                                                                    *" << G4endl  
             << "*   Energy in    [ "<< energyMin/GeV   << " - " << energyMax/GeV   << "] (GeV) "<< G4endl
             << "*   Phi angle in [ "<< phiMin          << " - " << phiMax          << "] (rad) "<< G4endl
             << "*                [ "<< phiMin/degree   << " - " << phiMax/degree   << "] (deg) "<< G4endl 
             << "*   Eta in       [ "<< etaMin          << " - " << etaMax          << "]       "<< G4endl
             << "*                                                                    *" << G4endl  
             << "*                                                                    *" << G4endl  
                   << "**********************************************************************" << G4endl;
    } else if (generatorInput == singleScan) {
      G4cout << G4endl
             << "**********************************************************************" << G4endl
             << "*                                                                    *" << G4endl  
             << "* CCalPrimaryGeneratorAction DEFAULT Scan Direction settings :       *" << G4endl
             << "*                                                                    *" << G4endl  
             << "*                                                                    *" << G4endl  
             << "*   Phi angle in [ " << phiMin/degree   << " - " << phiMax/degree << "] (deg) " << G4endl
             << "*   Eta in       [ " << etaMin          << " - " << etaMax        << "]       " << G4endl
             << "*   Steps along eta " << etaSteps << " and along phi " << phiSteps << G4endl
             << "*                                                                    *" << G4endl  
             << "*                                                                    *" << G4endl  
                   << "**********************************************************************" << G4endl;
    } else  if (generatorInput == singleFixed) {
      G4cout << G4endl
             << "*******************************************************************" << G4endl
             << "*                                                                 *" << G4endl  
             << "* CCalPrimaryGeneratorAction: Current settings :                  *" << G4endl
             << "*                                                                 *" << G4endl  
             << "* " << particleGun->GetNumberOfParticles() 
             << "  " << particleGun->GetParticleDefinition()->GetParticleName() 
             << " of " << particleGun->GetParticleEnergy()/GeV << " GeV" << G4endl
             << "* will be shot from " << particleGun->GetParticlePosition() << G4endl;
      G4cout << "* in direction " << particleGun->GetParticleMomentumDirection() << G4endl;
      G4cout << "*                                                                 *" << G4endl  
             << "*                                                                 *" << G4endl  
             << "*******************************************************************" << G4endl;
    }
  } 
}
