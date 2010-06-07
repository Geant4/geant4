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
// $Id: Tst50RunAction.cc,v 1.27 2010-06-07 10:08:39 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 17 May  2003   S. Guatelli   1st implementation
//
// -------------------------------------------------------------------
 
#include "G4ios.hh"
#include <cmath>
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "Tst50RunAction.hh"
#include "Tst50RunMessenger.hh"
#include "Tst50PrimaryGeneratorAction.hh"
#include "Tst50DetectorConstruction.hh"
#include "Tst50AnalysisManager.hh"

Tst50RunAction::Tst50RunAction()
{
  messenger = new Tst50RunMessenger(this);
  flag = false;
}

Tst50RunAction::~Tst50RunAction()
{
  delete messenger;
}

void Tst50RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  runID = aRun->GetRunID();
    
  if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager* UI = G4UImanager::GetUIpointer();
      UI->ApplyCommand("/vis/scene/notifyHandlers");
    } 

  gammaTransmitted = 0;
  numberEvents = 0;
  particleTransmitted = 0;
  particleBackscattered = 0;
  backscatteredEnergy = 0;
}

G4int Tst50RunAction::GetRunID()
{
  return runID;
}

G4bool Tst50RunAction::GetFlag()
{
  //transmission test or CSDA Range and Stopping Power test for massive 
  //particles 
  return flag;
}

void Tst50RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4RunManager* runManager = G4RunManager::GetRunManager();
  primary =
    (Tst50PrimaryGeneratorAction*)(runManager->GetUserPrimaryGeneratorAction());
  G4double primaryParticleEnergy = primary->GetInitialEnergy();
  G4String primaryParticleName = primary -> GetParticle();
  detector =
    (Tst50DetectorConstruction*)(runManager->GetUserDetectorConstruction());
  G4String absorberMaterialName = detector->GetMaterialName(); 
  G4double absorberMaterialDensity = detector->GetDensity();
  G4double targetThickness = detector -> GetTargetThickness();

  if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
    }

  numberEvents = aRun->GetNumberOfEvent();

#ifdef G4ANALYSIS_USE
  if (primaryParticleName == "gamma")
    {
      G4double gammaTransmittedFraction = (gammaTransmitted/numberEvents);
      G4double gammaTransmittedFractionError = 1/(gammaTransmitted*(targetThickness*absorberMaterialDensity)); 
      G4double gammaAttenuationCoefficient = -(std::log(gammaTransmittedFraction))/(targetThickness*absorberMaterialDensity);
      Tst50AnalysisManager* analysis = Tst50AnalysisManager::getInstance();
      analysis -> AttenuationGammaCoeffiecient(runID,primaryParticleEnergy/MeV,gammaAttenuationCoefficient/(cm2/g),gammaTransmittedFractionError/(cm2/g));
    }
 
  if (flag)
    {
     if (primaryParticleName == "e-" || primaryParticleName == "e+")
       {
         
	G4double particleTransmittedFraction = (particleTransmitted/numberEvents) ;
	G4double particleTransmittedFractionError = (std::sqrt(particleTransmitted))/numberEvents;
	G4double particleBackscatteredFraction = (particleBackscattered/numberEvents);
	G4double particleBackscatteredFractionError= (std::sqrt(particleBackscattered))/numberEvents;
	Tst50AnalysisManager* analysis = Tst50AnalysisManager::getInstance();
	analysis -> ParticleTransmission(runID,primaryParticleEnergy/MeV,particleTransmittedFraction,particleBackscatteredFraction,particleTransmittedFractionError,particleBackscatteredFractionError);

// Test: fraction of backscattered energy in respect to total initial
// energy of the beamOn
	G4int numberOfEvents = aRun -> GetNumberOfEvent();
        G4double energyFraction = backscatteredEnergy/( numberOfEvents *primaryParticleEnergy);
        analysis -> TransmittedEnergy(runID,primaryParticleEnergy/MeV, energyFraction);
	}
    }
#endif
}

void  Tst50RunAction::SetTransmissionTest(G4String newValue)
{
  if (newValue == "on"){flag = true;} else {flag = false;};
} 

void  Tst50RunAction::TransmittedGammaNumber()
{
  gammaTransmitted += 1;
}

void  Tst50RunAction::TransmittedParticleNumber()
{
  particleTransmitted +=1;
}

void  Tst50RunAction::BackscatteredParticleNumber()
{
  particleBackscattered +=1;
}
void  Tst50RunAction::BackscatteredEnergy(G4double energy)
{
  backscatteredEnergy += energy;
}
