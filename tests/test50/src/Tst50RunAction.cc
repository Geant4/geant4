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
// $Id: Tst50RunAction.cc,v 1.23 2003-07-28 15:05:52 guatelli Exp $
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
#include <math.h>
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
  
  //Tst50AnalysisManager* analysis = Tst50AnalysisManager::getInstance();
  //analysis->bookHistograms();
  
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
  
  if (primaryParticleName == "gamma")
    {
      G4double gammaTransmittedFraction = (gammaTransmitted/numberEvents);
      G4double gammaTransmittedFractionError = 1/(gammaTransmitted*(targetThickness*absorberMaterialDensity)); 
      G4double gammaAttenuationCoefficient = -(log(gammaTransmittedFraction))/(targetThickness*absorberMaterialDensity);
      Tst50AnalysisManager* analysis = Tst50AnalysisManager::getInstance();
      analysis -> AttenuationGammaCoeffiecient(runID,primaryParticleEnergy/MeV,gammaAttenuationCoefficient/(cm2/g),gammaTransmittedFractionError/(cm2/g));
    }
 
  if (flag)
    {
     if (primaryParticleName == "e-" || primaryParticleName == "e+")
       {
         
	G4double particleTransmittedFraction = (particleTransmitted/numberEvents) ;
	G4double particleTransmittedFractionError = (sqrt(particleTransmitted))/numberEvents;
	G4double particleBackscatteredFraction = (particleBackscattered/numberEvents);
	G4double particleBackscatteredFractionError= (sqrt(particleBackscattered))/numberEvents;
	Tst50AnalysisManager* analysis = Tst50AnalysisManager::getInstance();
	analysis -> ParticleTransmission(runID,primaryParticleEnergy/MeV,particleTransmittedFraction,particleBackscatteredFraction,particleTransmittedFractionError,particleBackscatteredFractionError);
        analysis -> TransmittedEnergy(runID,primaryParticleEnergy/MeV,backscatteredEnergy/MeV);
	}
    }
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
