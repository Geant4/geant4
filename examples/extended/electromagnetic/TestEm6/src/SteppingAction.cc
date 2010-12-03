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
// $Id: SteppingAction.cc,v 1.11 2010-12-03 14:54:55 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "RunAction.hh"
#include "G4SteppingManager.hh"
#include "G4VProcess.hh"
#include "G4ParticleTypes.hh"

#ifdef G4ANALYSIS_USE
 #include "AIDA/IHistogram1D.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(RunAction* RuAct)
:runAction(RuAct)
{ 
 muonMass = G4MuonPlus::MuonPlus()->GetPDGMass();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
 const G4VProcess* process = aStep->GetPostStepPoint()->GetProcessDefinedStep();
 if (process == 0) return;  
 G4String processName = process->GetProcessName();
 runAction->CountProcesses(processName); //count processes
  
 if (processName != "GammaToMuPair") return;
 
#ifdef G4ANALYSIS_USE
 G4StepPoint* PrePoint = aStep->GetPreStepPoint();  
 G4double      EGamma  = PrePoint->GetTotalEnergy();
 G4ThreeVector PGamma  = PrePoint->GetMomentum();
    
 G4double      Eplus(0), Eminus(0);
 G4ThreeVector Pplus   , Pminus;
 const G4TrackVector* secondary = fpSteppingManager->GetSecondary();
 for (size_t lp=0; lp<(*secondary).size(); lp++) {
   if ((*secondary)[lp]->GetDefinition()==G4MuonPlus::MuonPlusDefinition()) {
     Eplus  = (*secondary)[lp]->GetTotalEnergy();
     Pplus  = (*secondary)[lp]->GetMomentum();
   } else {
     Eminus = (*secondary)[lp]->GetTotalEnergy();
     Pminus = (*secondary)[lp]->GetMomentum();      	   
   }
 }
    	   
 G4double xPlus = Eplus/EGamma, xMinus = Eminus/EGamma;
 G4double thetaPlus = PGamma.angle(Pplus), thetaMinus = PGamma.angle(Pminus);
 G4double GammaPlus=EGamma*xPlus/muonMass;
 G4double GammaMinus=EGamma*xMinus/muonMass;

 runAction->GetHisto(0)->fill(1./(1.+std::pow(thetaPlus*GammaPlus,2)));
 runAction->GetHisto(1)->fill(std::log10(thetaPlus*GammaPlus));

 runAction->GetHisto(2)->fill(std::log10(thetaMinus*GammaMinus));
 runAction->GetHisto(3)->fill(std::log10(std::fabs(thetaPlus *GammaPlus
                                              -thetaMinus*GammaMinus)));
 
 runAction->GetHisto(4)->fill(xPlus);
 runAction->GetHisto(5)->fill(xMinus); 
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


