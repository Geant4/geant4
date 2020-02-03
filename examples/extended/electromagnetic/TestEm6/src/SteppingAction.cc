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
/// \file electromagnetic/TestEm6/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "RunAction.hh"
#include "G4SteppingManager.hh"
#include "G4VProcess.hh"
#include "G4ParticleTypes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(RunAction* RuAct)
:G4UserSteppingAction(),fRunAction(RuAct)
{ 
 fMuonMass = G4MuonPlus::MuonPlus()->GetPDGMass();
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
 fRunAction->CountProcesses(processName); //count processes
  
 if (processName != "GammaToMuPair") return;
 
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
 G4double GammaPlus = Eplus/fMuonMass;
 G4double GammaMinus= Eminus/fMuonMass;
 
 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

 if(0.0 == thetaPlus || 0.0 == thetaMinus) {
   G4cout << "SteppingAction: "
          << "thetaPlus= " << thetaPlus << " thetaMinus= " << thetaMinus
          << " gamPlus= " << GammaPlus << " gamMinus= " <<  GammaMinus
          << "  " << thetaPlus *GammaPlus - thetaMinus*GammaMinus << G4endl;
   return;
 }
 analysisManager->FillH1(1,1./(1.+std::pow(thetaPlus*GammaPlus,2)));
 analysisManager->FillH1(2,std::log10(thetaPlus*GammaPlus));

 analysisManager->FillH1(3,std::log10(thetaMinus*GammaMinus));
 analysisManager->FillH1(4,std::log10(std::fabs(thetaPlus *GammaPlus
                                              -thetaMinus*GammaMinus)));
 
 analysisManager->FillH1(5,xPlus);
 analysisManager->FillH1(6,xMinus);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


