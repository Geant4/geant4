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
//
// $Id: SteppingAction.cc,v 1.1 2002-05-23 13:30:44 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "RunAction.hh"
#include "G4SteppingManager.hh"
#include "G4VProcess.hh"
#include "G4ParticleTypes.hh"

#ifndef G4NOHIST
 #include "CLHEP/Hist/HBookFile.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(RunAction* RuAct)
:runAction(RuAct)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
 const G4VProcess* process = aStep->GetPostStepPoint()->GetProcessDefinedStep();
 if (process == 0) return;  
 G4String processName = process->GetProcessName();
 if (processName != "GammaToMuPair") return;
 
 G4StepPoint* PrePoint = aStep->GetPreStepPoint();  
 G4double      EGamma  = PrePoint->GetTotalEnergy();
 G4ThreeVector PGamma  = PrePoint->GetMomentum();
    
 G4double      Eplus(0), Eminus(0);
 G4ThreeVector Pplus   , Pminus;
 G4TrackVector* secondary = fpSteppingManager->GetSecondary();
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
	   
 #ifndef G4NOHIST
 static const G4double muonMass=G4MuonPlus::MuonPlus()->GetPDGMass();
  
 G4double GammaPlus=EGamma*xPlus/muonMass;
 runAction->GetHisto(0)->accumulate(1./(1.+pow(thetaPlus*GammaPlus,2)));
 runAction->GetHisto(1)->accumulate(log10(thetaPlus*GammaPlus));

 G4double GammaMinus=EGamma*xMinus/muonMass;
 runAction->GetHisto(2)->accumulate(log10(thetaMinus*GammaMinus));
 runAction->GetHisto(3)->accumulate(log10(fabs(thetaPlus *GammaPlus
                                              -thetaMinus*GammaMinus)));
 
 runAction->GetHisto(4)->accumulate(xPlus);
 runAction->GetHisto(5)->accumulate(xMinus); 
 #endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


