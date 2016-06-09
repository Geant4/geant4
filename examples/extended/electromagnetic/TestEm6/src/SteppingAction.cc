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
// $Id: SteppingAction.cc,v 1.5 2004/03/31 16:33:36 maire Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "RunAction.hh"
#include "G4SteppingManager.hh"
#include "G4VProcess.hh"
#include "G4ParticleTypes.hh"

#ifdef USE_AIDA
 #include "AIDA/IHistogram1D.h"
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
 
 static const G4double muonMass=G4MuonPlus::MuonPlus()->GetPDGMass();
 G4double GammaPlus=EGamma*xPlus/muonMass;
 G4double GammaMinus=EGamma*xMinus/muonMass;
   	   
#ifdef USE_AIDA
 runAction->GetHisto(0)->fill(1./(1.+pow(thetaPlus*GammaPlus,2)));
 runAction->GetHisto(1)->fill(log10(thetaPlus*GammaPlus));

 runAction->GetHisto(2)->fill(log10(thetaMinus*GammaMinus));
 runAction->GetHisto(3)->fill(log10(fabs(thetaPlus *GammaPlus
                                              -thetaMinus*GammaMinus)));
 
 runAction->GetHisto(4)->fill(xPlus);
 runAction->GetHisto(5)->fill(xMinus); 
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


