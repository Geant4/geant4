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
// $Id: SteppingAction.cc,v 1.3 2006-05-24 12:58:49 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "RunAction.hh"
#include "G4SteppingManager.hh"
#include "G4VProcess.hh"
#include "G4ParticleTypes.hh"
#include "G4UnitsTable.hh"

#ifdef G4ANALYSIS_USE
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
 static G4int iCalled=0;
 const  G4int nprint=0; // set to 10 to get debug print for first 10 calls

 if (processName == "SynRad")
 {
   iCalled++;
   G4StepPoint* PrePoint = aStep->GetPreStepPoint();
   G4TrackVector* secondary = fpSteppingManager->GetSecondary();

   //(*secondary)[lp-1]  points to the last photon generated
   //   
   size_t lp=(*secondary).size();
   if (lp)
   {
     G4double  Egamma =  (*secondary)[lp-1]->GetTotalEnergy();
     runAction->n_gam_sync++;
     runAction->e_gam_sync += Egamma;
     runAction->e_gam_sync2 += Egamma*Egamma;
     if (Egamma > runAction->e_gam_sync_max) runAction->e_gam_sync_max = Egamma;
     runAction->lam_gam_sync += aStep->GetStepLength();
     if (iCalled<nprint)
     {
       G4double      Eelec  = PrePoint->GetTotalEnergy();
       G4ThreeVector Pelec  = PrePoint->GetMomentum();
       G4ThreeVector PGamma = (*secondary)[lp-1]->GetMomentum();
       G4bool IsGamma = 
            ((*secondary)[lp-1]->GetDefinition() == G4Gamma::GammaDefinition());
	    
       G4cout << "UserSteppingAction processName=" << process->GetProcessName()
         << " Step Length=" << std::setw(6) 
	                    << G4BestUnit(aStep->GetStepLength(),"Length")
         << " Eelec=" << G4BestUnit(Eelec,"Energy")
         << " Pelec=" << G4BestUnit(Pelec,"Energy")
         << " IsGamma=" << IsGamma
         << " Egamma=" << G4BestUnit(Egamma,"Energy")
         << " PGamma=" << G4BestUnit(PGamma,"Energy")
         << " #secondaries lp=" << lp
         << '\n';
      }

#ifdef G4ANALYSIS_USE
     // fill histos
     if( runAction->GetHisto(0) ) // check the histos exist
     {
       runAction->GetHisto(0)->fill(Egamma/keV);
       runAction->GetHisto(1)->fill(Egamma/keV,Egamma/keV);
       runAction->GetHisto(2)->fill(aStep->GetStepLength()/m);
     }
#endif
   }
 }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
