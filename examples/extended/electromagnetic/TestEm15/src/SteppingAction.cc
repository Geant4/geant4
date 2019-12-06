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
/// \file electromagnetic/TestEm15/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"
#include "G4ParticleTypes.hh"

#include "G4RunManager.hh"

#include <G4ThreeVector.hh>
#include <G4RotationMatrix.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det,
                               RunAction* RuAct)
:G4UserSteppingAction(),fDetector(det), fRunAction(RuAct)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4StepPoint* prePoint = aStep->GetPreStepPoint();
  
  // if World --> return
  if(prePoint->GetTouchableHandle()->GetVolume()==fDetector->GetWorld()) return;
  
  // here we enter in the absorber Box
  // tag the event to be killed anyway after this step
  //
  G4RunManager::GetRunManager()->AbortEvent();
  
  //count processes and keep only Multiple Scattering or gamma converion
  //  
  G4StepPoint* endPoint = aStep->GetPostStepPoint();
  G4String procName = endPoint->GetProcessDefinedStep()->GetProcessName();
  fRunAction->CountProcesses(procName);
      
  if (procName == "msc" || procName == "muMsc" || procName == "stepMax") {
    
    //below, only multiple Scattering happens
    //
    G4ThreeVector position  = endPoint->GetPosition();
    G4ThreeVector direction = endPoint->GetMomentumDirection();
    
    G4double truePathLength = aStep->GetStepLength();      
    G4double geomPathLength = position.x() + 0.5*fDetector->GetBoxSize();
    G4double ratio = geomPathLength/truePathLength;
    fRunAction->SumPathLength(truePathLength,geomPathLength);
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();  
    analysisManager->FillH1(1,truePathLength);
    analysisManager->FillH1(2,geomPathLength);
    analysisManager->FillH1(3,ratio);
    
    G4double yend = position.y(), zend = position.z();
    G4double lateralDisplacement = std::sqrt(yend*yend + zend*zend);
    fRunAction->SumLateralDisplacement(lateralDisplacement);
    analysisManager->FillH1(4,lateralDisplacement);
    
    G4double psi = std::atan(lateralDisplacement/geomPathLength); 
    fRunAction->SumPsi(psi);
    analysisManager->FillH1(5,psi);
    
    G4double xdir = direction.x(),  ydir = direction.y(), zdir = direction.z();
    G4double tetaPlane = std::atan2(ydir, xdir); 
    fRunAction->SumTetaPlane(tetaPlane);
    analysisManager->FillH1(6,tetaPlane);
    tetaPlane = std::atan2(zdir, xdir); 
    fRunAction->SumTetaPlane(tetaPlane);
    analysisManager->FillH1(6,tetaPlane);
    
    G4double phiPos = std::atan2(zend, yend); 
    analysisManager->FillH1(7,phiPos);
    G4double phiDir = std::atan2(zdir, ydir); 
    analysisManager->FillH1(8,phiDir);

    G4double phiCorrel = 0.;
    if (lateralDisplacement > 0.)  
      phiCorrel = (yend*ydir + zend*zdir)/lateralDisplacement;
    fRunAction->SumPhiCorrel(phiCorrel);
    analysisManager->FillH1(9,phiCorrel);
  } else if (procName == "conv" || procName == "GammaToMuPair" ) {

    // gamma conversion
    
    G4StepPoint* PrePoint = aStep->GetPreStepPoint();
    G4double      EGamma  = PrePoint->GetTotalEnergy();
    G4ThreeVector PGamma  = PrePoint->GetMomentum();
    G4ThreeVector PolaGamma  = PrePoint->GetPolarization();

    G4double Eplus=-1;
    G4ThreeVector Pplus, Pminus, Precoil;

    const G4TrackVector* secondary = fpSteppingManager->GetSecondary();

    const size_t Nsecondaries = (*secondary).size();

    //No conversion , E < threshold
    if (Nsecondaries == 0) return;
  
    for (size_t lp=0; lp< std::min(Nsecondaries,size_t(2) ); lp++) {
      if  (((*secondary)[lp]->GetDefinition()==G4Electron::Definition())
           || ((*secondary)[lp]->GetDefinition()==G4MuonMinus::Definition()) )
	{
	  Pminus = (*secondary)[lp]->GetMomentum();
	}
      if (((*secondary)[lp]->GetDefinition()==G4Positron::Definition())
          || ((*secondary)[lp]->GetDefinition()==G4MuonPlus::Definition()) )
	{
	  Eplus  = (*secondary)[lp]->GetTotalEnergy();
	  Pplus  = (*secondary)[lp]->GetMomentum();
	}
    }

    if ( Nsecondaries >= 3 ) {
      Precoil  = (*secondary)[2]->GetMomentum();
    }

    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
    // Fill Histograms
    
    G4ThreeVector z = PGamma.unit(); // gamma direction
    G4ThreeVector x(1.,0.,0.);
    
    // pola perpendicular to direction

    if ( PolaGamma.mag() != 0.0 ) {
      x = PolaGamma.unit();
    } else { // Pola = 0 case
      x = z.orthogonal().unit();
    }

    G4ThreeVector y = z;
    y = y.cross(x);

    G4RotationMatrix GtoW(x,y,z); // from  gamma ref. sys. to World
    G4RotationMatrix WtoG = inverseOf(GtoW); // from World to gamma ref. sys.


    G4double angleE = Pplus.angle(Pminus) * EGamma;
    analysisManager->FillH1(10,angleE);
 
    if ( Nsecondaries >= 3 ) {
      // recoil returned
      analysisManager->FillH1(11,std::log10(Precoil.mag()));
      analysisManager->FillH1(12,Precoil.transform(WtoG).phi());
    }
    G4double phiPlus =  Pplus.transform(WtoG).phi();
    G4double phiMinus =  Pminus.transform(WtoG).phi();
    analysisManager->FillH1(13,phiPlus);
    analysisManager->FillH1(14,std::cos(phiPlus + phiMinus) * -2.0);
    analysisManager->FillH1(15,Eplus/EGamma);
    
    G4double phiPola =  PolaGamma.transform(WtoG).phi();
    analysisManager->FillH1(16, phiPola);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
