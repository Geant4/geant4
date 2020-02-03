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
/// \file runAndEvent/RE02/src/RE02RunAction.cc
/// \brief Implementation of the RE02RunAction class
//
//
// 
#include "RE02RunAction.hh"
#include "RE02Run.hh"

//-- In order to obtain detector information.
#include "G4RunManager.hh"
#include "RE02DetectorConstruction.hh"
#include "G4THitsMap.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"    

#include <fstream>
//=======================================================================
// RE02RunAction
//  
//
//
//=======================================================================

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Constructor
RE02RunAction::RE02RunAction()
  : G4UserRunAction(),
    fNx(0), fNy(0), fNz(0)
{
  // - Prepare data member for RE02Run.
  //   vector represents a list of MultiFunctionalDetector names.
  fSDName.push_back(G4String("PhantomSD"));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Destructor.
RE02RunAction::~RE02RunAction()
{
  fSDName.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//== 
G4Run* RE02RunAction::GenerateRun()
{
  // Generate new RUN object, which is specially
  // dedicated for MultiFunctionalDetector scheme.
  //  Detail description can be found in RE02Run.hh/cc.
  return new RE02Run(fSDName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//==
void RE02RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//== 
void RE02RunAction::EndOfRunAction(const G4Run* aRun)
{
  if(!IsMaster()) return;

  //- RE02Run object.
  RE02Run* re02Run = (RE02Run*)aRun;
  //--- Dump all socred quantities involved in RE02Run.
  re02Run->DumpAllScorer();
  //---

  //
  //- water phantom (Detector) Information.
  //-- Number of segments in the water phantom.
  const RE02DetectorConstruction* detector =
      (const RE02DetectorConstruction*)
      (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  detector->GetNumberOfSegmentsInPhantom(fNx,fNy,fNz); //Fill fNx,y,z.

  //---------------------------------------------
  // Dump accumulated quantities for this RUN.
  //  (Display only central region of x-y plane)
  //---------------------------------------------
  G4THitsMap<G4double>* totEdep  = re02Run->GetHitsMap("PhantomSD/totalEDep");
  G4THitsMap<G4double>* proEdep  = re02Run->GetHitsMap("PhantomSD/protonEDep");
  G4THitsMap<G4double>* proNstep = re02Run->GetHitsMap("PhantomSD/protonNStep");
  G4THitsMap<G4double>* passCFx  =
    re02Run->GetHitsMap("PhantomSD/chargedPassCellFlux");
  G4THitsMap<G4double>*     cFx  =
    re02Run->GetHitsMap("PhantomSD/chargedCellFlux");
  G4THitsMap<G4double>*  surfFx  =
    re02Run->GetHitsMap("PhantomSD/chargedSurfFlux");
  G4THitsMap<G4double>* gCurr00  =
    re02Run->GetHitsMap("PhantomSD/gammaSurfCurr000");
  G4THitsMap<G4double>* gCurr01  =
    re02Run->GetHitsMap("PhantomSD/gammaSurfCurr001");
  G4THitsMap<G4double>* gCurr02  =
    re02Run->GetHitsMap("PhantomSD/gammaSurfCurr002");
  G4THitsMap<G4double>* gCurr03  =
    re02Run->GetHitsMap("PhantomSD/gammaSurfCurr003");

  G4cout << "============================================================="
         << G4endl;
  G4cout << " Number of event processed : "<< aRun->GetNumberOfEvent() << G4endl;
  G4cout << "============================================================="
         << G4endl;
  G4cout << std::setw( 8) << "#Z Cell#";
  G4cout << std::setw(16) << totEdep->GetName();
  G4cout << std::setw(16) << proEdep->GetName();
  G4cout << std::setw(12) << proNstep->GetName();
  G4cout << std::setw(21) << passCFx->GetName();
  G4cout << std::setw(20) << cFx->GetName();
  G4cout << std::setw(20) << surfFx->GetName();
  G4cout << std::setw(20) << gCurr00->GetName();
  G4cout << std::setw(20) << gCurr01->GetName();
  G4cout << std::setw(20) << gCurr02->GetName();
  G4cout << std::setw(20) << gCurr03->GetName()
         << G4endl;
  G4int ix = fNx/2;  
  G4int iy = fNy/2;
  G4int iz;
  //G4double totE, proE, proN,pasCF,CF,surfF,gCr0,gCr1,gCr2,gCr3;
  for ( iz = 0; iz < fNz; iz++){ 
      G4double* totED = (*totEdep)[CopyNo(ix,iy,iz)];
      G4double* proED = (*proEdep)[CopyNo(ix,iy,iz)];
      G4double* proNS = (*proNstep)[CopyNo(ix,iy,iz)];
      G4double* pasCF = (*passCFx)[CopyNo(ix,iy,iz)];
      G4double* cF    = (*cFx)[CopyNo(ix,iy,iz)];
      G4double* sfx   = (*surfFx)[CopyNo(ix,iy,iz)];
      G4double* gcur0 = (*gCurr00)[CopyNo(ix,iy,iz)];
      G4double* gcur1 = (*gCurr01)[CopyNo(ix,iy,iz)];
      G4double* gcur2 = (*gCurr02)[CopyNo(ix,iy,iz)];
      G4double* gcur3 = (*gCurr03)[CopyNo(ix,iy,iz)];
      if ( !totED ) totED = new G4double(0.0);
      if ( !proED ) proED = new G4double(0.0);
      if ( !proNS ) proNS = new G4double(0.0);
      if ( !pasCF ) pasCF = new G4double(0.0);
      if ( !cF    )    cF = new G4double(0.0);
      if ( !sfx   )   sfx = new G4double(0.0);
      if ( !gcur0 ) gcur0 = new G4double(0.0);
      if ( !gcur1 ) gcur1 = new G4double(0.0);
      if ( !gcur2 ) gcur2 = new G4double(0.0);
      if ( !gcur3 ) gcur3 = new G4double(0.0);
      G4cout << std::setw( 6) << iz << "  "
             << std::setw(12) << G4BestUnit(*totED,"Energy")
             << std::setw(12) << G4BestUnit(*proED,"Energy")
             << std::setw(12) << (*proNS)        << "   "
             << std::setw(13) << (*pasCF)*cm*cm  <<" /cm2"
             << std::setw(15) << (*cF)*cm*cm     <<" /cm2"
             << std::setw(15) << (*sfx)*cm*cm    <<" /cm2"
             << std::setw(15) << (*gcur0)*cm*cm  <<" /cm2"
             << std::setw(15) << (*gcur1)*cm*cm  <<" /cm2"
             << std::setw(15) << (*gcur2)*cm*cm  <<" /cm2"
             << std::setw(15) << (*gcur3)*cm*cm  <<" /cm2"
             << G4endl;
  }
  G4cout << "============================================="<<G4endl;

  std::ofstream  file("totED.txt");
  for ( iz = 0; iz < fNz; iz++){   
    for ( iy = 0; iy < fNy; iy++){ 
      for ( ix = 0; ix < fNx; ix++){ 
        G4double* totED = (*totEdep)[CopyNo(ix,iy,iz)];
        if ( !totED ) totED = new G4double(0.0);
        file << ix << " "<<iy<<" "<<iz<<" "<< *totED/MeV << G4endl;

      }
    }
  }
  file.close();
}
//
// --
