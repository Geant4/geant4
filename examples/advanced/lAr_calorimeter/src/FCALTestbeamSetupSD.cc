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
// $Id: FCALTestbeamSetupSD.cc,v 1.11 2006-06-29 16:03:25 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "FCALTestbeamSetupSD.hh"

#include "FCALCalorHit.hh"

#include "FCALTestbeamSetup.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

#include "G4ios.hh"
#include <iostream>
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALTestbeamSetupSD::FCALTestbeamSetupSD(G4String name) : G4VSensitiveDetector(name)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALTestbeamSetupSD::~FCALTestbeamSetupSD()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALTestbeamSetupSD::Initialize(G4HCofThisEvent*)
{
  EBeamS1 = EBeamS2 = EBeamS3 = 0.;
  EHoleScint = EBeamHole = 0.;
  EBeamDead = 0;
  G4int j;
  for (j =0 ; j<8 ; j++) { ETailVis[j] = 0.;};
  for (j =0 ; j<7 ; j++) { ETailDep[j] = 0.;};
  TailCatcherID = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool FCALTestbeamSetupSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  
  G4double edep = aStep->GetTotalEnergyDeposit();
  if (edep==0.) return false;

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());  
  G4VPhysicalVolume* physVol = theTouchable->GetVolume();

  
  if(strcmp(physVol->GetName(),"ScintS1Physical")==0) { 
    EBeamS1 = EBeamS1 + edep;}
  
  if(strcmp(physVol->GetName(),"ScintS2Physical")==0) { 
    EBeamS2 = EBeamS2 + edep;}
  
  if(strcmp(physVol->GetName(),"ScintS3Physical")==0) {  
    EBeamS3 = EBeamS3 + edep;}

  if(strcmp(physVol->GetName(),"HoleScintPhysical")==0){ EHoleScint =+ edep;}
  if(strcmp(physVol->GetName(),"HoleCntrScintPhysical")==0){
    EBeamHole = EBeamHole + edep;}

  if(strcmp(physVol->GetName(),"MWPCPhysical")==0) { EBeamDead =+ edep;};
  if(strcmp(physVol->GetName(),"HoleCntrPbPhysical")==0) { EBeamDead =+ edep;};
  if(strcmp(physVol->GetName(),"HoleCntrAlPhysical")==0) { EBeamDead =+ edep;};
  if(strcmp(physVol->GetName(),"LeadWallPhysical")==0) { EBeamDead =+ edep;};
  if(strcmp(physVol->GetName(),"IronWallPhysical")==0) { EBeamDead =+ edep;};

  if(strcmp(physVol->GetName(),"BigScintPhysical")==0) {
    TailCatcherID = physVol->GetCopyNo();
    if(TailCatcherID > 0)  ETailVis[TailCatcherID] =+ edep;
  }
  if(strcmp(physVol->GetName(),"SmallScintPhysical")==0) {
    TailCatcherID = physVol->GetCopyNo();
    if(TailCatcherID > 0)  ETailVis[TailCatcherID + 3] =+ edep;
  }
  
  if(strcmp(physVol->GetName(),"BigIronPhysical")==0) {
    TailCatcherID = physVol->GetCopyNo();
    if(TailCatcherID > 0)  ETailDep[TailCatcherID] =+ edep;
  }
  if(strcmp(physVol->GetName(),"SmallIronPhysical")==0) {
    TailCatcherID = physVol->GetCopyNo();
    if(TailCatcherID > 0)  ETailDep[TailCatcherID+2] =+ edep;
  }

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALTestbeamSetupSD::EndOfEvent(G4HCofThisEvent*)
{
  G4cout << " Visisble Energy in S1 , S2 , S3 in (MeV)" << G4endl;
  G4cout << EBeamS1/MeV << " " << EBeamS2/MeV << " " << EBeamS3/MeV << " " << G4endl;

  G4cout << " Visible Energy in Hole Counter  (MeV) " << G4endl;
  G4cout << EHoleScint/MeV << " " << EBeamHole/MeV << G4endl;

  G4cout << " Visible Energy in Upstream Dead Materials " << G4endl;
  G4cout << EBeamDead/MeV << G4endl;

  G4cout << " Visible Energy in Tail Catcher Scintillator" << G4endl;
  G4int j;
  for (j=1; j<8 ; j++) {G4cout <<  ETailVis[j]/MeV << " " ;};  G4cout << G4endl;
 
  G4cout << " Visible Energy in Tail Catcher Absorber" << G4endl;
  for (j=1; j<7 ; j++) {G4cout <<  ETailDep[j]/MeV << " " ;};  G4cout << G4endl;

// Write data in File
  //-------------------
  G4String FileName = "Beam_802_1mm.dat";
  std::ios::openmode iostemp;
  if(InitBeam == 0) {
    iostemp = std::ios::out;
    InitBeam++;
  } else {
    iostemp = std::ios::out|std::ios::app;   
  };
  
  std::ofstream BeamDatafile(FileName, iostemp);
  // BeamDatafile.precision(5);

  BeamDatafile << EBeamS1/MeV << " " << EBeamS2/MeV << " " << EBeamS3/MeV << std::endl;
  BeamDatafile << EBeamHole/MeV << std::endl;
  BeamDatafile << EBeamDead/MeV << std::endl;
  for (j=1; j<8 ; j++) { BeamDatafile << ETailVis[j]/MeV << " " ;} ; BeamDatafile << std::endl;
  for (j=1; j<7 ; j++) { BeamDatafile << ETailDep[j]/MeV << " " ;} ; BeamDatafile << std::endl;
 
  BeamDatafile.close();

       
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALTestbeamSetupSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALTestbeamSetupSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALTestbeamSetupSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

