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
//
// $Id: ExDivRunAction.cc,v 1.2 2006-06-29 18:20:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4Timer.hh"
#include "ExDivRunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4ios.hh"
#include "ExVDivTester.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExDivRunAction::ExDivRunAction()
  : theNavigator(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExDivRunAction::~ExDivRunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExDivRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4int runID = aRun->GetRunID();
  G4cout << "### Run " << runID << " start." << G4endl;

  if (!runID)
  {
    theNavigator = G4TransportationManager::
      GetTransportationManager()->GetNavigatorForTracking();
    TestNavigator1();
    TestNavigator2();
  }

  if (G4VVisManager::GetConcreteInstance())
  {
    G4UImanager* UI = G4UImanager::GetUIpointer();
    UI->ApplyCommand("/vis/scene/notifyHandlers");
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExDivRunAction::EndOfRunAction(const G4Run*)
{
  if (G4VVisManager::GetConcreteInstance())
    {
     G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Test LocateGlobalPointAndSetup
//
void ExDivRunAction::TestNavigator1()
{
  G4Timer timer1;
  timer1.Start();
  G4cout << G4endl
         << " @@@@@@@@@@@@@@@@@@@@@  START testNavigator1  @@@@@@@@@@@@ "
         << G4endl;
  G4VPhysicalVolume *located=0;
   
  std::ofstream fout("output.points");
  std::ifstream fin("points.lis");
  G4double posX, posY, posZ;
  G4int ii = 0;
  for( ;;ii++ )
  {
    fin >> posX >> posY >> posZ;
    if( fin.eof() ) break;
    G4ThreeVector pos(posX,posY,posZ);
    located=theNavigator->LocateGlobalPointAndSetup(pos,0,false);
    //    G4cout << ii+1 << ". LOCATED POINT " <<  pos << " "
    //        << located->GetName() << " " << located->GetCopyNo() << G4endl;
    if( !ExVDivTester::bDivCylindrical ){
      fout << ii+1 << ". " << pos;
    } else {
      G4double phi = pos.phi()/deg;
      if( phi < 0. ) phi += 360.;
      fout << ii+1 << ". (" << pos.perp() << "," << phi << "," << pos.z() << ")";
    }
    fout << "  " << located->GetName() << " "
         << located->GetCopyNo() << G4endl;
  }
  timer1.Stop();
  G4cout << " TIME in TestNavigator1 " << timer1 << G4endl;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Test Stepping
//
void ExDivRunAction::TestNavigator2()
{
  G4Timer timer2;
  timer2.Start();
  G4cout << G4endl
         << " @@@@@@@@@@@@@@@@@@@@@  START testNavigator2  @@@@@@@@@@@@ "
         << G4endl;
  G4VPhysicalVolume *located=0;
  G4double Step,physStep,safety;
  G4ThreeVector* Hat = new G4ThreeVector[3];
  Hat[0] = G4ThreeVector(1,0,0);
  Hat[1] = G4ThreeVector(0,1,0);
  Hat[2] = G4ThreeVector(0,0,1);
   
  std::ofstream fout("output.step");
  std::ifstream fin("points.lis");
  G4double posX, posY, posZ;
  G4int ii = 0;
  for( ;;ii++ )
  {
    if( fin.eof() ) break;
    fin >> posX >> posY >> posZ;
    G4ThreeVector pos(posX,posY,posZ);
    located=theNavigator->LocateGlobalPointAndSetup(pos);
    physStep=kInfinity;
    if( !ExVDivTester::bDivCylindrical ){
      fout << ii+1 << ". " << pos;
    } else {
      G4double phi = pos.phi()/deg;
      if( phi < 0. ) phi += 360.;
      fout << ii+1 << ". R(" << pos.perp() << "," << phi << "," << pos.z() << ")";
    }
    for( G4int jj = 0; jj < 3; jj++ )
    {
      Step=theNavigator->ComputeStep(G4ThreeVector(posX,posY,posZ),
                             Hat[jj],physStep,safety);
      fout << "  " << Step << " ";
    }
    fout << G4endl;
  }
  timer2.Stop();
  G4cout << " TIME in TestNavigator2 " << timer2 << G4endl;
}

