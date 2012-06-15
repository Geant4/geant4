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
// $Id: ExN03SteppingAction.cc,v 1.15 2006/06/29 17:49:13 gunter Exp $
// GEANT4 tag $Name: geant4-08-01-patch-01-ref $
//
//
                                                                                
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......                                                                                
#include "ExN01SteppingAction.hh"
                                                                                
#include "ExN01DetectorConstruction.hh"
#include "ExN01EventAction.hh"
                                                                                
#include "G4Step.hh"
 #include "G4UnitsTable.hh"
                                                                                
////#include "G4RunManager.hh"
                                                                                
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
ExN01SteppingAction::ExN01SteppingAction(ExN01DetectorConstruction* det,
                                          ExN01EventAction* evt)
:detector(det), event_action(evt)                                         
{ }                                                                               
                                                      
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......                                                                                
ExN01SteppingAction::~ExN01SteppingAction()
{ }
                                                                                
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......                                                                                
void ExN01SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  // get volume of the current step
  // G4VPhysicalVolume* volume =
  // aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  // G4double EndStepX=0.;
  // G4double EndStepY=0.;
  // G4double EndStepZ=0.;
    G4String particleName = (aStep -> GetTrack() -> GetDynamicParticle()
                           -> GetDefinition() -> GetParticleName());
    if( particleName =="e-")
    {
      //      G4cout<<particleName<<G4endl;
      // EndStepX =
      // aStep->GetPreStepPoint()->GetPosition().x();
      // EndStepY =
      // aStep->GetPreStepPoint()->GetPosition().y();
      // EndStepZ =
      // aStep->GetPreStepPoint()->GetPosition().z();
        G4cout.precision(15);
	// G4cout<<particleName<<"  "<<"PreStep"<<" "
        //    << std::setw(6) << G4BestUnit(aStep->GetPreStepPoint()->GetPosition().x(),"Length")<<"  "
        //    << std::setw(6) << G4BestUnit(aStep->GetPreStepPoint()->GetPosition().y(),"Length")<<"  "
        //    << std::setw(6) << G4BestUnit(aStep->GetPreStepPoint()->GetPosition().z(),"Length")<<"  "
        //      << std::setw(6) << G4BestUnit(std::sqrt(EndStepX*EndStepX+(EndStepY)*(EndStepY)),"Length")<<"  " 
        //    
        //      << std::setw(6) << G4BestUnit(aStep->GetStepLength(),"Length")<<"  "
        //    << std::setw(10) <<aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName()<<"  "
        //    <<G4endl;
	//	G4cout<<particleName<<"  "<<"PostStep  "
	//  << std::setw(6) << G4BestUnit(aStep->GetPostStepPoint()->GetPosition().x(),"Length")<<"  "
	//  << std::setw(6) << G4BestUnit(aStep->GetPostStepPoint()->GetPosition().y(),"Length")<<"  "
	//  << std::setw(6) << G4BestUnit(aStep->GetPostStepPoint()->GetPosition().z(),"Length")<<"  "
	//    << std::setw(6) << G4BestUnit(aStep->GetStepLength(),"Length")<<"  "
	//   
	//  <<G4endl;
       //
       // const G4int verboseLevel=1;
       //G4cout.precision(6);
       // G4cout.setf(ios_base::fixed,ios_base::floatfield);
       //  G4cout << std::setw( 9) << "X(mm)" << " "
       //     << std::setw( 9) << "Y(mm)" << " "
       //     << std::setw( 9) << "Z(mm)" << " "
       //     << std::setw( 9) << " N_x " << " "
       //     << std::setw( 9) << " N_y " << " "
       //     << std::setw( 9) << " N_z " << " "
       //	 //  << std::setw( 9) << " Delta|N|" << " "
       //  //   << std::setw( 9) << " Delta(N_z) " << " "
       //  //   << std::setw( 9) << "StepLen" << " "
       //  //   << std::setw( 9) << "PhsStep" << " "
       //     << std::setw( 9) << "Safety" << " "
       //     << std::setw(18) << "NextVolume" << " "
       //     << G4endl;
    
 
   
       //G4cout.precision(6);
       //G4cout <<std::setw( 9) <<aStep->GetPreStepPoint()->GetPosition().x() << " "
       //     << std::setw( 9) <<aStep->GetPreStepPoint()->GetPosition().y() << " "
       //     << std::setw( 9) <<aStep->GetPreStepPoint()->GetPosition().z() << " "
       //	 //   << std::setw( 9) <<aStep->GetPreStepPoint()->GetVelocity.x() << " "
       //  //   << std::setw( 9) << aStep->GetPreStepPoint()->GetVelocity.y() << " "
       //  //    << std::setw( 9) <<aStep->GetPreStepPoint()->GetVelocity.z() << " ";
       //	 // G4cout.precision(2);
	 
       ///	 // << std::setw( 9) <<aStep->GetPreStepPoint()->GetSafety << " "
       //    << std::setw(12) << aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume() << " "
       //     <<G4endl;                                                                         



     }
      if( particleName =="chargedgeantino")
    {
      //      G4cout<<particleName<<G4endl;
      // EndStepX =
      // aStep->GetPreStepPoint()->GetPosition().x();
      // EndStepY =
      //   aStep->GetPreStepPoint()->GetPosition().y();
      // EndStepZ =
      //   aStep->GetPreStepPoint()->GetPosition().z();
       G4cout<<particleName<<"  "
            << std::setw(6) << G4BestUnit(aStep->GetPreStepPoint()->GetPosition().x(),"Length")<<"  "
            << std::setw(6) << G4BestUnit(aStep->GetPreStepPoint()->GetPosition().y(),"Length")<<"  "
            << std::setw(6) << G4BestUnit(aStep->GetPreStepPoint()->GetPosition().z(),"Length")<<"  "
            <<G4endl;
      
     }
      if( particleName =="mu+")
    {
      //      G4cout<<particleName<<G4endl;
      // EndStepX =
      // aStep->GetPreStepPoint()->GetPosition().x();
      // EndStepY =
      //   aStep->GetPreStepPoint()->GetPosition().y();
      // EndStepZ =
      //   aStep->GetPreStepPoint()->GetPosition().z();
       G4cout<<particleName<<"  "
            << std::setw(6) << G4BestUnit(aStep->GetPreStepPoint()->GetPosition().x(),"Length")<<"  "
            << std::setw(6) << G4BestUnit(aStep->GetPreStepPoint()->GetPosition().y(),"Length")<<"  "
            << std::setw(6) << G4BestUnit(aStep->GetPreStepPoint()->GetPosition().z(),"Length")<<"  "
            <<G4endl;
      
     }

  //end tanja
      //G4Exception("G4PropagatorInField::ComputeStep()", "Stop after first step",
      //        FatalException, 
      //          "check for propagation");

  // collect energy and track length step by step

  //end tanja
  //example of saving random number seed of this event, under condition
  //// if (condition) G4RunManager::GetRunManager()->rndmSaveThisEvent();
}
                                                                                
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
