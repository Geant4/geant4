//
// $Id: ThyroidSteppingVerboseTest.cc,v 1.1 2003-05-23 11:56:01 francy Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#include "ThyroidSteppingVerboseTest.hh"

#include "G4SteppingManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ThyroidSteppingVerbose::ThyroidSteppingVerbose()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ThyroidSteppingVerbose::~ThyroidSteppingVerbose()
{}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ThyroidSteppingVerbose::StepInfo()
{
  CopyState();
  G4int prec = G4cout.precision(3);
      G4cout << G4endl;    
      //    G4cout << G4std::setw( 5) << fTrack->GetCurrentStepNumber() << " "
      //  << G4std::setw( 6) << G4BestUnit(fTrack->GetPosition().x(),"Length")
      //   << G4std::setw( 6) << G4BestUnit(fTrack->GetPosition().y(),"Length")
      //   << G4std::setw( 6) << G4BestUnit(fTrack->GetPosition().z(),"Length")
      //   << G4std::setw( 6) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
      //   << G4std::setw( 6) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
      //   << G4std::setw( 6) << G4BestUnit(fStep->GetStepLength(),"Length")
      //   << G4std::setw( 6) << G4BestUnit(fTrack->GetTrackLength(),"Length"); 
    G4cout << G4std::setw( 5) << fTrack->GetCurrentStepNumber() << " "
       << G4std::setw( 6) << fTrack-> GetDefinition()->GetParticleName()<< " "
           << G4std::setw( 6) << fTrack->GetPosition().x() << " "
	   << G4std::setw( 6) << fTrack->GetPosition().y() << " "
	   << G4std::setw( 6) << fTrack->GetPosition().z() << " "
	   << G4std::setw( 6) << fTrack->GetKineticEnergy() << " "
	   << G4std::setw( 6) << fStep->GetTotalEnergyDeposit() << " "
	   << G4std::setw( 6) << fStep->GetStepLength() << " "
	   << G4std::setw( 6) << fTrack->GetTrackLength()  << " ";
 



if( fTrack->GetNextVolume() != 0 ) { 
      G4cout << G4std::setw(10) << fTrack->GetVolume()->GetName();} else {
	G4cout << G4std::setw(10) << "OutOfWorld"; }
    if(fStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL){
      G4cout << "  " 
             << G4std::setw(10) << fStep->GetPostStepPoint()->GetProcessDefinedStep()
	                                ->GetProcessName();  }
      else { G4cout << "   UserLimit"; }
    //    G4cout << G4endl;
}
void ThyroidSteppingVerbose::TrackingStarted()
{}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....



























