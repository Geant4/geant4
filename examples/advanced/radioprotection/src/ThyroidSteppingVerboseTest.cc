//
// $Id: ThyroidSteppingVerboseTest.cc,v 1.2 2003-06-19 14:38:28 gunter Exp $
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
      //    G4cout << std::setw( 5) << fTrack->GetCurrentStepNumber() << " "
      //  << std::setw( 6) << G4BestUnit(fTrack->GetPosition().x(),"Length")
      //   << std::setw( 6) << G4BestUnit(fTrack->GetPosition().y(),"Length")
      //   << std::setw( 6) << G4BestUnit(fTrack->GetPosition().z(),"Length")
      //   << std::setw( 6) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
      //   << std::setw( 6) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
      //   << std::setw( 6) << G4BestUnit(fStep->GetStepLength(),"Length")
      //   << std::setw( 6) << G4BestUnit(fTrack->GetTrackLength(),"Length"); 
    G4cout << std::setw( 5) << fTrack->GetCurrentStepNumber() << " "
       << std::setw( 6) << fTrack-> GetDefinition()->GetParticleName()<< " "
           << std::setw( 6) << fTrack->GetPosition().x() << " "
	   << std::setw( 6) << fTrack->GetPosition().y() << " "
	   << std::setw( 6) << fTrack->GetPosition().z() << " "
	   << std::setw( 6) << fTrack->GetKineticEnergy() << " "
	   << std::setw( 6) << fStep->GetTotalEnergyDeposit() << " "
	   << std::setw( 6) << fStep->GetStepLength() << " "
	   << std::setw( 6) << fTrack->GetTrackLength()  << " ";
 



if( fTrack->GetNextVolume() != 0 ) { 
      G4cout << std::setw(10) << fTrack->GetVolume()->GetName();} else {
	G4cout << std::setw(10) << "OutOfWorld"; }
    if(fStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL){
      G4cout << "  " 
             << std::setw(10) << fStep->GetPostStepPoint()->GetProcessDefinedStep()
	                                ->GetProcessName();  }
      else { G4cout << "   UserLimit"; }
    //    G4cout << G4endl;
}
void ThyroidSteppingVerbose::TrackingStarted()
{}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....



























