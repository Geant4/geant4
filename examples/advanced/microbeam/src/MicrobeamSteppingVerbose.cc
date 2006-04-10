// -------------------------------------------------------------------
// $Id: MicrobeamSteppingVerbose.cc,v 1.2 2006-04-10 14:47:32 sincerti Exp $
// -------------------------------------------------------------------

#include "G4UnitsTable.hh"

#include "MicrobeamSteppingVerbose.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MicrobeamSteppingVerbose::MicrobeamSteppingVerbose()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MicrobeamSteppingVerbose::~MicrobeamSteppingVerbose()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicrobeamSteppingVerbose::StepInfo()
{
  CopyState();
  
  G4int prec = G4cout.precision(6);
  
  if( verboseLevel >= 1 ){
    if( verboseLevel >= 4 ) VerboseTrack();
    if( verboseLevel >= 3 ){
      G4cout << G4endl;    
      G4cout << std::setw( 5) << "#Step#"     << " "
	     << std::setw( 6) << "X"          << "    "
	     << std::setw( 6) << "Y"          << "    "  
	     << std::setw( 6) << "Z"          << "    "
	     << std::setw( 9) << "KineE"      << " "
	     << std::setw( 9) << "dEStep"     << " "  
	     << std::setw(10) << "StepLeng"     
	     << std::setw(10) << "TrakLeng" 
	     << std::setw(10) << "NextVolu" 
	     << std::setw(10) << "Process"   << G4endl;	          
    }
   
    G4cout << std::setw( 5) << fTrack->GetCurrentStepNumber() << " "
	   << std::setw( 7) << G4BestUnit(fTrack->GetPosition().x(),"Length") << " "
	   << std::setw( 7) << G4BestUnit(fTrack->GetPosition().y(),"Length") << " "
	   << std::setw( 7) << G4BestUnit(fTrack->GetPosition().z(),"Length") << " "
	   << std::setw( 7) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy") << " "
	   << std::setw( 7) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy") << " "
	   << std::setw( 7) << G4BestUnit(fStep->GetStepLength(),"Length") << " "
	   << std::setw( 7) << G4BestUnit(fTrack->GetTrackLength(),"Length") << " ";
	  
    
    if( fTrack->GetNextVolume() != 0 ) { 
      G4cout << std::setw(11) << fTrack->GetNextVolume()->GetName() <<" ";
    } else {
      G4cout << std::setw(11) << "OutOfWorld";
    }
    
    if(fStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL){
      G4cout << std::setw(10) << fStep->GetPostStepPoint()->GetProcessDefinedStep()
	->GetProcessName();
    } else {
      G4cout << "User Limit";
    }
    
    G4cout << G4endl;
    
    if( verboseLevel == 2 ){
      G4int tN2ndariesTot = fN2ndariesAtRestDoIt +
	                    fN2ndariesAlongStepDoIt +
	                    fN2ndariesPostStepDoIt;
      if(tN2ndariesTot>0){
	G4cout << "    :----- List of 2ndaries - "
	       << "#SpawnInStep=" << std::setw(3) << tN2ndariesTot 
	       << "(Rest="  << std::setw(2) << fN2ndariesAtRestDoIt
	       << ",Along=" << std::setw(2) << fN2ndariesAlongStepDoIt
	       << ",Post="  << std::setw(2) << fN2ndariesPostStepDoIt
	       << "), "
	       << "#SpawnTotal=" << std::setw(3) << (*fSecondary).size()
	       << " ---------------"
	       << G4endl;
	
	for(size_t lp1=(*fSecondary).size()-tN2ndariesTot; 
	    lp1<(*fSecondary).size(); lp1++){
	  G4cout << "    : "
		 << std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length")
		 << std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length")
		 << std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length")
		 << std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy")
		 << std::setw(10)
		 << (*fSecondary)[lp1]->GetDefinition()->GetParticleName();
	  G4cout << G4endl;
	}
              
	G4cout << "    :-----------------------------"
	       << "----------------------------------"
	       << "-- EndOf2ndaries Info ---------------"
	       << G4endl;
      }
    }
    
  }
  G4cout.precision(prec);
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicrobeamSteppingVerbose::TrackingStarted()
{
  CopyState();
G4int prec = G4cout.precision(3);
  if( verboseLevel > 0 ){

    G4cout << std::setw( 5) << "Step#"      << " "
           << std::setw( 6) << "X"          << "    "
	   << std::setw( 6) << "Y"          << "    "  
	   << std::setw( 6) << "Z"          << "    "
	   << std::setw( 9) << "KineE"      << " "
	   << std::setw( 9) << "dEStep"     << " "  
	   << std::setw(10) << "StepLeng"  
	   << std::setw(10) << "TrakLeng"
	   << std::setw(10) << "NextVolu"
	   << std::setw(10) << "Process"    << G4endl;	     

    G4cout << std::setw( 5) << fTrack->GetCurrentStepNumber() << " "
	   << std::setw( 6) << G4BestUnit(fTrack->GetPosition().x(),"Length")
	   << std::setw( 6) << G4BestUnit(fTrack->GetPosition().y(),"Length")
	   << std::setw( 6) << G4BestUnit(fTrack->GetPosition().z(),"Length")
	   << std::setw( 6) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
	   << std::setw( 6) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
	   << std::setw( 6) << G4BestUnit(fStep->GetStepLength(),"Length")
	   << std::setw( 6) << G4BestUnit(fTrack->GetTrackLength(),"Length");

    if(fTrack->GetNextVolume()){
      G4cout << std::setw(11) << fTrack->GetNextVolume()->GetName() << " ";
    } else {
      G4cout << std::setw(11) << "OutOfWorld" << " ";
    }
    G4cout << std::setw(10) << "initStep" << G4endl;
  }
  G4cout.precision(prec);
  G4cout<< "exit MicrobeamSteppingVerbose::TrackingStarted()   " <<G4endl;
}
