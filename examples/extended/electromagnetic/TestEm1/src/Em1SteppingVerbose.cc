// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em1SteppingVerbose.cc,v 1.1 1999-11-10 17:06:43 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
// Em1SteppingVerbose.cc
//
// Description:
//    Implementation of  the Em1SteppingVerbose class
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

#include "Em1SteppingVerbose.hh"
#include "G4SteppingManager.hh"

#include "G4UnitsTable.hh"

////////////////////////////////////////////////
Em1SteppingVerbose::Em1SteppingVerbose()
////////////////////////////////////////////////
{
}

//////////////////////////////////////////////////
Em1SteppingVerbose::~Em1SteppingVerbose()
//////////////////////////////////////////////////
{
}

/////////////////////////////////////////
void Em1SteppingVerbose::StepInfo()
/////////////////////////////////////////
{
  CopyState();
  
  G4int prec = G4cout.precision(3);

  if( verboseLevel >= 1 ){
    if( verboseLevel >= 4 ) VerboseTrack();
    if( verboseLevel >= 3 ){
      G4cout << endl;    
      G4cout << setw( 5) << "#Step#"     << " "
	     << setw( 6) << "X"          << "    "
	     << setw( 6) << "Y"          << "    "  
	     << setw( 6) << "Z"          << "    "
	     << setw( 9) << "KineE"      << " "
	     << setw( 9) << "dEStep"     << " "  
	     << setw(10) << "StepLeng"     
	     << setw(10) << "TrakLeng" 
	     << setw(10) << "NextVolu" 
	     << setw(10) << "Process"   << endl;	          
    }

    G4cout << setw( 5) << fTrack->GetCurrentStepNumber() << " "
	   << setw( 6) << G4BestUnit(fTrack->GetPosition().x(),"Length")
	   << setw( 6) << G4BestUnit(fTrack->GetPosition().y(),"Length")
	   << setw( 6) << G4BestUnit(fTrack->GetPosition().z(),"Length")
	   << setw( 6) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
	   << setw( 6) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
	   << setw( 6) << G4BestUnit(fStep->GetStepLength(),"Length")
	   << setw( 6) << G4BestUnit(fTrack->GetTrackLength(),"Length");

    // if( fStepStatus != fWorldBoundary){ 
    if( fTrack->GetNextVolume() != 0 ) { 
      G4cout << setw(10) << fTrack->GetNextVolume()->GetName();
    } else {
      G4cout << setw(10) << "OutOfWorld";
    }

    if(fStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL){
      G4cout << setw(10) << fStep->GetPostStepPoint()->GetProcessDefinedStep()
	->GetProcessName();
    } else {
      G4cout << "User Limit";
    }

    G4cout << endl;

    if( verboseLevel == 2 ){
      G4int tN2ndariesTot = fN2ndariesAtRestDoIt +
	                    fN2ndariesAlongStepDoIt +
	                    fN2ndariesPostStepDoIt;
      if(tN2ndariesTot>0){
	G4cout << "    :----- List of 2ndaries - "
	       << "#SpawnInStep=" << setw(3) << tN2ndariesTot 
	       << "(Rest="  << setw(2) << fN2ndariesAtRestDoIt
	       << ",Along=" << setw(2) << fN2ndariesAlongStepDoIt
	       << ",Post="  << setw(2) << fN2ndariesPostStepDoIt
	       << "), "
	       << "#SpawnTotal=" << setw(3) << (*fSecondary).entries()
	       << " ---------------"
	       << endl;

	for(G4int lp1=(*fSecondary).entries()-tN2ndariesTot; 
                        lp1<(*fSecondary).entries(); lp1++){
	  G4cout << "    : "
		 << setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length")
		 << setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length")
		 << setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length")
		 << setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy")
		 << setw(10)
		 << (*fSecondary)[lp1]->GetDefinition()->GetParticleName();
	  G4cout << endl;
	}
              
	G4cout << "    :-----------------------------"
	       << "----------------------------------"
	       << "-- EndOf2ndaries Info ---------------"
	       << endl;
      }
    }
    
  }
  G4cout.precision(prec);
}

////////////////////////////////////////////////
void Em1SteppingVerbose::TrackingStarted()
////////////////////////////////////////////////
{

  CopyState();
G4int prec = G4cout.precision(3);
  if( verboseLevel > 0 ){

    G4cout << setw( 5) << "Step#"      << " "
           << setw( 6) << "X"          << "    "
	   << setw( 6) << "Y"          << "    "  
	   << setw( 6) << "Z"          << "    "
	   << setw( 9) << "KineE"      << " "
	   << setw( 9) << "dEStep"     << " "  
	   << setw(10) << "StepLeng"  
	   << setw(10) << "TrakLeng"
	   << setw(10) << "NextVolu"
	   << setw(10) << "Process"    << endl;	     

    G4cout << setw( 5) << fTrack->GetCurrentStepNumber() << " "
	   << setw( 6) << G4BestUnit(fTrack->GetPosition().x(),"Length")
	   << setw( 6) << G4BestUnit(fTrack->GetPosition().y(),"Length")
	   << setw( 6) << G4BestUnit(fTrack->GetPosition().z(),"Length")
	   << setw( 6) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
	   << setw( 6) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
	   << setw( 6) << G4BestUnit(fStep->GetStepLength(),"Length")
	   << setw( 6) << G4BestUnit(fTrack->GetTrackLength(),"Length");

    if(fTrack->GetNextVolume()){
      G4cout << setw(10) << fTrack->GetNextVolume()->GetName() << " ";
    } else {
      G4cout << setw(10) << "OutOfWorld" << " ";
    }
    G4cout << setw(10) << "initStep" << endl;
  }
  G4cout.precision(prec);
}
