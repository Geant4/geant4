// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02SteppingVerbose.cc,v 1.2 1999-12-15 14:49:22 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
// ExN02SteppingVerbose.cc
//
// Description:
//    Implementation of  the ExN02SteppingVerbose class
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

#include "ExN02SteppingVerbose.hh"
#include "G4SteppingManager.hh"

#include "G4UnitsTable.hh"

////////////////////////////////////////////////
ExN02SteppingVerbose::ExN02SteppingVerbose()
////////////////////////////////////////////////
{
}

//////////////////////////////////////////////////
ExN02SteppingVerbose::~ExN02SteppingVerbose()
//////////////////////////////////////////////////
{
}

/////////////////////////////////////////
void ExN02SteppingVerbose::StepInfo()
/////////////////////////////////////////
{
  CopyState();
  
  G4int prec = G4cout.precision(3);

  if( verboseLevel >= 1 ){
    if( verboseLevel >= 4 ) VerboseTrack();
    if( verboseLevel >= 3 ){
      G4cout << G4endl;    
      G4cout << G4std::setw( 5) << "#Step#"     << " "
	     << G4std::setw( 6) << "X"          << "    "
	     << G4std::setw( 6) << "Y"          << "    "  
	     << G4std::setw( 6) << "Z"          << "    "
	     << G4std::setw( 9) << "KineE"      << " "
	     << G4std::setw( 9) << "dEStep"     << " "  
	     << G4std::setw(10) << "StepLeng"     
	     << G4std::setw(10) << "TrakLeng" 
	     << G4std::setw(10) << "NextVolu" 
	     << G4std::setw(10) << "Process"   << G4endl;	          
    }

    G4cout << G4std::setw( 5) << fTrack->GetCurrentStepNumber() << " "
	   << G4std::setw( 6) << G4BestUnit(fTrack->GetPosition().x(),"Length")
	   << G4std::setw( 6) << G4BestUnit(fTrack->GetPosition().y(),"Length")
	   << G4std::setw( 6) << G4BestUnit(fTrack->GetPosition().z(),"Length")
	   << G4std::setw( 6) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
	   << G4std::setw( 6) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
	   << G4std::setw( 6) << G4BestUnit(fStep->GetStepLength(),"Length")
	   << G4std::setw( 6) << G4BestUnit(fTrack->GetTrackLength(),"Length");

    // if( fStepStatus != fWorldBoundary){ 
    if( fTrack->GetNextVolume() != 0 ) { 
      G4cout << G4std::setw(10) << fTrack->GetNextVolume()->GetName();
    } else {
      G4cout << G4std::setw(10) << "OutOfWorld";
    }

    if(fStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL){
      G4cout << G4std::setw(10) << fStep->GetPostStepPoint()->GetProcessDefinedStep()
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
	       << "#SpawnInStep=" << G4std::setw(3) << tN2ndariesTot 
	       << "(Rest="  << G4std::setw(2) << fN2ndariesAtRestDoIt
	       << ",Along=" << G4std::setw(2) << fN2ndariesAlongStepDoIt
	       << ",Post="  << G4std::setw(2) << fN2ndariesPostStepDoIt
	       << "), "
	       << "#SpawnTotal=" << G4std::setw(3) << (*fSecondary).entries()
	       << " ---------------"
	       << G4endl;

	for(G4int lp1=(*fSecondary).entries()-tN2ndariesTot; 
                        lp1<(*fSecondary).entries(); lp1++){
	  G4cout << "    : "
		 << G4std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length")
		 << G4std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length")
		 << G4std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length")
		 << G4std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy")
		 << G4std::setw(10)
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

////////////////////////////////////////////////
void ExN02SteppingVerbose::TrackingStarted()
////////////////////////////////////////////////
{

  CopyState();
G4int prec = G4cout.precision(3);
  if( verboseLevel > 0 ){

    G4cout << G4std::setw( 5) << "Step#"      << " "
           << G4std::setw( 6) << "X"          << "    "
	   << G4std::setw( 6) << "Y"          << "    "  
	   << G4std::setw( 6) << "Z"          << "    "
	   << G4std::setw( 9) << "KineE"      << " "
	   << G4std::setw( 9) << "dEStep"     << " "  
	   << G4std::setw(10) << "StepLeng"  
	   << G4std::setw(10) << "TrakLeng"
	   << G4std::setw(10) << "NextVolu"
	   << G4std::setw(10) << "Process"    << G4endl;	     

    G4cout << G4std::setw( 5) << fTrack->GetCurrentStepNumber() << " "
	   << G4std::setw( 6) << G4BestUnit(fTrack->GetPosition().x(),"Length")
	   << G4std::setw( 6) << G4BestUnit(fTrack->GetPosition().y(),"Length")
	   << G4std::setw( 6) << G4BestUnit(fTrack->GetPosition().z(),"Length")
	   << G4std::setw( 6) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
	   << G4std::setw( 6) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
	   << G4std::setw( 6) << G4BestUnit(fStep->GetStepLength(),"Length")
	   << G4std::setw( 6) << G4BestUnit(fTrack->GetTrackLength(),"Length");

    if(fTrack->GetNextVolume()){
      G4cout << G4std::setw(10) << fTrack->GetNextVolume()->GetName() << " ";
    } else {
      G4cout << G4std::setw(10) << "OutOfWorld" << " ";
    }
    G4cout << G4std::setw(10) << "initStep" << G4endl;
  }
  G4cout.precision(prec);
}
