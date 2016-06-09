//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: AnaEx01SteppingVerbose.cc,v 1.6 2003/06/20 14:55:45 gbarrand Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
//
//---------------------------------------------------------------
//
// AnaEx01SteppingVerbose.cc
//
// Description:
//    Implementation of  the AnaEx01SteppingVerbose class
//
//---------------------------------------------------------------

#include "AnaEx01SteppingVerbose.hh"
#include "G4SteppingManager.hh"

#include "G4UnitsTable.hh"

////////////////////////////////////////////////
AnaEx01SteppingVerbose::AnaEx01SteppingVerbose()
////////////////////////////////////////////////
{
}

//////////////////////////////////////////////////
AnaEx01SteppingVerbose::~AnaEx01SteppingVerbose()
//////////////////////////////////////////////////
{
}

/////////////////////////////////////////
void AnaEx01SteppingVerbose::StepInfo()
/////////////////////////////////////////
{
  CopyState();
  
  G4int prec = G4cout.precision(3);

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
	   << std::setw( 6) << G4BestUnit(fTrack->GetPosition().x(),"Length")
	   << std::setw( 6) << G4BestUnit(fTrack->GetPosition().y(),"Length")
	   << std::setw( 6) << G4BestUnit(fTrack->GetPosition().z(),"Length")
	   << std::setw( 6) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
	   << std::setw( 6) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
	   << std::setw( 6) << G4BestUnit(fStep->GetStepLength(),"Length")
	   << std::setw( 6) << G4BestUnit(fTrack->GetTrackLength(),"Length");

    // if( fStepStatus != fWorldBoundary){ 
    if( fTrack->GetNextVolume() != 0 ) { 
      G4cout << std::setw(10) << fTrack->GetNextVolume()->GetName();
    } else {
      G4cout << std::setw(10) << "OutOfWorld";
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

	for(G4int lp1=(*fSecondary).size()-tN2ndariesTot; 
                        lp1<(int)(*fSecondary).size(); lp1++){
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

////////////////////////////////////////////////
void AnaEx01SteppingVerbose::TrackingStarted()
////////////////////////////////////////////////
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
      G4cout << std::setw(10) << fTrack->GetNextVolume()->GetName() << " ";
    } else {
      G4cout << std::setw(10) << "OutOfWorld" << " ";
    }
    G4cout << std::setw(10) << "initStep" << G4endl;
  }
  G4cout.precision(prec);
}
