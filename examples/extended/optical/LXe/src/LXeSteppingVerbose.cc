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
#include "LXeSteppingVerbose.hh"

#include "G4SteppingManager.hh"
#include "G4UnitsTable.hh"

#include<strstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeSteppingVerbose::LXeSteppingVerbose()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeSteppingVerbose::~LXeSteppingVerbose()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeSteppingVerbose::StepInfo()
{
  CopyState();
  
  G4int prec = G4cout.precision(3);

  if( verboseLevel >= 1 ){
   G4int lengthLen=+G4UnitDefinition::GetUnitsTable()[G4BestUnit(0,"Length").GetIndexOfCategory()]->GetSymbMxLen();
   G4int energyLen=+G4UnitDefinition::GetUnitsTable()[G4BestUnit(0,"Energy").GetIndexOfCategory()]->GetSymbMxLen();

    if( verboseLevel >= 4 ) VerboseTrack();
    if( verboseLevel >= 3 ){
      G4cout << G4endl;    
      G4cout << std::setw( 5) << "Step#"      
           << std::setw( 6) << "X" << std::setw(lengthLen+1)<<" " 
	   << std::setw( 6) << "Y" << std::setw(lengthLen+1)<<" "
	   << std::setw( 6) << "Z" << std::setw(lengthLen+1)<<" "
	   << std::setw( 7) << "KineE" << std::setw(energyLen+1)<<" "
	   << std::setw( 7) << "dEStep"<< std::setw(energyLen+1)<<" "
	   << std::setw(10) << "StepLeng"<< std::setw(lengthLen+1)<<" "
	   << std::setw(10) << "TrakLeng"<< std::setw(lengthLen+1)<<" "
	   << std::setw(10) << "Volume"     
	   << std::setw(10) << "Process"    << G4endl;	     
    }

   //This is a less clean way of outputing the data than normal but it
   //produces nicer colums since G4BestUnit doesnt cooperate too well with
   //field widths. 
    G4cout << std::setw(5) << fTrack->GetCurrentStepNumber();
    strstream out;
    out.precision(3);
    out.rdbuf()->seekpos(0);
    out<<G4BestUnit(fTrack->GetPosition().x(),"Length")<<'\0';
    G4cout<<std::setw(7+lengthLen)<<out.str();
    out.rdbuf()->seekpos(0);
    out<<G4BestUnit(fTrack->GetPosition().y(),"Length")<<'\0';
    G4cout<<std::setw(7+lengthLen)<<out.str();
    out.rdbuf()->seekpos(0);
    out<<G4BestUnit(fTrack->GetPosition().z(),"Length")<<'\0';
    G4cout<<std::setw(7+lengthLen)<<out.str();
    out.rdbuf()->seekpos(0);
    out<<G4BestUnit(fTrack->GetKineticEnergy(),"Energy")<<'\0';
    G4cout<<std::setw(7+energyLen)<<out.str();
    out.rdbuf()->seekpos(0);
    out<<G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")<<'\0';
    G4cout<<std::setw(7+energyLen)<<out.str();
    out.rdbuf()->seekpos(0);
    out<<G4BestUnit(fStep->GetStepLength(),"Length")<<'\0';
    G4cout<<std::setw(10+lengthLen)<<out.str();
    out.rdbuf()->seekpos(0);
    out<<G4BestUnit(fTrack->GetTrackLength(),"Length")<<'\0';
    G4cout<<std::setw(10+lengthLen)<<out.str();
    
    // if( fStepStatus != fWorldBoundary){ 
    if( fTrack->GetNextVolume() != 0 ) { 
      G4cout << "   "<<fTrack->GetVolume()->GetName();
    } else { 
      G4cout << "   OutOfWorld";
    }

    if(fStep->GetPostStepPoint()->GetProcessDefinedStep() != 0){
      G4cout << " "
             << std::setw(10)
	     << fStep->GetPostStepPoint()->GetProcessDefinedStep()
	                                 ->GetProcessName();
    } else {
      G4cout << "   UserLimit";
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeSteppingVerbose::TrackingStarted()
{

  CopyState();
G4int prec = G4cout.precision(3);
 if( verboseLevel > 0 ){
   G4int lengthLen=+G4UnitDefinition::GetUnitsTable()[G4BestUnit(0,"Length").GetIndexOfCategory()]->GetSymbMxLen();
   G4int energyLen=+G4UnitDefinition::GetUnitsTable()[G4BestUnit(0,"Energy").GetIndexOfCategory()]->GetSymbMxLen();

   G4cout << std::setw( 5) << "Step#"      
           << std::setw( 6) << "X" << std::setw(lengthLen+1)<<" " 
	   << std::setw( 6) << "Y" << std::setw(lengthLen+1)<<" "
	   << std::setw( 6) << "Z" << std::setw(lengthLen+1)<<" "
	   << std::setw( 7) << "KineE" << std::setw(energyLen+1)<<" "
	   << std::setw( 7) << "dEStep"<< std::setw(energyLen+1)<<" "
	   << std::setw(10) << "StepLeng"<< std::setw(lengthLen+1)<<" "
	   << std::setw(10) << "TrakLeng"<< std::setw(lengthLen+1)<<" "
	   << std::setw(10) << "Volume"     
	   << std::setw(10) << "Process"    << G4endl;	     

   //This is a less clean way of outputing the data than normal but it
   //produces nicer colums since G4BestUnit doesnt cooperate too well with
   //field widths. 
   G4cout << std::setw( 5) << fTrack->GetCurrentStepNumber();
   strstream out;
   out.precision(3);
   out.rdbuf()->seekpos(0);
   out<<G4BestUnit(fTrack->GetPosition().x(),"Length")<<'\0';
   G4cout<<std::setw(7+lengthLen)<<out.str();
   out.rdbuf()->seekpos(0);
   out<<G4BestUnit(fTrack->GetPosition().y(),"Length")<<'\0';
   G4cout<<std::setw(7+lengthLen)<<out.str();
   out.rdbuf()->seekpos(0);
   out<<G4BestUnit(fTrack->GetPosition().z(),"Length")<<'\0';
   G4cout<<std::setw(7+lengthLen)<<out.str();
   out.rdbuf()->seekpos(0);
   out<<G4BestUnit(fTrack->GetKineticEnergy(),"Energy")<<'\0';
   G4cout<<std::setw(7+energyLen)<<out.str();
   out.rdbuf()->seekpos(0);
   out<<G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")<<'\0';
   G4cout<<std::setw(7+energyLen)<<out.str();
   out.rdbuf()->seekpos(0);
   out<<G4BestUnit(fStep->GetStepLength(),"Length")<<'\0';
   G4cout<<std::setw(10+lengthLen)<<out.str();
   out.rdbuf()->seekpos(0);
   out<<G4BestUnit(fTrack->GetTrackLength(),"Length")<<'\0';
   G4cout<<std::setw(10+lengthLen)<<out.str();


    if(fTrack->GetNextVolume()){
      G4cout << "   "<< fTrack->GetVolume()->GetName();
    } else {
      G4cout << std::setw(11+lengthLen)<< "OutOfWorld";
    }
    G4cout << "  initStep" << G4endl;
  }
  G4cout.precision(prec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

