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
// $Id: Mars01SteppingAction.cc,v 1.1 2001-12-13 14:58:48 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 

#include "Mars01SteppingAction.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4TrackVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ios.hh"
#include "g4std/iomanip"

Mars01SteppingAction::Mars01SteppingAction()
{
  nevent = 0;
}

Mars01SteppingAction::~Mars01SteppingAction()
{;}

void Mars01SteppingAction::UserSteppingAction(const G4Step* theStep)
{
  G4SteppingManager * SM = fpSteppingManager;
  G4Track * theTrack = theStep->GetTrack();
  
  G4cout << "==============================================" << G4endl;
  G4cout << " Event " << nevent << G4endl;
  G4cout << "==============================================" << G4endl;
  // check if it is alive
  if(theTrack->GetTrackStatus()==fAlive) { return; }
  
  G4cout << G4std::setw( 5) << "#Step#" << " "
	 << G4std::setw( 9) << "X(mm)" << " "
	 << G4std::setw( 9) << "Y(mm)" << " "
	 << G4std::setw( 9) << "Z(mm)" << " "
	 << G4std::setw( 9) << "KineE(MeV)" << " "
	 << G4std::setw( 9) << "dE(MeV)" << " "
	 << G4std::setw( 9) << "StepLeng" << " "
	 << G4std::setw( 9) << "TrackLeng" << " "
	 << G4std::setw(10) << "ProcName" << G4endl;
  G4cout.precision(3);
  G4cout << G4std::setw( 5) << theTrack->GetCurrentStepNumber() << " "
	 << G4std::setw( 9) << theTrack->GetPosition().x() / mm << " "
	 << G4std::setw( 9) << theTrack->GetPosition().y() / mm << " "
	 << G4std::setw( 9) << theTrack->GetPosition().z() / mm << " "
	 << G4std::setw( 9) << theTrack->GetKineticEnergy() / MeV << " "
	 << G4std::setw( 9) << theStep->GetTotalEnergyDeposit() /MeV << " "
	 << G4std::setw( 9) << theStep->GetStepLength() / mm << " "
	 << G4std::setw( 9) << theTrack->GetTrackLength() / mm << " ";
  if(theStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL){
    G4cout << theStep->GetPostStepPoint()->GetProcessDefinedStep()
      ->GetProcessName();
  } else {
    G4cout << "User Limit";
  }
  G4cout << "----------------------------------------------------" << G4endl;
  G4TrackVector* fSecondary = SM->GetfSecondary();
  G4cout << "   -- List of secondaries generated : "
	 << "(x,y,z,kE,t,PID) --" << G4endl;
  for( G4int lp1=0;lp1<(*fSecondary).size(); lp1++){
    G4cout << "      "
	   << G4std::setw( 9)
	   << (*fSecondary)[lp1]->GetPosition().x() / mm << " "
	   << G4std::setw( 9)
	   << (*fSecondary)[lp1]->GetPosition().y() / mm << " "
	   << G4std::setw( 9)
	   << (*fSecondary)[lp1]->GetPosition().z() / mm << " "
	   << G4std::setw( 9)
	   << (*fSecondary)[lp1]->GetKineticEnergy() / MeV << " "
	   << G4std::setw( 9)
	   << (*fSecondary)[lp1]->GetGlobalTime() / ns << " "
	   << G4std::setw(18)
	   << (*fSecondary)[lp1]->GetDefinition()
      ->GetParticleName();
    G4cout << G4endl;
  }

}

