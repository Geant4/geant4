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

#include "PhotInSteppingAction.hh"
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
#include <iomanip>

PhotInSteppingAction::PhotInSteppingAction() {;}

PhotInSteppingAction::~PhotInSteppingAction() {;}

void PhotInSteppingAction::UserSteppingAction(const G4Step* theStep)
{
  G4SteppingManager* SM = fpSteppingManager;
  G4Track* theTrack = theStep->GetTrack();
  G4cout<<"PhotInSteppingAction::UserSteppingAction: *** Material = "
        <<theTrack->GetMaterial()->GetName()<<" ***"<<G4endl;
  G4TrackVector* fSecondary = SM->GetfSecondary();
  G4int nSec = (*fSecondary).size(); // #of secondaries
  G4cout << std::setw( 5) << "#Step#"
         << std::setw( 9) << "X(mm)" << " "
         << std::setw( 9) << "Y(mm)" << " "
         << std::setw( 9) << "Z(mm)" << " "
         << std::setw( 9) << "KineE(MeV)"
         << std::setw( 9) << "dE(MeV)" << " "
         << std::setw( 9) << "StepLeng" << " "
         << std::setw( 9) << "TrackLeng" << " "
         << std::setw( 9) << "Particle" << "  "
         << std::setw( 9) << "ProcName" << G4endl;
  G4cout.precision(3);
  G4cout << std::setw( 5) << theTrack->GetCurrentStepNumber() << " "
         << std::setw( 9) << theTrack->GetPosition().x() / mm << " "
         << std::setw( 9) << theTrack->GetPosition().y() / mm << " "
         << std::setw( 9) << theTrack->GetPosition().z() / mm << " "
         << std::setw( 9) << theTrack->GetKineticEnergy() / MeV << " "
         << std::setw( 9) << theStep->GetTotalEnergyDeposit() /MeV << " "
         << std::setw( 9) << theStep->GetStepLength() / mm << " "
         << std::setw( 9) << theTrack->GetTrackLength() / mm << " "
         << std::setw( 9) << theTrack->GetDefinition()->GetParticleName()<< "   ";
  if(theStep->GetPostStepPoint()->GetProcessDefinedStep() != 0)
            G4cout<<theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  else G4cout<<"User Limit";
  G4cout<<G4endl;

  // check if it is alive and quit if no secondaries
  if(theTrack->GetTrackStatus()==fAlive)
  {
    G4cout<<"PhotInSteppingAction::UserSteppingAction:-TRACK IS ALIVE->, N="<<nSec<<G4endl;
    if(!nSec) return;
  }

  G4cout<<"PhotInSteppingAction::UserSteppingAction:Secondaries, N="<<nSec<<" ***"<<G4endl;
  G4cout<< "  "<<std::setw( 9)<<"X(mm)"     
        << ", "<<std::setw( 9)<<"Y(mm)"     
        << ", "<<std::setw( 9)<<"Z(mm)"     
        << ", "<<std::setw( 9)<<"KineE(MeV)"
        << ", "<<std::setw( 9)<<"Time(ns)"  
        << ", "<<std::setw(18)<<"Particle"  <<G4endl;
  if(nSec) for(G4int lp1=0; lp1<nSec; lp1++)
  {
    G4cout<<"  "<<std::setw( 9)<<(*fSecondary)[lp1]->GetPosition().x()  / mm
          <<", "<<std::setw( 9)<<(*fSecondary)[lp1]->GetPosition().y()  / mm
          <<", "<<std::setw( 9)<<(*fSecondary)[lp1]->GetPosition().z()  / mm
          <<", "<<std::setw( 9)<<(*fSecondary)[lp1]->GetKineticEnergy() / MeV
          <<", "<<std::setw( 9)<<(*fSecondary)[lp1]->GetGlobalTime()    / ns
          <<", "<<std::setw(18)<<(*fSecondary)[lp1]->GetDefinition()->GetParticleName();
    G4cout<<G4endl;
  }
}

