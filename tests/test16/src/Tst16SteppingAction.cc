
#include "Tst16SteppingAction.hh"
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
#include <iomanip.h>

Tst16SteppingAction::Tst16SteppingAction()
{;}

Tst16SteppingAction::~Tst16SteppingAction()
{;}

void Tst16SteppingAction::UserSteppingAction(const G4Step* theStep)
{
  G4SteppingManager * SM = fpSteppingManager;
  G4Track * theTrack = theStep->GetTrack();

  // check if it is alive
  if(theTrack->GetTrackStatus()==fAlive) { return; }

      G4cout << setw( 5) << "#Step#" << " "
        << setw( 9) << "X(mm)" << " "
          << setw( 9) << "Y(mm)" << " "
            << setw( 9) << "Z(mm)" << " "
              << setw( 9) << "KineE(MeV)" << " "
                << setw( 9) << "dE(MeV)" << " "
                  << setw( 9) << "StepLeng" << " "
                    << setw( 9) << "TrackLeng" << " "
                        << setw(10) << "ProcName" << endl;
    G4cout.precision(3);
    G4cout << setw( 5) << theTrack->GetCurrentStepNumber() << " "
      << setw( 9) << theTrack->GetPosition().x() / mm << " "
        << setw( 9) << theTrack->GetPosition().y() / mm << " "
          << setw( 9) << theTrack->GetPosition().z() / mm << " "
             << setw( 9) << theTrack->GetKineticEnergy() / MeV << " "
              << setw( 9) << theStep->GetTotalEnergyDeposit() /MeV << " "
                << setw( 9) << theStep->GetStepLength() / mm << " "
                  << setw( 9) << theTrack->GetTrackLength() / mm << " ";
    if(theStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL){
      G4cout << theStep->GetPostStepPoint()->GetProcessDefinedStep()
        ->GetProcessName();
    } else {
      G4cout << "User Limit";
    }
    G4cout << endl;

   G4TrackVector* fSecondary = SM->GetfSecondary();
       G4cout << "   -- List of secondaries generated : "
         << "(x,y,z,kE,t,PID) --" << endl;
       for( G4int lp1=0;lp1<(*fSecondary).entries(); lp1++){
         G4cout << "      "
           << setw( 9)
             << (*fSecondary)[lp1]->GetPosition().x() / mm << " "
               << setw( 9)
                 << (*fSecondary)[lp1]->GetPosition().y() / mm << " "
                   << setw( 9)
                     << (*fSecondary)[lp1]->GetPosition().z() / mm << " "
                       << setw( 9)
                         << (*fSecondary)[lp1]->GetKineticEnergy() / MeV << " "
                           << setw( 9)
                             << (*fSecondary)[lp1]->GetGlobalTime() / ns << " "
                               << setw(18)
                                 << (*fSecondary)[lp1]->GetDefinition()
                                   ->GetParticleName();
         G4cout << endl;
       }

}

