
#include "Tst06SteppingAction.hh"
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

Tst06SteppingAction::Tst06SteppingAction()
{;}

Tst06SteppingAction::~Tst06SteppingAction()
{;}

void Tst06SteppingAction::UserSteppingAction(const G4Step * theStep)
{
  G4SteppingManager * SM = fpSteppingManager;
  G4Track * theTrack = theStep->GetTrack();

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
    G4cout << G4endl;

   G4TrackVector* fSecondary = SM->GetfSecondary();
       G4cout << "   -- List of secondaries generated : "
         << "(x,y,z,kE,t,PID) --" << G4endl;
       for( G4int lp1=0;lp1<(*fSecondary).entries(); lp1++){
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



















