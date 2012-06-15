#include "G4ITTrackingInteractivity.hh"
#include "G4Track.hh"
#include "G4String.hh"

void G4ITTrackingInteractivity::TrackBanner(G4Track* track, const G4String& message)
{
    G4cout << G4endl;
    G4cout << "*******************************************************"
           << "**************************************************"
           << G4endl;
    if(message != "")
        G4cout << message ;
    G4cout << " * G4Track Information: "
           << "   Particle : " << track->GetDefinition()->GetParticleName()
           << ","
           << "   Track ID : " << track->GetTrackID()
           << ","
           << "   Parent ID : " << track->GetParentID()
           << G4endl;
    G4cout << "*******************************************************"
           << "**************************************************"
           << G4endl;
    G4cout << G4endl;
}

