#ifndef A01TRIGGERS_HH
#define A01TRIGGERS_HH

namespace A01Triggers {

  G4bool PrimaryTrackTrigger(G4Track* track)
  {
    G4cout<<"jane primary track "<<track->GetTrackID()<<G4endl;
    return (track->GetTrackID() == 1);
  }

  G4bool CalorimeterTrigger(const G4Track& track, const G4Step& step)
  {
    G4cout<<"jane cal trigger "<<track.GetVolume()->GetName()<<G4endl;
    return track.GetVolume()->GetName() == "cellPhysical";
  }

  G4bool Hadronic_LeadingParticleBiasing_Trigger(const G4Track& track, const G4Step& step)
  {
    G4cout<<"jane had lead particle trigger "<<track.GetDefinition()->GetParticleName()<<" "<<track.GetKineticEnergy()<<G4endl;
    return (track.GetKineticEnergy() < 5*GeV);
  }
}
#endif
