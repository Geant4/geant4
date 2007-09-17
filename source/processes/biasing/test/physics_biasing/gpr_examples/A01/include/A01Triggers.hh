#ifndef A01TRIGGERS_HH
#define A01TRIGGERS_HH

namespace A01Triggers {

  // Return true if track is the daughter of a primary
  G4bool DaughterOfPrimaryTrigger(G4Track* track)
  {
    return (track->GetParentID() == 1);
  }

  // Return true if in calorimeter volume
  G4bool CalorimeterTrigger(const G4Track& track, const G4Step& step)
  {
    return track.GetVolume()->GetName() == "cellPhysical";
  }

  // Return true if track energy is less than 5 GeV
  G4bool Hadronic_LeadingParticleBiasing_Trigger(const G4Track& track, const G4Step& step)
  {
    return (track.GetKineticEnergy() < 5*GeV);
  }
}

#endif
