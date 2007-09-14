#ifndef A01TRIGGERS_HH
#define A01TRIGGERS_HH

namespace A01Triggers {

  G4bool CalorimeterTrigger(const G4Track& track, const G4Step& step)
  {
    G4cout<<"jane cal trigger "<<track.GetVolume()->GetName()<<G4endl;
    return track.GetVolume()->GetName() == "cellPhysical";
  }
}
#endif
