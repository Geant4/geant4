#ifndef G4PionMinusNuclearAtRestChips_h
#define G4PionMinusNuclearAtRestChips_h

#include "globals.hh"
#include "G4VRestProcess.hh"
#include "G4StopElementSelector.hh"
#include "G4PionMinus.hh"
#include "G4ChiralInvariantPhaseSpace.hh"

class G4PionMinusNuclearAtRestChips : public G4VRestProcess
{
  private:
  // hide assignment operator as private 
      G4PionMinusNuclearAtRestChips& operator=(const G4PionMinusNuclearAtRestChips &right);
      G4PionMinusNuclearAtRestChips(const G4PionMinusNuclearAtRestChips& );
   
  public:
 
     G4PionMinusNuclearAtRestChips(const G4String& processName ="PionMinusAnnihilationAtRest")
      : G4VRestProcess (processName) {}
 
    ~G4PionMinusNuclearAtRestChips() {}

     G4bool IsApplicable(const G4ParticleDefinition& aParticle)
     {
       return ( &aParticle == G4PionMinus::PionMinusDefinition() );
     }

  // null physics table
     void BuildPhysicsTable(const G4ParticleDefinition&){}

     G4double AtRestGetPhysicalInteractionLength(const G4Track&track,
						 G4ForceCondition*condition);

  // zero mean lifetime
     G4double GetMeanLifeTime(const G4Track& aTrack,
			      G4ForceCondition* condition) {return 0.0;}

     G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&); 

  private:
    G4ChiralInvariantPhaseSpace theModel;
    G4StopElementSelector theSelector; // Assume identical laws as for muons
};

inline
G4VParticleChange * G4PionMinusNuclearAtRestChips::
AtRestDoIt(const G4Track& aTrack, const G4Step&aStep)
{
  if(aTrack.GetDynamicParticle()->GetDefinition() != G4PionMinus::PionMinus())
  {
    G4Exception("Calling G4PionMinusNuclearAtRestChips with particle other than pi-!!!");
  }
  
  // Create target
  G4Element * theTarget = theSelector.GetElement(aTrack.GetMaterial());
  G4Nucleus aTargetNucleus(theTarget->GetN() ,theTarget->GetZ());
  
  // Call chips
  return theModel.ApplyYourself(aTrack, aTargetNucleus);
}

G4double G4PionMinusNuclearAtRestChips::
AtRestGetPhysicalInteractionLength(const G4Track&track,
				   G4ForceCondition*condition)
{
  ResetNumberOfInteractionLengthLeft();
  *condition = NotForced;
  currentInteractionLength = GetMeanLifeTime(track, condition);
  if ((currentInteractionLength <0.0) || (verboseLevel>2))
  {
    G4cout << "G4PionMinusNuclearAtRestChips::AtRestGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" <<G4endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " << track.GetMaterial()->GetName() <<G4endl;
    G4cout << "MeanLifeTime = " << currentInteractionLength/ns << "[ns]" <<G4endl;
  }
  return theNumberOfInteractionLengthLeft * currentInteractionLength;
}
#endif
