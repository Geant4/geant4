#ifndef G4ProtonAntiProtonAtRestChips_h
#define G4ProtonAntiProtonAtRestChips_h

#include "globals.hh"
#include "G4VRestProcess.hh"
#include "G4ParticleTable.hh"
#include "G4Quasmon.hh"
#include "G4QHadronVector.hh"
#include "G4ParticleChange.hh"
#include "G4LorentzVector.hh"
#include "G4DynamicParticle.hh"
#include "G4IonTable.hh"
#include "G4Neutron.hh"
#include "G4StopElementSelector.hh"

class G4ProtonAntiProtonAtRestChips : public G4VRestProcess
{
  private:
  // hide assignment operator as private 
      G4AntiProtonAnnihilationAtRest& operator=(const G4AntiProtonAnnihilationAtRest &right);
      G4AntiProtonAnnihilationAtRest(const G4AntiProtonAnnihilationAtRest& );
   
  public:
 
     G4AntiProtonAnnihilationAtRest(const G4String& processName ="AntiProtonAnnihilationAtRest")
      : G4VRestProcess (processName) {}
 
    ~G4AntiProtonAnnihilationAtRest() {}

     G4bool IsApplicable(const G4ParticleDefinition& aParticle)
     {
       return ( &aParticle == G4AntiProton::AntiProtonDefinition() );
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
G4VParticleChange * G4ProtonAntiProtonAtRestChips::
AtRestDoIt(const G4Track& aTrack, const G4Step&aStep)
{
  if(aTrack.GetDynamicParticle()->GetDefinition() != G4AntiProton::AntiProton())
  {
    G4Exception("Calling G4ProtonAntiProtonAtRestChips with particle other than p-bar!!!");
  }
  if(aTargetNucleus.GetZ() != 1)
  {
    G4Exception("Calling G4ProtonAntiProtonAtRestChips for target other than Hydrogen!!!");
  }
    
  // Create target
  G4Element * theTarget = theSelector.GetElement(aTrack.GetMaterial());
  G4Nucleus aTargetNucleus(theTarget.GetN() ,theTarget.GetZ());
  
  // Call chips
  return theModel.ApplyYourself(aTrack, aTargetNucleus);
}

G4double G4ProtonAntiProtonAtRestChips::
AtRestGetPhysicalInteractionLength(const G4Track&track,
				   G4ForceCondition*condition)
{
  ResetNumberOfInteractionLengthLeft();
  *condition = NotForced;
  currentInteractionLength = GetMeanLifeTime(track, condition);
  if ((currentInteractionLength <0.0) || (verboseLevel>2))
  {
    G4cout << "G4ProtonAntiProtonAtRestChips::AtRestGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" <<G4endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " << track.GetMaterial()->GetName() <<G4endl;
    G4cout << "MeanLifeTime = " << currentInteractionLength/ns << "[ns]" <<G4endl;
  }
  return theNumberOfInteractionLengthLeft * currentInteractionLength;
}
#endif
