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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
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
#include "G4ChiralInvariantPhaseSpace.hh"

class G4ProtonAntiProtonAtRestChips : public G4VRestProcess
{
  private:
  // hide assignment operator as private 
      G4ProtonAntiProtonAtRestChips& operator=(const G4ProtonAntiProtonAtRestChips &right);
      G4ProtonAntiProtonAtRestChips(const G4ProtonAntiProtonAtRestChips& );
   
  public:
 
     G4ProtonAntiProtonAtRestChips(const G4String& processName ="AntiProtonAnnihilationAtRest")
      : G4VRestProcess (processName) {}
 
    ~G4ProtonAntiProtonAtRestChips() {}

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
  // Create target
  G4Element * theTarget = theSelector.GetElement(aTrack.GetMaterial());
  G4Nucleus aTargetNucleus(theTarget->GetN() ,theTarget->GetZ());

  // Check model validity - note this will be a sub-branch in the ordinary stopping @@@@@@
  // in the long haul. @@@@@@
  if(aTrack.GetDynamicParticle()->GetDefinition() != G4AntiProton::AntiProton())
  {
    G4Exception("Calling G4ProtonAntiProtonAtRestChips with particle other than p-bar!!!");
  }
  if(aTargetNucleus.GetZ() != 1)
  {
    G4Exception("Calling G4ProtonAntiProtonAtRestChips for target other than Hydrogen!!!");
  }
  
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
