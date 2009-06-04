//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// (Why only antiproton-proton, when the antiproton-nucleus is made? - M.K.)
// 17.02.2009 M.Kossov, now it is recommended to use the G4QCaptureAtRest process
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
#include "G4HadronicProcessType.hh"

class G4ProtonAntiProtonAtRestChips : public G4VRestProcess
{
  private:
  // hide assignment operator as private 
      G4ProtonAntiProtonAtRestChips& operator=(const G4ProtonAntiProtonAtRestChips &right);
      G4ProtonAntiProtonAtRestChips(const G4ProtonAntiProtonAtRestChips& );
   
  public:
 
     G4ProtonAntiProtonAtRestChips(const G4String& processName=
                                                            "AntiProtonAnnihilationAtRest")
      : G4VRestProcess (processName, fHadronic) 
     {
       SetProcessSubType(fHadronAtRest);
     }
 
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
    throw G4HadronicException(__FILE__, __LINE__,
                "Calling G4ProtonAntiProtonAtRestChips with particle other than p-bar!!!");
  }
  if(aTargetNucleus.GetZ() != 1)
  {
    throw G4HadronicException(__FILE__, __LINE__,
                "Calling G4ProtonAntiProtonAtRestChips for target other than Hydrogen!!!");
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
#ifdef CHIPSdebug
  if ((currentInteractionLength <0.0) || (verboseLevel>2))
  {
    G4cout << "G4ProtonAntiProtonAtRestChips::AtRestGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" <<G4endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " << track.GetMaterial()->GetName() <<G4endl;
    G4cout << "MeanLifeTime = " << currentInteractionLength/ns << "[ns]" <<G4endl;
  }
#endif
  return theNumberOfInteractionLengthLeft * currentInteractionLength;
}
#endif
