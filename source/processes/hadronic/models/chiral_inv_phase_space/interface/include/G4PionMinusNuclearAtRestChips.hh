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
// 17.02.2009 M.Kossov, now it is recommended to use the G4QCaptureAtRest process
//
#ifndef G4PionMinusNuclearAtRestChips_h
#define G4PionMinusNuclearAtRestChips_h

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4VRestProcess.hh"
#include "G4StopElementSelector.hh"
#include "G4PionMinus.hh"
#include "G4ChiralInvariantPhaseSpace.hh"
#include "G4HadronicProcessType.hh"
#include "G4HadronicDeprecate.hh"


class G4PionMinusNuclearAtRestChips : public G4VRestProcess
{
  private:
  // hide assignment operator as private 
      G4PionMinusNuclearAtRestChips& operator=(const G4PionMinusNuclearAtRestChips &right);
      G4PionMinusNuclearAtRestChips(const G4PionMinusNuclearAtRestChips& );
   
  public:
 
     G4PionMinusNuclearAtRestChips(const G4String& processName ="PionMinusCaptureAtRest")
      : G4VRestProcess (processName, fHadronic) 
     {
       G4HadronicDeprecate("G4PionMinusNuclearAtRestChips");
       SetProcessSubType(fHadronAtRest);
     }
 
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
    throw G4HadronicException(__FILE__, __LINE__,
                  "Calling G4PionMinusNuclearAtRestChips with particle other than pi-!!!");
  }
  
  // Create target
  G4Element* theTarget = theSelector.GetElement(aTrack.GetMaterial());
  G4Nucleus aTargetNucleus(theTarget->GetA_asInt(), theTarget->GetZ_asInt());
  
  // Call chips
  return theModel.ApplyYourself(aTrack, aTargetNucleus);
}

G4double G4PionMinusNuclearAtRestChips::
AtRestGetPhysicalInteractionLength(const G4Track&track, G4ForceCondition*condition)
{
  ResetNumberOfInteractionLengthLeft();
  *condition = NotForced;
  currentInteractionLength = GetMeanLifeTime(track, condition);
#ifdef CHIPSdebug
  if ((currentInteractionLength <0.0) || (verboseLevel>2))
  {
    G4cout << "G4PionMinusNuclearAtRestChips::AtRestGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" <<G4endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " << track.GetMaterial()->GetName() <<G4endl;
    G4cout << "MeanLifeTime = " << currentInteractionLength/CLHEP::ns << "[ns]" <<G4endl;
  }
#endif
  return theNumberOfInteractionLengthLeft * currentInteractionLength;
}
#endif
