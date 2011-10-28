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
// Identical to G4VRestProcess with dependency from G4VITProcess
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// Author: Mathieu Karamitros (kara@cenbg.in2p3.fr)
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------


#ifndef G4VITRestProcess_h
#define G4VITRestProcess_h 1

#include "G4VITProcess.hh"

class G4VITRestProcess : public G4VITProcess
{
  //  Abstract class which defines the public behavior of
  //  physics interactions at rest.

  public:
      G4VITRestProcess(const G4String&  ,
		     G4ProcessType   aType = fNotDefined );
      G4VITRestProcess(const G4VITRestProcess& );

      virtual ~G4VITRestProcess();

  public:   //  with description
      virtual G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4ForceCondition* condition
			    );

      virtual G4VParticleChange* AtRestDoIt(
			     const G4Track& ,
			     const G4Step&
			    );

     //  no operation in  PostStepDoIt and  AlongStepDoIt
      virtual G4double AlongStepGetPhysicalInteractionLength(
                             const G4Track&,
			     G4double  ,
			     G4double  ,
			     G4double& ,
	                     G4GPILSelection*
                           ){ return -1.0; }

      virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& ,
			     G4double   ,
			     G4ForceCondition*
                            ) { return -1.0; }

     //  no operation in  PostStepDoIt and  AlongStepDoIt
      virtual G4VParticleChange* PostStepDoIt(
			     const G4Track& ,
			     const G4Step&
                            ) {return 0;}

      virtual G4VParticleChange* AlongStepDoIt(
			     const G4Track& ,
			     const G4Step&
                            ) {return 0;}

  protected: //  with description

      virtual G4double GetMeanLifeTime(const G4Track& aTrack,G4ForceCondition* condition)=0;
      //  Calculates the mean life-time (i.e. for decays) of the
      //  particle at rest due to the occurence of the given process,
      //  or converts the probability of interaction (i.e. for
      //  annihilation) into the life-time of the particle for the
      //  occurence of the given process.

 protected:
  // hide default constructor and assignment operator as private
      G4VITRestProcess();
      G4VITRestProcess & operator=(const G4VITRestProcess &right);
};

// -----------------------------------------
//  inlined function members implementation
// -----------------------------------------
inline G4double G4VITRestProcess::AtRestGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4ForceCondition* condition
			    )
{
  // beggining of tracking
  ResetNumberOfInteractionLengthLeft();

  // condition is set to "Not Forced"
  *condition = NotForced;

  // get mean life time
  fState->currentInteractionLength = GetMeanLifeTime(track, condition);

#ifdef G4VERBOSE
 if ((fState->currentInteractionLength <0.0) || (verboseLevel>2)){
    G4cout << "G4VITRestProcess::AtRestGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" <<G4endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " << track.GetMaterial()->GetName() <<G4endl;
    G4cout << "MeanLifeTime = " << fState->currentInteractionLength/ns << "[ns]" <<G4endl;
  }
#endif

  return (fState->theNumberOfInteractionLengthLeft) * (fState->currentInteractionLength);
}


inline G4VParticleChange* G4VITRestProcess::AtRestDoIt(
			     const G4Track&,
			     const G4Step&
			    )
{
//  clear NumberOfInteractionLengthLeft
    ClearNumberOfInteractionLengthLeft();

    return pParticleChange;
}


#endif




