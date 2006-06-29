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
//
// $Id: G4VRestContinuousProcess.hh,v 1.6 2006-06-29 21:07:54 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// 
// Class Description 
//  Abstract class which defines the public behavior of
//  discrete physics interactions.
//
// ------------------------------------------------------------
//   New Physics scheme            8  Mar. 1997  H.Kurahige
// ------------------------------------------------------------
//   modified                     26  Mar. 1997  H.Kurahige
//   modified                     16  Apr. 1997  L.Urban    
//   modified                     18  Sep. 1997  H.Kurashige    
//   modified AlongStepGPIL etc.  17 Dec. 1997 H.Kurashige
//   fix bugs in GetGPILSelection() 24 Jan. 1998 H.Kurashige
//   modified for new ParticleChange 12 Mar. 1998  H.Kurashige

#ifndef G4VRestContinuousProcess_h
#define G4VRestContinuousProcess_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VProcess.hh"

class G4VRestContinuousProcess : public G4VProcess 
{
  //  Abstract class which defines the public behavior of
  //  discrete physics interactions.
  public:     

      G4VRestContinuousProcess(const G4String& ,
			       G4ProcessType   aType = fNotDefined );
      G4VRestContinuousProcess(G4VRestContinuousProcess &);

      virtual ~G4VRestContinuousProcess();

  public:   //  with description

      virtual G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& ,
			     G4ForceCondition* 
			    );

      virtual G4VParticleChange* AtRestDoIt(
			     const G4Track& ,
			     const G4Step&
			    );

      virtual G4double AlongStepGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4double previousStepSize,
			     G4double currentMinimumStep,
			     G4double& currentSafety,
                             G4GPILSelection* selection
			    );

      virtual G4VParticleChange* AlongStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    );
 
     //  no operation in  PostStepDoIt
      virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& ,
			     G4double ,
			     G4ForceCondition*
			    ){ return -1.0; };

     //  no operation in PostStepDoIt
      virtual G4VParticleChange* PostStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    ) {return 0;};

  protected: //  with description
    virtual G4double GetContinuousStepLimit(const G4Track& aTrack,
                             G4double  previousStepSize,
                             G4double currentMinimumStep,
			     G4double& currentSafety
                                                             )=0;
   // This pure virtual function is used to calculate step limit
    // for AlongStep in the derived processes  

  private:
    // this is the returnd value of  G4GPILSelection in 
    // the arguments of AlongStepGPIL()
    G4GPILSelection  valueGPILSelection;

  protected://  with description
    // these two methods are set/get methods for valueGPILSelection
    void SetGPILSelection(G4GPILSelection selection)
    { valueGPILSelection = selection;};

    G4GPILSelection GetGPILSelection() const{return valueGPILSelection;};

  protected: //  with description

      virtual G4double GetMeanLifeTime(const G4Track& aTrack,G4ForceCondition* condition)=0;
      //  Calculates the mean life-time (i.e. for decays) of the
      //  particle at rest due to the occurence of the given process,
      //  or converts the probability of interaction (i.e. for
      //  annihilation) into the life-time of the particle for the
      //  occurence of the given process.

  private:
  // hide default constructor and assignment operator as private 
      G4VRestContinuousProcess();
      G4VRestContinuousProcess & operator=(const G4VRestContinuousProcess &right);

};

// -----------------------------------------
//  inlined function members implementation
// -----------------------------------------
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4MaterialTable.hh"
#include "G4VParticleChange.hh"
inline G4double G4VRestContinuousProcess::AlongStepGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4double previousStepSize,
			     G4double currentMinimumStep,
			     G4double& currentSafety,
                             G4GPILSelection* selection
			    )
{
  // GPILSelection is set to defaule value of CandidateForSelection
  valueGPILSelection = CandidateForSelection;

  // get Step limit proposed by the process
  G4double steplength = GetContinuousStepLimit(track,previousStepSize,currentMinimumStep, currentSafety);

  // set return value for G4GPILSelection
  *selection = valueGPILSelection;
#ifdef G4VERBOSE
   if (verboseLevel>1){
    G4cout << "G4VRestContinuousProcess::AlongStepGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" <<G4endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " <<  track.GetMaterial()->GetName() <<G4endl;
    G4cout << "IntractionLength= " << steplength/cm <<"[cm] " <<G4endl;
  }
#endif
   return  steplength ;
}

inline G4double G4VRestContinuousProcess::AtRestGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4ForceCondition* condition
			    )
{
  // beggining of tracking 
  ResetNumberOfInteractionLengthLeft();

  // condition is set to "Not Forced"
  *condition = NotForced;

  // get mean life time
  currentInteractionLength = GetMeanLifeTime(track, condition);

#ifdef G4VERBOSE
   if ((currentInteractionLength <0.0) || (verboseLevel>2)){
    G4cout << "G4VRestContinuousProcess::AtRestGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" <<G4endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " << track.GetMaterial()->GetName() <<G4endl;
    G4cout << "MeanLifeTime = " << currentInteractionLength/ns << "[ns]" <<G4endl;
  }
#endif

  return theNumberOfInteractionLengthLeft * currentInteractionLength;
}


inline G4VParticleChange* G4VRestContinuousProcess::AtRestDoIt( 
			     const G4Track&,
			     const G4Step& 
			    )
{
//  clear NumberOfInteractionLengthLeft
    ClearNumberOfInteractionLengthLeft();

    return pParticleChange;
}

inline G4VParticleChange* G4VRestContinuousProcess::AlongStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    )
{ 
    return pParticleChange;
}

#endif






