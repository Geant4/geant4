// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VContinuousProcess.hh,v 1.3 1999-11-07 17:11:46 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD group
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
//	add G4VContinuousProcess(const G4String&) 24 Jul 1996, Hisaya kurashige
//
// Class Description
//  Abstract class which defines the public behavior of
//  Continuous physics interactions.
//
// ------------------------------------------------------------
//   New Physics scheme           18 Dec. 1996  H.Kurahige
// ------------------------------------------------------------
//   modified                     25 Feb. 1997  H.Kurahige
//   modified                     8 Mar.  1997 H.Kurashige
//   modified                     22 Mar.  1997 H.Kurashige
//   modified                     26 Mar.  1997 H.Kurashige
//   modified AlongStepGPIL etc.  17 Dec. 1997 H.Kurashige
//   fix bugs in GetGPILSelection() 24 Jan. 1998 H.Kurashige
//   modified for new ParticleChange 12 Mar. 1998  H.Kurashige

#ifndef G4VContinuousProcess_h
#define G4VContinuousProcess_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VProcess.hh"

class G4VContinuousProcess : public G4VProcess 
{
  //  Abstract class which defines the public behavior of
  //  Continuous physics interactions.
  public:     

      G4VContinuousProcess(const G4String& ,
			   G4ProcessType   aType = fNotDefined );
      G4VContinuousProcess(G4VContinuousProcess &);

      virtual ~G4VContinuousProcess();

  public:  // with description   
      virtual G4double AlongStepGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4double  previousStepSize,
			     G4double  currentMinimumStep,
			     G4double& proposedSafety,
                             G4GPILSelection* selection
			    );

      virtual G4VParticleChange* AlongStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    );

     //  no operation in  AtRestDoIt and  PostStepDoIt
      virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track&,
			     G4double,
			     G4ForceCondition* 
			    ){ return -1.0; };

      virtual G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& ,
			     G4ForceCondition* 
			    ) { return -1.0; };

     //  no operation in  AtRestDoIt and PostStepDoIt
      virtual G4VParticleChange* AtRestDoIt(
			     const G4Track& ,
			     const G4Step&
			    ) {return 0;};

      virtual G4VParticleChange* PostStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    ) {return 0;};
 
  protected: // with description   
    virtual G4double GetContinuousStepLimit(const G4Track& aTrack,
                             G4double  previousStepSize,
                             G4double  currentMinimumStep,
			     G4double& currentSafety
                                                             )=0;
    // This pure virtual function is used to calculate step limit
    // for AlongStep in the derived processes  

  private:
    // this is the returnd value of  G4GPILSelection in 
    // the arguments of AlongStepGPIL()
    G4GPILSelection  valueGPILSelection;

  protected: // with description 
    // these two methods are set/get methods for valueGPILSelection
    void SetGPILSelection(G4GPILSelection selection)
    { valueGPILSelection = selection;};

    G4GPILSelection GetGPILSelection() const{return valueGPILSelection;};


  private:
  // hide default constructor and assignment operator as private 
      G4VContinuousProcess();
      G4VContinuousProcess & operator=(const G4VContinuousProcess &right);

};
// -----------------------------------------
//  inlined function members implementation
// -----------------------------------------
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4MaterialTable.hh"
#include "G4VParticleChange.hh"

inline G4double G4VContinuousProcess::AlongStepGetPhysicalInteractionLength(
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
    G4cout << "G4VContinuousProcess::AlongStepGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" <<endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " <<  track.GetMaterial()->GetName() <<endl;
    G4cout << "IntractionLength= " << steplength/cm <<"[cm] " <<endl;
  }
#endif

  return  steplength ;
}

inline G4VParticleChange* G4VContinuousProcess::AlongStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    )
{ 
    return pParticleChange;
}

#endif




