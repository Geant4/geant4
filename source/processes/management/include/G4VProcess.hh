// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VProcess.hh,v 1.1 1999-01-07 16:13:56 gunter Exp $
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
// ------------------------------------------------------------
//   New Physics scheme           18 Dec. 1996  H.Kurahige
// ------------------------------------------------------------
//   change DoIt/GetPIL arguments type 20 Mar. 1997 H.Kurashige
//   modified AlongStepGPIL       17 Dec. 1997 H.Kurashige
//   modified for new ParticleChange 12 Mar. 1998  H.Kurashige
//   Add process trype            27 Mar. 1998  H.Kurashige
//   Remove thePhysicsTable       2 Aug. 1998   H.Kurashige

#ifndef G4VProcess_h 
#define G4VProcess_h 1

#include "globals.hh"
#include "G4ios.hh"

class G4ParticleDefinition;
class G4DynamicParticle;
class G4Track;
class G4Step;

#include "G4PhysicsTable.hh"
#include "G4VParticleChange.hh"
#include "G4ForceCondition.hh"
#include "G4GPILSelection.hh"
#include "G4ParticleChange.hh"
#include "G4ProcessType.hh"

class G4VProcess 
{
  //  A virtual class for physics process objects. It defines
  //  public methods which describe the behavior of a
  //  physics process.

  private:
  // hide default constructor and assignment operator as private 
  //  do not hide default constructor for alpha version 
  //  G4VProcess G4VProcess();  
      G4VProcess & operator=(const G4VProcess &right);

  public:
  //  constructor requires the process name and type
      G4VProcess(const G4String& aName =  "NoName",
		 G4ProcessType   aType = fNotDefined );

  //  copy constructor copys the name but does not copy the 
  //  physics table (NULL pointer is assigned)
      G4VProcess(G4VProcess &right);

      virtual ~G4VProcess();

      G4int operator==(const G4VProcess &right) const;
      G4int operator!=(const G4VProcess &right) const;

      virtual G4VParticleChange* PostStepDoIt(
			     const G4Track& track,
			     const G4Step&  stepData
			    ) = 0;

      virtual G4VParticleChange* AlongStepDoIt(
			     const G4Track& track,
			     const G4Step& stepData
			    ) = 0;
      virtual G4VParticleChange* AtRestDoIt(
			     const G4Track& track,
			     const G4Step& stepData
			    ) = 0;
      //  A virtual base class function that has to be overridden
      //  by any subclass. The DoIt method actually performs the
      //  physics process and determines either momentum change
      //  of the production of secondaries etc.
      //    arguments
      //      const G4Track&    track:
      //        reference to the current G4Track information
      //      const G4Step&     stepData:
      //        reference to the current G4Step information

      virtual G4double AlongStepGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4double  previousStepSize,
			     G4double  currentMinimumStep,
			     G4double& proposedSafety,
                             G4GPILSelection* selection) = 0;

      virtual G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4ForceCondition* condition
			    ) = 0;

      virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4double   previousStepSize,
			     G4ForceCondition* condition
			    ) = 0;
  
      //  Returns the Step-size (actual length) which is allowed 
      //  by "this" process. (for AtRestGetPhysicalInteractionLength,
      //  return value is Step-time) The NumberOfInteractionLengthLeft is
      //  recalculated by using previousStepSize and the Step-size is 
      //  calucalted accoding to the resultant NumberOfInteractionLengthLeft.
      //  using NumberOfInteractionLengthLeft, which is recalculated at 
      //    arguments
      //      const G4Track&    track:
      //        reference to the current G4Track information
      //      G4double*          previousStepSize: 
      //        the Step-size (actual length) of the previous Step 
      //        of this track. Negative calue indicates that
      //        NumberOfInteractionLengthLeft must be reset.
      //        the current physical interaction legth of this process
      //      G4ForceCondition* condition:
      //        the flag indicates DoIt of this process is forced 
      //        to be called
      //         Forced:    Corresponding DoIt is forced
      //         NotForced: Corresponding DoIt is called 
      //                    if the Step size of this Step is determined 
      //                    by this process
      //        !! AlongStepDoIt is always called !! 
      //      G4double& currentMinimumStep:
      //        this value is used for transformation of
      //        true path length to geometrical path length

      virtual G4bool IsApplicable(const G4ParticleDefinition&){return true;};
      // Returns true if this process object is applicable to
      // the particle type

      virtual void BuildPhysicsTable(const G4ParticleDefinition&){};
      // Messaged by the Particle definition (via the Process manager)
      // whenever cross section tables have to be rebuilt (i.e. if new
      // materials have been defined). 
      // It is overloaded by individual processes when they need physics
      // tables. 

      // Processes which Build (for example in their
      // constructors) physics tables independent of cuts
      // should preferably use a
      // private void BuildThePhysicsTable()
      // function. Not another BuildPhysicsTable, please.

      G4String GetProcessName() const;
      //  Returns the name of the process.

      G4ProcessType GetProcessType() const;
      //  Returns the process type.

      void SetProcessType(G4ProcessType );
      //  Set the process type.

      static G4String GetProcessTypeName(G4ProcessType );
      //  Returns the process type name

      virtual void StartTracking();
      virtual void EndTracking();
      // inform Start/End of tracking for each track to the physics process 

  protected:
       //--- Removed    ----//
      //  G4PhysicsTable* thePhysicsTable;
      //  A Physics Table can be either a cross-sections table or
      //  an energy table (or can be used for other specific
      //  purposes).

      G4VParticleChange* pParticleChange;
      //  The pointer to G4VParticleChange object 
      //  which is modified and returned by address by the DoIt() method.
      //  This pointer should be set in each physics process
      //  after construction of derived class object.  

      G4ParticleChange aParticleChange;
      //  This object is kept for compatibility with old scheme
      //  This will be removed in future

      G4double          theNumberOfInteractionLengthLeft;
     // The flight length left for the current tracking particle
     // in unit of "Interaction length".

      G4double          currentInteractionLength;
     // The InteractionLength in the current material

 public:
      virtual void      ResetNumberOfInteractionLengthLeft();
     // reset (determine the value of)NumberOfInteractionLengthLeft
 
 protected: 
     virtual void      SubtractNumberOfInteractionLengthLeft(
				  G4double previousStepSize
                                );
     // subtract NumberOfInteractionLengthLeft by the value corresponding to 
     // previousStepSize      
 
     virtual void      ClearNumberOfInteractionLengthLeft();
     // clear NumberOfInteractionLengthLeft 
     // !!! This method should be at the end of PostStepDoIt()
     // !!! and AtRestDoIt

 private: 
      G4String theProcessName;
      //  The name of the process

      G4ProcessType theProcessType;
      //  The type of the process

 public:
   virtual void  DumpInfo() const;
   
 public:
   void  SetVerboseLevel(G4int value);
   G4int GetVerboseLevel() const;

 protected:
   G4int verboseLevel;
   // controle flag for output message
   //  0: Silent
   //  1: Warning message
   //  2: More

};

// -----------------------------------------
//  inlined function members implementation
// -----------------------------------------
#include "Randomize.hh"              

inline G4String G4VProcess::GetProcessName() const
{
  return theProcessName;
}

inline      
 G4ProcessType G4VProcess::GetProcessType() const
{
  return theProcessType;
}

inline
 void G4VProcess::SetProcessType(G4ProcessType aType)
{
  theProcessType = aType;
}

inline  void G4VProcess::SetVerboseLevel(G4int value)
{
  verboseLevel = value;
}

inline  G4int G4VProcess::GetVerboseLevel() const
{
  return  verboseLevel;
}

inline void G4VProcess::ResetNumberOfInteractionLengthLeft()
{
  theNumberOfInteractionLengthLeft =  -log( G4UniformRand() );
}

inline void G4VProcess::ClearNumberOfInteractionLengthLeft()
{
  theNumberOfInteractionLengthLeft =  -1.0;
}

#endif
