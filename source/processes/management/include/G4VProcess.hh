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
// G4VProcess
//
// Class description:
//
// This class is the virtual class for physics process objects. 
// It defines public methods which describe the behavior of 
// a physics process.

// Authors:
// - 2 December 1995, G.Cosmo - First implementation, based on object model
// - 18 December 1996, H.Kurashige - New Physics scheme
// --------------------------------------------------------------------
#ifndef G4VProcess_hh 
#define G4VProcess_hh 1

#include <cmath>

#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"              

#include "G4PhysicsTable.hh"
#include "G4VParticleChange.hh"
#include "G4ForceCondition.hh"
#include "G4GPILSelection.hh"
#include "G4ParticleChange.hh"
#include "G4ProcessType.hh"

class G4ParticleDefinition;
class G4DynamicParticle;
class G4Track;
class G4Step;
class G4ProcessTable;

class G4VProcess 
{

  public:

    G4VProcess(const G4String& aName = "NoName",
               G4ProcessType aType = fNotDefined);
      // Constructor requires the process name and type

    G4VProcess(const G4VProcess& right);
      // Copy constructor copies the name but does not copy the 
      // physics table (null pointer is assigned instead)

    virtual ~G4VProcess();
      // Destructor 

    G4VProcess& operator=(const G4VProcess&) = delete;

    G4bool operator==(const G4VProcess& right) const;
    G4bool operator!=(const G4VProcess& right) const;
      // Equality operators

    ////////////////////////////
    // DoIt    /////////////////
    ////////////////////////////

    virtual G4VParticleChange* PostStepDoIt(
                             const G4Track& track,
                             const G4Step& stepData
                            ) = 0;

    virtual G4VParticleChange* AlongStepDoIt(
                             const G4Track& track,
                             const G4Step& stepData
                            ) = 0;
    virtual G4VParticleChange* AtRestDoIt(
                             const G4Track& track,
                             const G4Step& stepData
                            ) = 0;
      // A virtual base class function that has to be overridden
      // by any subclass. The DoIt() method actually performs the
      // physics process and determines either momentum change
      // of the production of secondaries etc.
      //    Arguments
      //      const G4Track& track:
      //        reference to the current G4Track information
      //      const G4Step&  stepData:
      //        reference to the current G4Step information

    //////////////////////////
    // GPIL    ///////////////
    //////////////////////////

    virtual G4double AlongStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double previousStepSize,
                             G4double currentMinimumStep,
                             G4double& proposedSafety,
                             G4GPILSelection* selection) = 0;

    virtual G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4ForceCondition* condition ) = 0;

    virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double previousStepSize,
                             G4ForceCondition* condition ) = 0;
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

    inline G4double GetCurrentInteractionLength() const;
      // Returns currentInteractionLength

    ////////// PIL factor ////////
    //
    inline void SetPILfactor(G4double value);
    inline G4double GetPILfactor() const;
      // Set/Get factor for PhysicsInteractionLength 
      // which is passed to G4SteppingManager for both AtRest and PostStep

    // These three GPIL methods are used by Stepping Manager.
    // They invoke virtual GPIL methods listed above.
    // As for AtRest and PostStep the returned value is multipled by
    // thePILfactor 
    // 
    inline G4double AlongStepGPIL( const G4Track& track,
                                   G4double previousStepSize,
                                   G4double currentMinimumStep,
                                   G4double& proposedSafety,
                                   G4GPILSelection* selection );

    inline G4double AtRestGPIL( const G4Track& track,
                                G4ForceCondition* condition );

    inline G4double PostStepGPIL( const G4Track& track,
                                  G4double previousStepSize,
                                  G4ForceCondition* condition );

    virtual G4bool IsApplicable(const G4ParticleDefinition&) { return true; }
      // Returns true if this process object is applicable to
      // the particle type. Process will not be registered to a
      // particle if IsApplicable is false   

    virtual void BuildPhysicsTable(const G4ParticleDefinition&) {}
      // Messaged by the Particle definition (via the Process manager)
      // whenever cross-section tables have to be rebuilt (i.e. if new
      // materials have been defined). 
      // It is overloaded by individual processes when they need physics
      // tables

    virtual void PreparePhysicsTable(const G4ParticleDefinition&) {}
      // Messaged by the Particle definition (via the Process manager)
      // whenever cross-section tables have to be prepared for rebuild
      // (i.e. if new materials have been defined). 
      // It is overloaded by individual processes when they need physics
      // tables

    // Processes which Build physics tables independent of cuts
    // (for example in their constructors) should preferably use private 
    // void BuildThePhysicsTable() and void PreparePhysicsTable().
    // *Not* another BuildPhysicsTable

    virtual G4bool StorePhysicsTable(const G4ParticleDefinition* ,
                                     const G4String&, G4bool) { return true; }
      // Store PhysicsTable in a file.
      // Return false in case of failure at I/O

    virtual G4bool RetrievePhysicsTable(const G4ParticleDefinition* ,
                                      const G4String&, G4bool) { return false; }
      // Retrieve Physics from a file.
      // Return true if the Physics Table can be built by using file.
      // Return false if the process has no functionality or in case
      // of failure. File name should be defined by each process and the
      // file should be placed under the directory specified by the argument

    const G4String& GetPhysicsTableFileName(const G4ParticleDefinition* ,
                                            const G4String& directory,
                                            const G4String& tableName,
                                            G4bool ascii = false);
      // This method is utility for Store/RetreivePhysicsTable

    inline const G4String& GetProcessName() const;
      // Returns the name of the process

    inline G4ProcessType GetProcessType() const;
      // Returns the process type

    inline void SetProcessType(G4ProcessType);
      // Sets the process type

    inline G4int GetProcessSubType() const;
      // Returns the process sub type

    inline void SetProcessSubType(G4int);
      // Sets the process sub type

    static const G4String& GetProcessTypeName(G4ProcessType);
      // Returns the process type name

    virtual const G4VProcess* GetCreatorProcess() const;
      // Returns the process to be used as CreatorProcess for secondaries
      // coming from this process

    virtual void StartTracking(G4Track*);
    virtual void EndTracking();
      // Inform Start/End of tracking for each track to the physics process 

    virtual void SetProcessManager(const G4ProcessManager*); 
      // A process manager sets its own pointer when the process
      // is registered in the process Manager
    virtual const G4ProcessManager* GetProcessManager(); 
      // Get the process manager which the process belongs to
  
    virtual void ResetNumberOfInteractionLengthLeft();
      // Reset (determine the value of) NumberOfInteractionLengthLeft

    inline G4double GetNumberOfInteractionLengthLeft() const;
      // Get NumberOfInteractionLengthLeft

    inline G4double GetTotalNumberOfInteractionLengthTraversed() const;
      // Get NumberOfInteractionLength after
      //  ResetNumberOfInteractionLengthLeft() is invoked

    inline G4bool isAtRestDoItIsEnabled() const;
    inline G4bool isAlongStepDoItIsEnabled() const;
    inline G4bool isPostStepDoItIsEnabled() const;
      // These methods indicate which DoIt is enabled.
      // They are used by G4ProcessManager to check
      // that ordering parameters are properly set
  
    virtual void  DumpInfo() const;
      // Dump out process information    

    virtual void ProcessDescription(std::ostream& outfile) const;
      // Write out to html file for automatic documentation

    inline void  SetVerboseLevel(G4int value);
    inline G4int GetVerboseLevel() const;
      // set/get control flag for output message
      //  0: Silent
      //  1: Warning message
      //  2: More

    virtual void SetMasterProcess(G4VProcess* masterP);
      // Sets the master thread process instance
    inline const G4VProcess* GetMasterProcess() const;
      // Returns the master thread process instance.
      // Can be used to initialise worker type processes
      // instances from master one (e.g. to share a read-only table)
      // if ( this != GetMasterProcess() ) { /*worker*/ }
      // else { /* master or sequential */ }

    virtual void BuildWorkerPhysicsTable(const G4ParticleDefinition& part);
      // Messaged by the Particle definition (via the Process manager)
      // in worker threads. See BuildWorkerPhyiscsTable() method.
      // Can be used to share among threads physics tables.
      // Use GetMasterProcess() to get pointer of master process from
      // worker thread.
      // By default this method makes a forward call to BuildPhysicsTable()
    
    virtual void PrepareWorkerPhysicsTable(const G4ParticleDefinition&);
      // Messaged by the Particle definition (via the Process manager)
      // in worker threads. See PreparephysicsTable().
      // Can be used to share among threads physics tables.
      // Use GetMasterProcess() to get pointer of master process from
      // worker thread
      // By default this method makes a forward call to PreparePhysicsTable()

  protected:

    inline void SubtractNumberOfInteractionLengthLeft(G4double prevStepSize);
      // Subtract NumberOfInteractionLengthLeft by the value corresponding
      // to previousStepSize      
 
    inline void ClearNumberOfInteractionLengthLeft();
      // This method should be at the end of PostStepDoIt() and AtRestDoIt()!

  protected:

    const G4ProcessManager* aProcessManager = nullptr; 
 
    G4VParticleChange* pParticleChange = nullptr;
      // The pointer to G4VParticleChange object 
      // which is modified and returned by address by the DoIt() method.
      // This pointer should be set in each physics process
      // after construction of derived class object

    G4ParticleChange aParticleChange;
      // This object is kept for compatibility with old scheme.
      // May be removed in future

    G4double theNumberOfInteractionLengthLeft = -1.0;
      // The flight length left for the current tracking particle
      // in unit of "Interaction length"

    G4double currentInteractionLength = -1.0;
      // The InteractionLength in the current material

    G4double theInitialNumberOfInteractionLength = -1.0;
      // The initial value when ResetNumberOfInteractionLengthLeft() is invoked

    G4String theProcessName;
      // The name of the process

    G4String thePhysicsTableFileName;

    G4ProcessType theProcessType = fNotDefined;
      // The type of the process

    G4int theProcessSubType = -1;
      // The sub type of the process

    G4double thePILfactor = 1.0;
      // Factor for PhysicsInteractionLength
      // which is passed to G4SteppingManager
 
    G4int verboseLevel = 0;
      // Controle flag for output message

    G4bool enableAtRestDoIt = true;
    G4bool enableAlongStepDoIt = true;
    G4bool enablePostStepDoIt = true;

  private:
 
    G4VProcess();  
      // Hidden default constructor

  private:

    G4VProcess* masterProcessShadow = nullptr;
      // For multi-threaded: pointer to the instance of this process
      // for the master thread

    G4ProcessTable* fProcessTable = nullptr;
};

// -----------------------------------------
//  inlined function members implementation
// -----------------------------------------

inline 
const G4String& G4VProcess::GetProcessName() const
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

inline
 G4int G4VProcess::GetProcessSubType() const
{
  return theProcessSubType;
}

inline
void G4VProcess::SetProcessSubType(G4int value)
{
  theProcessSubType = value;
}

inline
void G4VProcess::SetVerboseLevel(G4int value)
{
  verboseLevel = value;
}

inline
G4int G4VProcess::GetVerboseLevel() const
{
  return  verboseLevel;
}

inline
void G4VProcess::ClearNumberOfInteractionLengthLeft()
{
  theInitialNumberOfInteractionLength = -1.0; 
  theNumberOfInteractionLengthLeft =  -1.0;
}

inline
G4double G4VProcess::GetNumberOfInteractionLengthLeft() const
{
  return theNumberOfInteractionLengthLeft;
}

inline
G4double G4VProcess::GetTotalNumberOfInteractionLengthTraversed() const
{
  return theInitialNumberOfInteractionLength - theNumberOfInteractionLengthLeft;
}

inline
G4double G4VProcess::GetCurrentInteractionLength() const
{
  return currentInteractionLength;
}

inline
void G4VProcess::SetPILfactor(G4double value)
{
  if (value>0.) { thePILfactor = value; }
}

inline
G4double G4VProcess::GetPILfactor() const
{
  return thePILfactor;
}

inline
G4double G4VProcess::AlongStepGPIL( const G4Track& track,
                                    G4double  previousStepSize,
                                    G4double  currentMinimumStep,
                                    G4double& proposedSafety,
                                    G4GPILSelection* selection )
{
  return AlongStepGetPhysicalInteractionLength(track, previousStepSize,
                             currentMinimumStep, proposedSafety, selection);
}

inline
G4double G4VProcess::AtRestGPIL( const G4Track& track,
                                 G4ForceCondition* condition )
{
  return thePILfactor * AtRestGetPhysicalInteractionLength(track, condition);
}

inline
G4double G4VProcess::PostStepGPIL( const G4Track& track,
                                   G4double previousStepSize,
                                   G4ForceCondition* condition )
{
  return thePILfactor *
      PostStepGetPhysicalInteractionLength(track, previousStepSize, condition);
}
      
inline 
void G4VProcess::SetProcessManager(const G4ProcessManager* procMan)
{
  aProcessManager = procMan; 
}

inline
const G4ProcessManager* G4VProcess::GetProcessManager()
{
  return aProcessManager; 
}

inline
G4bool G4VProcess::isAtRestDoItIsEnabled() const
{
  return enableAtRestDoIt;
}

inline
G4bool G4VProcess::isAlongStepDoItIsEnabled() const
{
  return enableAlongStepDoIt;
}

inline
G4bool G4VProcess::isPostStepDoItIsEnabled() const
{
  return enablePostStepDoIt;
}

inline
const G4VProcess* G4VProcess::GetMasterProcess() const
{
  return masterProcessShadow;
}

inline
void G4VProcess::SubtractNumberOfInteractionLengthLeft( G4double prevStepSize )
{
  if (currentInteractionLength>0.0)
  {
    theNumberOfInteractionLengthLeft -= prevStepSize/currentInteractionLength;
    if(theNumberOfInteractionLengthLeft<0.)
    {
       theNumberOfInteractionLengthLeft=CLHEP::perMillion;
    }
  }
  else
  {
#ifdef G4VERBOSE
    if (verboseLevel>0)
    {
      G4cerr << "G4VProcess::SubtractNumberOfInteractionLengthLeft()";
      G4cerr << " [" << theProcessName << "]" <<G4endl;
      G4cerr << " currentInteractionLength = "
             << currentInteractionLength << " [mm]";
      G4cerr << " previousStepSize = " << prevStepSize << " [mm]";
      G4cerr << G4endl;
    }
#endif
    G4String msg = "Negative currentInteractionLength for ";
    msg += theProcessName;
    G4Exception("G4VProcess::SubtractNumberOfInteractionLengthLeft()",
                "ProcMan201", EventMustBeAborted, msg);
  }
}

#endif
