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
// $Id: G4VProcess.hh 105271 2017-07-18 07:35:12Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file
//
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
//
// Class Description
//  This class is the virtual class for physics process objects. 
//   It defines public methods which describe the behavior of 
//   a physics process.
//
// ------------------------------------------------------------
//   New Physics scheme           18 Dec. 1996  H.Kurahige
// ------------------------------------------------------------
//   change DoIt/GetPIL arguments type 20 Mar. 1997 H.Kurashige
//   modified AlongStepGPIL       17 Dec. 1997 H.Kurashige
//   modified for new ParticleChange 12 Mar. 1998  H.Kurashige
//   Add process trype            27 Mar. 1998  H.Kurashige
//   Remove thePhysicsTable       2 Aug. 1998   H.Kurashige
//   Add PILfactor and GPIL       3 Nov. 2000   H.Kurashige
//   Add Store/RetrievePhysicsTable 8  Nov. 2000   H.Kurashige
//   Modify Store/RetrievePhysicsTable methods 9 Mar. 2001   H.Kurashige
//   Added PreparePhysicsTable  20 Aug. 2004 H.Kurashige
//   Added isXXXXDoItIsEnabled   2 Oct. 2007 H.Kurashige
//   Added ProcessSubType   15 Nov. 2007 H.Kurashige

#ifndef G4VProcess_h 
#define G4VProcess_h 1

#include "globals.hh"
#include <cmath>
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

  public: // with description
  //  constructor requires the process name and type
      G4VProcess(const G4String& aName =  "NoName",
		 G4ProcessType   aType = fNotDefined );

  //  copy constructor copys the name but does not copy the 
  //  physics table (0 pointer is assigned)
      G4VProcess(const G4VProcess &right);

  public: 
  //  destructor 
      virtual ~G4VProcess();

  // equal opperators
      G4int operator==(const G4VProcess &right) const;
      G4int operator!=(const G4VProcess &right) const;

  public: // with description
  ////////////////////////////
  // DoIt    /////////////////
  ///////////////////////////
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

  //////////////////////////
  // GPIL    //////////////
  /////////////////////////  
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

      G4double GetCurrentInteractionLength() const;
      // Returns currentInteractionLength

      ////////// PIL factor ////////
      void SetPILfactor(G4double value);
      G4double GetPILfactor() const;
      // Set/Get factor for PhysicsInteractionLength 
      // which is passed to G4SteppingManager for both AtRest and PostStep

      // These three GPIL methods are used by Stepping Manager.
      // They invoke virtual GPIL methods listed above.
      // As for AtRest and PostStep the returned value is multipled by thePILfactor 
      // 
      G4double AlongStepGPIL( const G4Track& track,
                              G4double  previousStepSize,
                              G4double  currentMinimumStep,
                              G4double& proposedSafety,
                              G4GPILSelection* selection     );

      G4double AtRestGPIL( const G4Track& track,
                           G4ForceCondition* condition );

      G4double PostStepGPIL( const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition );

  ////////////////////// 
      virtual G4bool IsApplicable(const G4ParticleDefinition&){return true;}
      // Returns true if this process object is applicable to
      // the particle type
      // Process will not be registered to a particle if IsApplicable is false   

      virtual void BuildPhysicsTable(const G4ParticleDefinition&){}
      // Messaged by the Particle definition (via the Process manager)
      // whenever cross section tables have to be rebuilt (i.e. if new
      // materials have been defined). 
      // It is overloaded by individual processes when they need physics
      // tables. 

     virtual void PreparePhysicsTable(const G4ParticleDefinition&){}
      // Messaged by the Particle definition (via the Process manager)
      // whenever cross section tables have to be prepare for rebuilt 
      // (i.e. if new materials have been defined). 
      // It is overloaded by individual processes when they need physics
      // tables. 

      // Processes which Build physics tables independent of cuts
      // (for example in their constructors)
      // should preferably use private 
      // void BuildThePhysicsTable() and void PreparePhysicsTable().
      // Not another BuildPhysicsTable, please.


      virtual G4bool StorePhysicsTable(const G4ParticleDefinition* ,
				       const G4String&, G4bool){return true;}
      // Store PhysicsTable in a file.
      // (return false in case of failure at I/O )

      virtual G4bool RetrievePhysicsTable( const G4ParticleDefinition* ,
					   const G4String&, G4bool){return false;}
      // Retrieve Physics from a file.
      // (return true if the Physics Table can be build by using file)
      // (return false if the process has no functionality or in case of failure)
      // File name should be defined by each process
      // and the file should be placed under the directory specifed by the argument.
      const G4String& GetPhysicsTableFileName(const G4ParticleDefinition* ,
					      const G4String& directory,
					      const G4String& tableName,
					      G4bool ascii =false);
      // this method is utility for Store/RetreivePhysicsTable

  ////////////////////////////
      const G4String& GetProcessName() const;
      //  Returns the name of the process.

      G4ProcessType GetProcessType() const;
      //  Returns the process type.

      void SetProcessType(G4ProcessType );
      //  Set the process type.

      G4int GetProcessSubType() const;
      //  Returns the process sub type.

      void SetProcessSubType(G4int );
      //  Set the process sub type.

      static const G4String& GetProcessTypeName(G4ProcessType );
      //  Returns the process type name

      virtual void StartTracking(G4Track*);
      virtual void EndTracking();
      // inform Start/End of tracking for each track to the physics process 

  public:
      virtual void SetProcessManager(const G4ProcessManager*); 
      // A process manager set its own pointer when the process is registered
      // the process Manager
      virtual  const G4ProcessManager* GetProcessManager(); 
      // Get the process manager which the process belongs to
  
  protected:
      const G4ProcessManager* aProcessManager; 
 
  protected:
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

      G4double          theInitialNumberOfInteractionLength;
     // The initial value when ResetNumberOfInteractionLengthLeft is invoked

 public: // with description
      virtual void      ResetNumberOfInteractionLengthLeft();
     // reset (determine the value of)NumberOfInteractionLengthLeft

      G4double GetNumberOfInteractionLengthLeft() const;
     // get NumberOfInteractionLengthLeft

      G4double GetTotalNumberOfInteractionLengthTraversed() const;
     // get NumberOfInteractionLength 
     //   after  ResetNumberOfInteractionLengthLeft is invoked

 protected:  // with description
     void      SubtractNumberOfInteractionLengthLeft(
				  G4double previousStepSize
                                );
     // subtract NumberOfInteractionLengthLeft by the value corresponding to 
     // previousStepSize      
 
     void      ClearNumberOfInteractionLengthLeft();
     // clear NumberOfInteractionLengthLeft 
     // !!! This method should be at the end of PostStepDoIt()
     // !!! and AtRestDoIt

 public: // with description
    // These methods indicate which DoIt is enabled
    // These methods are used by G4ProcessManager to check
    // that ordering parameters are set properly
    G4bool isAtRestDoItIsEnabled() const;
    G4bool isAlongStepDoItIsEnabled() const;
    G4bool isPostStepDoItIsEnabled() const;
  
 protected: 
      G4String theProcessName;
      //  The name of the process

      G4String thePhysicsTableFileName;

      G4ProcessType theProcessType;
      //  The type of the process

      G4int theProcessSubType;
      //  The sub type of the process

      G4double thePILfactor;
      // factor for PhysicsInteractionLength 
      // which is passed to G4SteppingManager
 
      G4bool enableAtRestDoIt;
      G4bool enableAlongStepDoIt;
      G4bool enablePostStepDoIt;
      
 public: // with description
   virtual void  DumpInfo() const;
   // dump out process information    

   virtual void ProcessDescription(std::ostream& outfile) const;
   // write out to html file for automatic documentation

 public: // with description
   void  SetVerboseLevel(G4int value);
   G4int GetVerboseLevel() const;
   // set/get controle flag for output message
   //  0: Silent
   //  1: Warning message
   //  2: More


 protected:
   G4int verboseLevel;
   // controle flag for output message
    
private:
    G4VProcess* masterProcessShadow;
    //For multi-threaded: poitner to the instance of this process
    // for the master thread
public:
    virtual void SetMasterProcess( G4VProcess* masterP);
    // Sets the master thread process instance
    const G4VProcess* GetMasterProcess() const;
    // Returns the master thread process instnace
    // Can be used to initialize worker type processes
    // instances from master one (e.g. to share a read-only table)
    // if ( this != GetMasterProcess() ) { /*worker*/ }
    // else { /* master or sequential */ }

    virtual void BuildWorkerPhysicsTable(const G4ParticleDefinition& part);
    // Messaged by the Particle definition (via the Process manager)
    // in worker threads. See BuildWorkerBhyiscsTable method.
    // Can be used to share among threads physics tables. Use GetMasterProcess
    // to get pointer of master process from worker thread.
    // By default this method makes a forward call to
    // BuildPhysicsTable
    
    virtual void PrepareWorkerPhysicsTable(const G4ParticleDefinition&);
    // Messaged by the Particle definition (via the Process manager)
    // in worker threads. See PreparephysicsTable
    // Can be used to share among threads physics tables. Use GetMasterProcess
    // to get pointer of master process from worker thread
    // By default this method makes a forward call
    // to PreparePhysicsTable
};

// -----------------------------------------
//  inlined function members implementation
// -----------------------------------------
#include "Randomize.hh"              

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

inline  void G4VProcess::SetVerboseLevel(G4int value)
{
  verboseLevel = value;
}

inline  G4int G4VProcess::GetVerboseLevel() const
{
  return  verboseLevel;
}

inline void G4VProcess::ClearNumberOfInteractionLengthLeft()
{
  theInitialNumberOfInteractionLength = -1.0; 
  theNumberOfInteractionLengthLeft =  -1.0;
}

inline G4double G4VProcess::GetNumberOfInteractionLengthLeft() const
{
  return theNumberOfInteractionLengthLeft;
}

inline G4double G4VProcess::GetTotalNumberOfInteractionLengthTraversed() const
{
  return theInitialNumberOfInteractionLength - theNumberOfInteractionLengthLeft;}

inline G4double G4VProcess::GetCurrentInteractionLength() const
{
  return currentInteractionLength;
}

inline void G4VProcess::SetPILfactor(G4double value)
{
  if (value>0.) {
    thePILfactor = value;
  }
}

inline G4double G4VProcess::GetPILfactor() const
{
  return thePILfactor;
}

inline G4double G4VProcess::AlongStepGPIL( const G4Track& track,
                                     G4double  previousStepSize,
                                     G4double  currentMinimumStep,
                                     G4double& proposedSafety,
                                     G4GPILSelection* selection     )
{
  G4double value
   =AlongStepGetPhysicalInteractionLength(track, previousStepSize, currentMinimumStep, proposedSafety, selection);
  return value;
}

inline G4double G4VProcess::AtRestGPIL( const G4Track& track,
                                 G4ForceCondition* condition )
{
  G4double value
   =AtRestGetPhysicalInteractionLength(track, condition);
  return thePILfactor*value;
}

inline G4double G4VProcess::PostStepGPIL( const G4Track& track,
                                   G4double   previousStepSize,
                                   G4ForceCondition* condition )
{
  G4double value
   =PostStepGetPhysicalInteractionLength(track, previousStepSize, condition);
  return thePILfactor*value;
}
      
inline 
 void G4VProcess::SetProcessManager(const G4ProcessManager* procMan)
{
   aProcessManager = procMan; 
}

inline
 const G4ProcessManager* G4VProcess::GetProcessManager()
{
  return  aProcessManager; 
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
void G4VProcess::SubtractNumberOfInteractionLengthLeft(
                                  G4double previousStepSize )
{
  if (currentInteractionLength>0.0) {
    theNumberOfInteractionLengthLeft -= previousStepSize/currentInteractionLength;
    if(theNumberOfInteractionLengthLeft<0.) {
       theNumberOfInteractionLengthLeft=CLHEP::perMillion;
    }

  } else {
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cerr << "G4VProcess::SubtractNumberOfInteractionLengthLeft()";
      G4cerr << " [" << theProcessName << "]" <<G4endl;
      G4cerr << " currentInteractionLength = " << currentInteractionLength << " [mm]";
      G4cerr << " previousStepSize = " << previousStepSize << " [mm]";
      G4cerr << G4endl;
    }
#endif
    G4String msg = "Negative currentInteractionLength for ";
    msg +=      theProcessName;
    G4Exception("G4VProcess::SubtractNumberOfInteractionLengthLeft()",
                "ProcMan201",EventMustBeAborted,
                msg);
  }
}

#endif
