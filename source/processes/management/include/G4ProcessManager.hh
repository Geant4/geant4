// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ProcessManager.hh,v 1.3 1999-05-03 01:52:36 kurasige Exp $
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
//   ----------------  G4ProcessManager  -----------------
// History:
// revised by G.Cosmo, 06 May 1996
//    Added vector of processes at rest, 06 May 1996
// ------------------------------------------------------------
//   New Physics scheme           8 Jan. 1997  H.Kurahige
//   Add SetProcessOrdering methods     27 Mar 1998  H.Kurahige
//   Add copy constructor (deep copy)   28 June 1998 H.Kurashige
//   Add GetProcessActivation     3 May. 1999 H.Kurashige
// ------------------------------------------------------------

#ifndef G4ProcessManager_h
#define G4ProcessManager_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <rw/tpordvec.h>

#include "G4VProcess.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleDefinition.hh"

class G4ProcessManagerMessenger;
class G4ProcessAttribute;

//  Indexes for ProcessVector
enum G4ProcessVectorTypeIndex
{ 
  	typeGPIL = 0,	// for GetPhysicalInteractionLength 
	typeDoIt =1		// for DoIt
};
enum G4ProcessVectorDoItIndex
{
  	idxAll = -1,		// for all DoIt/GPIL 
  	idxAtRest = 0, 		// for AtRestDoIt/GPIL
	idxAlongStep = 1, 	// for AlongStepDoIt/GPIL
	idxPostStep =2		// for AlongSTepDoIt/GPIL
};

//  enumeration for Ordering Parameter      
enum G4ProcessVectorOrdering
{ 
   	ordInActive = -1,			// ordering parameter to indicate InActive DoIt
   	ordDefault = INT_MAX/2,		// default ordering parameter
   	ordLast    = INT_MAX		// ordering parameter to indicate the last DoIt
};

class G4ProcessManager 
{
  //  It collects all physics a particle can undertake as seven vectors.
  //  These vectors are 
  //   one vector for all processes (called as "process List")
  //   two vectors for processes with AtRestGetPhysicalInteractionLength
  //                                    and AtRestDoIt
  //   two vectors for processes with AlongStepGetPhysicalInteractionLength
  //                                    and AlongStepDoIt
  //   two vectors for processes with PostStepGetPhysicalInteractionLength
  //                                    and PostStepDoIt
  //  The tracking will message three types of GetPhysicalInteractionLength
  //  in order to limit the Step and select the occurence of processes. 
  //  It will message the corresponding DoIt() to apply the selected 
  //  processes. In addition, the Tracking will limit the Step
  //  and select the occurence of the processes according to
  //  the shortest physical interaction length computed (except for
  //  processes at rest, for which the Tracking will select the
  //  occurence of the process which returns the shortest mean
  //  life-time from the GetPhysicalInteractionLength()).

  public: 
      // copy constructor
      G4ProcessManager(G4ProcessManager &right);

  private:
      // hide default constructor and assignment operator
      G4ProcessManager & operator=(G4ProcessManager &right);
      G4ProcessManager();


  public:
 
      G4ProcessManager(const G4ParticleDefinition* aParticleType);
      //  Constructor

      ~G4ProcessManager();
      //  Destructor

      G4int operator==(const G4ProcessManager &right) const;
      G4int operator!=(const G4ProcessManager &right) const;

      G4ProcessVector* GetProcessList() const;
      //  Returns the address of the vector of all processes 

      G4int  GetProcessListLength() const;
      //  Returns the number of process in the ProcessVector 

      G4int GetProcessIndex(G4VProcess *) const;
      //  Returns the index of the process in the process List

      // --------------------------------------

      G4ProcessVector* GetProcessVector( 
			       G4ProcessVectorDoItIndex idx,
			       G4ProcessVectorTypeIndex typ = typeGPIL
			      ) const;
      //  Returns the address of the vector of processes 

      G4ProcessVector* GetAtRestProcessVector(
			       G4ProcessVectorTypeIndex typ = typeGPIL
                              ) const; 
      //  Returns the address of the vector of processes for
      //    AtRestGetPhysicalInteractionLength      idx =0
      //    AtRestGetPhysicalDoIt                   idx =1
      G4ProcessVector* GetAlongStepProcessVector(
			       G4ProcessVectorTypeIndex typ = typeGPIL
                              ) const;
      //  Returns the address of the vector of processes for
      //    AlongStepGetPhysicalInteractionLength      idx =0
      //    AlongStepGetPhysicalDoIt                   idx =1

      G4ProcessVector* GetPostStepProcessVector(
			       G4ProcessVectorTypeIndex typ = typeGPIL
                              ) const;
      //  Returns the address of the vector of processes for
      //    PostStepGetPhysicalInteractionLength      idx =0
      //    PostStepGetPhysicalDoIt                   idx =1

      G4int GetProcessVectorIndex(
                           G4VProcess* aProcess,
			   G4ProcessVectorDoItIndex idx,
			   G4ProcessVectorTypeIndex typ  = typeGPIL
			   ) const;
      G4int GetAtRestIndex(
                           G4VProcess* aProcess,
			   G4ProcessVectorTypeIndex typ  = typeGPIL
			   ) const;
      G4int GetAlongStepIndex(
                           G4VProcess* aProcess,
			   G4ProcessVectorTypeIndex typ  = typeGPIL
			   ) const;
      G4int GetPostStepIndex(
			   G4VProcess* aProcess,
			   G4ProcessVectorTypeIndex typ = typeGPIL
			   ) const;
      //  Returns the index for GPIL/DoIt process vector of the process  

      G4int AddProcess(
             G4VProcess *aProcess,
             G4int      ordAtRestDoIt = ordInActive,
             G4int      ordAlongSteptDoIt = ordInActive,
             G4int      ordPostStepDoIt = ordInActive
            );
      //  Add a process to the process List
      //  return values are index to the List. Negative return value 
      //  indicates that the process has not be added due to some errors
      //  The first argument is a pointer to process.
      //  Following arguments are ordering parameters of the process in 
      //  process vectors. If value is negative, the process is
      //  not added to the corresponding process vector. 
   
      //  following methods are provided for simple processes  
      //   AtRestProcess has only AtRestDoIt
      //   ContinuousProcess has only AlongStepDoIt
      //   DiscreteProcess has only PostStepDoIt
      //  if ord is not specified, the process is
      //  added at the end of List of processvectors for 
      //  both DoIt and GetPhysicalInteractionLength 

      G4int AddRestProcess(G4VProcess *aProcess, G4int ord = ordDefault);
      G4int AddDiscreteProcess(G4VProcess *aProcess, G4int ord = ordDefault);
      G4int AddContinuousProcess(G4VProcess *aProcess, G4int ord = ordDefault);

      // Methods for setting ordering parameters
       // Altanative methods for setting ordering parameters 
      //   Note: AddProcess method should precede these methods

      G4int GetProcessOrdering(
			       G4VProcess *aProcess,
			       G4ProcessVectorDoItIndex idDoIt
                               );

      void SetProcessOrdering(
			       G4VProcess *aProcess,
			       G4ProcessVectorDoItIndex idDoIt,
			       G4int      ordDoIt = ordDefault
                               );
      // Set ordering parameter for DoIt specified by typeDoIt.
             
       void SetProcessOrderingToFirst(
			       G4VProcess *aProcess,
			       G4ProcessVectorDoItIndex idDoIt
			       );
      // Set ordering parameter to the first of all processes 
      // for DoIt specified by idDoIt.
      //  Note: If you use this method for two processes,
      //        a process called later will be first.

      void SetProcessOrderingToLast(
			       G4VProcess *aProcess,
			       G4ProcessVectorDoItIndex idDoIt
			       );
      // Set ordering parameter to the last of all processes 
      // for DoIt specified by idDoIt.
      //  Note: If you use this method for two processes,
      //        a process called later will be the last one.

      G4VProcess*  RemoveProcess(G4VProcess *aProcess);
      G4VProcess*  RemoveProcess(G4int      index);
      //  Removes a process from the process List.
      //  return value is pointer to the removed process.
      //  (0 value will be returned in case of errors)

      G4VProcess* SetProcessActivation(G4VProcess *aProcess, G4bool fActive);
      G4VProcess* SetProcessActivation(G4int      index, G4bool fActive);
      //  Set activation flag. 
      //  return value is pointer to the applied process.
      //  (0 value will be returned in case of errors)

      G4bool GetProcessActivation(G4VProcess *aProcess) const;
      G4bool GetProcessActivation(G4int      index) const;
      //  Get activation flag. 

      G4ParticleDefinition*  GetParticleType() const;
      // get the particle type 
      void SetParticleType(const G4ParticleDefinition*);
      // set the particle type 

      void StartTracking();
      void EndTracking();
      // these two methods are used by G4TrackingManager 
      // in order to inform Start/End of tracking for each track
      // to the process manager and all physics processes 

  private:     
      G4ProcessAttribute* GetAttribute(G4int      index) const;
      G4ProcessAttribute* GetAttribute(G4VProcess *aProcess) const;
      // get Pointer to ProcessAttribute

      G4VProcess* ActivateProcess(G4int   index);
      G4VProcess* InActivateProcess(G4int  index);
      // Activate/InActivateProcess   Process
      
  private:     
      const G4ParticleDefinition*   theParticleType;
      //  particle which has this process manager object     

      G4int             numberOfProcesses;
      G4ProcessVector*  theProcessList;
      // vector for all processes (called as "process List")

  public:
      enum {SizeOfProcVectorArray = 6};
  private:
      G4ProcessVector* theProcVector[SizeOfProcVectorArray];
      // vector for processes with GetPhysicalInteractionLength/DoIt

      typedef RWTPtrOrderedVector<G4ProcessAttribute> G4ProcessAttrVector; 
      G4ProcessAttrVector*  theAttrVector;
      // vector for process attribute  

  protected:
      G4int InsertAt(G4int position, G4VProcess* process, G4int ivec);
      // insert process at position in theProcVector[ivec]

      G4int RemoveAt(G4int position, G4VProcess* process, G4int ivec);
      // remove process at position in theProcVector[ivec]

      G4int FindInsertPosition(G4int ord, G4int ivec);
      // find insert position according to ordering parameter 
      // in theProcVector[ivec]

      G4int GetProcessVectorId(G4ProcessVectorDoItIndex idx,
			       G4ProcessVectorTypeIndex typ  = typeGPIL) const;
 
  private:
      G4bool  duringTracking;
      void    CreateGPILvectors();
      void    SetIndexToProcessVector(G4int ivec);

 public:
   void  DumpInfo();
   void  SetVerboseLevel(G4int value);
   G4int GetVerboseLevel() const;

 protected:
   G4int verboseLevel;
   // controle flag for output message
   //  0: Silent
   //  1: Warning message
   //  2: More
 
 private:
   static G4ProcessManagerMessenger* fProcessManagerMessenger;
   static G4int                      counterOfObjects;
};
#include "G4ProcessAttribute.hh"

// -----------------------------------------
//  inlined function members implementation
// -----------------------------------------
inline  
 void G4ProcessManager::SetParticleType(const G4ParticleDefinition* aParticle)
{
  theParticleType = aParticle;
}

inline 
 G4ProcessVector* G4ProcessManager::GetProcessList() const
{
  return theProcessList;
}

inline
 G4int  G4ProcessManager::GetProcessListLength() const
{
  return numberOfProcesses;
}

inline 
 G4int  G4ProcessManager::GetProcessIndex(G4VProcess* aProcess) const
{
  G4int idx = theProcessList->index(aProcess);
  if (idx>=numberOfProcesses) idx = -1;
  return idx;
}

inline 
 G4int G4ProcessManager::GetProcessVectorId(G4ProcessVectorDoItIndex idx,
					    G4ProcessVectorTypeIndex typ) const
{
  if ( idx == idxAtRest ) {
    if (typ == typeGPIL) { return 0; }
    else                 { return 1; }
  } else if ( idx == idxAlongStep ) {
    if (typ == typeGPIL) { return 2; }
    else                 { return 3; }
  } else if ( idx == idxPostStep ) {
    if (typ == typeGPIL) { return 4; }
    else                 { return 5; }
  } else {
    return -1;
  }
}
 
inline  
 G4ProcessVector* G4ProcessManager::GetProcessVector(
				       G4ProcessVectorDoItIndex idx,  
				       G4ProcessVectorTypeIndex typ
                                      ) const
{
  G4int ivec = GetProcessVectorId(idx, typ);
  if ( ivec >=0 ) {
    return theProcVector[ivec];
  } else {
    return 0;
  }
}

inline 
 G4ProcessVector* G4ProcessManager::GetAtRestProcessVector(G4ProcessVectorTypeIndex typ) const
{
  if (typ == typeGPIL) { return theProcVector[0]; }
  else                { return theProcVector[1]; }
}

inline 
 G4ProcessVector* G4ProcessManager::GetAlongStepProcessVector(G4ProcessVectorTypeIndex typ) const
{
  if (typ == typeGPIL) { return theProcVector[2]; }
  else                { return theProcVector[3]; }
}

inline 
 G4ProcessVector* G4ProcessManager::GetPostStepProcessVector(G4ProcessVectorTypeIndex typ) const
{
  if (typ == typeGPIL) { return theProcVector[4]; }
  else                { return theProcVector[5]; }
}

inline
 G4int G4ProcessManager::GetAtRestIndex(
                           G4VProcess* aProcess,
			   G4ProcessVectorTypeIndex typ 
			   ) const
{
  return GetProcessVectorIndex(aProcess, idxAtRest, typ);
}

inline 
 G4int G4ProcessManager::GetAlongStepIndex(
                           G4VProcess* aProcess,
			   G4ProcessVectorTypeIndex typ 
			   ) const
{
  return GetProcessVectorIndex(aProcess, idxAlongStep, typ);
}

inline
 G4int G4ProcessManager::GetPostStepIndex(
                           G4VProcess* aProcess,
			   G4ProcessVectorTypeIndex typ 
                         ) const
{
  return GetProcessVectorIndex(aProcess, idxPostStep, typ);
}

inline 
 G4int G4ProcessManager::AddRestProcess(G4VProcess *aProcess,G4int ord)
{
  return AddProcess(aProcess, ord, ordInActive, ordInActive);
}

inline 
 G4int G4ProcessManager::AddContinuousProcess(G4VProcess *aProcess,G4int ord)
{
  return AddProcess(aProcess, ordInActive, ord, ordInActive);
}

inline 
 G4int G4ProcessManager::AddDiscreteProcess(G4VProcess *aProcess,G4int ord)
{
  return AddProcess(aProcess, ordInActive, ordInActive, ord);
}

inline 
 G4ParticleDefinition* G4ProcessManager::GetParticleType() const
{ 
  return (G4ParticleDefinition* )theParticleType; 
}


inline 
 void G4ProcessManager::SetVerboseLevel(G4int value)
{
  verboseLevel = value;
}

inline  
 G4int G4ProcessManager::GetVerboseLevel() const
{
  return  verboseLevel;
}

#endif

