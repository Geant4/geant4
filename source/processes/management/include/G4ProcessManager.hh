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
// G4ProcessManager
//
// Class Description:
//
// G4ProcessManager collects all physics a particle can undertake as
// vectors. These vectors are:
// - one vector for all processes (called as "process List")
// - two vectors for processes with AtRestGetPhysicalInteractionLength()
//                              and AtRestDoIt()
// - two vectors for processes with AlongStepGetPhysicalInteractionLength()
//                              and AlongStepDoIt()
// - two vectors for processes with PostStepGetPhysicalInteractionLength()
//                              and PostStepDoIt()
// The tracking will message three types of GetPhysicalInteractionLength()
// in order to limit the Step and select the occurrence of processes. 
// It will message the corresponding DoIt() to apply the selected 
// processes. In addition, the Tracking will limit the Step
// and select the occurrence of the processes according to
// the shortest physical interaction length computed (except for
// processes at rest, for which the Tracking will select the
// occurrence of the process which returns the shortest mean
// life-time from the GetPhysicalInteractionLength()).

// Authors:
// - 02.12.1995, G.Cosmo - First implementation, based on object model
// - 06.05.1996, G.Cosmo - Revised; added vector of processes at rest
// - 08.01.1997, H.Kurashige - New Physics scheme
// --------------------------------------------------------------------
#ifndef G4ProcessManager_hh
#define G4ProcessManager_hh 1

#include <vector>

#include "globals.hh"
#include "G4ios.hh"

#include "G4VProcess.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleDefinition.hh"

class G4ProcessManagerMessenger;
class G4ProcessAttribute;

// Indexes for ProcessVector
//
enum G4ProcessVectorTypeIndex
{ 
  typeGPIL = 0,        // for GetPhysicalInteractionLength 
  typeDoIt =1          // for DoIt
};
enum G4ProcessVectorDoItIndex
{
  idxAll = -1,         // for all DoIt/GPIL 
  idxAtRest = 0,       // for AtRestDoIt/GPIL
  idxAlongStep = 1,    // for AlongStepDoIt/GPIL
  idxPostStep =2,      // for AlongSTepDoIt/GPIL
  NDoit =3
};

// Enumeration for Ordering Parameter
//
enum G4ProcessVectorOrdering
{ 
  ordInActive = -1,    // ordering parameter to indicate InActive DoIt
  ordDefault = 1000,   // default ordering parameter
  ordLast    = 9999    // ordering parameter to indicate the last DoIt
};

class G4ProcessManager 
{
  using G4ProcessAttrVector = std::vector<G4ProcessAttribute*>; 

  public: 

    G4ProcessManager(const G4ParticleDefinition* aParticleType);
      //  Constructor

    G4ProcessManager(G4ProcessManager& right);
      // copy constructor

    G4ProcessManager() = delete;
    G4ProcessManager& operator=(const G4ProcessManager&) = delete;
      // Default constructor and assignment operator not allowed

    ~G4ProcessManager();
      // Destructor

    G4bool operator==(const G4ProcessManager &right) const;
    G4bool operator!=(const G4ProcessManager &right) const;

    inline G4ProcessVector* GetProcessList() const;
      // Returns the address of the vector of all processes 

    inline G4int GetProcessListLength() const;
      // Returns the number of process in the ProcessVector 

    inline G4int GetProcessIndex(G4VProcess*) const;
      // Returns the index of the process in the process List

    inline G4ProcessVector* GetProcessVector( 
                               G4ProcessVectorDoItIndex idx,
                               G4ProcessVectorTypeIndex typ = typeGPIL
                              ) const;
      // Returns the address of the vector of processes 

    inline G4ProcessVector* GetAtRestProcessVector(
                               G4ProcessVectorTypeIndex typ = typeGPIL
                              ) const; 
      // Returns the address of the vector of processes for
      //    AtRestGetPhysicalInteractionLength      idx =0
      //    AtRestGetPhysicalDoIt                   idx =1

    inline G4ProcessVector* GetAlongStepProcessVector(
                               G4ProcessVectorTypeIndex typ = typeGPIL
                              ) const;
      // Returns the address of the vector of processes for
      //    AlongStepGetPhysicalInteractionLength      idx =0
      //    AlongStepGetPhysicalDoIt                   idx =1

    inline G4ProcessVector* GetPostStepProcessVector(
                               G4ProcessVectorTypeIndex typ = typeGPIL
                              ) const;
      // Returns the address of the vector of processes for
      //    PostStepGetPhysicalInteractionLength      idx =0
      //    PostStepGetPhysicalDoIt                   idx =1

    G4int GetProcessVectorIndex(
                           G4VProcess* aProcess,
                           G4ProcessVectorDoItIndex idx,
                           G4ProcessVectorTypeIndex typ  = typeGPIL
                           ) const;
    inline G4int GetAtRestIndex(
                           G4VProcess* aProcess,
                           G4ProcessVectorTypeIndex typ  = typeGPIL
                           ) const;
    inline G4int GetAlongStepIndex(
                           G4VProcess* aProcess,
                           G4ProcessVectorTypeIndex typ  = typeGPIL
                           ) const;
    inline G4int GetPostStepIndex(
                           G4VProcess* aProcess,
                           G4ProcessVectorTypeIndex typ = typeGPIL
                           ) const;
      // Returns the index for GPIL/DoIt process vector of the process  

    G4int AddProcess( G4VProcess* aProcess,
                      G4int ordAtRestDoIt = ordInActive,
                      G4int ordAlongSteptDoIt = ordInActive,
                      G4int ordPostStepDoIt = ordInActive );
      // Adds a process to the process List
      // Return values is the index to the List. Negative return value 
      // indicates that the process has not been added due to some errors
      // The first argument is a pointer to the process.
      // Successive arguments are ordering parameters of the process in 
      // process vectors. If value is negative, the process is
      // not added to the corresponding process vector
   
    /////////////////////////////////////////////// 
    // The following methods are provided for simple processes  
    //   AtRestProcess has only AtRestDoIt
    //   ContinuousProcess has only AlongStepDoIt
    //   DiscreteProcess has only PostStepDoIt
    // If the ordering parameter is not specified, the process is
    // added at the end of List of process vectors 
    // If a process with same ordering parameter exists, 
    // this new process will be added just after processes 
    // with same ordering parameter (except for processes assigned to LAST
    // explicitly) for both DoIt() and GetPhysicalInteractionLength()
    /////////////////////////////////////////////// 

    inline G4int AddRestProcess(G4VProcess* aProcess, G4int ord=ordDefault);
    inline G4int AddDiscreteProcess(G4VProcess* aProcess, G4int ord=ordDefault);
    inline G4int AddContinuousProcess(G4VProcess* aProcess, G4int ord=ordDefault);

    /////////////////////////////////////////////// 
    // Alternative methods for setting ordering parameters
    // Note: AddProcess() method should precede calls to these methods
    /////////////////////////////////////////////// 

    G4int GetProcessOrdering(
                            G4VProcess* aProcess,
                            G4ProcessVectorDoItIndex idDoIt
                            );

    void SetProcessOrdering(
                            G4VProcess* aProcess,
                            G4ProcessVectorDoItIndex idDoIt,
                            G4int ordDoIt = ordDefault
                           );
      // Set ordering parameter for DoIt() specified by typeDoIt.
      // If a process with same ordering parameter exists, 
      // this new process will be added just after processes 
      // with same ordering parameter  
      // Note: Ordering parameter will bet set to non-zero 
      //       even if you set ordDoIt = 0
            
    void SetProcessOrderingToFirst(
                            G4VProcess* aProcess,
                            G4ProcessVectorDoItIndex idDoIt
                           );
      // Set ordering parameter to the first of all processes 
      // for DoIt() specified by idDoIt.
      // Note: If you use this method for two processes,
      //       a process called later will be first

    void SetProcessOrderingToSecond(
                            G4VProcess* aProcess,
                            G4ProcessVectorDoItIndex idDoIt
                           );
      // Set ordering parameter to 1 for DoIt() specified by idDoIt
      // and the process will be added just after 
      // the processes with ordering parameter equal to zero
      // Note: If you use this method for two processes,
      //       a process called later will be first

    void SetProcessOrderingToLast(
                            G4VProcess* aProcess,
                            G4ProcessVectorDoItIndex idDoIt
                           );
      // Set ordering parameter to the last of all processes 
      // for DoIt() specified by idDoIt.
      // Note: If you use this method for two processes,
      //       a process called later will precede.

    /////////////////////////////////////////////// 

    G4VProcess* RemoveProcess(G4VProcess* aProcess);
    G4VProcess* RemoveProcess(G4int index);
      // Removes a process from the process List.
      // Returns pointer to the removed process.
      // (nullptr value will be returned in case of errors)

    G4VProcess* SetProcessActivation(G4VProcess* aProcess, G4bool fActive);
    G4VProcess* SetProcessActivation(G4int index, G4bool fActive);
      // Set activation flag. 
      // Returns pointer to the applied process.
      // (nullptr value will be returned in case of errors)

    G4bool GetProcessActivation(G4VProcess* aProcess) const;
    G4bool GetProcessActivation(G4int index) const;
      // Get activation flag. 

    inline G4ParticleDefinition* GetParticleType() const;
      // Get the particle type 
    inline void SetParticleType(const G4ParticleDefinition*);
      // Set the particle type

    G4VProcess* GetProcess (const G4String&) const;
      // Get process by process name

    void StartTracking(G4Track* aTrack = nullptr);
    void EndTracking();
      // These two methods are used by G4TrackingManager 
      // in order to inform Start/End of tracking for each track
      // to the process manager and all physics processes

    void DumpInfo();

    inline void SetVerboseLevel(G4int value);
    inline G4int GetVerboseLevel() const;
      // Control flag for output message
      //  0: Silent
      //  1: Warning message
      //  2: More

    enum {SizeOfProcVectorArray = 6};

  private:

    G4int InsertAt(G4int position, G4VProcess* process, G4int ivec);
      // Insert process at position in theProcVector[ivec]

    G4int RemoveAt(G4int position, G4VProcess* process, G4int ivec);
      // Remove process at position in theProcVector[ivec]

    G4int FindInsertPosition(G4int ord, G4int ivec);
      // Find insert position according to ordering parameter 
      // in theProcVector[ivec]

    inline G4int GetProcessVectorId(G4ProcessVectorDoItIndex idx,
                             G4ProcessVectorTypeIndex typ = typeGPIL) const;

    void CheckOrderingParameters(G4VProcess*) const;
      // Check consistencies between ordering parameters and 
      // validity of DoIt() of the Process 

    G4ProcessAttribute* GetAttribute(G4int index) const;
    G4ProcessAttribute* GetAttribute(G4VProcess* aProcess) const;
      // Get Pointer to ProcessAttribute

    void  CreateGPILvectors();
    void  SetIndexToProcessVector(G4int ivec);

    G4VProcess* ActivateProcess(G4int index);
    G4VProcess* InActivateProcess(G4int index);
      // Activate/InActivate process
      
  private:

    G4ProcessVector* theProcVector[SizeOfProcVectorArray];
      // Vector for processes with GetPhysicalInteractionLength()/DoIt()

    G4ProcessAttrVector* theAttrVector = nullptr;
      // Vector for process attribute  

    const G4ParticleDefinition* theParticleType = nullptr;
      // Particle which has this process manager object     

    G4int numberOfProcesses = 0;
    G4ProcessVector* theProcessList = nullptr;
      // Vector for all processes (called as "process List")

    G4bool duringTracking = false;

    G4bool isSetOrderingFirstInvoked[NDoit];
    G4bool isSetOrderingLastInvoked[NDoit];

    G4int verboseLevel = 1;
  
    static G4ThreadLocal G4ProcessManagerMessenger* fProcessManagerMessenger;
    static G4ThreadLocal G4int counterOfObjects;
};

#include "G4ProcessManager.icc"

#endif
