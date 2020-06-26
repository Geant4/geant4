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
// G4ProcessTable
//
// Class description:
//
// This class is used for "book keeping" of all processes 
// which are registered for all particles

// Author: H.Kurashige, 4 August 1998
// --------------------------------------------------------------------
#ifndef G4ProcessTable_hh
#define G4ProcessTable_hh 1

#include <vector>

#include "globals.hh"
#include "G4ProcTblElement.hh"
#include "G4ProcessVector.hh"
#include "G4ThreadLocalSingleton.hh"

class G4UImessenger;
class G4ProcessTableMessenger;

class G4ProcessTable
{
  friend class G4ThreadLocalSingleton<G4ProcessTable>;

  public:

    using G4ProcTableVector = std::vector<G4ProcTblElement*>;
    using G4ProcNameVector = std::vector<G4String>;

    ~G4ProcessTable();
      // Destructor
  
    G4ProcessTable(const G4ProcessTable&) = delete;
    G4ProcessTable& operator=(const G4ProcessTable&) = delete;
    G4bool operator==(const G4ProcessTable &right) const = delete;
    G4bool operator!=(const G4ProcessTable &right) const = delete;
      // Copy constructor and operators not allowed  
 
    static G4ProcessTable* GetProcessTable();
      // Return the pointer to the G4ProcessTable object
      // As "singleton" one can get the instance pointer by this function

    inline G4int Length() const;
      // Return the number of processes in the table

    G4int Insert(G4VProcess* aProcess, G4ProcessManager* aProcMgr);
    G4int Remove(G4VProcess* aProcess, G4ProcessManager* aProcMgr);  
      // Each process object is registered with information of process
      // managers that use it

    G4VProcess* FindProcess(const G4String& processName, 
                            const G4String& particleName) const;
    inline G4VProcess* FindProcess(const G4String& processName, 
                                   const G4ParticleDefinition* particle) const;
    G4VProcess* FindProcess(const G4String& processName, 
                            const G4ProcessManager* processManager) const;
    G4VProcess* FindProcess(G4ProcessType processType, 
                            const G4ParticleDefinition* particle) const;
    G4VProcess* FindProcess(G4int processSubType, 
                            const G4ParticleDefinition* particle) const;
      // Return the process pointer

    void RegisterProcess(G4VProcess*);
    void DeRegisterProcess(G4VProcess*);
      // Implementation of registration mechanism
  
    inline G4ProcessVector* FindProcesses();
    inline G4ProcessVector* FindProcesses( const G4ProcessManager* pManager );
    inline G4ProcessVector* FindProcesses( const G4String& processName );
    inline G4ProcessVector* FindProcesses( G4ProcessType processType );
      // Return pointer of a process vector which includes processes specified
      // Note: user is responsible to delete this process vector object  

    void SetProcessActivation( const G4String& processName, 
                               G4bool          fActive );
    void SetProcessActivation( const G4String& processName, 
                               const G4String& particleName, 
                               G4bool          fActive );
    inline void SetProcessActivation( const G4String& processName, 
                                      const G4ParticleDefinition* particle, 
                                      G4bool          fActive );
    void SetProcessActivation( const G4String& processName, 
                               G4ProcessManager* processManager, 
                               G4bool          fActive );
    void SetProcessActivation( G4ProcessType   processType, 
                               G4bool          fActive );
    void SetProcessActivation( G4ProcessType   processType,
                               const G4String& particleName, 
                               G4bool          fActive );
    inline void SetProcessActivation( G4ProcessType   processType,
                                      const G4ParticleDefinition* particle, 
                                      G4bool          fActive );
    void SetProcessActivation( G4ProcessType   processType,
                               G4ProcessManager* processManager, 
                               G4bool          fActive );
      // These methods are provided to activate or inactivate processes

    inline G4ProcNameVector* GetNameList();
      // Return pointer of the list of process name

    inline G4ProcTableVector* GetProcTableVector();
      // Return pointer of the vector of G4ProcTblElement

    void DumpInfo(G4VProcess* process, 
                  const G4ParticleDefinition* particle = nullptr);
      // Dump out information of the process table. The second argument
      // is used to specify processes designated by a particle 
    
    inline void  SetVerboseLevel(G4int value);
    inline G4int GetVerboseLevel() const;
      // Set/Get control flag for output message
      //  0: Silent
      //  1: Warning message
      //  2: More

  private:

    G4ProcessTable();
      // Private default constructor
  
    G4ProcTableVector* Find(const G4String& processName );
    G4ProcTableVector* Find(G4ProcessType   processType );
      // Return pointer of a ProcTableVector which includes
      // ProcTbleElement specified

    G4ProcessVector* ExtractProcesses(G4ProcTableVector*) const;
      // Extract all process objects from the process table 
 
  private:

    static G4ThreadLocal G4ProcessTable* fProcessTable;
    G4ProcessTableMessenger* fProcTblMessenger = nullptr;

    G4ProcTableVector* fProcTblVector = nullptr;
    G4ProcNameVector* fProcNameVector = nullptr;
  
    G4ProcTableVector* tmpTblVector = nullptr;
      // Used only internally as temporary buffer

    std::vector<G4VProcess*> fListProcesses;
      // Used for registration of process instances

    G4int verboseLevel = 1;
      // Control flag for output message
};

#include "G4ProcessTable.icc"

#endif
