
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ProcessTable.hh,v 1.5 1999-11-11 15:37:52 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, IT Division, ASD group
//	History: first implementation, based on object model of
//	4th Aug 1998, H.Kurashige
//
// Class Description
//  This class is used for "book keeping" of all processes 
//  which are registered in all particles
//
//  History:
//      Added G4ProcessTableMesseneger      16 Aug. 1998, H.Kurashige
//
// ------------------------------------------------------------

#ifndef G4ProcessTable_h
#define G4ProcessTable_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "g4rw/tpordvec.h"
#include "g4rw/tvordvec.h"

#include "G4ProcTblElement.hh"
#include "G4ProcessVector.hh"
class G4ProcessTableMessenger;

class G4ProcessTable
{
 public:
  G4ProcessTable();
  //  Constructors
  
  ~G4ProcessTable();
  //  Destructor
  
 private:
  G4ProcessTable(const G4ProcessTable &right);
  G4ProcessTable & operator=(const G4ProcessTable &right);
  // Assignment operator
  G4int operator==(const G4ProcessTable &right) const;
  G4int operator!=(const G4ProcessTable &right) const;
  // equal / unequal operator
  
 
 public: // with description
  static G4ProcessTable* GetProcessTable();
  // return the pointer to G4ProcessTable object
  //   G4ProcessTable is a "singleton" and can get its pointer by this function

  G4int  Length() const;
  // return the number of processes in the table

  G4int  Insert(G4VProcess* aProcess, G4ProcessManager* aProcMgr);
  G4int  Remove(G4VProcess* aProcess, G4ProcessManager* aProcMgr);  
  // insert and remove methods
  //  each process object is registered with information of process managers
  //  that use it.

  G4VProcess* FindProcess(const G4String& processName, 
			  const G4String& particleName) const;
  G4VProcess* FindProcess(const G4String& processName, 
			  const G4ParticleDefinition* particle) const;
  G4VProcess* FindProcess(const G4String& processName, 
			  const G4ProcessManager* processManager) const;
  // return the process pointer   
  
  G4ProcessVector* FindProcesses();
  G4ProcessVector* FindProcesses( const G4ProcessManager* processManager );
  G4ProcessVector* FindProcesses( const G4String& processName );
  G4ProcessVector* FindProcesses( G4ProcessType   processType );
  // return pointer of a process vector 
  //  which includes processes specified
  //  Note:: User is responsible to delete this process vector object  

  void SetProcessActivation( const G4String& processName, 
			     G4bool          fActive);
  void SetProcessActivation( const G4String& processName, 
		             const G4String& particleName, 
			     G4bool          fActive );
  void SetProcessActivation( const G4String& processName, 
		             G4ParticleDefinition* particle, 
			     G4bool          fActive );
  void SetProcessActivation( const G4String& processName, 
		             G4ProcessManager* processManager, 
			     G4bool          fActive  );
  void SetProcessActivation( G4ProcessType   processType, 
			     G4bool          fActive  );
  void SetProcessActivation( G4ProcessType   processType,
		             const G4String& particleName, 
			     G4bool          fActive  );
  void SetProcessActivation( G4ProcessType   processType,
		             G4ParticleDefinition* particle, 
			     G4bool          fActive );
  void SetProcessActivation( G4ProcessType   processType,
		             G4ProcessManager* processManager, 
			     G4bool          fActive  );
  // These methods are provided to activate or inactivate processes

 public:
  typedef G4RWTPtrOrderedVector<G4ProcTblElement>  G4ProcTableVector;
  typedef G4RWTValOrderedVector<G4String> G4ProcNameVector;

 public: // with description
  G4ProcNameVector*  GetNameList();
  // return pointer of the list of process name

  G4ProcTableVector* GetProcTableVector();
  // return pointer of the vector of G4ProcTblElement
  
 private:
  G4ProcTableVector* Find(  G4ProcTableVector* procTableVector,
			    const G4String& processName );
  G4ProcTableVector* Find(  G4ProcTableVector* procTableVector,
			    G4ProcessType   processType );
  // return pointer of a ProcTableVector 
  //  which includes ProcTbleElement specified

  G4ProcessVector*   ExtractProcesses( G4ProcTableVector* procTableVector);
  // extract all process objects from the process table 
 
 public: // with description
  void DumpInfo(G4VProcess* process, G4ParticleDefinition* particle=0);
  // dump out information of the process table
  //  second argument is used to specify processes designated by a particle 
  
 public: // with description
   G4UImessenger* CreateMessenger();
   void           DeleteMessenger();
   // These methods are used by RunManager to let the process table
   // know the timing of creation/destructuion of  messengers
  
 public: // with description
   void  SetVerboseLevel(G4int value);
   G4int GetVerboseLevel() const;
   // Set/Get controle flag for output message
   //  0: Silent
   //  1: Warning message
   //  2: More


 private:
  static G4ProcessTable*    fProcessTable;
  G4ProcessTableMessenger*   fProcTblMessenger;

 private:
  G4ProcTableVector*        fProcTblVector;
  G4ProcNameVector*         fProcNameVector;
  // list of G4ProcTblElement
  
  G4ProcTableVector*        tmpTblVector;
  // used only internaly for temporary buffer.

 private:
   G4int verboseLevel;
   // controle flag for output message

};

inline 
 void  G4ProcessTable::SetVerboseLevel(G4int value)
{
  verboseLevel = value;
}

inline 
 G4int G4ProcessTable::GetVerboseLevel() const
{
  return verboseLevel;
}

inline 
 G4int  G4ProcessTable::Length() const
{
  return fProcTblVector->length();
}

#include "G4ProcessTable.icc"
#endif
