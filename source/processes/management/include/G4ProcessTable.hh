
// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ProcessTable.hh,v 1.2 1999-04-13 09:45:02 kurasige Exp $
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
//  History:
//      Added G4ProcessTableMesseneger      16 Aug. 1998, H.Kurashige
//
// ------------------------------------------------------------

#ifndef G4ProcessTable_h
#define G4ProcessTable_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <rw/tpordvec.h>
#include <rw/tvordvec.h>

#include "G4ProcTblElement.hh"
#include "G4ProcessVector.hh"
class G4ProcessTableMessenger;

class G4ProcessTable
{
 // this class is used by G4ProcessTable ONLY for booking !!!
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
  
 public:
  static G4ProcessTable* GetProcessTable();
  // return the pointer to G4ProcessTable object
  //   G4ProcessTable is a "singleton" and can get its pointer by this function

  G4int  Length() const;
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

 public:
  typedef RWTPtrOrderedVector<G4ProcTblElement>  G4ProcTableVector;
  typedef RWTValOrderedVector<G4String> G4ProcNameVector;

  G4ProcNameVector*  GetNameList();
  G4ProcTableVector* GetProcTableVector();

 private:
  G4ProcTableVector* Find(  G4ProcTableVector* procTableVector,
			    const G4String& processName );
  G4ProcTableVector* Find(  G4ProcTableVector* procTableVector,
			    G4ProcessType   processType );
  // return pointer of a ProcTableVector 
  //  which includes ProcTbleElement specified

  G4ProcessVector*   ExtractProcesses( G4ProcTableVector* procTableVector);
   
 public: 
  void DumpInfo(G4VProcess* process, G4ParticleDefinition* particle=0);

 public:
   G4UImessenger* CreateMessenger();
   void           DeleteMessenger();

  
 public:
   void  SetVerboseLevel(G4int value);
   G4int GetVerboseLevel() const;


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
   //  0: Silent
   //  1: Warning message
   //  2: More
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
