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
// $Id: G4ProcessTable.hh 71231 2013-06-12 13:06:28Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	History: first implementation, based on object model of
//	4th Aug 1998, H.Kurashige
//
// Class Description
//  This class is used for "book keeping" of all processes 
//  which are registered in all particles
//
//  History:
//   Added G4ProcessTableMesseneger         16 Aug. 1998, H.Kurashige
//   Use STL vector instead of RW vector    1. Mar 00 H.Kurashige
//
// ------------------------------------------------------------

#ifndef G4ProcessTable_h
#define G4ProcessTable_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <vector>

#include "G4ProcTblElement.hh"
#include "G4ProcessVector.hh"
class G4UImessenger;
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
  typedef std::vector<G4ProcTblElement*>  G4ProcTableVector;
  typedef std::vector<G4String> G4ProcNameVector;

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
  static G4ThreadLocal G4ProcessTable*    fProcessTable;
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

#include "G4ProcessTable.icc"
#endif
