// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ProcTblElement.hh,v 1.1 1999-01-07 16:13:53 gunter Exp $
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
// ------------------------------------------------------------

#ifndef G4ProcTblElement_h
#define G4ProcTblElement_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <rw/tpordvec.h>

#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"
#include "G4ProcessManager.hh"

class G4ProcTblElement
{
  // this class is used by G4ProcessTable ONLY for booking !!!
  private:
    G4ProcTblElement();

  public:
    G4ProcTblElement(const G4ProcTblElement& right);
    G4ProcTblElement(G4VProcess* aProcess);
    //  Constructors

    ~G4ProcTblElement();
    //  Destructor

    G4ProcTblElement & operator=(G4ProcTblElement &right);
    // Assignment operator

    G4int operator==(const G4ProcTblElement &right) const;
    G4int operator!=(const G4ProcTblElement &right) const;
    // equal / unequal operator

    typedef RWTPtrOrderedVector<G4ProcessManager> G4ProcMgrVector;

    G4int Length() const ;
    void  Insert(G4ProcessManager* aProcMgr);
    void  Remove(G4ProcessManager* aProcMgr);

    G4VProcess*       GetProcess() const;
    G4String          GetProcessName() const;
 
    G4ProcessManager* GetProcessManager(G4int index) const;
    G4ProcMgrVector*  GetProcMgrVector() const;
  
    G4int             GetIndex(const G4ProcessManager* pManager) const ;
    G4bool            Contains(const G4ProcessManager* pManager) const ;

  private:
    G4VProcess*       pProcess;
    // pointer to G4VProcess

    G4ProcMgrVector*  pProcMgrVector;
};

#include "G4ProcTblElement.icc"
#endif




