// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ProcTblElement.hh,v 1.3 1999-06-01 14:52:35 gcosmo Exp $
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
  friend class G4ProcessTable;
  protected:
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

  protected:
    typedef RWTPtrOrderedVector<G4ProcessManager> G4ProcMgrVector;

    inline G4int Length() const ;
    inline void  Insert(G4ProcessManager* aProcMgr);
    inline void  Remove(G4ProcessManager* aProcMgr);

    inline G4VProcess*       GetProcess() const;
    inline const G4String&   GetProcessName() const;
 
    inline G4ProcessManager* GetProcessManager(G4int index) const;
    const G4ProcMgrVector*   GetProcMgrVector() const { return pProcMgrVector; }
  
    inline G4int      GetIndex(const G4ProcessManager* pManager) const ;
    inline G4bool     Contains(const G4ProcessManager* pManager) const ;

  private:
    G4VProcess*       pProcess;
    // pointer to G4VProcess

    G4ProcMgrVector*  pProcMgrVector;
};

#include "G4ProcTblElement.icc"
#endif




