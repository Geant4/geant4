// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ProcTblElement.hh,v 1.7 2000-03-02 01:45:10 kurasige Exp $
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
//   Use STL vector instead of RW vector    1. Mar 00 H.Kurashige
//
// Class Description  
//  This class is used by G4ProcessTable ONLY for booking !!!
//

#ifndef G4ProcTblElement_h
#define G4ProcTblElement_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/vector"

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
    // Use STL Vector 
    typedef G4std::vector<G4ProcessManager*> G4ProcMgrVector;

    G4int Length() const ;
    void  Insert(G4ProcessManager* aProcMgr);
    void  Remove(G4ProcessManager* aProcMgr);

    G4VProcess*       GetProcess() const;
    const G4String&   GetProcessName() const;
 
    G4ProcessManager* GetProcessManager(G4int index) const;

    inline 
     const G4ProcMgrVector*   GetProcMgrVector() const
      { return pProcMgrVector;}
  
    G4int      GetIndex(const G4ProcessManager* pManager) const ;
    G4bool     Contains(const G4ProcessManager* pManager) const ;

  private:
    G4VProcess*       pProcess;
    // pointer to G4VProcess

    G4ProcMgrVector*  pProcMgrVector;
};

#include "G4ProcTblElement.icc"
#endif




