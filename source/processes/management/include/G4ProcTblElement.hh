//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4ProcTblElement.hh,v 1.9 2001-07-11 10:08:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
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

    G4ProcTblElement & operator=(const G4ProcTblElement &right);
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




