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
// $Id: G4ProcTblElement.hh 71231 2013-06-12 13:06:28Z gcosmo $
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
#include <vector>

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
    typedef std::vector<G4ProcessManager*> G4ProcMgrVector;

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




