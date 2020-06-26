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
// G4ProcTblElement
//
// Class description:
//
// This class is exclusively used by G4ProcessTable for booking.

// Author: H.Kurashige, 4 August 1998
// --------------------------------------------------------------------
#ifndef G4ProcTblElement_hh
#define G4ProcTblElement_hh 1

#include <vector>

#include "globals.hh"
#include "G4ios.hh"

#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"
#include "G4ProcessManager.hh"

class G4ProcTblElement
{
  friend class G4ProcessTable;

  using G4ProcMgrVector = std::vector<G4ProcessManager*>;

  public:

    G4ProcTblElement(const G4ProcTblElement& right);
    G4ProcTblElement(G4VProcess* aProcess);
      // Constructors

    ~G4ProcTblElement();
      // Destructor

    G4ProcTblElement& operator=(const G4ProcTblElement& right);
      // Assignment operator

    G4bool operator==(const G4ProcTblElement& right) const;
    G4bool operator!=(const G4ProcTblElement& right) const;
      // Equality operators

  protected:

    G4ProcTblElement();

    inline G4int Length() const ;
    inline void Insert(G4ProcessManager* aProcMgr);
    inline void Remove(G4ProcessManager* aProcMgr);

    inline G4VProcess*       GetProcess() const;
    inline const G4String&   GetProcessName() const;
    inline G4ProcessManager* GetProcessManager(G4int index) const;

    inline const G4ProcMgrVector* GetProcMgrVector() const; 
  
    inline G4int  GetIndex(const G4ProcessManager* pManager) const;
    inline G4bool Contains(const G4ProcessManager* pManager) const;

  private:

    G4VProcess* pProcess = nullptr;
    G4ProcMgrVector* pProcMgrVector = nullptr;
};

#include "G4ProcTblElement.icc"

#endif
