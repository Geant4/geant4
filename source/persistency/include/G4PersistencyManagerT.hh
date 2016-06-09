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
// File: G4PersistencyManagerT.hh
//
// History:
//   '01.08.10  Youhei Morita  Initial creation (with "fadsclass3")

#ifndef PERSISTENCY_MANAGER_T_HH
#define PERSISTENCY_MANAGER_T_HH 1

#include "G4PersistencyCenter.hh"

// Class inherited:
#include "G4PersistencyManager.hh"

// Class Description:
//   Template class of G4PersistencyManager for late binding

template <class T> class G4PersistencyManagerT
 : public G4PersistencyManager
{
    public: // With description
      G4PersistencyManagerT(G4PersistencyCenter* pc, std::string n)
       : G4PersistencyManager(pc, n), pm(0)
      {
         if ( m_verbose > 2 ) {
           G4cout << "G4PersistencyManagerT: Registering G4PersistencyManager \""
                  << n << "\"" << G4endl;
         }
         G4PersistencyCenter::GetPersistencyCenter()->
                       RegisterPersistencyManager(this);
      }
      // Constructor

      ~G4PersistencyManagerT() {};
      // Destructor

    public: // With description
      G4PersistencyManager* Create()
      {
        pm = new T(f_pc, GetName());
        return pm;
      }
      // Create a new persistency manager

      void Delete() { if (pm!=0) delete pm; };
      // Delete a persistency mamanger

      G4VPEventIO* EventIO()
      {
        if (pm) return pm->EventIO();
        else    return 0;
      };
      // Returns the current event I/O handling manager

      G4VPHitIO* HitIO()
      {
        if (pm) return pm->HitIO();
        else    return 0;
      };
      // Returns the current hit I/O handling manager

      G4VPDigitIO* DigitIO()
      {
        if (pm) return pm->DigitIO();
        else    return 0;
      };
      // Returns the current digit I/O handling manager

      G4VHepMCIO* HepMCIO()
      {
        if (pm) return pm->HepMCIO();
        else    return 0;
      };
      // Returns the current digit I/O handling manager

      G4VMCTruthIO* MCTruthIO()
      {
        if (pm) return pm->MCTruthIO();
        else    return 0;
      };
      // Returns the current digit I/O handling manager

      G4VTransactionManager* TransactionManager()
      {
        if (pm) return pm->TransactionManager();
        else    return 0;
      };
      // Returns the current transaction manager

      void Initialize() {};
      // Initialize the persistency package.

      void SetVerboseLevel(int v)
      {
        if (pm) return pm->SetVerboseLevel();
        else    return 0;
      };
      // Sets the verbose level to the persistency manager

    private:
      G4PersistencyManager* pm;

}; // End of class G4PersistencyManagerT

#endif

