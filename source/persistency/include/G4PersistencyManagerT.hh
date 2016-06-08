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
      G4PersistencyManagerT(G4PersistencyCenter* pc, G4std::string n)
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

