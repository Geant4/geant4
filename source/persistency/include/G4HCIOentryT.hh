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
// File: G4HCIOentryT.hh
//
// History:
//   '01.09.12  Youhei Morita  Initial creation

#ifndef HCIO_ENTRY_T_HH
#define HCIO_ENTRY_T_HH 1

#include <string>
#include "G4Types.hh"
#include "G4VPHitsCollectionIO.hh"

// Class inherited:
#include "G4VHCIOentry.hh"

// Class Description:
//   Template class of HitsCollection I/O Manager for late binding

template <class T> class G4HCIOentryT
 : public G4VHCIOentry
{
    public: // With description
      G4HCIOentryT<T>(G4std::string n)
       : G4VHCIOentry(n), f_manager(0)
      {
         if ( m_verbose > 2 ) {
           G4cout << "G4HCIOentryT: Registering HitsCollection IO manager"
                  << " for \"" << n << "\"" <<  G4endl;
         }
      }
      // Constructor

      ~G4HCIOentryT() {};
      // Destructor

    public: // With description
      void CreateHCIOmanager(G4std::string detName, G4std::string colName)
      {
        if ( f_manager == 0 ) {
          f_manager = new T( detName, colName );
          if ( m_verbose > 2 ) {
            G4cout << "G4HCIOentryT: Constructing HitsCollection IO manager"
                   << " for \"" << detName << "\" " << f_manager <<  G4endl;
          }
          G4HCIOcatalog::GetHCIOcatalog()->RegisterHCIOmanager(f_manager);
          if ( m_verbose > 2 ) {
            G4HCIOcatalog::GetHCIOcatalog()->PrintHCIOmanager();
          }
        }
      }
      // Create a new hits collection I/O manager

      void DeleteHCIOmanager() { if (f_manager!=0) delete f_manager; };
      // Delete a hits collection I/O manager

    private:
      G4VPHitsCollectionIO* f_manager;

}; // End of class G4HCIOentryT

#endif

