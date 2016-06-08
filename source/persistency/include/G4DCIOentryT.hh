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
// File: G4DCIOentryT.hh
//
// History:
//   '01.09.12  Youhei Morita  Initial creation

#ifndef DCIO_ENTRY_T_HH
#define DCIO_ENTRY_T_HH 1

#include <string>
#include "G4Types.hh"
#include "G4VPDigitsCollectionIO.hh"

// Class inherited:
#include "G4VDCIOentry.hh"

// Class Description:
//   Template class of DigitsCollection I/O Manager for late binding

template <class T> class G4DCIOentryT
 : public G4VDCIOentry
{
    public: // With description
      G4DCIOentryT<T>(G4std::string n)
       : G4VDCIOentry(n), f_manager(0)
      {
         if ( m_verbose > 2 ) {
           G4cout << "G4DCIOentryT: Registering DigitsCollection IO manager"
                  << " for \"" << n << "\"" <<  G4endl;
         }
      }
      // Constructor

      ~G4DCIOentryT() {};
      // Destructor

    public: // With description
      void CreateDCIOmanager(G4std::string detName, G4std::string colName)
      {
        if ( f_manager == 0 ) {
          f_manager = new T( detName, colName );
          if ( m_verbose > 2 ) {
            G4cout << "G4DCIOentryT: Constructing DigitsCollection IO manager"
                   << " for \"" << detName << "\" " << f_manager <<  G4endl;
          }
          G4DCIOcatalog::GetDCIOcatalog()->RegisterDCIOmanager(f_manager);
          if ( m_verbose > 2 ) {
            G4DCIOcatalog::GetDCIOcatalog()->PrintDCIOmanager();
          }
        }
      }
      // Create a new digits collection I/O manager

      void DeleteDCIOmanager() { if (f_manager!=0) delete f_manager; };
      // Delete a digits collection I/O manager

    private:
      G4VPDigitsCollectionIO* f_manager;

}; // End of class G4DCIOentryT

#endif

