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
      G4HCIOentryT<T>(std::string n)
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
      void CreateHCIOmanager(std::string detName, std::string colName)
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

