// $Id: G4HCIOentryT.hh,v 1.1 2002-11-24 13:45:23 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4HCIOentryT.hh
//
// History:
//   '01.09.12  Youhei Morita  Initial creation

#ifndef HCIO_ENTRY_T_HH
#define HCIO_ENTRY_T_HH 1

#include <string>
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
           std::cout << "G4HCIOentryT: Registering HitsCollection IO manager"
                << " for \"" << n << "\"" <<  std::endl;
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
            std::cout << "G4HCIOentryT: Constructing HitsCollection IO manager"
                << " for \"" << detName << "\" " << f_manager <<  std::endl;
          }
          G4HCIOcatalog::GetG4HCIOcatalog()->RegisterHCIOmanager(f_manager);
          if ( m_verbose > 2 ) {
            G4HCIOcatalog::GetG4HCIOcatalog()->PrintHCIOmanager();
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

