// $Id: G4HCIOcatalog.hh,v 1.2 2002-12-04 10:25:49 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4HCIOcatalog.hh
//
// History:
//   '01.09.12  Youhei Morita  Initial creation

#ifndef HCIO_CATALOG_HH
#define HCIO_CATALOG_HH 1

#include "g4std/map"
#include "G4Types.hh"
#include "G4VPHitsCollectionIO.hh"

class G4VHCIOentry;

typedef G4std::map<G4std::string, G4VHCIOentry*, G4std::less<G4std::string> > HCIOmap;

typedef G4std::map<G4std::string, G4VPHitsCollectionIO*, G4std::less<G4std::string> > HCIOstore;

// Class Description:
//   Catalog for the I/O manager of hits collection for each detector.

class G4HCIOcatalog
{
    public: // With description
      G4HCIOcatalog();
      // Constructor

      virtual ~G4HCIOcatalog() {};
      // Destructor

    public: // With description
      static G4HCIOcatalog* GetG4HCIOcatalog();
      // Construct G4HCIOcatalog and returns the pointer

      void SetVerboseLevel(int v) { m_verbose = v; };
      // Set verbose level.

      void RegisterEntry(G4VHCIOentry* d);
      // Register I/O manager entry

      void RegisterHCIOmanager(G4VPHitsCollectionIO* d);
      // Register I/O manager

      G4VHCIOentry* GetEntry(G4std::string name);
      // Returns the I/O manager entry

      G4VPHitsCollectionIO* GetHCIOmanager(G4std::string name);
      // Returns the registered I/O manager entry

      void PrintEntries();
      // Prints the list of I/O manager entries

      G4std::string CurrentHCIOmanager();
      // Returns the list of I/O managers

      void PrintHCIOmanager();
      // Prints the list of I/O managers

      size_t NumberOfHCIOmanager() { return theStore.size(); };
      // Returns the number of registered I/O managers.

      G4VPHitsCollectionIO* GetHCIOmanager(int n);
      // Returns the n-th registered I/O manager entry

    private:
      int m_verbose;
      static G4HCIOcatalog* f_thePointer;
      HCIOmap theCatalog;
      HCIOstore theStore;

}; // End of class G4HCIOcatalog

#endif

