// $Id: G4DCIOcatalog.hh,v 1.2 2002-12-04 10:25:48 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4DCIOcatalog.hh
//
// History:
//   '01.09.12  Youhei Morita  Initial creation

#ifndef DCIO_CATALOG_HH
#define DCIO_CATALOG_HH 1

#include "g4std/map"
#include "G4Types.hh"
#include "G4VPDigitsCollectionIO.hh"

class G4VDCIOentry;

typedef G4std::map<G4std::string, G4VDCIOentry*, G4std::less<G4std::string> > DCIOmap;

typedef G4std::map<G4std::string, G4VPDigitsCollectionIO*, G4std::less<G4std::string> > DCIOstore;

// Class Description:
//   Catalog for the I/O manager of digits collection for each detector.

class G4DCIOcatalog
{
    public: // With description
      G4DCIOcatalog();
      // Constructor

      virtual ~G4DCIOcatalog() {};
      // Destructor

    public: // With description
      static G4DCIOcatalog* GetG4DCIOcatalog();
      // Construct G4DCIOcatalog and returns the pointer

      void SetVerboseLevel(int v) { m_verbose = v; };
      // Set verbose level.

      void RegisterEntry(G4VDCIOentry* d);
      // Register I/O manager entry

      void RegisterDCIOmanager(G4VPDigitsCollectionIO* d);
      // Register I/O manager

      G4VDCIOentry* GetEntry(G4std::string name);
      // Returns the I/O manager entry

      G4VPDigitsCollectionIO* GetDCIOmanager(G4std::string name);
      // Returns the registered I/O manager entry

      void PrintEntries();
      // Prints the list of I/O manager entries

      G4std::string CurrentDCIOmanager();
      // Returns the list of I/O managers

      void PrintDCIOmanager();
      // Prints the list of I/O managers

      size_t NumberOfDCIOmanager() { return theStore.size(); };
      // Returns the number of registered I/O managers.

      G4VPDigitsCollectionIO* GetDCIOmanager(int n);
      // Returns the n-th registered I/O manager entry

    private:
      int m_verbose;
      static G4DCIOcatalog* f_thePointer;
      DCIOmap theCatalog;
      DCIOstore theStore;

}; // End of class G4DCIOcatalog

#endif

