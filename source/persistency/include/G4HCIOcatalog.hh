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
//
// File: G4HCIOcatalog.hh
//
// History:
//   '01.09.12  Youhei Morita  Initial creation

#ifndef HCIO_CATALOG_HH
#define HCIO_CATALOG_HH 1

#include <map>
#include "G4Types.hh"
#include "G4VPHitsCollectionIO.hh"

class G4VHCIOentry;

typedef std::map<std::string, G4VHCIOentry*, std::less<std::string> > HCIOmap;

typedef std::map<std::string, G4VPHitsCollectionIO*, std::less<std::string> > HCIOstore;

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
      static G4HCIOcatalog* GetHCIOcatalog();
      // Construct G4HCIOcatalog and returns the pointer

      void SetVerboseLevel(int v) { m_verbose = v; };
      // Set verbose level.

      void RegisterEntry(G4VHCIOentry* d);
      // Register I/O manager entry

      void RegisterHCIOmanager(G4VPHitsCollectionIO* d);
      // Register I/O manager

      G4VHCIOentry* GetEntry(std::string name);
      // Returns the I/O manager entry

      G4VPHitsCollectionIO* GetHCIOmanager(std::string name);
      // Returns the registered I/O manager entry

      void PrintEntries();
      // Prints the list of I/O manager entries

      std::string CurrentHCIOmanager();
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

