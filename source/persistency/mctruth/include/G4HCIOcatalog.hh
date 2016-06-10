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
      static G4ThreadLocal G4HCIOcatalog* f_thePointer;
      HCIOmap theCatalog;
      HCIOstore theStore;

}; // End of class G4HCIOcatalog

#endif

