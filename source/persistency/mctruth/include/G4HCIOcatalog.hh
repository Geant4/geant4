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
// G4HCIOcatalog
//
// Class Description:
//
// Catalog for the I/O manager of hits collection for each detector.

// Author: Youhei Morita, 12.09.2001
// --------------------------------------------------------------------
#ifndef G4HCIOCATALOG_HH
#define G4HCIOCATALOG_HH 1

#include <map>
#include "G4Types.hh"
#include "G4String.hh"
#include "G4VPHitsCollectionIO.hh"

class G4VHCIOentry;

using HCIOmap = std::map<std::string, G4VHCIOentry*, std::less<std::string>>;

using HCIOstore = std::map<G4String, G4VPHitsCollectionIO*,
                           std::less<G4String>>;

class G4HCIOcatalog
{
  public:

    G4HCIOcatalog();
      // Constructor

    virtual ~G4HCIOcatalog() {}
      // Destructor

    static G4HCIOcatalog* GetHCIOcatalog();
      // Construct G4HCIOcatalog and returns the pointer

    void SetVerboseLevel(G4int v) { m_verbose = v; }
      // Set verbose level

    void RegisterEntry(G4VHCIOentry* d);
      // Register I/O manager entry

    void RegisterHCIOmanager(G4VPHitsCollectionIO* d);
      // Register I/O manager

    G4VHCIOentry* GetEntry(const G4String& name);
      // Returns the I/O manager entry

    G4VPHitsCollectionIO* GetHCIOmanager(const G4String& name);
      // Returns the registered I/O manager entry

    void PrintEntries();
      // Prints the list of I/O manager entries

    G4String CurrentHCIOmanager();
      // Returns the list of I/O managers

    void PrintHCIOmanager();
      // Prints the list of I/O managers

    std::size_t NumberOfHCIOmanager() { return theStore.size(); }
      // Returns the number of registered I/O managers.

    G4VPHitsCollectionIO* GetHCIOmanager(G4int n);
      // Returns the n-th registered I/O manager entry

  private:

    G4int m_verbose = 0;
    static G4ThreadLocal G4HCIOcatalog* f_thePointer;
    HCIOmap theCatalog;
    HCIOstore theStore;
};

#endif
