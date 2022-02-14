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
// G4DCIOcatalog
//
// Class Description:
//
// Catalog for the I/O manager of digits collection for each detector

// Author: Youhei Morita, 12.09.2001
// --------------------------------------------------------------------
#ifndef G4DCIOCATALOG_HH
#define G4DCIOCATALOG_HH 1

#include <map>
#include "G4Types.hh"
#include "G4VPDigitsCollectionIO.hh"

class G4VDCIOentry;

using DCIOmap = std::map<G4String, G4VDCIOentry*,
                         std::less<G4String>>;
using DCIOstore = std::map<G4String, G4VPDigitsCollectionIO*,
                           std::less<G4String>>;

class G4DCIOcatalog
{
  public:

    G4DCIOcatalog();
      // Constructor

    virtual ~G4DCIOcatalog() {}
      // Destructor

    static G4DCIOcatalog* GetDCIOcatalog();
      // Construct G4DCIOcatalog and returns the pointer

    void SetVerboseLevel(int v) { m_verbose = v; }
      // Set verbose level

    void RegisterEntry(G4VDCIOentry* d);
      // Register I/O manager entry

    void RegisterDCIOmanager(G4VPDigitsCollectionIO* d);
      // Register I/O manager

    G4VDCIOentry* GetEntry(const G4String& name);
      // Returns the I/O manager entry

    G4VPDigitsCollectionIO* GetDCIOmanager(const G4String& name);
      // Returns the registered I/O manager entry

    void PrintEntries();
      // Prints the list of I/O manager entries

    G4String CurrentDCIOmanager();
      // Returns the list of I/O managers

    void PrintDCIOmanager();
      // Prints the list of I/O managers

    std::size_t NumberOfDCIOmanager() { return theStore.size(); }
      // Returns the number of registered I/O managers

    G4VPDigitsCollectionIO* GetDCIOmanager(G4int n);
      // Returns the n-th registered I/O manager entry

  private:

    G4int m_verbose = 0;
    static G4ThreadLocal G4DCIOcatalog* f_thePointer;
    DCIOmap theCatalog;
    DCIOstore theStore;
};

#endif
