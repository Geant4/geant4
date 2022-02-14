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
// G4VPDigitsCollectionIO
//
// Class Description:
//
// Abstract base class for storing and retrieving digit collections.

// Author: Youhei Morita, 16.08.2001
// --------------------------------------------------------------------
#ifndef G4VPDIGITSCOLLECTIONIO_HH
#define G4VPDIGITSCOLLECTIONIO_HH 1

#include "G4VDigiCollection.hh"

class G4VPDigitsCollectionIO
{
  public:

    G4VPDigitsCollectionIO(const G4String& detName, const G4String& colName);
      // Constructor

    virtual ~G4VPDigitsCollectionIO() {}
      // Destructor

    virtual G4bool Store(const G4VDigiCollection*) = 0;
      // Pure virtual method for storing the digit collection.
      // Each persistency package should implement a concrete method
      // with this signature

    virtual G4bool Retrieve(G4VDigiCollection*&) = 0;
      // Pure virtual method for retrieving the digit collection.
      // Each persistency package should implement a concrete method
      // with this signature

    G4bool operator==(const G4VPDigitsCollectionIO& right) const;
      // Virtual operator for comparing digit collections with names

    const G4String& DMname() { return f_detName; }
      // Returns the digitizer module name

    const G4String& CollectionName() { return f_colName; }
      // Returns the digit collection name

    void SetVerboseLevel(G4int v) { m_verbose = v; }
      // Sets the verbose level

  protected:

    G4int m_verbose = 0;
    G4String f_detName;
    G4String f_colName;
};

#endif
