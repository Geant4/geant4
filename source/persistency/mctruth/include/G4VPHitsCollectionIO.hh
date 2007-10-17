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
// File: G4VPHitsCollectionIO.hh
//
// History:
//   '01.08.16  Youhei Morita  Initial creation

#ifndef V_P_HITS_COLLECTION_I_O_HH
#define V_P_HITS_COLLECTION_I_O_HH 1

#include "G4VHitsCollection.hh"

// Class Description:
//   Abstract base class for storing and retrieving hit collections

class G4VPHitsCollectionIO
{
    public: // With description
      G4VPHitsCollectionIO( std::string detName, std::string colName );
      // Constructor

      virtual ~G4VPHitsCollectionIO() {};
      // Destructor

    public: // With description
      virtual G4bool Store(const G4VHitsCollection*) =0;
      // Pure virtual method for storing the hit collection.
      // Each persistency package should implement a concrete method
      // with this signature.

      virtual G4bool Retrieve(G4VHitsCollection*&) =0;
      // Pure virtual method for retrieving the hit collection.
      // Each persistency package should implement a concrete method
      // with this signature.

      G4bool operator== (const G4VPHitsCollectionIO& right) const;
      // virtual operator for comparing hit collections with names.

      std::string SDname() { return f_detName; };
      // Returns the sensitive detector name.

      std::string CollectionName() { return f_colName; };
      // Returns the hit collection name.

      void SetVerboseLevel(int v) { m_verbose = v; };
      // Sets the verbose level

    protected:
      G4int m_verbose;
      std::string f_detName;
      std::string f_colName;

}; // End of class G4VPHitsCollectionIO

#endif

