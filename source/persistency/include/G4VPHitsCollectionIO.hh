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
      G4VPHitsCollectionIO( G4std::string detName, G4std::string colName );
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

      G4std::string SDname() { return f_detName; };
      // Returns the sensitive detector name.

      G4std::string CollectionName() { return f_colName; };
      // Returns the hit collection name.

      void SetVerboseLevel(int v) { m_verbose = v; };
      // Sets the verbose level

    protected:
      G4int m_verbose;
      G4std::string f_detName;
      G4std::string f_colName;

}; // End of class G4VPHitsCollectionIO

#endif

