// $Id: G4VPHitsCollectionIO.hh,v 1.2 2002-12-04 10:25:49 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

