// $Id: G4VPDigitsCollectionIO.hh,v 1.2 2002-12-04 10:25:49 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4VPDigitsCollectionIO.hh
//
// History:
//   '01.08.16  Youhei Morita  Initial creation

#ifndef V_P_DIGITS_COLLECTION_I_O_HH
#define V_P_DIGITS_COLLECTION_I_O_HH 1

#include "G4VDigiCollection.hh"

// Class Description:
//   Abstract base class for storing and retrieving digit collections

class G4VPDigitsCollectionIO
{
    public: // With description
      G4VPDigitsCollectionIO( G4std::string detName, G4std::string colName );
      // Constructor

      virtual ~G4VPDigitsCollectionIO() {};
      // Destructor

    public: // With description
      virtual G4bool Store(const G4VDigiCollection*) =0;
      // Pure virtual method for storing the digit collection.
      // Each persistency package should implement a concrete method
      // with this signature.

      virtual G4bool Retrieve(G4VDigiCollection*&) =0;
      // Pure virtual method for retrieving the digit collection.
      // Each persistency package should implement a concrete method
      // with this signature.

      G4bool operator== (const G4VPDigitsCollectionIO& right) const;
      // virtual operator for comparing digit collections with names.

      G4std::string DMname() { return f_detName; };
      // Returns the digitizer module name.

      G4std::string CollectionName() { return f_colName; };
      // Returns the digit collection name.

      void SetVerboseLevel(int v) { m_verbose = v; };
      // Sets the verbose level

    protected:
      G4int m_verbose;
      G4std::string f_detName;
      G4std::string f_colName;

}; // End of class G4VPDigitsCollectionIO

#endif

