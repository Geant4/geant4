// $Id: G4VPDigitsCollectionIO.hh,v 1.1 2002-11-24 13:45:23 morita Exp $
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
      G4VPDigitsCollectionIO( std::string detName, std::string colName );
      // Constructor

      virtual ~G4VPDigitsCollectionIO() {};
      // Destructor

    public: // With description
      virtual bool Store(const G4VDigiCollection*) =0;
      // Pure virtual method for storing the digit collection.
      // Each persistency package should implement a concrete method
      // with this signature.

      virtual bool Retrieve(G4VDigiCollection*&) =0;
      // Pure virtual method for retrieving the digit collection.
      // Each persistency package should implement a concrete method
      // with this signature.

      bool operator== (const G4VPDigitsCollectionIO& right) const;
      // virtual operator for comparing digit collections with names.

      std::string DMname() { return f_detName; };
      // Returns the digitizer module name.

      std::string CollectionName() { return f_colName; };
      // Returns the digit collection name.

      void SetVerboseLevel(int v) { m_verbose = v; };
      // Sets the verbose level

    protected:
      int m_verbose;
      std::string f_detName;
      std::string f_colName;

}; // End of class G4VPDigitsCollectionIO

#endif

