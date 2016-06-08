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

