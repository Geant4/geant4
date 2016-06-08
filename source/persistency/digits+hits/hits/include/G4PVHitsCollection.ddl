// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PVHitsCollection.ddl,v 1.14 2000/12/15 08:04:14 morita Exp $
// GEANT4 tag $Name: geant4-03-01 $
//

// Class Description:
//   This is a persistent version of container class which stores
// the associations to the hits of one sensitive detector.
// User should inherit his/her own persistent hits collection
// from this class to store user derived hits.
//   At the end of each event, user sensitive detector should
// create a new collection with a detector name and a collection
// name.  Comparison operator and GetName(), GetSDName() are
// provided.
//

#ifndef G4PVHitsCollection_h
#define G4PVHitsCollection_h 1

#include "G4Pglobals.hh"
#include "G4PersistentTypes.hh"
#include "G4PersistentSchema.hh"

#include "HepODBMS/odbms/HepODBMS.h"

class G4PVHitsCollection 
 : public HepPersObj
{
  public: // with description
      G4PVHitsCollection(G4String detName,G4String colNam);
      // Constructor.
      ~G4PVHitsCollection();
      // Destructor.
      int operator==(const G4PVHitsCollection &right) const;
      // Comparison with detName and colName.

  protected:

      // Persistent Collection name
      G4PString pcollectionName;
      G4PString pSDname;

  public: // with description
      inline G4String GetName()
      { return (G4String) pcollectionName; }
      // returns the collection name.
      inline G4String GetSDname()
      { return (G4String) pSDname; }
      // returns the detector name.

};

#endif

