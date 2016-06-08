// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SensitiveVolumeList.hh,v 1.2.2.1.2.1 1999/12/07 20:47:42 gunter Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
// ------------------------------------------------------------
//      GEANT 4 class header file --- Copyright CERN 1996
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on Hits+Digi domain
//      object model of April 1996, S.Piperov
//
//   ----------------  G4SensitiveVolumeList  -----------------

#ifndef G4SensitiveVolumeList_h
#define G4SensitiveVolumeList_h 1

#include "g4rw/tpordvec.h"
#include "g4rw/tvordvec.h"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

// class description:
//
//  This class object can have lists of logical and physical volumes.
// In case a sensitive detector is shared by several logical volumes and/or
// a logical volume is shared by several physical volumes, this class can be
// used by the veto list for individual logical/physical volumes.
//

class G4SensitiveVolumeList 
{

  public:
    //Constructors
      G4SensitiveVolumeList();
      G4SensitiveVolumeList(const G4SensitiveVolumeList &right);

    //Destructor
      ~G4SensitiveVolumeList();

    //Assignment Operation
      const G4SensitiveVolumeList & operator=(const G4SensitiveVolumeList &right
);

    //Equality Operations
      G4int operator==(const G4SensitiveVolumeList &right) const;
      G4int operator!=(const G4SensitiveVolumeList &right) const;


    //Other Operations
      G4bool CheckPV(const G4VPhysicalVolume *pvp) const;
      G4bool CheckLV(const G4LogicalVolume *lvp) const;

    //Get and Set Operations for Has Relationships
      const G4RWTPtrOrderedVector<G4VPhysicalVolume>& GetThePhysicalVolumeList() const;
      void SetThePhysicalVolumeList(const G4RWTPtrOrderedVector<G4VPhysicalVolume> value);

      const G4RWTPtrOrderedVector<G4LogicalVolume>& GetTheLogicalVolumeList() const;
      void SetTheLogicalVolumeList(const G4RWTPtrOrderedVector<G4LogicalVolume> value);



  private: 

    //Data Members for Has Relationships

      G4RWTPtrOrderedVector<G4VPhysicalVolume> thePhysicalVolumeList;
      G4RWTPtrOrderedVector<G4LogicalVolume> theLogicalVolumeList;

};


#endif

