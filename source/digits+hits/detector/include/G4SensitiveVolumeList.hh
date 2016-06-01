// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SensitiveVolumeList.hh,v 2.1 1998/07/12 02:53:20 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
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

#include <rw/tpordvec.h>
#include <rw/tvordvec.h>
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"


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
      const RWTPtrOrderedVector<G4VPhysicalVolume>& GetThePhysicalVolumeList() const;
      void SetThePhysicalVolumeList(const RWTPtrOrderedVector<G4VPhysicalVolume> value);

      const RWTPtrOrderedVector<G4LogicalVolume>& GetTheLogicalVolumeList() const;
      void SetTheLogicalVolumeList(const RWTPtrOrderedVector<G4LogicalVolume> value);



  private: 

    //Data Members for Has Relationships

      RWTPtrOrderedVector<G4VPhysicalVolume> thePhysicalVolumeList;
      RWTPtrOrderedVector<G4LogicalVolume> theLogicalVolumeList;

};


#endif

