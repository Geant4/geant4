// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SensitiveVolumeList.cc,v 2.1 1998/07/12 02:53:31 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD Group
//      History: first implementation June'96, based on Hits+Digi
//      domain model of April 1996, S.Piperov
//
// ----------------  G4SensitiveVolumeList  -----------------

#include "G4SensitiveVolumeList.hh"

 //Constructors
   G4SensitiveVolumeList::G4SensitiveVolumeList()
   {
   }

   G4SensitiveVolumeList::G4SensitiveVolumeList(const G4SensitiveVolumeList &right)
   {
     thePhysicalVolumeList = right.thePhysicalVolumeList;
     theLogicalVolumeList = right.theLogicalVolumeList;
   }


 //Destructor
   G4SensitiveVolumeList::~G4SensitiveVolumeList()
   {
   }


 //Assignment Operation
   const G4SensitiveVolumeList & G4SensitiveVolumeList::operator=(const G4SensitiveVolumeList &right)
   {
     thePhysicalVolumeList = right.thePhysicalVolumeList;
     theLogicalVolumeList = right.theLogicalVolumeList;
     return *this;
   }


 //Equality Operations
   G4int G4SensitiveVolumeList::operator==(const G4SensitiveVolumeList &right) const
   {
     return (this == (G4SensitiveVolumeList *) &right);
   }

   G4int G4SensitiveVolumeList::operator!=(const G4SensitiveVolumeList &right) const
   {
     return (this != (G4SensitiveVolumeList *) &right);
   }



 //Other Operations 
   G4bool G4SensitiveVolumeList::CheckPV(const G4VPhysicalVolume * pvp) const
   {
     if (thePhysicalVolumeList.entries()==0) return false;
     for(int i=0;i<thePhysicalVolumeList.entries();i++)
     { if(thePhysicalVolumeList(i)==pvp) return true; }
     return false;  
   }


   G4bool G4SensitiveVolumeList::CheckLV(const G4LogicalVolume * lvp) const
   {
    if (theLogicalVolumeList.entries()==0) return false;
      for(int i=0;i<theLogicalVolumeList.entries();i++)
      { if(theLogicalVolumeList(i)==lvp) return true; }
    return false;  
   }


