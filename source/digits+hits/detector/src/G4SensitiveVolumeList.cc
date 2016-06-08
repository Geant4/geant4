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
//
//
// $Id: G4SensitiveVolumeList.cc,v 1.5 2001/07/13 15:00:09 gcosmo Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//
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
     if (thePhysicalVolumeList.size()==0) return false;
     for(size_t i=0;i<thePhysicalVolumeList.size();i++)
     { if(thePhysicalVolumeList[i]==pvp) return true; }
     return false;  
   }


   G4bool G4SensitiveVolumeList::CheckLV(const G4LogicalVolume * lvp) const
   {
    if (theLogicalVolumeList.size()==0) return false;
      for(size_t i=0;i<theLogicalVolumeList.size();i++)
      { if(theLogicalVolumeList[i]==lvp) return true; }
    return false;  
   }


