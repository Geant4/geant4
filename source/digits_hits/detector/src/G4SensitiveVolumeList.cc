//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4SensitiveVolumeList.cc 94771 2015-12-09 09:44:05Z gcosmo $
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
    G4SensitiveVolumeList & G4SensitiveVolumeList::operator=(const G4SensitiveVolumeList &right)
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
     for(auto pv : thePhysicalVolumeList) 
     { if(pv==pvp) return true; }
     return false;  
   }


   G4bool G4SensitiveVolumeList::CheckLV(const G4LogicalVolume * lvp) const
   {
     for(auto lv : theLogicalVolumeList)
     { if(lv==lvp) return true; }
     return false;  
   }


