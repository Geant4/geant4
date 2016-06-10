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
// $Id: G4SensitiveVolumeList.hh 81087 2014-05-20 15:44:27Z gcosmo $
//
// ------------------------------------------------------------
//      GEANT 4 class header file --- Copyright CERN 1996
//      CERN Geneva Switzerland
//
//      History: first implementation, based on Hits+Digi domain
//      object model of April 1996, S.Piperov
//
//   ----------------  G4SensitiveVolumeList  -----------------

#ifndef G4SensitiveVolumeList_h
#define G4SensitiveVolumeList_h 1

//#include "g4rw/tpordvec.h"
//#include "g4rw/tvordvec.h"
#include <vector>
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
    G4SensitiveVolumeList & operator=(const G4SensitiveVolumeList &right
);

    //Equality Operations
      G4int operator==(const G4SensitiveVolumeList &right) const;
      G4int operator!=(const G4SensitiveVolumeList &right) const;


    //Other Operations
      G4bool CheckPV(const G4VPhysicalVolume *pvp) const;
      G4bool CheckLV(const G4LogicalVolume *lvp) const;

    //Get and Set Operations for Has Relationships
      const std::vector<G4VPhysicalVolume*>& GetThePhysicalVolumeList() const;
      void SetThePhysicalVolumeList(const std::vector<G4VPhysicalVolume*> value);

      const std::vector<G4LogicalVolume*>& GetTheLogicalVolumeList() const;
      void SetTheLogicalVolumeList(const std::vector<G4LogicalVolume*> value);



  private: 

    //Data Members for Has Relationships

      std::vector<G4VPhysicalVolume*> thePhysicalVolumeList;
      std::vector<G4LogicalVolume*> theLogicalVolumeList;

};


#endif

