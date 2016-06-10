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
// $Id: G3toG4BuildTree.cc 67982 2013-03-13 10:36:03Z gcosmo $
//
// modified by I. Hrivnacova, 2.8.99 

#include "globals.hh"
#include "G3toG4BuildTree.hh"
#include "G3RotTable.hh"
#include "G3MedTable.hh"
#include "G3VolTable.hh"
#include "G3SensVolVector.hh"
#include "G3Pos.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ReflectionFactory.hh"
#include "G4Transform3D.hh"

void G3toG4BuildTree(G3VolTableEntry* curVTE, G3VolTableEntry* motherVTE)
{
  G3toG4BuildLVTree(curVTE, motherVTE);
  G3toG4BuildPVTree(curVTE);
}  

void G3toG4BuildLVTree(G3VolTableEntry* curVTE, G3VolTableEntry* motherVTE)
{
  // check existence of the solid
  if (curVTE->GetSolid()) {
    G4LogicalVolume* curLog = curVTE->GetLV();
    if (!curLog) {
      // skip creating logical volume
      // in case it already exists
    
      // material
      G4Material* material = 0;
      G3MedTableEntry* mte = G3Med.get(curVTE->GetNmed());
      if (mte) material = mte->GetMaterial();
      if (!material) {
        G4String err_message = "VTE " + curVTE->GetName()
                             + " has not defined material!!";
        G4Exception("G3toG4BuildLVTree()", "G3toG40001",
                    FatalException, err_message);
        return;
      } 

      // logical volume
      curLog = 
        new G4LogicalVolume(curVTE->GetSolid(), material, curVTE->GetName());
      curVTE->SetLV(curLog);
      
      // insert logical volume to G3SensVol vector
      // in case it is sensitive
      if (mte->GetISVOL()) G3SensVol.push_back(curLog);
    }  
  }
  else {
    if ( !(curVTE->GetDivision() && motherVTE->GetMasterClone() == motherVTE &&
           motherVTE->GetNoClones()>1)) {
      // ignore dummy vte's 
      // (this should be the only case when the vte is dummy but
      // is present in mother <-> daughters tree
      G4String err_message = "VTE " + curVTE->GetName()
                           + " has not defined solid!!";
      G4Exception("G3toG4BuildLVTree()", "G3toG40002",
                  FatalException, err_message);
      return;
    }
  }  
  
  // process daughters
  G4int Ndau = curVTE->GetNoDaughters();
  for (int Idau=0; Idau<Ndau; Idau++){
    G3toG4BuildLVTree(curVTE->GetDaughter(Idau), curVTE);
  }
}

void G3toG4BuildPVTree(G3VolTableEntry* curVTE)
{
  // check existence of the solid
  if (curVTE->GetSolid()) {
    G4LogicalVolume* curLog = curVTE->GetLV();
  
    // positions in motherVTE
    for (G4int i=0; i<curVTE->NPCopies(); i++){

      G3Pos* theG3Pos = curVTE->GetG3PosCopy(i);
      if (theG3Pos) {

        // loop over all mothers 
        for (G4int im=0; im<curVTE->GetNoMothers(); im++) {

          G3VolTableEntry* motherVTE = curVTE->GetMother(im);
          if (theG3Pos->GetMotherName() == motherVTE->GetMasterClone()->GetName()) {
     
            // get mother logical volume
            G4LogicalVolume* mothLV=0;
            G4String motherName = motherVTE->GetName();    
            if (!curVTE->FindMother(motherName)) continue;
            if (curVTE->FindMother(motherName)->GetName() != motherName) {
              // check consistency - tbr
              G4String err_message =
                       "G3toG4BuildTree: Inconsistent mother <-> daughter !!";
              G4Exception("G3toG4BuildPVTree()", "G3toG40003",
                          FatalException, err_message);
              return;
            }
            mothLV = motherVTE->GetLV();
    
            // copy number
            // (in G3 numbering starts from 1 but in G4 from 0)
            G4int copyNo = theG3Pos->GetCopy() - 1;
      
            // position it if not top-level volume

	    if (mothLV != 0) {

              // transformation
              G4int irot = theG3Pos->GetIrot();
              G4RotationMatrix* theMatrix = 0;
              if (irot>0) theMatrix = G3Rot.Get(irot);
              G4Rotate3D rotation;
              if (theMatrix) {            
  	        rotation = G4Rotate3D(*theMatrix);
	      }

              #ifndef G3G4_NO_REFLECTION
              G4Translate3D translation(*(theG3Pos->GetPos()));
	      G4Transform3D transform3D = translation * (rotation.inverse());

              G4ReflectionFactory::Instance()
	        ->Place(transform3D,       // transformation
	                curVTE->GetName(), // PV name
			curLog,            // its logical volume 
			mothLV,            // mother logical volume
			false,             // only
			copyNo);           // copy
              #else
              new G4PVPlacement(theMatrix,            // rotation matrix
                              *(theG3Pos->GetPos()),  // its position
                              curLog,                 // its LogicalVolume 
                              curVTE->GetName(),      // PV name
                              mothLV,                 // Mother LV
                              0,                      // only
                              copyNo);                // copy
              #endif

              // verbose
  	      #ifdef G3G4DEBUG
	        G4cout << "PV: " << i << "th copy of " << curVTE->GetName()
	  	       << "  in " << motherVTE->GetName() << "  copyNo: " 
		       << copyNo << "  irot: " << irot << "  pos: " 
		       << *(theG3Pos->GetPos()) << G4endl;
	      #endif
            }	    
	  }
        }
        // clear this position
        curVTE->ClearG3PosCopy(i);
        i--;
      }
    }

    // divisions     
    if (curVTE->GetDivision()) {
      curVTE->GetDivision()->CreatePVReplica();
      // verbose
      #ifdef G3G4DEBUG
	G4cout << "CreatePVReplica: " << curVTE->GetName() 
	       << " in "  <<  curVTE->GetMother()->GetName() << G4endl;
      #endif

      // clear this divison
      curVTE->ClearDivision(); 
    }
  }
  
  // process daughters
  G4int Ndau = curVTE->GetNoDaughters();
  for (int Idau=0; Idau<Ndau; Idau++){
    G3toG4BuildPVTree(curVTE->GetDaughter(Idau));
  }
}

