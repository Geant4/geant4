// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3toG4BuildTree.cc,v 1.9 1999-12-09 00:04:59 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// modified by I. Hrivnacova, 2.8.99 

#include "G3toG4BuildTree.hh"
#include "G3RotTable.hh"
#include "G3MedTable.hh"
#include "G3VolTable.hh"
#include "G3SensVolVector.hh"
#include "G3Pos.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G3toG4Debug.hh"

void G3toG4BuildTree(G3VolTableEntry* curVTE, G3VolTableEntry* motherVTE)
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
        G4Exception("VTE " + curVTE->GetName() + " has not defined material!!");
      } 

      // logical volume
      curLog = 
        new G4LogicalVolume(curVTE->GetSolid(), material, curVTE->GetName());
      curVTE->SetLV(curLog);
      
      // insert logical volume to G3SensVol vector
      // in case it is sensitive
      if (mte->GetISVOL()) G3SensVol.insert(curLog);
    }  
  
    // get mother logical volume
    G4LogicalVolume* mothLV;
    if (motherVTE) {
      G4String motherName = motherVTE->GetName();    
      if (curVTE->FindMother(motherName)->GetName() != motherName) {
        // check consistency - tbr
        G4Exception("G3toG4BuildTree: Inconsistent mother <-> daughter !!");
      }
      mothLV = motherVTE->GetLV();
    }  
    else {  	    
       mothLV = 0;
    }  
    
    // positions in motherVTE
    for (G4int i=0; i<curVTE->NPCopies(); i++){

      G3Pos* theG3Pos = curVTE->GetG3PosCopy(i);
      if (theG3Pos->GetMotherName() == motherVTE->GetMasterClone()->GetName()) {

        // rotation matrix
        G4int irot = theG3Pos->GetIrot();
        G4RotationMatrix* theMatrix = 0;
        if (irot>0) theMatrix = G3Rot.Get(irot);

        // copy number
        // (in G3 numbering starts from 1 but in G4 from 0)
        G4int copyNo = theG3Pos->GetCopy() - 1;
      
        // position it if not top-level volume
	if (mothLV != 0) {
	  new G4PVPlacement(theMatrix,          // rotation matrix
			    *(theG3Pos->GetPos()),  // its position
			    curLog,                 // its LogicalVolume 
			    curVTE->GetName(),      // PV name
			    mothLV,                 // Mother LV
			    0,                      // only
			    copyNo);                // copy
	
        // verbose
	  
	  if (G3toG4Debug != 0) 
	    G4cout << "PV: " << i << "th copy of " << curVTE->GetName()
		   << "  in " << motherVTE->GetName() << "  copyNo: " 
		   << copyNo << "  irot: " << irot << "  pos: " 
		   << *(theG3Pos->GetPos()) << endl;
	}
      }
    }

    // divisions     
    if (curVTE->GetDivision()) {
      curVTE->GetDivision()->CreatePVReplica();
      // verbose
      if (G3toG4Debug != 0) {
	G4cout << "CreatePVReplica: " << curVTE->GetName() 
	       << " in "  << motherVTE->GetName() << endl;
      }
    }
  }
  else {
    if ( !(curVTE->GetDivision() && motherVTE->GetMasterClone() == motherVTE &&
           motherVTE->GetNoClones()>1)) {
      // ignore dummy vte's 
      // (this should be the only case when the vte is dummy but
      // is present in mother <-> daughters tree
      G4Exception("VTE " + curVTE->GetName() + " has not defined solid!!");
    }
  }  
  
  // process daughters
  G4int Ndau = curVTE->GetNoDaughters();
  for (int Idau=0; Idau<Ndau; Idau++){
    G3toG4BuildTree(curVTE->GetDaughter(Idau), curVTE);
  }
}

