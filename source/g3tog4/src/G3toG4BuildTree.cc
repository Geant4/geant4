#include "globals.hh"
#include "G3toG4BuildTree.hh"
#include "G3VolTable.hh"
#include "VolTableEntry.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G3Pos.hh"
#include "G3RotTable.hh"
void
G3toG4BuildTree(VolTableEntry* CurVTE){
  // create logical volume
  G4LogicalVolume* CurLog = 
    new G4LogicalVolume(CurVTE->GetSolid(),   // solids pointer
			CurVTE->GetMaterial(),// material pointer
			CurVTE->GetName());; // lv name
  CurVTE->SetLV(CurLog);
  
  for (int ICopy=0; ICopy < CurVTE->NPCopies(); ICopy++){
    G3Pos* TheG3Pos = CurVTE->GetG3PosCopy(ICopy); 
    // pointer to mother
    G4LogicalVolume* MothLV = CurVTE->GetMother()->GetLV();
    // rotation matrix
    G4int irot = TheG3Pos->GetIrot();
    //G4cout << "Positioning PV " << TheG3Pos->GetName() << " irot " << irot
    //   << " inside mother " << MothLV->GetName() << " copy "
    //   << TheG3Pos->GetCopy() << endl;
    // Position it
    new G4PVPlacement((G4RotationMatrix*) 
		      G3Rot.Get(TheG3Pos->GetIrot()),   // rotation matrix
		      *(TheG3Pos->GetPos()),            // its position
		      CurLog,                           // its LogicalVolume 
		      TheG3Pos->GetName(),              // PV name
		      MothLV,                           // Mother LV
		      0,                                // only
		      TheG3Pos->GetCopy());             // copy
  }
  // number of G3Pos daughters
  G4int Ndau = CurVTE->GetNoDaughters();
  for (int Idau=0; Idau<Ndau; Idau++){
    G3toG4BuildTree(CurVTE->GetDaughter(Idau));
  }
};

