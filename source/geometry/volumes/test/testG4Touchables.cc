// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testG4Touchables.cc,v 1.2 1999-12-15 14:50:29 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Test file for G4TouchableHistory/G4GRSVolume/G4GRSSolid
// Paul Kent August 1996

#include "ApproxEqual.hh"
#include "G4ThreeVector.hh"

#include "G4Navigator.hh"
#include "G4GRSVolume.hh"
#include "G4GRSSolid.hh"
#include "G4TouchableHistory.hh"

#include "G4GeometryManager.hh"

#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"

// Sample Parameterisation
class G4LinScale : public G4VPVParameterisation
{
  virtual void ComputeTransformation(const G4int n,
				     G4VPhysicalVolume *pRep) const
  {
    pRep->SetTranslation(G4ThreeVector(n*4-2,n*4-2,n*4-2));
  }
  
  virtual void ComputeDimensions(G4Box &pBox,
				 const G4int n,
				 const G4VPhysicalVolume *pRep) const
  {
    pBox.SetXHalfLength(n+1);
    pBox.SetYHalfLength(n+1);
    pBox.SetZHalfLength(n+1);
  }
};

G4LinScale myParam;

G4VPhysicalVolume* BuildGeometry()
{
  G4Box *worldBox=new G4Box("WorldBox",1000,1000,1000);
  G4Box *posBox=new G4Box("PosBox",10,10,10);
  G4Box *posBoxSlice=new G4Box("PosBoxSlice",10,10,2);
  G4Box *paramBox=new G4Box("ParamBox",1,1,1);

  G4LogicalVolume *worldLog=new G4LogicalVolume(worldBox,0,
						"WorldLog",0,0,0);
  G4LogicalVolume *posLog=new G4LogicalVolume(posBox,0,
					      "PosLog",0,0,0);
  G4LogicalVolume *pos2Log=new G4LogicalVolume(posBox,0,
					      "Pos2Log",0,0,0);
  G4LogicalVolume *pos3Log=new G4LogicalVolume(posBox,0,
					      "Pos3Log",0,0,0);
  G4LogicalVolume *pos4Log=new G4LogicalVolume(posBox,0,
					      "Pos4Log",0,0,0);
  G4LogicalVolume *posSliceLog=new G4LogicalVolume(posBoxSlice,0,
						   "PosSliceLog",0,0,0);
  G4LogicalVolume *paramLog=new G4LogicalVolume(paramBox,0,
						"ParamLog",0,0,0);

  G4PVPlacement *worldPhys=new G4PVPlacement(0,G4ThreeVector(0,0,0),
					    "WorldPhys",worldLog,0,false,0);
  G4PVPlacement *posPhys1=new G4PVPlacement(0,G4ThreeVector(10,11,12),
					    "PosPhys1",posLog,worldPhys,
					    false,0);
  G4PVPlacement *posPhys2=new G4PVPlacement(0,G4ThreeVector(0,0,0),
					    "PosPhys2",pos2Log,posPhys1,
					    false,0);
  G4PVPlacement *posPhys3=new G4PVPlacement(0,G4ThreeVector(-10,-11,-12),
					    "PosPhys3",pos3Log,worldPhys,
					    false,0);
  G4PVPlacement *posPhys4=new G4PVPlacement(0,G4ThreeVector(10,0,0),
					    "PosPhys4",pos4Log,worldPhys,
					    false,0);
  G4PVReplica *repPhys=new G4PVReplica("RepPhys",posSliceLog,posPhys3,
				      kZAxis,5,4);

  G4PVParameterised *paramPhys=new G4PVParameterised("ParamPhys",
						     paramLog,
						     posPhys4,
						     kXAxis,2,&myParam);
  return worldPhys;
}


G4bool testGRSVolume(G4Navigator& nav)
{
  G4ThreeVector pos;

  pos=G4ThreeVector(11,12,13);
  nav.LocateGlobalPointAndSetup(pos);
  G4GRSVolume *touch=nav.CreateGRSVolume();
  assert(touch);
  assert(touch->GetVolume()->GetName()=="PosPhys2");
  assert(touch->GetSolid()->GetName()=="PosBox");
  assert(ApproxEqual(touch->GetTranslation(),G4ThreeVector(10,11,12)));
  assert(!touch->GetRotation()||touch->GetRotation()->isIdentity());
  delete touch;
  return true;
}

G4bool testGRSSolid(G4Navigator& nav)
{
  G4ThreeVector pos;

  pos=G4ThreeVector(11,12,13);
  nav.LocateGlobalPointAndSetup(pos);
  G4GRSSolid *touch=nav.CreateGRSSolid();
  assert(touch);
  assert(touch->GetSolid()->GetName()=="PosBox");
  assert(ApproxEqual(touch->GetTranslation(),G4ThreeVector(10,11,12)));
  assert(!touch->GetRotation()||touch->GetRotation()->isIdentity());
  delete touch;
  return true;
}

G4bool testTouchableHistory(G4Navigator& nav)
{
  G4ThreeVector pos(11,12,13),pos2(9,-1,-1);
  G4VPhysicalVolume *pvol;

  pvol=nav.LocateGlobalPointAndSetup(pos);
  assert(pvol->GetName()=="PosPhys2");
  assert(pvol->GetMother()->GetName()=="PosPhys1");
  assert(pvol->GetMother()->GetMother()->GetName()=="WorldPhys");

  G4TouchableHistory *touch=nav.CreateTouchableHistory();
  assert(touch);
  assert(touch->GetVolume()->GetName()=="PosPhys2");
  assert(touch->GetSolid()->GetName()=="PosBox");
  assert(ApproxEqual(touch->GetTranslation(),G4ThreeVector(10,11,12)));
  assert(!touch->GetRotation()||touch->GetRotation()->isIdentity());
  assert(touch->GetHistory()->GetDepth()==2);

  pvol=nav.LocateGlobalPointAndSetup(-pos);
  assert(pvol->GetName()=="RepPhys");
  assert(pvol->GetMother()->GetName()=="PosPhys3");
  assert(pvol->GetMother()->GetMother()->GetName()=="WorldPhys");
  assert(ApproxEqual(nav.GetCurrentLocalCoordinate(),G4ThreeVector(-1,-1,-1)));

  G4TouchableHistory *touch2=nav.CreateTouchableHistory();
  assert(touch2);
  assert(touch2->GetVolume()->GetName()=="RepPhys");
  assert(touch2->GetSolid()->GetName()=="PosBoxSlice");
  assert(ApproxEqual(touch2->GetTranslation(),G4ThreeVector(-10,-11,-12)));
  assert(!touch2->GetRotation()||touch2->GetRotation()->isIdentity());
  assert(touch2->GetHistory()->GetDepth()==2);

  pvol=nav.LocateGlobalPointAndSetup(pos2);
  assert(pvol->GetName()=="ParamPhys");
  assert(pvol->GetMother()->GetName()=="PosPhys4");
  assert(pvol->GetMother()->GetMother()->GetName()=="WorldPhys");
  assert(ApproxEqual(nav.GetCurrentLocalCoordinate(),G4ThreeVector(1,1,1)));

  G4TouchableHistory *touch3=nav.CreateTouchableHistory();
  // Relocate to another parameterised volume causing modification of
  // physical volume + solid
  pvol=nav.LocateGlobalPointAndSetup(pos2+G4ThreeVector(2,2,2));
  assert(touch3);
  assert(touch3->GetVolume()->GetName()=="ParamPhys");
  assert(touch3->GetSolid()->GetName()=="ParamBox");
  assert(ApproxEqual(touch3->GetTranslation(),G4ThreeVector(8,-2,-2)));
  assert(!touch3->GetRotation()||touch3->GetRotation()->isIdentity());
  assert(touch3->GetHistory()->GetDepth()==2);
  
  pvol=nav.LocateGlobalPointAndSetup(pos,*touch);
  assert(ApproxEqual(nav.GetCurrentLocalCoordinate(),G4ThreeVector(1,1,1)));
  assert(pvol->GetName()=="PosPhys2");
  assert(pvol->GetMother()->GetName()=="PosPhys1");
  assert(pvol->GetMother()->GetMother()->GetName()=="WorldPhys");
  assert(ApproxEqual(nav.NetTranslation(),G4ThreeVector(10,11,12)));

  pvol=nav.LocateGlobalPointAndSetup(-pos,*touch2);
  assert(ApproxEqual(nav.GetCurrentLocalCoordinate(),G4ThreeVector(-1,-1,-1)));
  assert(pvol->GetName()=="RepPhys");
  assert(pvol->GetMother()->GetName()=="PosPhys3");
  assert(pvol->GetMother()->GetMother()->GetName()=="WorldPhys");
  assert(ApproxEqual(nav.NetTranslation(),G4ThreeVector(-10,-11,-12)));

  pvol=nav.LocateGlobalPointAndSetup(pos2,*touch3);
  assert(pvol->GetName()=="ParamPhys");
  assert(pvol->GetMother()->GetName()=="PosPhys4");
  assert(pvol->GetMother()->GetMother()->GetName()=="WorldPhys");
  assert(ApproxEqual(nav.NetTranslation(),G4ThreeVector(8,-2,-2)));

  delete touch;
  delete touch2;
  delete touch3;
  return true;
}

int main()
{
  G4VPhysicalVolume *myTopNode=BuildGeometry();
  G4GeometryManager::GetInstance()->CloseGeometry(false);
  G4Navigator nav;
  nav.SetWorldVolume(myTopNode);

  assert(testGRSVolume(nav));
  assert(testGRSSolid(nav));
  assert(testTouchableHistory(nav));

  G4GeometryManager::GetInstance()->OpenGeometry();
  return EXIT_SUCCESS;
}
