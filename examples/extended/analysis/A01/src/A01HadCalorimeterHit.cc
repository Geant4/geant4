// $Id: A01HadCalorimeterHit.cc,v 1.1 2002-11-13 07:23:26 duns Exp $
// --------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#include "A01HadCalorimeterHit.hh"
#include "A01DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4RotationMatrix.hh"
#include "G4Box.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4ios.hh"

G4Allocator<A01HadCalorimeterHit> A01HadCalorimeterHitAllocator;

A01HadCalorimeterHit::A01HadCalorimeterHit()
{
  columnID = -1;
  rowID = -1;
  edep = 0.;
}

A01HadCalorimeterHit::A01HadCalorimeterHit(G4int iCol,G4int iRow)
{
  columnID = iCol;
  rowID = iRow;
  edep = 0.;
}

A01HadCalorimeterHit::~A01HadCalorimeterHit()
{;}

A01HadCalorimeterHit::A01HadCalorimeterHit(const A01HadCalorimeterHit &right)
{
  columnID = right.columnID;
  rowID = right.rowID;
  edep = right.edep;
  pos = right.pos;
  rot = right.rot;
}

const A01HadCalorimeterHit& A01HadCalorimeterHit::operator=(const A01HadCalorimeterHit &right)
{
  columnID = right.columnID;
  rowID = right.rowID;
  edep = right.edep;
  pos = right.pos;
  rot = right.rot;
  return *this;
}

int A01HadCalorimeterHit::operator==(const A01HadCalorimeterHit &right) const
{
  return (columnID==right.columnID&&rowID==right.rowID);
}

void A01HadCalorimeterHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager&&(edep>0.))
  {
    // Draw a calorimeter cell with depth propotional to the energy deposition
    G4Transform3D trans(rot.inverse(),pos);
    G4VisAttributes attribs;
    G4Colour colour(1.,0.,0.);
    attribs.SetColour(colour);
    attribs.SetForceWireframe(false);
    attribs.SetForceSolid(true);
    G4Box box("dummy",15.*cm,15.*cm,1.*m*edep/(0.1*GeV));
    pVVisManager->Draw(box,attribs,trans);
  }
}

void A01HadCalorimeterHit::Print()
{
  G4cout << "  Cell[" << rowID << ", " << columnID << "] " << edep/MeV << " (MeV) " << pos << G4endl;
}


