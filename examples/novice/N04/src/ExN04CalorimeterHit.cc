
#include "ExN04CalorimeterHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"

G4Allocator<ExN04CalorimeterHit> ExN04CalorimeterHitAllocator;

ExN04CalorimeterHit::ExN04CalorimeterHit()
{pLogV=NULL;}

ExN04CalorimeterHit::ExN04CalorimeterHit(G4LogicalVolume* logVol,G4int z,G4int phi)
:pLogV(logVol),ZCellID(z),PhiCellID(phi)
{;}

ExN04CalorimeterHit::~ExN04CalorimeterHit()
{;}

ExN04CalorimeterHit::ExN04CalorimeterHit(const ExN04CalorimeterHit &right)
{
  ZCellID = right.ZCellID;
  PhiCellID = right.PhiCellID;
  edep = right.edep;
  pos = right.pos;
  rot = right.rot;
  pLogV = right.pLogV;
}

const ExN04CalorimeterHit& ExN04CalorimeterHit::operator=(const ExN04CalorimeterHit &right)
{
  ZCellID = right.ZCellID;
  PhiCellID = right.PhiCellID;
  edep = right.edep;
  pos = right.pos;
  rot = right.rot;
  pLogV = right.pLogV;
  return *this;
}

int ExN04CalorimeterHit::operator==(const ExN04CalorimeterHit &right) const
{
  return ((ZCellID==right.ZCellID)&&(PhiCellID==right.PhiCellID));
}

void ExN04CalorimeterHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Transform3D trans(rot,pos);
    G4VisAttributes attribs;
    const G4VisAttributes* pVA = pLogV->GetVisAttributes();
    if(pVA) attribs = *pVA;
    G4Colour colour(1.,0.,0.);
    attribs.SetColour(colour);
    attribs.SetForceWireframe(false);
    attribs.SetForceSolid(true);
    pVVisManager->Draw(*pLogV,attribs,trans);
  }
}

void ExN04CalorimeterHit::Print()
{;}


