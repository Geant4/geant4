//    ********************************
//    *                              *
//    *     BrachyWaterBoxHit.cc     *
//    *                              *
//    ********************************

#include "BrachyWaterBoxHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"

G4Allocator<BrachyWaterBoxHit> BrachyWaterBoxHitAllocator;

//....

BrachyWaterBoxHit::BrachyWaterBoxHit(G4LogicalVolume* logVol,G4int XID,G4int ZID)
:m_pLogV(logVol),m_XID(XID),m_ZID(ZID)
{
 m_Edep=0;
}

//....

BrachyWaterBoxHit::~BrachyWaterBoxHit()
{
}

//....

BrachyWaterBoxHit::BrachyWaterBoxHit(const BrachyWaterBoxHit &right)
{
 m_XID = right.m_XID;
 m_ZID = right.m_ZID;
 m_Edep = right.m_Edep;
 m_Pos = right.m_Pos;
 m_Rot = right.m_Rot;
 m_pLogV = right.m_pLogV;
}

//....

const BrachyWaterBoxHit& BrachyWaterBoxHit::operator=(const BrachyWaterBoxHit &right)
{
 m_XID = right.m_XID;
 m_ZID = right.m_ZID;
 m_Edep = right.m_Edep;
 m_Pos = right.m_Pos;
 m_Rot = right.m_Rot;
 m_pLogV = right.m_pLogV;
 return *this;
}

//....

int BrachyWaterBoxHit::operator==(const BrachyWaterBoxHit &right) const
{
 return((m_XID==right.m_XID)&&(m_ZID==right.m_ZID));
}

//....

void BrachyWaterBoxHit::Draw()
{
}

//....

void BrachyWaterBoxHit::Print()
{
}
