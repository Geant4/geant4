//
//
//    ********************************
//    *                              *
//    *         ThyroidHit.cc        *
//    *                              *
//    ********************************

#include "ThyroidHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"

G4Allocator<ThyroidHit> ThyroidHitAllocator;

//....

ThyroidHit::ThyroidHit(G4LogicalVolume*logVol):m_pLogV(logVol)
{
 m_Edep=0;
}

//....

ThyroidHit::~ThyroidHit()
{
}

//....

ThyroidHit::ThyroidHit(const ThyroidHit &right)
{
 m_XID = right.m_XID;
 m_YID = right.m_YID;
 m_ZID = right.m_ZID;
 m_Edep = right.m_Edep;
 m_Pos = right.m_Pos;
 m_Rot = right.m_Rot;
 m_pLogV = right.m_pLogV;
}

//....

const ThyroidHit& ThyroidHit::operator=(const ThyroidHit &right)
{
 m_XID = right.m_XID;
 m_YID = right.m_YID;
 m_ZID = right.m_ZID;
 m_Edep = right.m_Edep;
 m_Pos = right.m_Pos;
 m_Rot = right.m_Rot;
 m_pLogV = right.m_pLogV;
 return *this;
}

//....

int ThyroidHit::operator==(const ThyroidHit &right) const
{
 return((m_XID==right.m_XID)&&(m_YID==right.m_YID)&&(m_ZID==right.m_ZID));
}

//....

void ThyroidHit::Draw()
{
}

//....

void ThyroidHit::Print()
{
}








