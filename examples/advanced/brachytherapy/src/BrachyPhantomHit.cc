//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// Code developed by:
// S. Agostinelli, F. Foppiano, S. Garelli , M. Tropeano, S.Guatelli
//
//    ********************************
//    *                              *
//    *     BrachyPhantomHit.cc     *
//    *                              *
//    ********************************
//
// $Id: BrachyPhantomHit.cc,v 1.3 2002-11-18 15:18:38 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "BrachyPhantomHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"

G4Allocator<BrachyPhantomHit> BrachyPhantomHitAllocator;

//....

BrachyPhantomHit::BrachyPhantomHit(G4LogicalVolume* logVol,G4int XID,G4int YID,G4int ZID)
:m_pLogV(logVol),m_XID(XID),m_ZID(ZID),m_YID(YID)
{
 m_Edep=0;
}

//....

BrachyPhantomHit::~BrachyPhantomHit()
{
}

//....

BrachyPhantomHit::BrachyPhantomHit(const BrachyPhantomHit &right)
{
 m_XID = right.m_XID;
 m_ZID = right.m_ZID;
 m_YID = right.m_YID;
 m_Edep = right.m_Edep;
 m_Pos = right.m_Pos;
 m_Rot = right.m_Rot;
 m_pLogV = right.m_pLogV;
}

//....

const BrachyPhantomHit& BrachyPhantomHit::operator=(const BrachyPhantomHit &right)
{
 m_XID = right.m_XID;
 m_ZID = right.m_ZID;
 m_YID = right.m_YID;
 m_Edep = right.m_Edep;
 m_Pos = right.m_Pos;
 m_Rot = right.m_Rot;
 m_pLogV = right.m_pLogV;
 return *this;
}

//....

int BrachyPhantomHit::operator==(const BrachyPhantomHit &right) const
{
 return((m_XID==right.m_XID)&&(m_ZID==right.m_ZID)&&(m_YID==right.m_YID));
}

//....

void BrachyPhantomHit::Draw()
{
}

//....

void BrachyPhantomHit::Print()
{
}
