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
// $Id: G4SmartVoxelStat.cc,v 1.1 2001-10-22 16:08:05 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------
// GEANT 4 class source file
//
// class G4SmartVoxelStat
//
// Store information on the performance of the smart voxel algorithm
// for an individual logical volume
//
// Author: D.C.Williams, UCSC (davidw@scipp.ucsc.edu)
//
#include "G4SmartVoxelStat.hh"

#include "G4SmartVoxelHeader.hh"
#include "G4SmartVoxelNode.hh"
#include "G4SmartVoxelProxy.hh"

//
// Constructor
//
G4SmartVoxelStat::G4SmartVoxelStat( const G4LogicalVolume *theVolume,
                                    const G4SmartVoxelHeader *theVoxel,
                                          G4double theSysTime,
                                          G4double theUserTime )
    : volume(theVolume),
      voxel(theVoxel),
      sysTime(theSysTime),
      userTime(theUserTime),
      heads(1),
      nodes(0),
      pointers(0)
{
  CountHeadsAndNodes( voxel );
}


//
// Simple Accessors
//
const G4LogicalVolume *G4SmartVoxelStat::GetVolume() const
{
  return volume;
}

const G4SmartVoxelHeader *G4SmartVoxelStat::GetVoxel() const
{
  return voxel;
}

G4double G4SmartVoxelStat::GetSysTime() const
{
  return sysTime;
}

G4double G4SmartVoxelStat::GetUserTime() const
{
  return userTime;
}

G4double G4SmartVoxelStat::GetTotalTime() const
{
  return sysTime + userTime;
}

G4long G4SmartVoxelStat::GetNumberHeads() const
{
  return heads;
}

G4long G4SmartVoxelStat::GetNumberNodes() const
{
  return nodes;
}

G4long G4SmartVoxelStat::GetNumberPointers() const
{
  return pointers;
}


//
// Return approximate memory use
//
G4long G4SmartVoxelStat::GetMemoryUse() const
{
  static const G4long headSize = sizeof(G4SmartVoxelHeader)
                               + sizeof(G4SmartVoxelProxy);

  static const G4long nodeSize = sizeof(G4SmartVoxelNode)
                               + sizeof(G4SmartVoxelProxy);

  static const G4long pointerSize = sizeof(G4SmartVoxelProxy*);
  
  return nodes*nodeSize + heads*headSize + pointers*pointerSize;
}



//
// CountHeadsAndNodes
//
// Recursively count the number of voxel headers and nodes,
// updating class member variables heads and nodes on the way.
//
void G4SmartVoxelStat::CountHeadsAndNodes( const G4SmartVoxelHeader *head ) 
{
  G4int numSlices = head->GetNoSlices();
  
  pointers += numSlices;
  
  const G4SmartVoxelProxy *lastProxy = 0;
  
  for(G4int i=0;i<numSlices;++i) {
    const G4SmartVoxelProxy *proxy = head->GetSlice(i);
    if (proxy == lastProxy) continue;
    
    lastProxy = proxy;
    
    if (proxy->IsNode()) {
      nodes++;
    }
    else {
      heads++;
      CountHeadsAndNodes(proxy->GetHeader());
    }
  }
}
