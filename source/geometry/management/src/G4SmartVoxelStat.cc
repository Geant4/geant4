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
// G4SmartVoxelStat class implementation
//
// Stores information on the performance of the smart voxel algorithm
// for an individual logical volume.
//
// Author: D.C.Williams, UCSC (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------

#include "G4SmartVoxelStat.hh"
#include "G4SmartVoxelHeader.hh"
#include "G4SmartVoxelNode.hh"
#include "G4SmartVoxelProxy.hh"

//
// Constructor
//
G4SmartVoxelStat::G4SmartVoxelStat( const G4LogicalVolume* theVolume,
                                    const G4SmartVoxelHeader* theVoxel,
                                          G4double theSysTime,
                                          G4double theUserTime )
    : volume(theVolume),
      voxel(theVoxel),
      sysTime(theSysTime),
      userTime(theUserTime)
{
  CountHeadsAndNodes( voxel );
}

//
// Simple Accessors
//
const G4LogicalVolume* G4SmartVoxelStat::GetVolume() const
{
  return volume;
}

const G4SmartVoxelHeader* G4SmartVoxelStat::GetVoxel() const
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
void G4SmartVoxelStat::CountHeadsAndNodes( const G4SmartVoxelHeader* head ) 
{
  std::size_t numSlices = head->GetNoSlices();
  
  pointers += numSlices;
  
  const G4SmartVoxelProxy* lastProxy = nullptr;
  
  for(std::size_t i=0; i<numSlices; ++i)
  {
    const G4SmartVoxelProxy *proxy = head->GetSlice(i);
    if (proxy == lastProxy) continue;
    
    lastProxy = proxy;
    
    if (proxy->IsNode())
    {
      ++nodes;
    }
    else
    {
      ++heads;
      CountHeadsAndNodes(proxy->GetHeader());
    }
  }
}
