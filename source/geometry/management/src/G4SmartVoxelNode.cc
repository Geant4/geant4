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
// $Id: G4SmartVoxelNode.cc,v 1.4 2002/04/19 08:20:22 gcosmo Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
// Class G4SmartVoxelNode
//
// Implementation
//

#include "G4SmartVoxelNode.hh"

// Empty destructor
//
G4SmartVoxelNode::~G4SmartVoxelNode()
{
}

// Return true if contents equal
//
// Preconditions:
//
// Node contents were entered in the same order
G4bool G4SmartVoxelNode::operator == (const G4SmartVoxelNode& v) const
{
    G4int maxNode=GetNoContained();
    if (maxNode==v.GetNoContained())
	{
	    for (G4int node=0;node<maxNode;node++)
		{
		    if (GetVolume(node)!=v.GetVolume(node))
			{
			    return false;
			}
		}
	    return true;
	}
    return false;
}
