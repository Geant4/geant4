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
// The code was written by :
//	*Louis Archambault louis.archambault@phy.ulaval.ca,
//      *Luc Beaulieu beaulieu@phy.ulaval.ca
//      +Vincent Hubert-Tremblay at tigre.2@sympatico.ca
//
//
// *Centre Hospitalier Universitaire de Quebec (CHUQ),
// Hotel-Dieu de Quebec, departement de Radio-oncologie
// 11 cote du palais. Quebec, QC, Canada, G1R 2J6
// tel (418) 525-4444 #6720
// fax (418) 691 5268
//
// + Université Laval, Québec (QC) Canada
//*******************************************************


#include "DicomOctree.hh"
#include "DicomOctreeNode.hh"
#include "DicomOctreeMiddleNode.hh"
#include "DicomOctreeTerminalNode.hh"

#include "G4ios.hh"

DicomOctree::DicomOctree( G4int noLevels, G4double size )
{
  std::cout << "++++++ DicomOctree was instantiated now." << std::endl;
  mNoLevels = noLevels;
  mSize = size;
  mRoot = new DicomOctreeMiddleNode(0);
}

DicomOctree::~DicomOctree()
{
  delete mRoot;
  mRoot = 0; 
}

DicomOctreeNode* DicomOctree::CreateNode( G4double nodeX, G4double nodeY, G4double nodeZ, G4int level )
{
  DicomOctreeNode* current = mRoot;
  G4double currentX = 0;
  G4double currentY = 0;
  G4double currentZ = 0;
  indexChild = 8;   
  for ( G4int i = 0; i < indexChild; i++ ) // Make children
    {
      G4double childLevelResolution = (1 << (i+1) );
      G4double childSize = mSize / childLevelResolution;

      G4int dirX = G4int( nodeX >= currentX + childSize );
      G4int dirY = G4int( nodeY >= currentY + childSize );
      G4int dirZ = G4int( nodeZ >= currentZ + childSize );
      G4int direction = dirX + ( dirY << 1 ) + ( dirZ << 2 );

      if ( (*current)[direction] == 0 ) 
        {
	  if ( i < level - 1 )
            {
	      (*current)[direction] = new DicomOctreeMiddleNode( current );
            } else
	      {
                (*current)[direction] = new DicomOctreeTerminalNode( current );
	      }
        }

      current = (*current)[direction];
      currentX += dirX*childSize;
      currentY += dirY*childSize;
      currentZ += dirZ*childSize;
    }
  return current;
}

DicomOctreeNode* DicomOctree::operator()( G4double nodeX, 
                                G4double nodeY, 
                                G4double nodeZ, 
                                G4int level )
{
    DicomOctreeNode* current = mRoot;
    G4double currentX = 0;
    G4double currentY = 0;
    G4double currentZ = 0;
    // Make children
    for ( G4int i = 0; i < level; i++ )
    {
        G4double childLevelResolution = ( 1 << (i+1) );
        G4double childSize = mSize / childLevelResolution;

        G4int dirX = G4int( nodeX >= currentX + childSize );
        G4int dirY = G4int( nodeY >= currentY + childSize );
        G4int dirZ = G4int( nodeZ >= currentZ + childSize );
        G4int direction = dirX + ( dirY << 1 ) + ( dirZ << 2 );

        if ( (*current)[direction] == 0 ) return 0;

        current = (*current)[direction];
        currentX += dirX*childSize;
        currentY += dirY*childSize;
        currentZ += dirZ*childSize;
    }
    return current;
}

void DicomOctree::DeleteTree()
{
  delete mRoot;
  mRoot = 0;
}

void DicomOctree::CountRecursive( DicomOctreeNode* pNode, 
                             G4int rMiddle, 
                             G4int rTerminal )
{
  if ( pNode->Type() == MIDDLE_NODE )rMiddle++;
  else rTerminal++;
      
  for ( G4int i = 0; i < 8; i++ )
    {
      if ( (*pNode)[i] != 0 )
	CountRecursive( (*pNode)[i], rMiddle, rTerminal );       
    }
}

G4int DicomOctree::CountMemory( G4int rMiddle, G4int rTerminal )
{
  CountRecursive( mRoot, rMiddle, rTerminal );
  G4int total = rMiddle*sizeof(DicomOctreeMiddleNode) + 
                rTerminal*sizeof(DicomOctreeTerminalNode);
  return total;
}


