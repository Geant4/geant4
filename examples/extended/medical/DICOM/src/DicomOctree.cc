
//#include <stddef.h>
#include "DicomOctree.hh"
#include "DicomOctreeNode.hh"
#include "DicomOctreeMiddleNode.hh"
#include "DicomOctreeTerminalNode.hh"

DicomOctree::DicomOctree( G4int noLevels, G4double size )
{
  mNoLevels = noLevels;
  mSize = size;
  mRoot = new MiddleNode(0);
}

DicomOctree::~DicomOctree()
{
  delete mRoot;
  mRoot = 0; 
}

OctreeNode* DicomOctree::CreateNode( G4double nodeX, G4double nodeY, G4double nodeZ, G4int level )
{
  OctreeNode* current = mRoot;
  G4double currentX = 0;
  G4double currentY = 0;
  G4double currentZ = 0;
  
  for ( G4int i = 0; i < 8; i++ ) // Make children
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
	      (*current)[direction] = new MiddleNode( current );
            } else
	      {
                (*current)[direction] = new TerminalNode( current );
	      }
        }

      current = (*current)[direction];
      currentX += dirX*childSize;
      currentY += dirY*childSize;
      currentZ += dirZ*childSize;
    }
  return current;
}

OctreeNode* DicomOctree::operator()( G4double nodeX, 
                                G4double nodeY, 
                                G4double nodeZ, 
                                G4int level )
{
    OctreeNode* current = mRoot;
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

void DicomOctree::CountRecursive( OctreeNode* pNode, 
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
  G4int total = rMiddle*sizeof(MiddleNode) + rTerminal*sizeof(TerminalNode);
  return total;
}


