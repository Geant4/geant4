
#include "DicomOctreeNode.hh"
#include "DicomOctreeMiddleNode.hh"

DicomOctreeMiddleNode::DicomOctreeMiddleNode()
{
    ResetFamily();
}

DicomOctreeMiddleNode::~DicomOctreeMiddleNode()
{
  for ( G4int i = 0; i < 8; i++ )
    {
      if ( mChildren[i] != 0 ) delete mChildren[i];
    }
}
void DicomOctreeMiddleNode::ResetFamily()
{
  // ---- MGP ---- Remove explicit numbers in code
  for ( G4int i = 0; i < 8; i++ ) mChildren[i] = 0;
}

DicomOctreeMiddleNode::DicomOctreeMiddleNode( DicomOctreeNode* pParent ) : DicomOctreeNode( pParent )
{
  ResetFamily();
}

G4int DicomOctreeMiddleNode::FindChild( const DicomOctreeNode* pNode )
{
  for ( G4int i = 0; i < 8; i++ )
    {
      if ( mChildren[i] == pNode ) return i;
    }
  return -1;
}

G4int DicomOctreeMiddleNode::MemSize()
{
  return sizeof(DicomOctreeMiddleNode);
}


