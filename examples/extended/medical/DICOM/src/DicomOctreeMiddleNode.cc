
#include "DicomOctreeNode.hh"
#include "DicomOctreeMiddleNode.hh"

MiddleNode::MiddleNode()
{
    ResetFamily();
}

MiddleNode::~MiddleNode()
{
  for ( G4int i = 0; i < 8; i++ )
    {
      if ( mChildren[i] != 0 ) delete mChildren[i];
    }
}
void MiddleNode::ResetFamily()
{
  // ---- MGP ---- Remove explicit numbers in code
  for ( G4int i = 0; i < 8; i++ ) mChildren[i] = 0;
}

MiddleNode::MiddleNode( OctreeNode* pParent ) : OctreeNode( pParent )
{
  ResetFamily();
}

G4int MiddleNode::FindChild( const OctreeNode* pNode )
{
  for ( G4int i = 0; i < 8; i++ )
    {
      if ( mChildren[i] == pNode ) return i;
    }
  return -1;
}

G4int MiddleNode::MemSize()
{
  return sizeof(MiddleNode);
}


