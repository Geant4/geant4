
#include <stddef.h>
#include "DicomOctreeNode.hh"
using namespace std;

G4int OctreeNode::mInstanceCounter = 0;

OctreeNode::OctreeNode()
{
  mParent = 0;
  mInstanceCounter++;
}

OctreeNode::~OctreeNode()
{
  mInstanceCounter--;
}

OctreeNode::OctreeNode( OctreeNode* pParent )
{
  mParent = pParent;
  mInstanceCounter++;
}

OctreeNode* TerminalNode::mNull = NULL;

OctreeNode*& TerminalNode::operator []( G4int )
{
  return mNull;
}

OctreeNodeType TerminalNode::Type()
{
  return TERMINAL_NODE;
}

G4int TerminalNode::FindChild( const OctreeNode*  )
{
  return -1;
}

TerminalNode::~TerminalNode()
{

}

TerminalNode::TerminalNode( OctreeNode* pParent) : OctreeNode( pParent )
{

}

G4int TerminalNode::MemSize()
{
  return sizeof(TerminalNode);
}
























//---------------------------------------------------------------------------
void MiddleNode::ResetFamily()
{
    for ( unsigned int i = 0; i < 8; i++ ) mChildren[i] = NULL;
}
//---------------------------------------------------------------------------
MiddleNode::MiddleNode()
{
    ResetFamily();
}
//---------------------------------------------------------------------------
MiddleNode::~MiddleNode()
{
    for ( unsigned int i = 0; i < 8; i++ )
    {
        if ( mChildren[i] != NULL ) delete mChildren[i];
    }
}
//---------------------------------------------------------------------------
MiddleNode::MiddleNode( OctreeNode* pParent ) : OctreeNode( pParent )
{
    ResetFamily();
}
//---------------------------------------------------------------------------
OctreeNode*& MiddleNode::operator []( G4int index )
{
//    if ( index > 7 ) throw runtime_error("Index out of bounds.");
    return mChildren[index];
}
//---------------------------------------------------------------------------
OctreeNodeType MiddleNode::Type()
{
    return MIDDLE_NODE;
}
//---------------------------------------------------------------------------
G4int MiddleNode::FindChild( const OctreeNode* pNode )
{
    for ( unsigned int i = 0; i < 8; i++ )
    {
        if ( mChildren[i] == pNode ) return i;
    }
    return -1;
}
//---------------------------------------------------------------------------
G4int MiddleNode::MemSize()
{
    return sizeof(MiddleNode);
}
//---------------------------------------------------------------------------

