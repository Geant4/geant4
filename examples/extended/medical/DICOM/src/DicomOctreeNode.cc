
#include <stddef.h>
#include "DicomOctreeNode.hh"
using namespace std;

int OctreeNode::mInstanceCounter = 0;

OctreeNode::OctreeNode()
{
  mParent = NULL;
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

float& OctreeNode::Density()
{
  return mDensity;
}

OctreeNode* OctreeNode::Parent()
{
  return mParent;
}

int OctreeNode::InstanceCounter()
{
  return mInstanceCounter;
}

OctreeNode* TerminalNode::mNull = NULL;

OctreeNode*& TerminalNode::operator []( int )
{
  return mNull;
}

OctreeNodeType TerminalNode::Type()
{
  return TERMINAL_NODE;
}

int TerminalNode::FindChild( OctreeNode* pNode )
{
  return -1;
}

TerminalNode::~TerminalNode()
{

}

TerminalNode::TerminalNode( OctreeNode* pParent ) : OctreeNode( pParent )
{

}

int TerminalNode::MemSize()
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
OctreeNode*& MiddleNode::operator []( int index )
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
int MiddleNode::FindChild( OctreeNode* pNode )
{
    for ( unsigned int i = 0; i < 8; i++ )
    {
        if ( mChildren[i] == pNode ) return i;
    }
    return -1;
}
//---------------------------------------------------------------------------
int MiddleNode::MemSize()
{
    return sizeof(MiddleNode);
}
//---------------------------------------------------------------------------

