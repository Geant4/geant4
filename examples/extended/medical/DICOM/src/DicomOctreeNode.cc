
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


