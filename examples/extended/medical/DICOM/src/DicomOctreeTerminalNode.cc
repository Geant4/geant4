#include "DicomOctreeNode.hh"
#include "DicomOctreeTerminalNode.hh"


OctreeNode* TerminalNode::mNull = 0;

TerminalNode::TerminalNode( const OctreeNode* pParent) : OctreeNode( pParent )
{

}
TerminalNode::~TerminalNode()
{

}

OctreeNode*& TerminalNode::operator []( G4int )
{
  return mNull;
}

G4int TerminalNode::MemSize()
{
  return sizeof(TerminalNode);
}


