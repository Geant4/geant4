#ifndef DicomOctreeTerminalNode_h
#define DicomOctreeTerminalNode_h

class Octree;
class TerminalNode : public OctreeNode
{
public:

  TerminalNode(OctreeNode* pParent );
  ~TerminalNode();

public:
  OctreeNode*& operator []( G4int index );
  G4int MemSize();

  OctreeNodeType Type()                      {return TERMINAL_NODE;}
  G4int FindChild( const OctreeNode* pNode ) {return -1;}

private:
  static OctreeNode* mNull;
};
#endif
