#ifndef OCTREE_NODE_H
#define OCTREE_NODE_H
#include "globals.hh"

//---------------------------------------------------------------------------
class Octree;
//---------------------------------------------------------------------------
enum OctreeNodeType { MIDDLE_NODE, TERMINAL_NODE };
//---------------------------------------------------------------------------
/// \brief OctreeNode represents a single node in the octree.
///
/// Longer description.
class OctreeNode
{
    friend class Octree;

public:
    OctreeNode();
    OctreeNode( OctreeNode* pParent );
    virtual ~OctreeNode();

    float& Density();

    virtual OctreeNode*& operator []( G4int index ) = 0;

    virtual OctreeNodeType Type() = 0;

    const OctreeNode* Parent();

    virtual G4int FindChild( const OctreeNode* pNode ) = 0;

    static G4int InstanceCounter();

    virtual G4int MemSize() = 0;

private:
    static G4int mInstanceCounter;
    OctreeNode* mParent;
    float mDensity;
};

//---------------------------------------------------------------------------
class MiddleNode : public OctreeNode
{
    friend class Octree;

public:
    OctreeNode*  mChildren[8];

    void ResetFamily();
    MiddleNode();
    MiddleNode( OctreeNode* pParent );
    ~MiddleNode();

    OctreeNode*& operator []( G4int index );
    OctreeNodeType Type();
    G4int FindChild( const OctreeNode* pNode );
    G4int MemSize();
};
//---------------------------------------------------------------------------
class TerminalNode : public OctreeNode
{
    friend class Octree;

public:
    OctreeNode*& operator []( G4int index );

    OctreeNodeType Type();

    TerminalNode(OctreeNode* pParent );

    ~TerminalNode();

    G4int FindChild( const OctreeNode* pNode );

    G4int MemSize();

private:
    static OctreeNode* mNull;
};
//---------------------------------------------------------------------------

#endif // __OCTREE_NODE_H__
