
#include <stddef.h>
#include "DicomOctree.hh"

Octree::Octree( G4int noLevels, G4double size )
{
  mNoLevels = noLevels;
  mSize = size;
  mRoot = new MiddleNode(0);
}

Octree::~Octree()
{
  delete mRoot;
  mRoot = 0; 
}

OctreeNode* Octree::CreateNode( G4double nodeX, G4double nodeY, G4double nodeZ, G4int level )
{
    OctreeNode* current = mRoot;
    float curr_x = 0;
    float curr_y = 0;
    float curr_z = 0;
    // Make children
    for ( G4int i = 0; i < 8; i++ )
    {
        float child_level_resolution = (1 << (i+1) );
        float child_size = mSize / child_level_resolution;

        G4int dir_x = G4int( nodeX >= curr_x + child_size );
        G4int dir_y = G4int( nodeY >= curr_y + child_size );
        G4int dir_z = G4int( nodeZ >= curr_z + child_size );
        G4int direction = dir_x + ( dir_y << 1 ) + ( dir_z << 2 );

        if ( (*current)[direction] == 0 ) // NULL
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
        curr_x += dir_x*child_size;
        curr_y += dir_y*child_size;
        curr_z += dir_z*child_size;
    }
    return current;
}

OctreeNode* Octree::operator()( G4double nodeX, 
                                G4double nodeY, 
                                G4double nodeZ, 
                                G4int level )
{
    OctreeNode* current = mRoot;
    float curr_x = 0;
    float curr_y = 0;
    float curr_z = 0;
    // Make children
    for ( G4int i = 0; i < level; i++ )
    {
        float child_level_resolution = ( 1 << (i+1) );
        float child_size = mSize / child_level_resolution;

        G4int dir_x = G4int( nodeX >= curr_x + child_size );
        G4int dir_y = G4int( nodeY >= curr_y + child_size );
        G4int dir_z = G4int( nodeZ >= curr_z + child_size );
        G4int direction = dir_x + ( dir_y << 1 ) + ( dir_z << 2 );

        if ( (*current)[direction] == 0/*NULL*/ ) return 0/*NULL*/;

        current = (*current)[direction];
        curr_x += dir_x*child_size;
        curr_y += dir_y*child_size;
        curr_z += dir_z*child_size;
    }
    return current;
}



OctreeNode* Octree::Root()
{
    return mRoot;
}

G4int Octree::NoLevels()
{
    return mNoLevels;
}

G4int Octree::Resolution()
{
    return ( 1 << mNoLevels );
}

void Octree::DeleteTree()
{
    delete mRoot;
    mRoot = NULL;
}

void Octree::CountRecursive( OctreeNode* pNode, G4int rMiddle, G4int rTerminal )
{
    if ( pNode->Type() == MIDDLE_NODE )
    {
        rMiddle++;
    } else
    {
        rTerminal++;
    }
    for ( unsigned int i = 0; i < 8; i++ )
    {
        if ( (*pNode)[i] != NULL )
        {
            CountRecursive( (*pNode)[i], rMiddle, rTerminal );
        }
    }
}
//---------------------------------------------------------------------------
G4int Octree::CountMemory( G4int rMiddle, G4int rTerminal )
{
    CountRecursive( mRoot, rMiddle, rTerminal );
    G4int total = rMiddle*sizeof(MiddleNode) + rTerminal*sizeof(TerminalNode);
    return total;
}
//---------------------------------------------------------------------------

