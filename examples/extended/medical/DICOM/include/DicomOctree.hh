#ifndef OCTREE_HH
#define OCTREE_HH

#include "globals.hh"
#include "DicomOctreeNode.hh"

//---------------------------------------------------------------------------
///  Octree encapsulates a volumetric compressed representation of the
/// vector field used for the reconstruction of a 3D model from multiple
/// range images or curves.
/// 
/// It is assumed that the space occupied by octree is a cube whose
/// bottom-left-front corner is at origin (0,0,0) while the opposite corner
/// (top,right,back) is at ( mSize,mSize,mSize). Resolution of an octree
/// denotes the number of voxels along each axis such
/// that the resolution = 2^mNoLevels.


class Octree
{

public:
  
  Octree( G4int noLevels, G4double size );
  ~Octree();
  void DeleteTree();

  OctreeNode* CreateNode( G4double i, G4double j, G4double k, G4int level );
  OctreeNode* operator()( G4double nodeX, G4double nodeY, 
                          G4double nodeZ, G4int level );
    
  OctreeNode* Root(){return mRoot;} 
  G4double Size() {return mSize;}
  G4int CountMemory( G4int rMiddle, G4int rTerminal );
  
  G4int NoLevels(){return mNoLevels;}
  G4int Resolution(){return ( 1 << mNoLevels );}
  
private:


  void CountRecursive(OctreeNode* pNode, G4int rMiddle, G4int rTerminal );

  // Root node of the tree
  MiddleNode* mRoot;
  
  // In an octree, size denotes physical size of the cube,
  // ie length of its sides (which are assumed equal).
  G4double mSize;
  
  // Number of levels denotes the maximal number of
  // nodes in a single branch, starting from the root node.
  G4int mNoLevels;

};
#endif 

