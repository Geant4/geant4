//version: Thu Jul 22 09:52:44 CEST 1999

#ifndef G4DrawVoxels_HH
#define G4DrawVoxels_HH

//***********what I need to use (include and forward declarations) FOR DRAWING VOXELS****************
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4PlacedPolyhedron.hh" //#include "G4Polyhedron.hh" included

//#include <rw/tpordvec.h>
//#include <rw/tvordvec.h>

#include "G4SmartVoxelHeader.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4Vector3D.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"

#define voxel_width 0;

#define G4DrawVoxelsDebug

//***************************************************************************************************
class G4DrawVoxels{
  public: 
    //constructor. It initialises the members data to default colors
    G4DrawVoxels();
    ~G4DrawVoxels(){};
    //Copy constructor Assignment operator not supported (array fvoxelcolours ...)
    
    void DrawVoxels(const G4LogicalVolume* lv);
    G4PlacedPolyhedronList* CreatePlacedPolyhedra(const G4LogicalVolume*);

    void SetVoxelColours(G4Colour&,G4Colour&,G4Colour&);
    void SetBoundingBoxColour(G4Colour&);

  private:
    //Member data
    G4Colour fvoxelcolours[3];
    G4Colour fboundingboxcolour;
    
    void ComputeVoxelPolyhedra(const G4LogicalVolume*,const G4SmartVoxelHeader*,G4VoxelLimits&,G4PlacedPolyhedronList*);
    
    G4AffineTransform GetAbsoluteTransformation(const G4VPhysicalVolume*);
    
    //Copy constructor Assignment operator not supported (array fvoxelcolours ...)
    G4DrawVoxels(const G4DrawVoxels&);	
    G4DrawVoxels operator=(const G4DrawVoxels&);	
};

#endif
