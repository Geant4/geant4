// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DrawVoxels.hh,v 1.5 1999-08-03 09:51:57 graignac Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4DrawVoxels
//
// Implementation
//
// Define G4DrawVoxelsDebug for debugging information on G4cout
//
// History:
// 03/08/1999 The G4VisAttributes have been made member data for lifetime reasons / visualisation  L.G (see John Allison for further explanation) 
// 29.07.99 first comitted version L.G.


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

    void SetVoxelsVisAttributes(G4VisAttributes&,G4VisAttributes&,G4VisAttributes&);
    void SetBoundingBoxVisAttributes(G4VisAttributes&);

  private:
    //Member data
    G4VisAttributes fVoxelsVisAttributes[3];
    G4VisAttributes fBoundingBoxVisAttributes;
    
    void ComputeVoxelPolyhedra(const G4LogicalVolume*,const G4SmartVoxelHeader*,G4VoxelLimits&,G4PlacedPolyhedronList*);
    
    G4AffineTransform GetAbsoluteTransformation(const G4VPhysicalVolume*);
    
    //Copy constructor Assignment operator not supported (array fvoxelcolours ...)
    G4DrawVoxels(const G4DrawVoxels&);	
    G4DrawVoxels operator=(const G4DrawVoxels&);	
};

#endif
