//version: Tue Jul 13 09:16:46 MET DST 1999

#ifndef G4DrawVoxels_HH
#define G4DrawVoxels_HH

//***********what I need to use (include and forward declarations) FOR DRAWING VOXELS****************
#include "G4VVisManager.hh" //./../../../intercoms/include/
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4Polyhedron.hh"
#include <rw/tpordvec.h>

#include "G4SmartVoxelHeader.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4Vector3D.hh"
#include "G4LogicalVolume.hh"
//#include "G4VPhysicalVolume.hh" already included in G4LogicalVolume.hh
#include "G4VSolid.hh"

#define voxel_width 0;

#define G4DrawVoxelsDebug

typedef RWTPtrOrderedVector<const G4LogicalVolume> G4LogicalVolumeList;
//***************************************************************************************************


//This class is a singleton. You can instantiate it only once.
//It stores a list of logical volumes whose voxelisation will have to be drawn.

class G4DrawVoxels{
public:
  //In order to get the adress of the G4DrawVoxels's unique instance
  static G4DrawVoxels* G4DrawVoxels::GetInstance();

  //Functions that allow basic operations	#****************************#
  //on the list of logical volumes to be drawn	#****************************#
  
  //Returns the number of logical volumes stored here.
  G4int G4DrawVoxels::GetNoLogicalVolumes();

  //Returns a pointer to the ith logical volume stored here.
  //No check is provided in case of out of range
  const G4LogicalVolume* G4DrawVoxels::GetLogicalVolume(const G4int);

  //Stores lv as the last logical volume in the list
  //No check is provided for the case lv is already in the list.
  void G4DrawVoxels::AddLogicalVolume(const G4LogicalVolume* lv);

  //Removes lv from the list of logical volumes.
  //Does nothing if lv is not a member of th list
  void G4DrawVoxels::RemoveLogicalVolume(const G4LogicalVolume* lv);

  //Returns true if lv is a member of the list of logical volumes,
  //Returns false otherwise.
  G4bool G4DrawVoxels::IsContained(const G4LogicalVolume* lv);

  //if lv is a member of the listlogical volumes, draws lv's voxels using the specified colors
  //Otherwise does nothing. 
  void G4DrawVoxels::DrawVoxels(const G4LogicalVolume* lv) {
        G4VoxelLimits limits;  // Working object for recursive call.
	DrawVoxels(lv,lv->GetVoxelHeader(),limits);}
  //const G4Colour& boundingbox_color=G4Colour(0.,1.,0.),
  //const G4Colour& slices_color=G4Colour(1.,0.,0.)

  void G4DrawVoxels::SetVoxelColours(G4Colour&,G4Colour&,G4Colour&);
  void G4DrawVoxels::SetBoundingBoxColour(G4Colour&);

private:
  //Member data
  static G4DrawVoxels* fgInstance;
  G4LogicalVolumeList fLogicalVolumes;
  G4Colour fvoxelcolours[3];
  G4Colour fboundingboxcolour;
  
  //constructor. It initialises the members data
  G4DrawVoxels::G4DrawVoxels();
  //To compute the transformation of a given G4PhysicalVolume according to the
  //the world volume.
  G4AffineTransform G4DrawVoxels::GetAbsoluteTransformation(const G4VPhysicalVolume* pv);
  /*Has become redundant
  G4Transform3D G4DrawVoxels::GetAbsolute3DTransformation(const G4VPhysicalVolume* pv);
  */
  void G4DrawVoxels::DrawVoxels(const G4LogicalVolume*,const G4SmartVoxelHeader*,G4VoxelLimits&);
};

#endif
