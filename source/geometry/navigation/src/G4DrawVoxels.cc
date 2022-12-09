//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// class G4DrawVoxels implementation
//
// Define G4DrawVoxelsDebug for debugging information on G4cout
//
// 29/07/1999 first comitted version L.G.
// --------------------------------------------------------------------

#include "G4DrawVoxels.hh"
#include "G4AffineTransform.hh"
#include "G4SmartVoxelHeader.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4TransportationManager.hh"
#include "G4TouchableHistoryHandle.hh"

#define voxel_width 0

// Private Constructor
//
G4DrawVoxels::G4DrawVoxels()
{
  fVoxelsVisAttributes[0].SetColour(G4Colour(1.,0.,0.));
  fVoxelsVisAttributes[1].SetColour(G4Colour(0.,1.,0.));
  fVoxelsVisAttributes[2].SetColour(G4Colour(0.,0.,1.));
  fBoundingBoxVisAttributes.SetColour(G4Colour(.3,0.,.2));
}

// Destructor
//
G4DrawVoxels::~G4DrawVoxels()
{
}

// Methods that allow changing colors of the drawing
//
void G4DrawVoxels::SetVoxelsVisAttributes(G4VisAttributes& VA_voxelX,
                                          G4VisAttributes& VA_voxelY,
                                          G4VisAttributes& VA_voxelZ)
{
  fVoxelsVisAttributes[0] = VA_voxelX;
  fVoxelsVisAttributes[1] = VA_voxelY;
  fVoxelsVisAttributes[2] = VA_voxelZ;
}

void G4DrawVoxels::SetBoundingBoxVisAttributes(G4VisAttributes& VA_boundingbox)
{
  fBoundingBoxVisAttributes = VA_boundingbox;
}

// --------------------------------------------------------------------

void
G4DrawVoxels::ComputeVoxelPolyhedra(const G4LogicalVolume* lv,
                                    const G4SmartVoxelHeader* header,
                                          G4VoxelLimits& limit,
                                          G4PlacedPolyhedronList* ppl) const
{
  // Let's draw the selected voxelisation now !
 
   G4VSolid* solid = lv->GetSolid();
  
   G4double dx=kInfinity, dy=kInfinity, dz=kInfinity;
   G4double xmax=0, xmin=0, ymax=0, ymin=0, zmax=0, zmin=0;
   
   if (lv->GetNoDaughters()<=0)
   {
     return;
   }
   
   // Let's get the data for the voxelisation

   solid->CalculateExtent(kXAxis,limit,G4AffineTransform(),xmin,xmax);
     // G4AffineTransform() is identity
   solid->CalculateExtent(kYAxis,limit,G4AffineTransform(),ymin,ymax);
     // extents according to the axis of the local frame
   solid->CalculateExtent(kZAxis,limit,G4AffineTransform(),zmin,zmax);
   dx = xmax-xmin;
   dy = ymax-ymin;
   dz = zmax-zmin;

   // Preparing the colored bounding polyhedronBox for the pVolume
   //
   G4PolyhedronBox bounding_polyhedronBox(dx*0.5,dy*0.5,dz*0.5);
   bounding_polyhedronBox.SetVisAttributes(&fBoundingBoxVisAttributes);
   G4ThreeVector t_centerofBoundingBox((xmin+xmax)*0.5,
                                       (ymin+ymax)*0.5,
                                       (zmin+zmax)*0.5);
   
   ppl->push_back(G4PlacedPolyhedron(bounding_polyhedronBox,
                                     G4Translate3D(t_centerofBoundingBox)));
   
   G4ThreeVector t_FirstCenterofVoxelPlane;
   const G4VisAttributes* voxelsVisAttributes = nullptr;

   G4ThreeVector unit_translation_vector;
   G4ThreeVector current_translation_vector;
   
   switch(header->GetAxis())
   {
     case kXAxis:
       dx=voxel_width;
       unit_translation_vector=G4ThreeVector(1,0,0);
       t_FirstCenterofVoxelPlane=G4ThreeVector(xmin,(ymin+ymax)*0.5,
                                                    (zmin+zmax)*0.5);
       voxelsVisAttributes=&fVoxelsVisAttributes[0];
       break;
     case kYAxis:
       dy=voxel_width;
       t_FirstCenterofVoxelPlane=G4ThreeVector((xmin+xmax)*0.5,ymin,
                                               (zmin+zmax)*0.5);
       unit_translation_vector=G4ThreeVector(0,1,0);
       voxelsVisAttributes=&fVoxelsVisAttributes[1];
       break;
     case kZAxis:
       dz=voxel_width;
       t_FirstCenterofVoxelPlane=G4ThreeVector((xmin+xmax)*0.5,
                                               (ymin+ymax)*0.5,zmin);
       unit_translation_vector=G4ThreeVector(0,0,1);
       voxelsVisAttributes=&fVoxelsVisAttributes[2];
       break;
     default:
       break;
   };
     
   G4PolyhedronBox voxel_plane(dx*0.5,dy*0.5,dz*0.5);
   voxel_plane.SetVisAttributes(voxelsVisAttributes);
   
   G4SmartVoxelProxy* slice = header->GetSlice(0);
   std::size_t slice_no = 0, no_slices = header->GetNoSlices();
   G4double beginning = header->GetMinExtent(),
            step = (header->GetMaxExtent()-beginning)/no_slices;

   while (slice_no<no_slices)
   {    
     if (slice->IsHeader())
     {
       G4VoxelLimits newlimit(limit);
       newlimit.AddLimit(header->GetAxis(), beginning+step*slice_no,
         beginning+step*(slice->GetHeader()->GetMaxEquivalentSliceNo()+1));
       ComputeVoxelPolyhedra(lv,slice->GetHeader(), newlimit, ppl);
     }
     current_translation_vector = unit_translation_vector;
     current_translation_vector *= step*slice_no;
   
     ppl->push_back(G4PlacedPolyhedron(voxel_plane,
                    G4Translate3D(current_translation_vector
                                 + t_FirstCenterofVoxelPlane)));
     slice_no = (slice->IsHeader()
               ? slice->GetHeader()->GetMaxEquivalentSliceNo()+1
               : slice->GetNode()->GetMaxEquivalentSliceNo()+1);
     if (slice_no<no_slices) { slice=header->GetSlice(slice_no); }
   }
}

// --------------------------------------------------------------------

G4PlacedPolyhedronList*
G4DrawVoxels::CreatePlacedPolyhedra(const G4LogicalVolume* lv) const
{
  G4PlacedPolyhedronList* pplist = new G4PlacedPolyhedronList;
  G4VoxelLimits limits;  // Working object for recursive call.
  ComputeVoxelPolyhedra(lv,lv->GetVoxelHeader(),limits,pplist);
  return pplist; //it s up to the calling program to destroy it then!
}

// --------------------------------------------------------------------

void G4DrawVoxels::DrawVoxels(const G4LogicalVolume* lv) const
{   
   G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

   if (lv->GetNoDaughters()<=0)
   {
     return;
   }

   // Computing the transformation according to the world volume 
   // (the drawing is directly in the world volume while the axis
   // are relative to the mother volume of lv's daughter.)

   G4TouchableHistoryHandle aTouchable =
     G4TransportationManager::GetTransportationManager()->
     GetNavigatorForTracking()->CreateTouchableHistoryHandle();
   G4AffineTransform globTransform =
     aTouchable->GetHistory()->GetTopTransform().Inverse();
   G4Transform3D transf3D(globTransform.NetRotation(),
                          globTransform.NetTranslation());

   G4PlacedPolyhedronList* pplist = CreatePlacedPolyhedra(lv);
   if(pVVisManager != nullptr)
   {
     // Drawing the bounding and voxel polyhedra for the pVolume
     //
     for (size_t i=0; i<pplist->size(); ++i)
     {
       pVVisManager->Draw((*pplist)[i].GetPolyhedron(),
                          (*pplist)[i].GetTransform()*transf3D);
     }
   }
   else
   {
     G4Exception("G4DrawVoxels::DrawVoxels()",
                 "GeomNav1002", JustWarning,
                 "Pointer to visualization manager is null!");
   }
   delete pplist;
}
