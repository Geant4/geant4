// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DrawVoxels.cc,v 1.9 1999-12-15 14:49:52 gunter Exp $
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
// 29/07/1999 first comitted version L.G.


//Including the necessary data
#include "G4DrawVoxels.hh"
#include "G4VVisManager.hh"
/****************************/

//Methods that allow changing colors of the drawing
void G4DrawVoxels::SetVoxelsVisAttributes(G4VisAttributes& VA_voxelX,G4VisAttributes& VA_voxelY,G4VisAttributes& VA_voxelZ){
   fVoxelsVisAttributes[0]=VA_voxelX;	//operateur de copie ...
   fVoxelsVisAttributes[1]=VA_voxelY;
   fVoxelsVisAttributes[2]=VA_voxelZ;
}
void G4DrawVoxels::SetBoundingBoxVisAttributes(G4VisAttributes& VA_boundingbox){
   fBoundingBoxVisAttributes=VA_boundingbox;	//operateur de copie ...
}


//#********************************************************************************************************************#
//#********************************************************************************************************************#
void  G4DrawVoxels::ComputeVoxelPolyhedra(const G4LogicalVolume* lv,
					 const G4SmartVoxelHeader* header,
					 G4VoxelLimits& limit,
					 G4PlacedPolyhedronList* ppl) const
{
//######################################################################################################################
// Let's draw the selected voxelisation now !
 
   #ifdef G4DrawVoxelsDebug
 	G4cout << "**** !!! G4DrawVoxels::Draw_Voxels : a query for drawing voxelisation";
 	G4cout << "	Logical Volume:" << lv->GetName() <<G4endl;
   #endif
 

   G4VSolid* solid=lv->GetSolid();
  
   G4double dx=kInfinity,dy=kInfinity,dz=kInfinity;
   G4double xmax=0,xmin=0,ymax=0,ymin=0,zmax=0,zmin=0;
   G4Transform3D translation_for_voxel_plane;	// According to the solid's own frame
   
   if (lv->GetNoDaughters()<=0) {
     #ifdef G4DrawVoxelsDebug
        G4cout << "**** !!! G4DrawVoxels::Draw_Voxels : a query for drawing voxelisation" << G4endl;
        G4cout << "	of a logical volume that has no (or <=0) daughter has been detected !!! ****" <<G4endl;
     #endif
     return;
   }
   
      	
   //Let's get the data for the voxelisation
   solid->CalculateExtent(kXAxis,limit,G4AffineTransform(),xmin,xmax); //G4AffineTransform() is identity
   solid->CalculateExtent(kYAxis,limit,G4AffineTransform(),ymin,ymax); //extents according to the axis of the local frame
   solid->CalculateExtent(kZAxis,limit,G4AffineTransform(),zmin,zmax);
   dx=xmax-xmin;
   dy=ymax-ymin;
   dz=zmax-zmin;

   //Preparing the colored bounding polyhedronBox for the pVolume
   G4PolyhedronBox bounding_polyhedronBox(dx*0.5,dy*0.5,dz*0.5);
   bounding_polyhedronBox.SetVisAttributes(&fBoundingBoxVisAttributes);
   G4ThreeVector t_centerofBoundingBox((xmin+xmax)*0.5,(ymin+ymax)*0.5,(zmin+zmax)*0.5);
   
   //ppl->resize(ppl->entries()+1);	//manual resize to avoid Rogue grabbing RW_DEFAULT
   ppl->insert(G4PlacedPolyhedron(bounding_polyhedronBox,G4Translate3D(t_centerofBoundingBox)));
   
   G4ThreeVector t_FirstCenterofVoxelPlane;
   const G4VisAttributes* voxelsVisAttributes;

   G4ThreeVector unit_translation_vector;
   G4ThreeVector current_translation_vector;
   
   switch(header->GetAxis())
	{
		case kXAxis:
		    dx=voxel_width;
		    unit_translation_vector=G4ThreeVector(1,0,0);
		    t_FirstCenterofVoxelPlane=G4ThreeVector(xmin,(ymin+ymax)*0.5,(zmin+zmax)*0.5);
		    voxelsVisAttributes=&fVoxelsVisAttributes[0];
		    break;
		case kYAxis:
		    dy=voxel_width;
		    t_FirstCenterofVoxelPlane=G4ThreeVector((xmin+xmax)*0.5,ymin,(zmin+zmax)*0.5);
		    unit_translation_vector=G4ThreeVector(0,1,0);
		    voxelsVisAttributes=&fVoxelsVisAttributes[1];
		    break;
		case kZAxis:
		    dz=voxel_width;
		    t_FirstCenterofVoxelPlane=G4ThreeVector((xmin+xmax)*0.5,(ymin+ymax)*0.5,zmin);
		    unit_translation_vector=G4ThreeVector(0,0,1);
		    voxelsVisAttributes=&fVoxelsVisAttributes[2];
		    break;
		default:
		//erreur interne
		/************************/
		G4cout << "PANIC: in Draw_Voxels:DrawVoxels(...) header is decayed." <<G4endl;
		G4cout << "header->GetAxis returns an invalid axis." <<G4endl; 
		/************************/
		break;
	};
     
   G4PolyhedronBox voxel_plane(dx*0.5,dy*0.5,dz*0.5);
   voxel_plane.SetVisAttributes(voxelsVisAttributes);
   
   G4SmartVoxelProxy* slice=header->GetSlice(0);
   G4int slice_no=0,no_slices=header->GetNoSlices();
   G4double beginning=header->GetMinExtent(),step=(header->GetMaxExtent()-beginning)/no_slices;
		
   #ifdef G4DrawVoxelsDebug
      G4cout << "	Axis of the voxelisation:" << header->GetAxis();
      G4cout << "	No slices of the voxelisation:" << header->GetNoSlices();
      G4cout << "	step:" << step <<G4endl;	
   #endif

   while (slice_no<no_slices){
		  
	#ifdef G4DrawVoxelsDebug
	   G4cout << "		slice_no:" << slice_no;
	   G4cout << "		header/node:" << slice->IsHeader() <<G4endl;					
  	#endif
  		  
	if (slice->IsHeader()){
	   G4VoxelLimits newlimit(limit);
	   newlimit.AddLimit(header->GetAxis(),beginning+step*slice_no,
	   beginning+step*(slice->GetHeader()->GetMaxEquivalentSliceNo()+1));
	   ComputeVoxelPolyhedra(lv,slice->GetHeader(),newlimit,ppl);
	}
	current_translation_vector=unit_translation_vector;
	current_translation_vector*=step*slice_no;
		   
	//ppl->resize(ppl->entries()+1);	//manual resize to avoid Rogue grabbing RW_DEFAULT
	ppl->insert(G4PlacedPolyhedron(voxel_plane,G4Translate3D(current_translation_vector+t_FirstCenterofVoxelPlane)));
		   
	slice_no=(slice->IsHeader()?slice->GetHeader()->GetMaxEquivalentSliceNo()+1
		  			     :slice->GetNode()->GetMaxEquivalentSliceNo()+1);
	slice=header->GetSlice(slice_no);	
   }
}//end of ComputeVoxelPolyhedra...
//######################################################################################################################

//Private Constructor
G4DrawVoxels::G4DrawVoxels(){
	fVoxelsVisAttributes[0].SetColour(G4Colour(1.,0.,0.));
	fVoxelsVisAttributes[1].SetColour(G4Colour(0.,1.,0.));
	fVoxelsVisAttributes[2].SetColour(G4Colour(0.,0.,1.));
	fBoundingBoxVisAttributes.SetColour(G4Colour(.3,0.,.2));
}

G4AffineTransform G4DrawVoxels::GetAbsoluteTransformation(const G4VPhysicalVolume* pv) const{
  //Initialisation of the transformation to be computed
  G4AffineTransform transf; //default constructor ie Id

  const G4VPhysicalVolume* current= pv;
  const G4VPhysicalVolume* mother = pv->GetMother();

  //Now we loop up to the world physical volume
  while (mother!=NULL){
    //the following has to be verified since it could be the inverse transformations to be used: 
    //the explanation is not clear in G4VPhysicalVolume.hh.
    transf=G4AffineTransform(current->GetObjectRotation(),current->GetObjectTranslation())*transf;
    current=mother;
    mother=current->GetMother();
  }
  //current is the world volume
  return transf;
}



G4PlacedPolyhedronList* G4DrawVoxels::CreatePlacedPolyhedra(const G4LogicalVolume* lv) const{
   G4PlacedPolyhedronList* pplist=new G4PlacedPolyhedronList;
   G4VoxelLimits limits;  // Working object for recursive call.
   ComputeVoxelPolyhedra(lv,lv->GetVoxelHeader(),limits,pplist);
   return pplist; //it s up to the calling program to destroy it then!
}


void G4DrawVoxels::DrawVoxels(const G4LogicalVolume* lv) const{
   
   G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

   if (lv->GetNoDaughters()<=0) {
   #ifdef G4DrawVoxelsDebug
        G4cout << "**** !!! G4DrawVoxels::Draw_Voxels : a query for drawing voxelisation" << G4endl;
        G4cout << "	of a logical volume that has no (or <=0) daughter has been detected !!! ****" <<G4endl;
   #endif
   return;
   }
   //Computing the transformation according to the world volume 
   //(the drawing is directly in the world volume while the axis are relative to the mother volume of lv's daughter.)
   G4AffineTransform transf=GetAbsoluteTransformation(lv->GetDaughter(0)->GetMother());
	//in the hope that all daughters have same mother
	//and that the 0th is indeed the first to be filled.
   //G4Transform3D transf3D=GetAbsolute3DTransformation(lv->GetDaughter(0)->GetMother());
   G4Transform3D transf3D(transf.NetRotation(),transf.NetTranslation());

   G4PlacedPolyhedronList* pplist=CreatePlacedPolyhedra(lv);
   if(pVVisManager) {
	//Drawing the bounding and voxel polyhedra for the pVolume
	for (G4int i=0;i<pplist->entries();i++)
	   pVVisManager->Draw((*pplist)(i).GetPolyhedron(),(*pplist)(i).GetTransform()*transf3D);	
   }
   else G4cout << "@@@@ void G4DrawVoxels::DrawVoxels(const G4LogicalVolume*) pVVisManager is null! @@@@" <<G4endl; 	
   delete pplist;
}












