//version: Tue Jul 13 09:16:46 MET DST 1999

//Including the necessary data
#include "G4DrawVoxels.hh"
//****************************

// Static class variable: ptr to single instance of class
G4DrawVoxels* G4DrawVoxels::fgInstance = 0;

//To get a reference to the singleton 	        #****************************#
G4DrawVoxels* G4DrawVoxels::GetInstance(){
    static G4DrawVoxels draw_voxels;
    if (!fgInstance) fgInstance = &draw_voxels;
    return fgInstance;    
}

//Functions that allow basic operations		#****************************#
//on the list of logical volumes to be drawn	#****************************#
G4int G4DrawVoxels::GetNoLogicalVolumes(){
	return fLogicalVolumes.entries();}

const G4LogicalVolume* G4DrawVoxels::GetLogicalVolume(const G4int i){
	return fLogicalVolumes(i);}

void G4DrawVoxels::AddLogicalVolume(const G4LogicalVolume* lv){
	//manual resize to avoid Rogue grabbing RW_DEFAULT
	fLogicalVolumes.resize(fLogicalVolumes.entries()+1);
	fLogicalVolumes.insert(lv);
}
void G4DrawVoxels::RemoveLogicalVolume(const G4LogicalVolume* lv){
	fLogicalVolumes.remove(lv);
	fLogicalVolumes.resize(fLogicalVolumes.entries());
}
G4bool G4DrawVoxels::IsContained(const G4LogicalVolume* lv){
	return fLogicalVolumes.contains(lv);}


//Methods that allow changing colors of the drawing
void G4DrawVoxels::SetVoxelColours(G4Colour& voxelX,G4Colour& voxelY,G4Colour& voxelZ){
  fvoxelcolours[0]=voxelX;
  fvoxelcolours[1]=voxelY;
  fvoxelcolours[2]=voxelZ;
}
void G4DrawVoxels::SetBoundingBoxColour(G4Colour& boundingboxcolour){
  fboundingboxcolour=boundingboxcolour;
}



//#********************************************************************************************************************#
//#********************************************************************************************************************#
void G4DrawVoxels::DrawVoxels(const G4LogicalVolume* lv,const G4SmartVoxelHeader* header,G4VoxelLimits& limit)
{
//######################################################################################################################
// Let's draw the selected voxelisation now !
//CAUTION: it is assumed that the voxelisation is done for the WORLD VOLUME => After, the transformation/World will be needed
// CAUTION: the drawing is not yet limited to the mother volume
//	CAUTION: G4Translate3D must change OK but also the new center for the voxel_plane.
//Notice that the voxel_width should perhaps( ??) be parameterisable
//Notice that it would be better not to draw equivalent slices.
//Notice that nothing is coded to check if step and voxels' thin voxels_width allow intercrosses of slices.
 
 #ifdef G4DrawVoxelsDebug
 	G4cout << "**** !!! G4DrawVoxels::Draw_Voxels : a query for drawing voxelisation";
 	G4cout << "	Logical Volume:" << lv->GetName() <<endl;
 #endif
 
 if (fLogicalVolumes.contains(lv)){
   G4VSolid* solid=lv->GetSolid();
   //G4SmartVoxelHeader* header=lv->GetVoxelHeader();
   /*G4Translate3D t_centeroflimits((limit.GetMaxXExtent()+limit.GetMinXExtent())*0.5,
   				  (limit.GetMaxYExtent()+limit.GetMinYExtent())*0.5,
   				  (limit.GetMaxZExtent()+limit.GetMinZExtent())*0.5);*/
   
   G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
   
   G4double dx=kInfinity,dy=kInfinity,dz=kInfinity;
   G4double xmax=0,xmin=0,ymax=0,ymin=0,zmax=0,zmin=0;
   G4Transform3D translation_for_voxel_plane;	// According to the world volume
   
   if (lv->GetNoDaughters()<=0) {
 #ifdef G4DrawVoxelsDebug
      G4cout << "**** !!! G4DrawVoxels::Draw_Voxels : a query for drawing voxelisation" << endl;
      G4cout << "	of a logical volume that has no (or <=0) daughter has been detected !!! ****" <<endl;
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
      	
   //Let's get the data for the voxelisation
   solid->CalculateExtent(kXAxis,limit,transf,xmin,xmax);
   solid->CalculateExtent(kYAxis,limit,transf,ymin,ymax);
   solid->CalculateExtent(kZAxis,limit,transf,zmin,zmax);
   dx=xmax-xmin;
   dy=ymax-ymin;
   dz=zmax-zmin;

   //Preparing the colored bounding polyhedronBox for the pVolume
   G4PolyhedronBox bounding_polyhedronBox(dx*0.5,dy*0.5,dz*0.5);
   bounding_polyhedronBox.SetVisAttributes(G4VisAttributes(fboundingboxcolour));
   G4Vector3D t_centerofBoundingBox((xmin+xmax)*0.5,(ymin+ymax)*0.5,(zmin+zmax)*0.5);
   G4Vector3D t_FirstCenterofVoxelPlane;
   G4Colour voxelcolour;

   G4Vector3D unit_translation_vector;
   G4Vector3D current_translation_vector;
   
   switch(header->GetAxis())
	{
		case kXAxis:
		    dx=voxel_width;
		    unit_translation_vector=transf3D*G4Vector3D(1,0,0);
		    t_FirstCenterofVoxelPlane=G4Vector3D(xmin,(ymin+ymax)*0.5,(zmin+zmax)*0.5);
		    voxelcolour=fvoxelcolours[0];
		    break;
		case kYAxis:
		    dy=voxel_width;
		    t_FirstCenterofVoxelPlane=G4Vector3D((xmin+xmax)*0.5,ymin,(zmin+zmax)*0.5);
		    unit_translation_vector=transf3D*G4Vector3D(0,1,0);
		    voxelcolour=fvoxelcolours[1];
		    break;
		case kZAxis:
		    dz=voxel_width;
		    t_FirstCenterofVoxelPlane=G4Vector3D((xmin+xmax)*0.5,(ymin+ymax)*0.5,zmin);
		    unit_translation_vector=transf3D*G4Vector3D(0,0,1);
		    voxelcolour=fvoxelcolours[2];
		    break;
		default:
		//erreur interne
		/************************
		PANIC: in Draw_Voxels: G4SmartVoxelHeader header is decayed.
		header->GetAxis returns an invalid axis. 
		************************/
		break;
	};
     
   G4PolyhedronBox voxel_plane(dx*0.5,dy*0.5,dz*0.5);
   voxel_plane.SetVisAttributes(G4VisAttributes(voxelcolour));
   
   if(pVVisManager) {
  #ifdef G4DrawVoxelsDebug
 	//G4cout << "  ** G4DrawVoxels::DrawVoxels(const G4LogicalVolume*,...) -DrawingVoxels- pVVisManager is OK! **" <<endl;
	G4cout << "	Axis of the voxelisation:" << header->GetAxis();
	G4cout << "	No slices of the voxelisation:" << header->GetNoSlices();
	G4cout << "	Colour of the slices:" << voxelcolour <<endl;	
  #endif
		//Drawing the green bounding polyhedronBox for the pVolume
		pVVisManager->Draw(bounding_polyhedronBox,G4Transform3D(transf3D.getRotation(),t_centerofBoundingBox));  //modified
		
		G4SmartVoxelProxy* slice=header->GetSlice(0);
		G4int slice_no=0,no_slices=header->GetNoSlices();
		G4double beginning=header->GetMinExtent(),/*end=header->GetMaxExtent()*/
            		 step=(header->GetMaxExtent()-beginning)/no_slices;
		current_translation_vector=unit_translation_vector;
		
		while (slice_no<no_slices){
		  #ifdef G4DrawVoxelsDebug
			G4cout << "		slice_no:" << slice_no;
			G4cout << "				header/node:" << slice->IsHeader() <<endl;					
  		  #endif
  		if (slice->IsHeader()){
		     G4VoxelLimits newlimit(limit);
		     newlimit.AddLimit(header->GetAxis(),beginning+step*slice_no,
		     	beginning+step*(slice->GetHeader()->GetMaxEquivalentSliceNo()+1));
		     DrawVoxels(lv,slice->GetHeader(),newlimit);
		     }
		  current_translation_vector*=step*slice_no;
		  pVVisManager->Draw(voxel_plane,G4Transform3D(transf3D.getRotation(),current_translation_vector+t_FirstCenterofVoxelPlane));
		  current_translation_vector=unit_translation_vector;
		  slice_no=(slice->IsHeader()?slice->GetHeader()->GetMaxEquivalentSliceNo()+1
		  			     :slice->GetNode()->GetMaxEquivalentSliceNo()+1);
		  slice=header->GetSlice(slice_no);	
		}
	}
   else G4cout << "@@@@ G4DrawVoxels::DrawVoxels(const G4LogicalVolume*,...) pVVisManager is null! @@@@"; 
 } //endif (fLogicalVolumes.contains(lv))
 else G4cout << "@@@@ G4DrawVoxels::DrawVoxels(const G4LogicalVolume*,...) lv is not inscribed I cannot draw it @@@@";
}//end of Draw_Voxels...
//######################################################################################################################

//Private Constructor
G4DrawVoxels::G4DrawVoxels(){
	fgInstance=0;
	fLogicalVolumes(0);
	fvoxelcolours[0]=G4Colour(1.,0.,0.);
	fvoxelcolours[1]=G4Colour(0.,1.,0.);
	fvoxelcolours[2]=G4Colour(0.,0.,1.);
	fboundingboxcolour=G4Colour(.3,0.,.3);
}


//## Computation of the transformation of a given G4VPhysicalVolume	************************************************
//## according to the world volume					************************************************	
G4AffineTransform G4DrawVoxels::GetAbsoluteTransformation(const G4VPhysicalVolume* pv)
{
  //Initialisation of the transformation to be computed
  G4AffineTransform transf; //default constructor ie Id

  const G4VPhysicalVolume* current(pv);
  const G4VPhysicalVolume* mother(pv->GetMother());

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

/*
//## Computation of the transformation of a given G4VPhysicalVolume	************************************************
//## according to the world volume					************************************************	
G4Transform3D G4DrawVoxels::GetAbsolute3DTransformation(const G4VPhysicalVolume* pv)
{
  //Initialisation of the transformation to be computed
  G4Transform3D transf3D; //default constructor ie Id

  const G4VPhysicalVolume* current(pv);
  const G4VPhysicalVolume* mother(pv->GetMother());

  //Now we loop up to the world physical volume
  while (mother!=NULL){
    //the following has to be verified since it could be the inverse transformations to be used: 
    //the explanation is not clear in G4VPhysicalVolume.hh.
    transf3D=G4Transform3D(*(current->GetObjectRotation()),current->GetObjectTranslation())*transf3D;
    current=mother;
    mother=current->GetMother();
  }
  //current is the world volume
  return transf3D;
}
*/















