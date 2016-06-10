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
// $Id: G4AdjointCrossSurfChecker.cc 66872 2013-01-15 01:25:57Z japost $
//
/////////////////////////////////////////////////////////////////////////////
//      Class Name:	G4AdjointCrossSurfChecker
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////

#include "G4AdjointCrossSurfChecker.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4PhysicalVolumeStore.hh" 
#include "G4VSolid.hh"
#include "G4AffineTransform.hh"

//////////////////////////////////////////////////////////////////////////////
// 
G4ThreadLocal G4AdjointCrossSurfChecker* G4AdjointCrossSurfChecker::instance = 0;

//////////////////////////////////////////////////////////////////////////////
//
G4AdjointCrossSurfChecker::G4AdjointCrossSurfChecker()
{;
}
///////////////////////////////////////////////////////////////////////////////
//
G4AdjointCrossSurfChecker::~G4AdjointCrossSurfChecker()
{
  delete instance;
}
////////////////////////////////////////////////////////////////////////////////
//
G4AdjointCrossSurfChecker* G4AdjointCrossSurfChecker::GetInstance()
{
  if (!instance) instance = new G4AdjointCrossSurfChecker();
  return instance;
}		  
/////////////////////////////////////////////////////////////////////////////////
//
G4bool G4AdjointCrossSurfChecker::CrossingASphere(const G4Step* aStep,G4double sphere_radius, G4ThreeVector sphere_center,G4ThreeVector& crossing_pos, G4double& cos_th , G4bool& GoingIn)
{
  G4ThreeVector pos1=  aStep->GetPreStepPoint()->GetPosition() - sphere_center;
  G4ThreeVector pos2=  aStep->GetPostStepPoint()->GetPosition() - sphere_center;
  G4double r1= pos1.mag();
  G4double r2= pos2.mag();
  G4bool did_cross =false; 
  
  if (r1<=sphere_radius && r2>sphere_radius){
 	did_cross=true;
	GoingIn=false;
  } 
  else if (r2<=sphere_radius && r1>sphere_radius){
  	did_cross=true;
	GoingIn=true;
  }

  if (did_cross) { 
  	
	G4ThreeVector dr=pos2-pos1;
	G4double r12 = r1*r1;
	G4double rdr = dr.mag();
	G4double a,b,c,d;
	a = rdr*rdr;
	b = 2.*pos1.dot(dr);
	c = r12-sphere_radius*sphere_radius;
	d=std::sqrt(b*b-4.*a*c);
	G4double l= (-b+d)/2./a;
	if (l > 1.) l=(-b-d)/2./a;
	crossing_pos=pos1+l*dr;
	cos_th = std::abs(dr.cosTheta(crossing_pos));
	
  
  }
  return did_cross;
}
/////////////////////////////////////////////////////////////////////////////////
//
G4bool G4AdjointCrossSurfChecker::GoingInOrOutOfaVolume(const G4Step* aStep,const G4String& volume_name, G4double& , G4bool& GoingIn) //from external surface
{
  G4bool step_at_boundary = (aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary);
  G4bool did_cross =false;
  if (step_at_boundary){
  	const G4VTouchable* postStepTouchable = aStep->GetPostStepPoint()->GetTouchable();
	const G4VTouchable* preStepTouchable = aStep->GetPreStepPoint()->GetTouchable();
	if (preStepTouchable && postStepTouchable && postStepTouchable->GetVolume() && preStepTouchable->GetVolume()){
		G4String post_vol_name = postStepTouchable->GetVolume()->GetName();
		G4String pre_vol_name = preStepTouchable->GetVolume()->GetName();
		
		if (post_vol_name == volume_name ){
			GoingIn=true;
			did_cross=true;
		}
		else if (pre_vol_name == volume_name){
			GoingIn=false;
			did_cross=true;
			
		}
		
	} 
  }
  return did_cross;
  //still need to compute the cosine of the direction
}
/////////////////////////////////////////////////////////////////////////////////
//
G4bool G4AdjointCrossSurfChecker::GoingInOrOutOfaVolumeByExtSurface(const G4Step* aStep,const G4String& volume_name, const G4String& mother_logical_vol_name, G4double& , G4bool& GoingIn) //from external surface
{
  G4bool step_at_boundary = (aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary);
  G4bool did_cross =false;
  if (step_at_boundary){
  	const G4VTouchable* postStepTouchable = aStep->GetPostStepPoint()->GetTouchable();
	const G4VTouchable* preStepTouchable = aStep->GetPreStepPoint()->GetTouchable();
	if (preStepTouchable && postStepTouchable && postStepTouchable->GetVolume() && preStepTouchable->GetVolume()){
		G4String post_vol_name = postStepTouchable->GetVolume()->GetName();
		G4String post_log_vol_name = postStepTouchable->GetVolume()->GetLogicalVolume()->GetName();
		G4String pre_vol_name = preStepTouchable->GetVolume()->GetName();
		G4String pre_log_vol_name = preStepTouchable->GetVolume()->GetLogicalVolume()->GetName();
		if (post_vol_name == volume_name && pre_log_vol_name ==  mother_logical_vol_name){
			GoingIn=true;
			did_cross=true;
		}
		else if (pre_vol_name == volume_name && post_log_vol_name ==  mother_logical_vol_name ){
			GoingIn=false;
			did_cross=true;
			
		}
		
	} 
  }
  return did_cross;
  //still need to compute the cosine of the direction
}
/////////////////////////////////////////////////////////////////////////////////
//
G4bool G4AdjointCrossSurfChecker::CrossingAGivenRegisteredSurface(const G4Step* aStep,const G4String& surface_name,G4ThreeVector& crossing_pos, G4double& cos_to_surface, G4bool& GoingIn)
{
  G4int ind = FindRegisteredSurface(surface_name);
  G4bool did_cross = false;
  if (ind >=0){
  	did_cross = CrossingAGivenRegisteredSurface(aStep, ind, crossing_pos,cos_to_surface, GoingIn);
  }
  return did_cross;
}
/////////////////////////////////////////////////////////////////////////////////
//
G4bool G4AdjointCrossSurfChecker::CrossingAGivenRegisteredSurface(const G4Step* aStep, int ind,G4ThreeVector& crossing_pos,   G4double& cos_to_surface, G4bool& GoingIn)
{
   G4String surf_type = ListOfSurfaceType[ind];
   G4double radius = ListOfSphereRadius[ind];
   G4ThreeVector  center = ListOfSphereCenter[ind];
   G4String vol1 = ListOfVol1Name[ind];
   G4String vol2 = ListOfVol2Name[ind];
	
   G4bool did_cross = false;	
   if (surf_type == "Sphere"){
	did_cross = CrossingASphere(aStep, radius, center,crossing_pos, cos_to_surface, GoingIn);
   }
   else if (surf_type == "ExternalSurfaceOfAVolume"){
  
	did_cross = GoingInOrOutOfaVolumeByExtSurface(aStep, vol1, vol2, cos_to_surface, GoingIn);
	crossing_pos= aStep->GetPostStepPoint()->GetPosition();
		
   }
   else if (surf_type == "BoundaryBetweenTwoVolumes"){
	did_cross = CrossingAnInterfaceBetweenTwoVolumes(aStep, vol1, vol2,crossing_pos, cos_to_surface, GoingIn);
   }
   return did_cross;
	
  
}
/////////////////////////////////////////////////////////////////////////////////
//
//
G4bool G4AdjointCrossSurfChecker::CrossingOneOfTheRegisteredSurface(const G4Step* aStep,G4String& surface_name,G4ThreeVector& crossing_pos, G4double& cos_to_surface, G4bool& GoingIn)
{
 for (size_t i=0;i <ListOfSurfaceName.size();i++){
  	if (CrossingAGivenRegisteredSurface(aStep, int(i),crossing_pos,   cos_to_surface, GoingIn)){
		surface_name = ListOfSurfaceName[i];
		return true;
	}
 }
 return false;
}
/////////////////////////////////////////////////////////////////////////////////
//
G4bool G4AdjointCrossSurfChecker::CrossingAnInterfaceBetweenTwoVolumes(const G4Step* aStep,const G4String& vol1_name,const G4String& vol2_name,G4ThreeVector& , G4double& , G4bool& GoingIn)
{
  G4bool step_at_boundary = (aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary);
  G4bool did_cross =false;
  if (step_at_boundary){
  	const G4VTouchable* postStepTouchable = aStep->GetPostStepPoint()->GetTouchable();
	const G4VTouchable* preStepTouchable = aStep->GetPreStepPoint()->GetTouchable();
	if (preStepTouchable && postStepTouchable){
		
		G4String post_vol_name = postStepTouchable->GetVolume()->GetName();
		if (post_vol_name =="") post_vol_name = postStepTouchable->GetVolume()->GetLogicalVolume()->GetName();
		G4String pre_vol_name = preStepTouchable->GetVolume()->GetName();
		if (pre_vol_name =="") pre_vol_name = preStepTouchable->GetVolume()->GetLogicalVolume()->GetName();
		
		
		if ( pre_vol_name == vol1_name && post_vol_name == vol2_name){
			GoingIn=true;
			did_cross=true;
		}
		else if (pre_vol_name == vol2_name && post_vol_name == vol1_name){
			GoingIn=false;
			did_cross=true;
		}
		
	} 
  }
  return did_cross;
  //still need to compute the cosine of the direction
}

/////////////////////////////////////////////////////////////////////////////////
//   
G4bool G4AdjointCrossSurfChecker::AddaSphericalSurface(const G4String& SurfaceName, G4double radius, G4ThreeVector pos, G4double& Area)
{
  G4int ind = FindRegisteredSurface(SurfaceName);
  Area= 4.*pi*radius*radius;
  if (ind>=0) {
  	ListOfSurfaceType[ind]="Sphere";
	ListOfSphereRadius[ind]=radius;
	ListOfSphereCenter[ind]=pos;
   	ListOfVol1Name[ind]="";
   	ListOfVol2Name[ind]="";
	AreaOfSurface[ind]=Area;
  }
  else {
  	ListOfSurfaceName.push_back(SurfaceName);
	ListOfSurfaceType.push_back("Sphere");
	ListOfSphereRadius.push_back(radius);
	ListOfSphereCenter.push_back(pos);
   	ListOfVol1Name.push_back("");
   	ListOfVol2Name.push_back("");
	AreaOfSurface.push_back(Area);
  } 
  return true;
}
/////////////////////////////////////////////////////////////////////////////////
//   
G4bool G4AdjointCrossSurfChecker::AddaSphericalSurfaceWithCenterAtTheCenterOfAVolume(const G4String& SurfaceName, G4double radius, const G4String& volume_name, G4ThreeVector & center, G4double& area)
{ 
 
  G4VPhysicalVolume*  thePhysicalVolume = 0;
  G4PhysicalVolumeStore* thePhysVolStore =G4PhysicalVolumeStore::GetInstance();
  for ( unsigned int i=0; i< thePhysVolStore->size();i++){
	if ((*thePhysVolStore)[i]->GetName() == volume_name){
		thePhysicalVolume = (*thePhysVolStore)[i];
	};
	
  }
  if (thePhysicalVolume){
  	G4VPhysicalVolume* daughter =thePhysicalVolume;
 	G4LogicalVolume* mother = thePhysicalVolume->GetMotherLogical();
 	G4AffineTransform theTransformationFromPhysVolToWorld = G4AffineTransform();
	while (mother){
 		theTransformationFromPhysVolToWorld *=
		G4AffineTransform(daughter->GetFrameRotation(),daughter->GetObjectTranslation());
 		/*G4cout<<"Mother "<<mother->GetName()<<std::endl;
		G4cout<<"Daughter "<<daughter->GetName()<<std::endl;
		G4cout<<daughter->GetObjectTranslation()<<std::endl;
		G4cout<<theTransformationFromPhysVolToWorld.NetTranslation()<<std::endl;*/
		for ( unsigned int i=0; i< thePhysVolStore->size();i++){
			if ((*thePhysVolStore)[i]->GetLogicalVolume() == mother){
				daughter = (*thePhysVolStore)[i];
				mother =daughter->GetMotherLogical();
				break;
			};		
		}
	
 	}
	center = theTransformationFromPhysVolToWorld.NetTranslation();
  	G4cout<<"Center of the spherical surface is at the position: "<<center/cm<<" cm"<<std::endl;
	
  }
  else {
  	G4cout<<"The physical volume with name "<<volume_name<<" does not exist!!"<<std::endl;
	return false;
  }
  return AddaSphericalSurface(SurfaceName, radius, center, area);

}
/////////////////////////////////////////////////////////////////////////////////
//
G4bool G4AdjointCrossSurfChecker::AddanExtSurfaceOfAvolume(const G4String& SurfaceName, const G4String& volume_name, G4double& Area)
{
  G4int ind = FindRegisteredSurface(SurfaceName);

  G4VPhysicalVolume*  thePhysicalVolume = 0;
  G4PhysicalVolumeStore* thePhysVolStore =G4PhysicalVolumeStore::GetInstance();
  for ( unsigned int i=0; i< thePhysVolStore->size();i++){
	if ((*thePhysVolStore)[i]->GetName() == volume_name){
		thePhysicalVolume = (*thePhysVolStore)[i];
	};
	
  }
  if (!thePhysicalVolume){
  	G4cout<<"The physical volume with name "<<volume_name<<" does not exist!!"<<std::endl;
	return false;
  }	
  Area = thePhysicalVolume->GetLogicalVolume()->GetSolid()->GetSurfaceArea();
  G4String mother_vol_name = ""; 
  G4LogicalVolume* theMother = thePhysicalVolume->GetMotherLogical();
 
  if (theMother) mother_vol_name= theMother->GetName();
  if (ind>=0) {
  	ListOfSurfaceType[ind]="ExternalSurfaceOfAVolume";
	ListOfSphereRadius[ind]=0.;
	ListOfSphereCenter[ind]=G4ThreeVector(0.,0.,0.);
   	ListOfVol1Name[ind]=volume_name;
   	ListOfVol2Name[ind]=mother_vol_name;
	AreaOfSurface[ind]=Area;
  }
  else {
  	ListOfSurfaceName.push_back(SurfaceName);
	ListOfSurfaceType.push_back("ExternalSurfaceOfAVolume");
	ListOfSphereRadius.push_back(0.);
	ListOfSphereCenter.push_back(G4ThreeVector(0.,0.,0.));
   	ListOfVol1Name.push_back(volume_name);
   	ListOfVol2Name.push_back(mother_vol_name);
	AreaOfSurface.push_back(Area);
  }
  return true;
}
/////////////////////////////////////////////////////////////////////////////////
//
G4bool G4AdjointCrossSurfChecker::AddanInterfaceBetweenTwoVolumes(const G4String& SurfaceName, const G4String& volume_name1, const G4String& volume_name2,G4double& Area)
{
  G4int ind = FindRegisteredSurface(SurfaceName);
  Area=-1.; //the way to compute the surface is not known yet
  if (ind>=0) {
  	ListOfSurfaceType[ind]="BoundaryBetweenTwoVolumes";
	ListOfSphereRadius[ind]=0.;
	ListOfSphereCenter[ind]=G4ThreeVector(0.,0.,0.);
   	ListOfVol1Name[ind]=volume_name1;
   	ListOfVol2Name[ind]=volume_name2;
	AreaOfSurface[ind]=Area;
	
  }
  else {
  	ListOfSurfaceName.push_back(SurfaceName);
	ListOfSurfaceType.push_back("BoundaryBetweenTwoVolumes");
	ListOfSphereRadius.push_back(0.);
	ListOfSphereCenter.push_back(G4ThreeVector(0.,0.,0.));
   	ListOfVol1Name.push_back(volume_name1);
   	ListOfVol2Name.push_back(volume_name2);
	AreaOfSurface.push_back(Area);
  }
  return true;
}
/////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointCrossSurfChecker::ClearListOfSelectedSurface()
{
  ListOfSurfaceName.clear();
  ListOfSurfaceType.clear();
  ListOfSphereRadius.clear();
  ListOfSphereCenter.clear();
  ListOfVol1Name.clear();
  ListOfVol2Name.clear();
}
/////////////////////////////////////////////////////////////////////////////////
//
G4int G4AdjointCrossSurfChecker::FindRegisteredSurface(const G4String& name)
{
 G4int ind=-1;
 for (size_t i = 0; i<ListOfSurfaceName.size();i++){
 	if (name == ListOfSurfaceName[i]) {
		ind = int (i);
		return ind;
	}
 }
 return ind;
}

