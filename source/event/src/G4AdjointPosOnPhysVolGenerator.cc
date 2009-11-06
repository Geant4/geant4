#include "G4AdjointPosOnPhysVolGenerator.hh"
#include "G4VSolid.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "Randomize.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"

G4AdjointPosOnPhysVolGenerator* G4AdjointPosOnPhysVolGenerator::theInstance = 0;

////////////////////////////////////////////////////
//
G4AdjointPosOnPhysVolGenerator* G4AdjointPosOnPhysVolGenerator::GetInstance()
{
  if(theInstance == 0) {
    static G4AdjointPosOnPhysVolGenerator manager;
    theInstance = &manager;
  }
  return theInstance;
}


////////////////////////////////////////////////////
//
G4AdjointPosOnPhysVolGenerator::~G4AdjointPosOnPhysVolGenerator()
{ 
}

////////////////////////////////////////////////////
//
G4AdjointPosOnPhysVolGenerator::G4AdjointPosOnPhysVolGenerator()
{ 
  theSolid=0;
  NStat =1000000;
  epsilon=0.001;
  ModelOfSurfaceSource = "OnSolid"; //OnSolid, ExternalSphere, ExternalBox
  thePhysicalVolume = 0;
  theTransformationFromPhysVolToWorld = G4AffineTransform();
  UseSphere =true;
}

/////////////////////////////////////////////////////////////////////////////////////////
//
G4VPhysicalVolume* G4AdjointPosOnPhysVolGenerator::DefinePhysicalVolume(G4String aName)
{
  thePhysicalVolume = 0;
  theSolid =0;
  G4PhysicalVolumeStore* thePhysVolStore =G4PhysicalVolumeStore::GetInstance();
  for ( unsigned int i=0; i< thePhysVolStore->size();i++){
  	G4String vol_name =(*thePhysVolStore)[i]->GetName();
	if (vol_name == ""){
		vol_name = (*thePhysVolStore)[i]->GetLogicalVolume()->GetName();
	}
	if (vol_name == aName){
		thePhysicalVolume = (*thePhysVolStore)[i];
	};
	
  }
  if (thePhysicalVolume){
  	theSolid = thePhysicalVolume->GetLogicalVolume()->GetSolid();
	ComputeTransformationFromPhysVolToWorld();
	/*AreaOfExternalSurfaceOfThePhysicalVolume=ComputeAreaOfExternalSurface(1.e-3);
	G4cout<<"Monte Carlo  Estimate of the  area of the external surface :"<<AreaOfExternalSurfaceOfThePhysicalVolume/m/m<<" m2"<<std::endl;*/
	
  }
  else {
  	G4cout<<"The physical volume with name "<<aName<<" does not exist!!"<<std::endl;
	G4cout<<"Before generating a source on an external surface of a volume you should select another physical volume"<<std::endl; 
  }
  return thePhysicalVolume;
}
/////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPosOnPhysVolGenerator::DefinePhysicalVolume1(G4String aName)
{  thePhysicalVolume = DefinePhysicalVolume(aName);
  
}
////////////////////////////////////////////////////
//
double G4AdjointPosOnPhysVolGenerator::ComputeAreaOfExternalSurface()
{ 

   return ComputeAreaOfExternalSurface(theSolid); 
}
////////////////////////////////////////////////////
//
double G4AdjointPosOnPhysVolGenerator::ComputeAreaOfExternalSurface(G4int NStat)
{  return ComputeAreaOfExternalSurface(theSolid,NStat); 
}

////////////////////////////////////////////////////
//
double G4AdjointPosOnPhysVolGenerator::ComputeAreaOfExternalSurface(double epsilon)
{   return ComputeAreaOfExternalSurface(theSolid,epsilon); 
}
////////////////////////////////////////////////////
//
double G4AdjointPosOnPhysVolGenerator::ComputeAreaOfExternalSurface(G4VSolid* aSolid)
{  return ComputeAreaOfExternalSurface(aSolid,1.e-3); 
}
////////////////////////////////////////////////////
//
double G4AdjointPosOnPhysVolGenerator::ComputeAreaOfExternalSurface(G4VSolid* aSolid,G4int NStat)
{ 


  if (ModelOfSurfaceSource == "OnSolid" ){
	if (UseSphere){
		return ComputeAreaOfExternalSurfaceStartingFromSphere(aSolid,NStat);
	
	}
	else {
		return ComputeAreaOfExternalSurfaceStartingFromBox(aSolid,NStat);
		
	}
  }
  else {
  	G4ThreeVector p,dir;
	if (ModelOfSurfaceSource == "ExternalSphere" ) return GenerateAPositionOnASphereBoundary(aSolid, p,dir);
  	return GenerateAPositionOnABoxBoundary(aSolid, p,dir);
 }
  
}
////////////////////////////////////////////////////
//
double G4AdjointPosOnPhysVolGenerator::ComputeAreaOfExternalSurface(G4VSolid* aSolid,double epsilon)
{ int Nstat = int (1./(epsilon*epsilon));
  return ComputeAreaOfExternalSurface(aSolid,Nstat);
}

////////////////////////////////////////////////////
void G4AdjointPosOnPhysVolGenerator::GenerateAPositionOnTheExternalSurfaceOfASolid(G4VSolid* aSolid,G4ThreeVector& p, G4ThreeVector& direction)
{ G4double area;
  area =1.;
  if (ModelOfSurfaceSource == "OnSolid" ){
	return GenerateAPositionOnASolidBoundary(aSolid, p,direction);
  }
  if (ModelOfSurfaceSource == "ExternalSphere" ) {
  	area = GenerateAPositionOnASphereBoundary(aSolid, p, direction);
  	return;
  }	
  	area = GenerateAPositionOnABoxBoundary(aSolid, p, direction);
	return;

}
////////////////////////////////////////////////////
void G4AdjointPosOnPhysVolGenerator::GenerateAPositionOnTheExternalSurfaceOfTheSolid(G4ThreeVector& p, G4ThreeVector& direction)
{ GenerateAPositionOnTheExternalSurfaceOfASolid(theSolid,p,direction);

}
////////////////////////////////////////////////////
//
double G4AdjointPosOnPhysVolGenerator::ComputeAreaOfExternalSurfaceStartingFromBox(G4VSolid* aSolid,G4int Nstat)
{ G4double area=1.;
  int i=0;
  int j=0;
  while (i<Nstat){
  	G4ThreeVector p, direction;
  	area = GenerateAPositionOnABoxBoundary( aSolid,p, direction);
	G4double dist_to_in = aSolid->DistanceToIn(p,direction);
	if (dist_to_in<kInfinity/2.) i++;
	j++;
 }
 area=area*double(i)/double(j);
 return area;
}
/////////////////////////////////////////////////////////////////////////////////////////
//
double G4AdjointPosOnPhysVolGenerator::ComputeAreaOfExternalSurfaceStartingFromSphere(G4VSolid* aSolid,G4int Nstat)
{ G4double area=1.;
  int i=0;
  int j=0;
  while (i<Nstat){
  	G4ThreeVector p, direction;
  	area = GenerateAPositionOnASphereBoundary( aSolid,p, direction);
	G4double dist_to_in = aSolid->DistanceToIn(p,direction);
	if (dist_to_in<kInfinity/2.) i++;
	j++;
 }
 area=area*double(i)/double(j);
 
 return area;
}
/////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPosOnPhysVolGenerator::GenerateAPositionOnASolidBoundary(G4VSolid* aSolid,G4ThreeVector& p, G4ThreeVector&  direction)
{ 
  bool find_pos =false;
  double area=1.;
  while (!find_pos){
  	 if (UseSphere) area = GenerateAPositionOnASphereBoundary( aSolid,p, direction);
  	 else  area = GenerateAPositionOnABoxBoundary( aSolid,p, direction);
	 G4double dist_to_in = aSolid->DistanceToIn(p,direction);
	 if (dist_to_in<kInfinity/2.) {
		find_pos =true;
		G4ThreeVector p1=p+ 0.99999*direction*dist_to_in;
		G4ThreeVector norm =aSolid->SurfaceNormal(p1);
		p+= 0.999999*direction*dist_to_in;
		CosThDirComparedToNormal=direction.dot(-norm);
		//std::cout<<CosThDirComparedToNormal<<std::endl;
		
		return;
	} 

	
  }
}
/////////////////////////////////////////////////////////////////////////////////////////
//
double G4AdjointPosOnPhysVolGenerator::GenerateAPositionOnASphereBoundary(G4VSolid* aSolid,G4ThreeVector& p, G4ThreeVector&  direction)
{ 

  G4double minX,maxX,minY,maxY,minZ,maxZ;
  G4bool yesno;
 

  // values needed for CalculateExtent signature

  G4VoxelLimits limit;                // Unlimited
  G4AffineTransform origin;

  // min max extents of pSolid along X,Y,Z

  yesno = aSolid->CalculateExtent(kXAxis,limit,origin,minX,maxX);
  yesno = aSolid->CalculateExtent(kYAxis,limit,origin,minY,maxY);
  yesno = aSolid->CalculateExtent(kZAxis,limit,origin,minZ,maxZ);
  

  
  G4ThreeVector center = G4ThreeVector((minX+maxX)/2.,(minY+maxY)/2.,(minZ+maxZ)/2.);
  
  double dX=(maxX-minX)/2.;
  double dY=(maxY-minY)/2.;
  double dZ=(maxZ-minZ)/2.;
  double scale=1.01;
  double r=scale*std::sqrt(dX*dX+dY*dY+dZ*dZ);
  
 
  G4double cos_th2 = G4UniformRand();
  G4double theta = std::acos(std::sqrt(cos_th2));
  G4double phi=G4UniformRand()*3.1415926*2;
  direction.setRThetaPhi(1.,theta,phi);
  direction=-direction;
  G4double cos_th = (1.-2.*G4UniformRand());
  theta = std::acos(cos_th);
  if (G4UniformRand() <0.5) theta=3.1415926-theta;
  phi=G4UniformRand()*3.1415926*2;
  p.setRThetaPhi(r,theta,phi);
  p+=center;
  direction.rotateY(theta);
  direction.rotateZ(phi);
  return 4.*3.1415926*r*r;;
}
/////////////////////////////////////////////////////////////////////////////////////////
//
double G4AdjointPosOnPhysVolGenerator::GenerateAPositionOnABoxBoundary(G4VSolid* aSolid,G4ThreeVector& p, G4ThreeVector&  direction)
{

  G4double ran_var,px,py,pz,minX,maxX,minY,maxY,minZ,maxZ;
  G4bool yesno;
  
  // values needed for CalculateExtent signature

  G4VoxelLimits limit;                // Unlimited
  G4AffineTransform origin;

  // min max extents of pSolid along X,Y,Z

  yesno = aSolid->CalculateExtent(kXAxis,limit,origin,minX,maxX);
  yesno = aSolid->CalculateExtent(kYAxis,limit,origin,minY,maxY);
  yesno = aSolid->CalculateExtent(kZAxis,limit,origin,minZ,maxZ);
  
  double scale=.1;
  minX-=scale*std::abs(minX);
  minY-=scale*std::abs(minY);
  minZ-=scale*std::abs(minZ);
  maxX+=scale*std::abs(maxX);
  maxY+=scale*std::abs(maxY);
  maxZ+=scale*std::abs(maxZ);
  
  double dX=(maxX-minX);
  double dY=(maxY-minY);
  double dZ=(maxZ-minZ);
  
  
  double XY_prob=2.*dX*dY;
  double YZ_prob=2.*dY*dZ;
  double ZX_prob=2.*dZ*dX;
  double area=XY_prob+YZ_prob+ZX_prob;
  XY_prob/=area;
  YZ_prob/=area;
  ZX_prob/=area;
  
  ran_var=G4UniformRand();
  G4double cos_th2 = G4UniformRand();
  G4double sth = std::sqrt(1.-cos_th2);
  G4double cth = std::sqrt(cos_th2);
  G4double phi=G4UniformRand()*3.1415926*2;
  G4double dirX = sth*std::cos(phi);
  G4double dirY = sth*std::sin(phi);
  G4double dirZ = cth;
  if (ran_var <=XY_prob){ //on the XY faces
	double ran_var1=ran_var/XY_prob;
	double ranX=ran_var1;
	if (ran_var1<=0.5){
		pz=minZ;
		direction=G4ThreeVector(dirX,dirY,dirZ);
		ranX=ran_var1*2.;
	} 
	else{
		pz=maxZ;
		direction=-G4ThreeVector(dirX,dirY,dirZ);
		ranX=(ran_var1-0.5)*2.;
	}
	double ranY=G4UniformRand();
	px=minX+(maxX-minX)*ranX;
	py=minY+(maxY-minY)*ranY;
  }
  else if (ran_var <=(XY_prob+YZ_prob)){ //on the YZ faces 
	double ran_var1=(ran_var-XY_prob)/YZ_prob;
	double ranY=ran_var1;
	if (ran_var1<=0.5){
		px=minX;
		direction=G4ThreeVector(dirZ,dirX,dirY);
		ranY=ran_var1*2.;
	} 
  	else{
		px=maxX;
		direction=-G4ThreeVector(dirZ,dirX,dirY);
		ranY=(ran_var1-0.5)*2.;
  	}
	double ranZ=G4UniformRand();
	py=minY+(maxY-minY)*ranY;
	pz=minZ+(maxZ-minZ)*ranZ;
		
  }
  else{ //on the ZX faces 
	double ran_var1=(ran_var-XY_prob-YZ_prob)/ZX_prob;
	double ranZ=ran_var1;
	if (ran_var1<=0.5){
		py=minY;
		direction=G4ThreeVector(dirY,dirZ,dirX);
		ranZ=ran_var1*2.;
	} 
	else{
		py=maxY;
		direction=-G4ThreeVector(dirY,dirZ,dirX);
		ranZ=(ran_var1-0.5)*2.;
	}
	double ranX=G4UniformRand();
	px=minX+(maxX-minX)*ranX;
	pz=minZ+(maxZ-minZ)*ranZ;
  }
  
  p=G4ThreeVector(px,py,pz);
  return area;
  
  
}
/////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPosOnPhysVolGenerator::GenerateAPositionOnTheExternalSurfaceOfThePhysicalVolume(G4ThreeVector& p, G4ThreeVector&  direction)
{ if (!thePhysicalVolume) {
  	G4cout<<"Before generating a source on an external surface of volume you should select a physical volume"<<std::endl; 
  	return;
  };
  GenerateAPositionOnTheExternalSurfaceOfTheSolid(p,direction);
  p = theTransformationFromPhysVolToWorld.TransformPoint(p);
  direction = theTransformationFromPhysVolToWorld.TransformAxis(direction);
}
/////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPosOnPhysVolGenerator::GenerateAPositionOnTheExternalSurfaceOfThePhysicalVolume(G4ThreeVector& p, G4ThreeVector&  direction,
  										G4double& costh_to_normal)
{GenerateAPositionOnTheExternalSurfaceOfThePhysicalVolume(p, direction);
 costh_to_normal = CosThDirComparedToNormal;
}										
/////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPosOnPhysVolGenerator::ComputeTransformationFromPhysVolToWorld()
{G4VPhysicalVolume* daughter =thePhysicalVolume;
 G4LogicalVolume* mother = thePhysicalVolume->GetMotherLogical();
 theTransformationFromPhysVolToWorld = G4AffineTransform();
 G4PhysicalVolumeStore* thePhysVolStore =G4PhysicalVolumeStore::GetInstance();
 while (mother){
 	theTransformationFromPhysVolToWorld *=
	G4AffineTransform(daughter->GetFrameRotation(),daughter->GetObjectTranslation());
 	for ( unsigned int i=0; i< thePhysVolStore->size();i++){
		if ((*thePhysVolStore)[i]->GetLogicalVolume() == mother){
			daughter = (*thePhysVolStore)[i];
			mother =daughter->GetMotherLogical();
			break;
		};		
	}
	
 }
}
