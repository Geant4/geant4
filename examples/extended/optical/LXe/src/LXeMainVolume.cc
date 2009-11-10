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
#include "LXeMainVolume.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "LXePMTSD.hh"
#include "LXeScintSD.hh"

LXeScintSD* LXeMainVolume::scint_SD;
LXePMTSD* LXeMainVolume::pmt_SD;

G4LogicalVolume* LXeMainVolume::housing_log=NULL;

LXeMainVolume::LXeMainVolume(G4RotationMatrix *pRot,
			     const G4ThreeVector &tlate,
			     G4LogicalVolume *pMotherLogical,
			     G4bool pMany,
			     G4int pCopyNo,
			     LXeDetectorConstruction* c)
  //Pass info to the G4PVPlacement constructor
  :G4PVPlacement(pRot,tlate,
		 //Temp logical volume must be created here
		 new G4LogicalVolume(new G4Box("temp",1,1,1),
				     G4Material::GetMaterial("Vacuum"),
				     "temp",0,0,0),
		 "housing",pMotherLogical,pMany,pCopyNo),constructor(c)
{
  CopyValues();

  if(!housing_log || updated){
    
    G4double housing_x=scint_x+d_mtl;
    G4double housing_y=scint_y+d_mtl;
    G4double housing_z=scint_z+d_mtl;
    
    //*************************** housing and scintillator
    scint_box = new G4Box("scint_box",scint_x/2.,scint_y/2.,scint_z/2.);
    housing_box = new G4Box("housing_box",housing_x/2.,housing_y/2.,
			    housing_z/2.);
    
    scint_log = new G4LogicalVolume(scint_box,G4Material::GetMaterial("LXe"),
				    "scint_log",0,0,0);
    housing_log = new G4LogicalVolume(housing_box,
				      G4Material::GetMaterial("Al"),
				      "housing_log",0,0,0);
    

    scint_phys = new G4PVPlacement(0,G4ThreeVector(),scint_log,"scintillator",
				   housing_log,false,0);  
    
    //*************** Miscellaneous sphere to demonstrate skin surfaces
    sphere = new G4Sphere("sphere",0.*mm,2.*cm,0.*deg,360.*deg,0.*deg,
			  360.*deg);
    sphere_log = new G4LogicalVolume(sphere,G4Material::GetMaterial("Al"),
				     "sphere_log");
    if(sphereOn)
      sphere_phys = new G4PVPlacement(0,G4ThreeVector(5.*cm,5.*cm,5.*cm),
				      sphere_log,"sphere",scint_log,false,0);
    
        
    //****************** Build PMTs
    G4double innerRadius_pmt = 0.*cm;
    G4double height_pmt = d_mtl/2.;
    G4double startAngle_pmt = 0.*deg;
    G4double spanningAngle_pmt = 360.*deg;
    
    pmt = new G4Tubs("pmt_tube",innerRadius_pmt,outerRadius_pmt,
		     height_pmt,startAngle_pmt,spanningAngle_pmt);
    
    //the "photocathode" is a metal slab at the back of the glass that
    //is only a very rough approximation of the real thing since it only
    //absorbs or detects the photons based on the efficiency set below
    photocath = new G4Tubs("photocath_tube",innerRadius_pmt,outerRadius_pmt,
			   height_pmt/2,startAngle_pmt,spanningAngle_pmt);
    
    pmt_log = new G4LogicalVolume(pmt,G4Material::GetMaterial("Glass"),
				  "pmt_log");
    photocath_log = new G4LogicalVolume(photocath,
					G4Material::GetMaterial("Al"),
					"photocath_log");
    
    photocath_phys = new G4PVPlacement(0,G4ThreeVector(0,0,-height_pmt/2),
				       photocath_log,"photocath",
				       pmt_log,false,0);
    
    
    
    //***********Arrange pmts around the outside of housing**********
    //---pmt sensitive detector
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    
    if(!pmt_SD){
      pmt_SD = new LXePMTSD("/LXeDet/pmtSD");
      SDman->AddNewDetector(pmt_SD);
      //Created here so it exists as pmts are being placed
    }
    pmt_SD->InitPMTs((nx*ny+nx*nz+ny*nz)*2); //let pmtSD know # of pmts
    //-------
    
    G4double dx = scint_x/nx;
    G4double dy = scint_y/ny;
    G4double dz = scint_z/nz;
    
    G4double x,y,z;
    G4double xmin = -scint_x/2. - dx/2.;
    G4double ymin = -scint_y/2. - dy/2.;
    G4double zmin = -scint_z/2. - dz/2.;
    G4int k=0;
    
    z = -scint_z/2. - height_pmt;      //front
    PlacePMTs(pmt_log,0,x,y,dx,dy,xmin,ymin,nx,ny,x,y,z,k,pmt_SD);
    G4RotationMatrix* rm_z = new G4RotationMatrix();
    rm_z->rotateY(180*deg);
    z = scint_z/2. + height_pmt;       //back
    PlacePMTs(pmt_log,rm_z,x,y,dx,dy,xmin,ymin,nx,ny,x,y,z,k,pmt_SD);
    
    G4RotationMatrix* rm_y1 = new G4RotationMatrix();
    rm_y1->rotateY(-90*deg);
    x = -scint_x/2. - height_pmt;      //left
    PlacePMTs(pmt_log,rm_y1,y,z,dy,dz,ymin,zmin,ny,nz,x,y,z,k,pmt_SD);
    G4RotationMatrix* rm_y2 = new G4RotationMatrix();
    rm_y2->rotateY(90*deg);
    x = scint_x/2. + height_pmt;      //right
    PlacePMTs(pmt_log,rm_y2,y,z,dy,dz,ymin,zmin,ny,nz,x,y,z,k,pmt_SD);
    
    G4RotationMatrix* rm_x1 = new G4RotationMatrix();
    rm_x1->rotateX(90*deg);
    y = -scint_y/2. - height_pmt;     //bottom
    PlacePMTs(pmt_log,rm_x1,x,z,dx,dz,xmin,zmin,nx,nz,x,y,z,k,pmt_SD);
    G4RotationMatrix* rm_x2 = new G4RotationMatrix();
    rm_x2->rotateX(-90*deg);
    y = scint_y/2. + height_pmt;      //top
    PlacePMTs(pmt_log,rm_x2,x,z,dx,dz,xmin,zmin,nx,nz,x,y,z,k,pmt_SD);
    
    //**********Setup Sensitive Detectors***************
    if(!scint_SD){//determine if it has already been created
      scint_SD = new LXeScintSD("/LXeDet/scintSD");
      SDman->AddNewDetector(scint_SD);    
    }
    scint_log->SetSensitiveDetector(scint_SD);
    
    //sensitive detector is not actually on the photocathode.
    //processHits gets done manually by the stepping action.
    //It is used to detect when photons hit and get absorbed&detected at the
    //boundary to the photocathode (which doesnt get done by attaching it to a
    //logical volume.
    //It does however need to be attached to something or else it doesnt get
    //reset at the begining of events
    photocath_log->SetSensitiveDetector(pmt_SD);

    VisAttributes();
    SurfaceProperties();
  }

  SetLogicalVolume(housing_log);
}

void LXeMainVolume::CopyValues(){
  updated=constructor->GetUpdated();

  scint_x=constructor->GetScintX();
  scint_y=constructor->GetScintY();
  scint_z=constructor->GetScintZ();
  d_mtl=constructor->GetHousingThickness();
  nx=constructor->GetNX();
  ny=constructor->GetNY();
  nz=constructor->GetNZ();
  outerRadius_pmt=constructor->GetPMTRadius();
  sphereOn=constructor->GetSphereOn();
  refl=constructor->GetHousingReflectivity();
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void LXeMainVolume::PlacePMTs(G4LogicalVolume* pmt_log,
			      G4RotationMatrix *rot,
			      G4double &a, G4double &b, G4double da,
			      G4double db, G4double amin,
			      G4double bmin, G4int na, G4int nb,
			      G4double &x, G4double &y, G4double &z,
			      G4int &k,LXePMTSD* sd){
/*PlacePMTs : a different way to parameterize placement that does not depend on
  calculating the position from the copy number
  
  pmt_log = logical volume for pmts to be placed
  rot = rotation matrix to apply
  a,b = coordinates to vary(ie. if varying in the xy plane then pass x,y)
  da,db = value to increment a,b by
  amin,bmin = start values for a,b
  na,nb = number of repitions in a and b
  x,y,z = just pass x,y, and z by reference (the same ones passed for a,b)
  k = copy number to start with
  sd = sensitive detector for pmts
*/
  a=amin;
  for(G4int j=1;j<=na;j++){
    a+=da;
    b=bmin;
    for(G4int i=1;i<=nb;i++){
      b+=db;
      new G4PVPlacement(rot,G4ThreeVector(x,y,z),pmt_log,"pmt",
			housing_log,false,k);
      sd->SetPMTPos(k,x,y,z);
      k++;
    }
  }
}

void LXeMainVolume::VisAttributes(){
  G4VisAttributes* housing_va = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
  housing_log->SetVisAttributes(housing_va);

  G4VisAttributes* sphere_va = new G4VisAttributes();
  sphere_va->SetForceSolid(true);
  sphere_log->SetVisAttributes(sphere_va);
}

void LXeMainVolume::SurfaceProperties(){    
  const G4int num = 2;
  G4double Ephoton[num] = {7.0*eV, 7.14*eV};  

  //**Scintillator housing properties 
  G4double Reflectivity[num] = {refl, refl};
  G4double Efficiency[num] = {0.0, 0.0}; 
  G4MaterialPropertiesTable* scintHsngPT = new G4MaterialPropertiesTable(); 
  scintHsngPT->AddProperty("REFLECTIVITY", Ephoton, Reflectivity, num);
  scintHsngPT->AddProperty("EFFICIENCY", Ephoton, Efficiency, num);
  G4OpticalSurface* OpScintHousingSurface =
    new G4OpticalSurface("HousingSurface",unified,polished,dielectric_metal);
  OpScintHousingSurface->SetMaterialPropertiesTable(scintHsngPT);
  
  //**Sphere surface properties
  G4double SphereReflectivity[num] = {1.0, 1.0};
  G4double SphereEfficiency[num] = {0.0, 0.0}; 
  G4MaterialPropertiesTable* spherePT = new G4MaterialPropertiesTable();
  spherePT->AddProperty("REFLECTIVITY", Ephoton, SphereReflectivity, num);
  spherePT->AddProperty("EFFICIENCY", Ephoton, SphereEfficiency, num);
  G4OpticalSurface* OpSphereSurface =
    new G4OpticalSurface("SphereSurface",unified,polished,dielectric_metal);
  OpSphereSurface->SetMaterialPropertiesTable(spherePT);
  
  //**Photocathode surface properties
  G4double photocath_EFF[num]={1.,1.}; //Enables 'detection' of photons
  G4double photocath_ReR[num]={1.92,1.92};
  G4double photocath_ImR[num]={1.69,1.69};
  G4MaterialPropertiesTable* photocath_mt = new G4MaterialPropertiesTable();
  photocath_mt->AddProperty("EFFICIENCY",Ephoton,photocath_EFF,num);
  photocath_mt->AddProperty("REALRINDEX",Ephoton,photocath_ReR,num);
  photocath_mt->AddProperty("IMAGINARYRINDEX",Ephoton,photocath_ImR,num);
  G4OpticalSurface* photocath_opsurf= 
    new G4OpticalSurface("photocath_opsurf",glisur,polished,
			 dielectric_metal);
  photocath_opsurf->SetMaterialPropertiesTable(photocath_mt);


  //**Create logical skin surfaces
  new G4LogicalSkinSurface("photocath_surf",housing_log,
			   OpScintHousingSurface);
  new G4LogicalSkinSurface("sphere_surface",sphere_log,OpSphereSurface);  
  new G4LogicalSkinSurface("photocath_surf",photocath_log,photocath_opsurf);
}



