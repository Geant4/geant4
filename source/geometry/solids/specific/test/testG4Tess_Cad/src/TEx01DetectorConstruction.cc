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
//
// $Id: TEx01DetectorConstruction.cc,v 1.2 2006-06-29 18:50:11 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#include "TEx01DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4SDManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4Tokenizer.hh"

#include "G4ios.hh"

#include "G4TriangularFacet.hh"
#include "G4VFacet.hh"

#include <fstream>
#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
TEx01DetectorConstruction::TEx01DetectorConstruction()
:solidWorld(0),  logicWorld(0),  physiWorld(0),
 solidTarget(0), logicTarget(0), physiTarget(0)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
TEx01DetectorConstruction::~TEx01DetectorConstruction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* TEx01DetectorConstruction::Construct()
{
//
//--------- Material definition ---------

  G4double a, z;
  G4double density, temperature, pressure;
  G4int nel;

  //Air
  G4Element* N = new G4Element("Nitrogen", "N", z=7., a= 14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8., a= 16.00*g/mole);
   
  G4Material* Air = new G4Material("Air", density= 1.29*mg/cm3, nel=2);
  Air->AddElement(N, 70*perCent);
  Air->AddElement(O, 30*perCent);

  //Lead
  G4Material* Pb = 
  new G4Material("Lead", z=82., a= 207.19*g/mole, density= 11.35*g/cm3);
    
  // Print all the materials defined.
  //
  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  //------------------------------ 
  // Tesselated Solid
  //------------------------------
  
  G4String solidName, logicName, physiName;
  std::vector <G4TessellatedSolid*> solidStore;
    
  //
  //  G4String geomFile ;
  
  //  char val[100];
  //  std::ostrstream os(val,100);
  //  os <<dirName <<"/z" <<Z <<".a" <<A <<'\0';
  //  G4String file(val);
  std::ifstream GeomFile("./test.geom");
  //solidTarget = new G4Box("target",targetSize,targetSize,targetSize);

  if (!GeomFile){
    G4cerr << "Geometry file not opened !!!" << G4endl;
  }else {
    G4cout << " Tessellated Geometry opened "<<G4endl;
    char inputChars[150]={' '};
    G4String inputLine;
    G4String recordType("");
    G4double a1,a2,a3,a4,a5,a6,a7,a8,a9;
    G4TriangularFacet *facet;
    G4bool firstSolid (true);
    //
    while (!GeomFile.getline(inputChars,150).eof()) {
      inputLine = inputChars;
      inputLine = inputLine.strip(1);
      //      G4cout << inputLine << G4endl;
      if (inputChars[0] == 'f') {
	if (firstSolid) {
	  firstSolid = false;
	}else {
	  solidTarget->SetSolidClosed(true);
	  solidStore.push_back(solidTarget);
	  G4cout << solidName <<" has been created " <<G4endl;
	}
        std::istringstream tmpstream(inputLine);
        tmpstream >>recordType >>solidName ;	
	solidTarget = new G4TessellatedSolid(solidName);		      
      } else if (inputChars[0] == 'p' ) {
        std::istringstream tmpstream(inputLine);
        tmpstream >>recordType >>a1 >> a1 >> a2 >> a3 >> a4 >> a5 >> a6
		  >> a7 >> a8  >> a9  ;
	facet = new
	  G4TriangularFacet (G4ThreeVector(a1,a2,a3),
			     G4ThreeVector(a4,a5,a6),
			     G4ThreeVector(a7,a8,a9),
			     ABSOLUTE);
	solidTarget->AddFacet((G4VFacet*) facet);
      }
    }
    // need to close and store the last solid 
    solidTarget->SetSolidClosed(true);
    solidStore.push_back(solidTarget);
    G4cout << solidName <<" has been created " <<G4endl;
    GeomFile.close();
    G4cout << " Tessellated Geometry file completed "<<G4endl;

  }

//--------- Sizes of the principal geometrical components (solids)  ---------
  G4double fWorld_XMin, fWorld_YMin, fWorld_ZMin;
  fWorld_XMin=fWorld_YMin=fWorld_ZMin = 0.;
  G4double fWorld_XMax, fWorld_YMax, fWorld_ZMax;
  fWorld_XMax=fWorld_YMax=fWorld_ZMax = 0.;
  std::vector<G4TessellatedSolid*>::iterator itr;
  for (itr = solidStore.begin(); itr != solidStore.end(); itr++) {
    if ((*itr)->GetMinXExtent() < fWorld_XMin ) fWorld_XMin = (*itr)->GetMinXExtent();
    if ((*itr)->GetMinYExtent() < fWorld_YMin ) fWorld_YMin = (*itr)->GetMinYExtent();
    if ((*itr)->GetMinZExtent() < fWorld_ZMin ) fWorld_ZMin = (*itr)->GetMinZExtent();
    if ((*itr)->GetMaxXExtent() > fWorld_XMax ) fWorld_XMax = (*itr)->GetMaxXExtent();
    if ((*itr)->GetMaxYExtent() > fWorld_YMax ) fWorld_YMax = (*itr)->GetMaxYExtent();
    if ((*itr)->GetMaxZExtent() > fWorld_ZMin ) fWorld_ZMax = (*itr)->GetMaxZExtent();
  }
  G4ThreeVector Total_Translate = G4ThreeVector( (fWorld_XMin+fWorld_XMax)/2.,
						 (fWorld_YMin+fWorld_YMax)/2.,
						 (fWorld_ZMin+fWorld_ZMax)/2.);
  G4double fWorldHX = (-fWorld_XMin+fWorld_XMax)/2.;
  G4double fWorldHY = (-fWorld_YMin+fWorld_YMax)/2.;
  G4double fWorldHZ = (-fWorld_ZMin+fWorld_ZMax)/2.;
   
//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------
  
  //------------------------------ 
  // World
  //------------------------------ 
  
  solidWorld= new G4Box("world",fWorldHX,fWorldHY,fWorldHZ);
  logicWorld= new G4LogicalVolume( solidWorld, Air, "World", 0, 0, 0);
  
  //  Must place the World Physical volume unrotated at (0,0,0).
  // 
  physiWorld = new G4PVPlacement(0,               // no rotation
                                 G4ThreeVector(), // at (0,0,0)
                                 logicWorld,      // its logical volume
				 "World",         // its name
                                 0,               // its mother  volume
                                 false,           // no boolean operations
                                 0);              // no field specific to volume
  //
  //
  // now the logical and physical volumes
  //
  std::ifstream TreeFile("./test.tree");
  if (!TreeFile){
    G4cerr << "Geometry tree file not opened !!!" << G4endl;
  }else {
    G4cout << " Tessellated Geometry tree file opened "<<G4endl;
    G4String inputLine;
    G4String lsName, copy;
    G4String recordType("");
    G4int level;
    char *sname, *tokenPtr;
    G4double a[17];
    G4RotationMatrix rM;
    //
    while (!std::getline(TreeFile,inputLine).eof()) {
      //      G4cout << inputLine <<G4endl;
      if (inputLine.substr(0,1) == "g") {
	std::istringstream tmpstream(inputLine);
        tmpstream >>recordType >> level >> lsName >>a[0] >>a[1] >> a[2] >>a[3] >>a[4] >>
	  a[5] >>a[6] >>a[7] >>a[8] >>a[9] >>a[10] >>a[11] >>a[12] >>a[13] >>a[14] >>
	  a[15] >> a[16];
	G4cout << inputLine <<G4endl;
	sname       = new char[strlen(lsName)+1];
	strcpy(sname,lsName);
	tokenPtr = strtok(sname,"_");
	lsName = G4String(tokenPtr);
	tokenPtr = strtok( NULL, "_");
	copy = G4String(tokenPtr);
	tokenPtr = strtok( NULL, "_");
	while (tokenPtr != NULL) {
	  lsName += "_";
	  lsName += copy;
	  copy = G4String(tokenPtr);
	  tokenPtr = strtok( NULL, "_");
	}
	// the name
	logicName = "L_"+lsName;
	physiName = "P_"+lsName;    
	// the copy number
	std::istringstream is2(copy);
	is2 >> level;
	//	G4cout << logicName << " " << physiName << " " << level << G4endl;
	//
	std::vector<G4TessellatedSolid*>::iterator itr=solidStore.begin();
	while ((*itr)->GetName() != lsName && itr != solidStore.end()) 
	  itr++ ;
	if (itr != solidStore.end()) {
	  solidTarget = (*itr);
	} else {
	  G4cerr << " the .geom and .tree files don't match!" << G4endl;
	} 
	logicTarget = new G4LogicalVolume(solidTarget,Pb,logicName,0,0,0);

	rM = G4RotationMatrix (G4ThreeVector(a[0],a[1],a[2]),
			       G4ThreeVector(a[4],a[5],a[6]),
			       G4ThreeVector(a[8],a[9],a[10]));
	physiTarget = new G4PVPlacement( G4Transform3D(rM,   //rotation 
			      (G4ThreeVector(a[12], a[13], a[14]) - Total_Translate)),// position
					 logicTarget,     // its logical volume				  
					 physiName,        // its name
					 logicWorld,      // its mother  volume
					 false,           // no boolean operations
					 level-1);        // copy no. 
      }
    }
    TreeFile.close();
    G4cout << " Tessellated Geometry tree file completed "<<G4endl;
  }
//--------- Visualization attributes -------------------------------
  
  G4VisAttributes* BoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  BoxVisAtt->SetVisibility(false);
  G4VisAttributes* FacetVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  FacetVisAtt->SetVisibility(true);
  
  logicWorld  ->SetVisAttributes(BoxVisAtt);  
  G4int i = 1;
  G4double r, g, b;
  for (itr = solidStore.begin(); itr != solidStore.end(); itr++) {
    r = (256-i)/256;
    g= 1.- r;
    b = 1- r*r;
    FacetVisAtt= new G4VisAttributes(G4Colour(r, g, b));
    // works only on logic volume
    //    (*itr) ->SetVisAttributes(FacetVisAtt);
  }  
//--------- example of User Limits -------------------------------

  // below is an example of how to set tracking constraints in a given
  // logical volume(see also in N02PhysicsList how to setup the process
  // G4UserSpecialCuts).  
  // Sets a max Step length in the tracker region
  // G4double maxStep = 0.5*ChamberWidth, maxLength = 2*fTrackerLength;
  // G4double maxTime = 0.1*ns, minEkin = 10*MeV;
  // logicTracker->SetUserLimits(new G4UserLimits(maxStep,maxLength,maxTime,
  //                                               minEkin));
  
  return physiWorld;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
