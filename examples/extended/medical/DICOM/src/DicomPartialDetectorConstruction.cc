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
#include "globals.hh"

#include "DicomPartialDetectorConstruction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "G4PartialPhantomParameterisation.hh"
#include "G4NistManager.hh"
#include "G4UIcommand.hh"
#include "G4Tubs.hh"

#include "G4tgbVolumeMgr.hh"
#include "G4tgbMaterialMgr.hh"
#include "G4tgrMessenger.hh"
#include "G4tgrMessenger.hh"

//---------------------------------------------------------------------------
DicomPartialDetectorConstruction::DicomPartialDetectorConstruction()
{
}

//---------------------------------------------------------------------------
DicomPartialDetectorConstruction::~DicomPartialDetectorConstruction()
{
}


//-------------------------------------------------------------
G4VPhysicalVolume* DicomPartialDetectorConstruction::Construct()
{
  InitialisationOfMaterials();

  //----- Build world
  G4double worldXDimension = 1.*m;
  G4double worldYDimension = 1.*m;
  G4double worldZDimension = 1.*m;

  G4Box* world_solid = new G4Box( "WorldSolid",
                          worldXDimension,
                          worldYDimension,
                          worldZDimension );

  world_logic = new G4LogicalVolume( world_solid, 
                                    air, 
                                    "WorldLogical", 
                                    0, 0, 0 );

  G4VPhysicalVolume* world_phys = new G4PVPlacement( 0,
                                  G4ThreeVector(0,0,0),
                                  "World",
                                  world_logic,
                                  0,
                                  false,
                                  0 );

  ReadPhantomData();

  ConstructPhantomContainer();
  ConstructPhantom();

  return world_phys;
}


//-------------------------------------------------------------
void DicomPartialDetectorConstruction::ReadPhantomData()
{
  std::ifstream finDF("Data.dat");
  G4String fname;
  if(finDF.good() != 1 ) {
    G4Exception(" DicomPartialDetectorConstruction::ReadPhantomData.  Problem reading data file: Data.dat");
  }

  G4int compression;
  finDF >> compression; // not used here

  G4int fNoFiles;
  finDF >> fNoFiles;  // only 1 file supported for the moment
  if( fNoFiles != 1 ) {
    G4Exception("DicomPartialDetectorConstruction::ReadPhantomData",
		"More than 1 DICOM file is not supported",
		FatalErrorInArgument,
		"");
  }
  for(G4int i = 0; i < fNoFiles; i++ ) {    
    finDF >> fname;
    //--- Read one data file
    fname += ".g4pdcm";
    ReadPhantomDataFile(fname);
  }

  finDF.close();

}


//---------------------------------------------------------------------------
void DicomPartialDetectorConstruction::ReadPhantomDataFile(const G4String& fname)
{
  //  G4String filename = "phantom.g4pdcm";
#ifdef G4VERBOSE
  G4cout << " DicomDetectorConstruction::ReadPhantomDataFile opening file " << fname << G4endl;
#endif 
  std::ifstream fin(fname.c_str(), std::ios_base::in);
  if( !fin.is_open() ) {
    G4Exception("DicomPartialDetectorConstruction::ReadPhantomDataFile. File not found " + fname );
  }
 G4int nMaterials;
  fin >> nMaterials;
  G4String stemp;
  G4int nmate;
  for( G4int ii = 0; ii < nMaterials; ii++ ){
    fin >> nmate >> stemp;
    G4cout << "DicomPartialDetectorConstruction::ReadPhantomData reading nmate " << ii << " = " << nmate << " mate " << stemp << G4endl;
    if( ii != nmate ) G4Exception("DicomPartialDetectorConstruction::ReadPhantomData material number should be in increasing order: wrong material number ");
  }

  fin >> nVoxelX >> nVoxelY >> nVoxelZ;
  G4cout << "DicomPartialDetectorConstruction::ReadPhantomData nVoxel X/Y/Z " << nVoxelX << " " << nVoxelY << " " << nVoxelZ << G4endl;
  fin >> offsetX >> dimX;
  dimX = (dimX - offsetX)/nVoxelX;
  fin >> offsetY >> dimY;
  dimY = (dimY - offsetY)/nVoxelY;
  fin >> offsetZ >> dimZ;
  dimZ = (dimZ - offsetZ)/nVoxelZ;
  G4cout << "DicomPartialDetectorConstruction::ReadPhantomData voxelDimX " << dimX << " offsetX " << offsetX << G4endl;
  G4cout << "DicomPartialDetectorConstruction::ReadPhantomData voxelDimY " << dimY << " offsetY " << offsetY << G4endl;
  G4cout << "DicomPartialDetectorConstruction::ReadPhantomData voxelDimZ " << dimZ << " offsetZ " << offsetZ << G4endl;

  //--- Read voxels that are filled
  theNVoxels = 0;
  //  G4bool* isFilled = new G4bool[nVoxelX*nVoxelY*nVoxelZ];
  //  theFilledIDs = new size_t[nVoxelZ*nVoxelY+1];
  //?  theFilledIDs.insert(0);
  int ifxmin1, ifxmax1;
  for( G4int iz = 0; iz < nVoxelZ; iz++ ) {
    std::map< G4int, G4int > ifmin, ifmax;
    for( G4int iy = 0; iy < nVoxelY; iy++ ) {
      fin >> ifxmin1 >> ifxmax1;      
      // check coherence ...

      ifmin[iy] = ifxmin1;
      ifmax[iy] = ifxmax1;
      G4int ifxdiff = ifxmax1-ifxmin1+1;
      if( ifxmax1 == -1 && ifxmin1 == -1 ) ifxdiff = 0;
      //      theFilledIDs.insert(std::multimap<G4int,size_t>::value_type(ifxdiff+theNVoxels, ifxmin1));
      theFilledIDs.insert(std::pair<size_t,size_t>(ifxdiff+theNVoxels-1, ifxmin1));
      G4cout << "DicomPartialDetectorConstruction::ReadPhantomData insert FilledIDs " << ifxdiff+theNVoxels-1 << " min " << ifxmin1 << " N= " << theFilledIDs.size() << G4endl;
      //filledIDs[iz*nVoxelY+iy+1] = ifxmax1-ifxmin1+1;
      for( G4int ix = 0; ix < nVoxelX; ix++ ) {
	//	G4cout << "REGNAV DicomPartialDetectorConstruction::ReadPhantomData checking to add voxel iz " << iz << " iy " << iy << " ix " << ix << "  min= " << ifxmin1 << " max= " << ifxmax1 << G4endl;
	//	G4int copyNo = ix + (iy)*nVoxelX + (iz)*nVoxelX*nVoxelY;
	if( ix >= G4int(ifxmin1) && ix <= G4int(ifxmax1) ) {
	  //	  isFilled[copyNo] = true;
	  theNVoxels++;
	} else {
	  //	  isFilled[copyNo] = false;
	}
      }
    }
    theFilledMins[iz] = ifmin;
    theFilledMaxs[iz] = ifmax;
  }

  //--- Read material IDs
  G4int mateID1;
  mateIDs = new size_t[nVoxelX*nVoxelY*nVoxelZ];
  G4int copyNo = 0;
  for( G4int iz = 0; iz < nVoxelZ; iz++ ) {
    std::map< G4int, G4int > ifmin = theFilledMins[iz];
    std::map< G4int, G4int > ifmax = theFilledMaxs[iz];
    for( G4int iy = 0; iy < nVoxelY; iy++ ) {
      ifxmin1 = ifmin[iy];
      ifxmax1 = ifmax[iy];
      for( G4int ix = 0; ix < nVoxelX; ix++ ) {
	if( ix >= G4int(ifxmin1) && ix <= G4int(ifxmax1) ) {
	  fin >> mateID1;
	  if( mateID1 < 0 || mateID1 >= nMaterials ) {
	    G4Exception("DicomPartialDetectorConstruction::ReadPhantomData",
			"Wrong index in phantom file",
			FatalException,
			G4String("It should be between 0 and "
				 + G4UIcommand::ConvertToString(nMaterials-1) 
				 + ", while it is " 
				 + G4UIcommand::ConvertToString(mateID1)).c_str());
	  }
	  mateIDs[copyNo] = mateID1;
	  copyNo++;
	}
      }
    }
  }

  ReadVoxelDensitiesPartial( fin );

  fin.close();
}


//------------------------------------------------------------------------
void DicomPartialDetectorConstruction::ReadVoxelDensitiesPartial( std::ifstream& fin )
{
  std::map<G4int, std::pair<G4double,G4double> > densiMinMax;
  std::map<G4int, std::pair<G4double,G4double> >::iterator mpite;
  for( size_t ii = 0; ii < fOriginalMaterials.size(); ii++ ){
    densiMinMax[ii] = std::pair<G4double,G4double>(DBL_MAX,-DBL_MAX);
  }

  //----- Define density differences (maximum density difference to create a new material)
  char* part = getenv( "DICOM_CHANGE_MATERIAL_DENSITY" );
  G4double densityDiff = -1.;
  if( part ) densityDiff = G4UIcommand::ConvertToDouble(part);

  std::map<G4int,G4double> densityDiffs;
  if( densityDiff != -1. ) {
    for( size_t ii = 0; ii < fOriginalMaterials.size(); ii++ ){
      densityDiffs[ii] = densityDiff; //currently all materials with same step
    }
  }else {
    if( fMaterials.size() == 0 ) { // do it only for first slice
      for( unsigned int ii = 0; ii < fOriginalMaterials.size(); ii++ ){
	fMaterials.push_back( fOriginalMaterials[ii] );
      }
    }
  }

  //  densityDiffs[0] = 0.0001; //air

  //--- Calculate the average material density for each material/density bin
  std::map< std::pair<G4Material*,G4int>, matInfo* > newMateDens;
 
  G4double dens1;
  size_t ifxmin1, ifxmax1;

  //---- Read the material densities
  G4int copyNo = 0;
  for( G4int iz = 0; iz < nVoxelZ; iz++ ) {
    std::map< G4int, G4int > ifmin = theFilledMins[iz];
    std::map< G4int, G4int > ifmax = theFilledMaxs[iz];
    for( G4int iy = 0; iy < nVoxelY; iy++ ) {
      ifxmin1 = ifmin[iy];
      ifxmax1 = ifmax[iy];
      for( G4int ix = 0; ix < nVoxelX; ix++ ) {
	if( ix >= G4int(ifxmin1) && ix <= G4int(ifxmax1) ) {
	  fin >> dens1;
	  //	G4cout << ix << " " << iy << " " << iz << " filling mateIDs " << copyNo << " = " <<  atoi(stemp.c_str())-1 << " " << stemp << G4endl;
	  
	//--- store the minimum and maximum density for each material (just for printing)
	  mpite = densiMinMax.find( mateIDs[copyNo] );
	  if( dens1 < (*mpite).second.first ) (*mpite).second.first = dens1;
	  if( dens1 > (*mpite).second.second ) (*mpite).second.second = dens1;
	  
	  //--- Get material from original list of material in file
	  int mateID = mateIDs[copyNo];
	  G4Material* mate = fOriginalMaterials[mateID];
	  //	G4cout << copyNo << " mateID " << mateID << G4endl;
	  //--- Check if density is equal to the original material density
	  if( std::fabs(dens1 - mate->GetDensity()/g*cm3 ) < 1.e-9 ) continue;
	  
	  //--- Build material name with fOriginalMaterials name + density
	  //	float densityBin = densityDiffs[mateID] * (G4int(dens1/densityDiffs[mateID])+0.5);
	  G4String newMateName = mate->GetName();

	  G4int densityBin = 0;
	  //	G4cout << " densityBin " << densityBin << " " << dens1 << " " <<densityDiffs[mateID] << G4endl; 	  

	  if( densityDiff != -1.) {
	    densityBin = (G4int(dens1/densityDiffs[mateID]));
	    newMateName = mate->GetName()+G4UIcommand::ConvertToString(densityBin);
	  }

	  //--- Look if it is the first voxel with this material/densityBin
	  std::pair<G4Material*,G4int> matdens(mate, densityBin );
	  
	  std::map< std::pair<G4Material*,G4int>, matInfo* >::iterator mppite = newMateDens.find( matdens );
	  if( mppite != newMateDens.end() ){
	    matInfo* mi = (*mppite).second;
	    mi->sumdens += dens1;
	    mi->nvoxels++;
	    mateIDs[copyNo] = fOriginalMaterials.size()-1 + mi->id;
	    //	  G4cout << copyNo << " mat new again " << fOriginalMaterials.size()-1 + mi->id << " " << mi->id << G4endl;
	  } else {
	    matInfo* mi = new matInfo;
	    mi->sumdens = dens1;
	    mi->nvoxels = 1;
	    mi->id = newMateDens.size()+1;
	    newMateDens[matdens] = mi;
	    mateIDs[copyNo] = fOriginalMaterials.size()-1 + mi->id;
	    //	  G4cout << copyNo << " mat new first " << fOriginalMaterials.size()-1 + mi->id << G4endl;
	  }
	  copyNo++;
	//	G4cout << ix << " " << iy << " " << iz << " filling mateIDs " << copyNo << " = " << atoi(cid)-1 << G4endl;
				      //	mateIDs[copyNo] = atoi(cid)-1;
	}
      }
    }
  }

  //----- Build the list of phantom materials that go to Parameterisation
  //--- Add original materials
  std::vector<G4Material*>::const_iterator mimite;
  for( mimite = fOriginalMaterials.begin(); mimite != fOriginalMaterials.end(); mimite++ ){
    thePhantomMaterials.push_back( (*mimite) );
  }
  // 
  //---- Build and add new materials
  std::map< std::pair<G4Material*,G4int>, matInfo* >::iterator mppite;
  for( mppite= newMateDens.begin(); mppite != newMateDens.end(); mppite++ ){
    G4double averdens = (*mppite).second->sumdens/(*mppite).second->nvoxels;
    G4double saverdens = G4int(1000.001*averdens)/1000.;
    G4cout << "GmReadPhantomGeometry::ReadVoxelDensities AVER DENS " << averdens << " -> " << saverdens << " -> " << G4int(1000*averdens) << " " << G4int(1000*averdens)/1000 << " " <<  G4int(1000*averdens)/1000. << G4endl;

    G4cout << "GmReadPhantomGeometry::ReadVoxelDensities AVER DENS " << averdens << " -> " << saverdens << " -> " << G4UIcommand::ConvertToString(saverdens) << " Nvoxels " <<(*mppite).second->nvoxels << G4endl;
    G4String mateName = ((*mppite).first).first->GetName() + "_" + G4UIcommand::ConvertToString(saverdens);
    thePhantomMaterials.push_back( BuildMaterialWithChangingDensity( (*mppite).first.first, averdens, mateName ) );
  }

}


//-----------------------------------------------------------------------
void DicomPartialDetectorConstruction::ConstructPhantomContainer()
{
  //build a clinder that encompass the partial phantom built in the example
  // /dicom/intersectWithUserVolume 0. 0. 0. 0.*deg 0. 0. TUBE 0. 200. 100.
  //---- Extract number of voxels and voxel dimensions

  G4String contname= "PhantomContainer";

  //----- Define the volume that contains all the voxels
  G4Tubs* container_solid = new G4Tubs(contname,0., 200., 100., 0., 360*deg);

  container_logic = 
    new G4LogicalVolume( container_solid, 
			 air,  
			 contname, 
			 0, 0, 0 );

  G4ThreeVector posCentreVoxels(0.,0.,0.);
  //G4cout << " PhantomContainer posCentreVoxels " << posCentreVoxels << G4endl;
  G4RotationMatrix* rotm = new G4RotationMatrix;

  container_phys = 
    new G4PVPlacement(rotm,  // rotation
		      posCentreVoxels,
		      container_logic,     // The logic volume
		      contname,  // Name
		      world_logic,  // Mother
		      false,           // No op. bool.
		      1);              // Copy number

}


//---------------------------------------------------------------------------
void DicomPartialDetectorConstruction::ConstructPhantom()
{
  G4String OptimAxis = "kXAxis";
  G4bool bSkipEqualMaterials = 0;
  G4int regStructureID = 2;

  G4ThreeVector posCentreVoxels(0.,0.,0.);

  //----- Build phantom
  G4String voxelName = "phantom";
  G4Box* phantom_solid = new G4Box(voxelName, dimX/2., dimY/2., dimZ/2. );
  G4LogicalVolume* phantom_logic = 
    new G4LogicalVolume( phantom_solid, 
			 air,
			 voxelName, 
			 0, 0, 0 );
  G4bool pVis = 0;
  if( !pVis ) {
    G4VisAttributes* visAtt = new G4VisAttributes( FALSE );
    phantom_logic->SetVisAttributes( visAtt );
  }

  G4double theSmartless = 0.5;
  //  container_logic->SetSmartless( theSmartless );
  phantom_logic->SetSmartless( theSmartless );

  thePartialPhantomParam = new G4PartialPhantomParameterisation();

  thePartialPhantomParam->SetMaterials( thePhantomMaterials );
  thePartialPhantomParam->SetVoxelDimensions( dimX/2., dimY/2., dimZ/2. );
  thePartialPhantomParam->SetNoVoxel( nVoxelX, nVoxelY, nVoxelZ );
  thePartialPhantomParam->SetMaterialIndices( mateIDs );

  thePartialPhantomParam->SetFilledIDs(theFilledIDs);

  thePartialPhantomParam->SetFilledMins(theFilledMins);

  thePartialPhantomParam->BuildContainerWalls();

  //  G4cout << " Number of Materials " << thePhantomMaterials.size() << G4endl;
  //  G4cout << " SetMaterialIndices(0) " << mateIDs[0] << G4endl;

  G4PVParameterised * phantom_phys = 0;
  if( OptimAxis == "kUndefined" ) {
    phantom_phys = new G4PVParameterised(voxelName,phantom_logic,container_logic,
					 kUndefined, theNVoxels, thePartialPhantomParam);
  } else   if( OptimAxis == "kXAxis" ) {
    //    G4cout << " optim kX " << G4endl;
    phantom_phys = new G4PVParameterised(voxelName,phantom_logic,container_logic,
					 //					  kXAxis, nVoxelX*nVoxelY*nVoxelZ, thePartialPhantomParam);
					  kXAxis, theNVoxels, thePartialPhantomParam);
    thePartialPhantomParam->SetSkipEqualMaterials(bSkipEqualMaterials);
  } else {
    G4Exception("GmReadPhantomGeometry::ConstructPhantom","Wrong argument to parameter GmReadPhantomGeometry:Phantom:OptimAxis",FatalErrorInArgument,(G4String("Only allowed 'kUndefined' or 'kXAxis', it is: "+OptimAxis).c_str()));
  }

  phantom_phys->SetRegularStructureId(regStructureID);  
  
    //    G4cout << " GmReadPhantomGeometry  phantom_phys " << phantom_phys << " trans " << phantom_phys->GetTranslation() << G4endl;
  //-  thePartialPhantomParam->BuildContainerSolid(cont_phys);
  
  //  G4cout << " GmReadPhantomGeometry::constructPhantom ended " << G4endl;
  
}
